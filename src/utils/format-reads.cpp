/* format_reads: a program to ensure SAM and BAM format reads are
 * conforming to expectations of dnmtools software
 *
 * Copyright (C) 2020-2023 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew Smith and Guilherme Sena
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

/* The output of this program should include mapped reads that are
 * T-rich, which might require reverse-complementing sequences, and
 * switching their strand, along with any tag that indicates T-rich
 * vs. A-rich.
 *
 * Focusing only on single end reads, for each supported read mapping
 * tool, we require a means of determining whether or not the read is
 * A-rich and then changing that format to indicate T-rich.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <tuple>

#include <config.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "htslib_wrapper.hpp"
#include "sam_record.hpp"
#include "cigar_utils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::runtime_error;
using std::unordered_map;
using std::swap;
using std::to_string;
using std::ostringstream;

// ADS: do we need to check that both mates are on the same strand? Or
// that they are on opposite strands?

static bool
abismal_is_a_rich(const sam_rec &aln) {
  auto the_cv_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "CV:") == 0;
                            });
  if (the_cv_tag != end(aln.tags))
    return the_cv_tag->back() == 'A';
  return false;
}

static bool
is_a_rich(const sam_rec &aln) {
  return abismal_is_a_rich(aln);
}

bool
is_rc(const sam_rec &aln) {
  return check_flag(aln, samflags::read_rc);
}

static void
flip_strand(sam_rec &aln) {
  if (check_flag(aln, samflags::read_rc))
    unset_flag(aln, samflags::read_rc);
  else
    set_flag(aln, samflags::read_rc);
}

static void
flip_conversion(sam_rec &aln) {
  flip_strand(aln); // set strand to opposite of current value
  revcomp_inplace(aln.seq); // reverse complement sequence
  std::reverse(begin(aln.qual), end(aln.qual)); // and quality scores

  // ADS: assuming abismal here
  auto the_cv_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "CV:") == 0;
                            });
  if (the_cv_tag != end(aln.tags)) {
    if (abismal_is_a_rich(aln))
      the_cv_tag->back() = 'T';
    else
      the_cv_tag->back() = 'A';
  }
}


static bool
are_opposite_strands(const sam_rec &one, const sam_rec &two) {
  return (check_flag(one, samflags::template_first) !=
          check_flag(two, samflags::template_first)) &&
         (check_flag(one, samflags::template_last) !=
          check_flag(two, samflags::template_last));
}

static bool
are_mates(const sam_rec &one, const sam_rec &two) {
  return ((one.rnext == "=" && two.rnext == "=") ||
          (one.rnext == two.qname)) &&
           one.pnext == two.pos &&
           two.pnext == one.pos &&
           are_opposite_strands(one, two);
}

size_t
get_edit_distance(const sam_rec &aln) {
  auto the_nm_tag = find_if(begin(aln.tags), end(aln.tags),
                            [](const string &t) {
                              return t.compare (0, 3, "NM:") == 0;
                            });

  // GS: NM:i:xxxx, get xxxx
  assert(the_nm_tag->size() > 5);
  return stoi(the_nm_tag->substr(5));
}


static size_t
merge_mates(const size_t range,
            const sam_rec &one, const sam_rec &two, sam_rec &merged) {

  if (!are_mates(one, two)) {
    return -std::numeric_limits<int>::max();
  }


  // ADS: not sure this can be consistent across mappers
  // GS: true after standardization

  merged = one;
  merged.rnext = "*";

  // arithmetic easier using base 0 so subtracting 1 from pos
  const int one_s = one.pos - 1;
  const int one_e = one_s + cigar_rseq_ops(one.cigar);
  const int two_s = two.pos - 1;
  const int two_e = two_s + cigar_rseq_ops(two.cigar);
  assert(one_s >= 0 && two_s >= 0);

  const int spacer = two_s - one_e;
  if (spacer >= 0) {
    /* fragments longer enough that there is space between them: this
     * size of the spacer ("_") is determined based on the reference
     * positions of the two ends, and here we assume "one" maps to
     * positive genome strand.
     *
     * left                                                         right
     * one_s                    one_e      two_s                    two_e
     * [------------end1------------]______[------------end2------------]
     */
    merged.seq += revcomp(two.seq);
    // ADS: need to take care of soft clipping in between;
    merged.cigar += to_string(spacer) + "N";
    merged.cigar += two.cigar;
    // merged.qual = (one.qual == "*" ? one.qual : one.qual + revcomp(two.qual));
  }
  else {
    const int head = two_s - one_s;
    if (head >= 0) {
      /* fragment longer than or equal to the read length, but shorter
       * than twice the read length: this is determined by obtaining
       * the size of the "head" in the diagram below: the portion of
       * end1 that is not within [=]. If the read maps to the positive
       * strand, this depends on the reference start of end2 minus the
       * reference start of end1. For negative strand, this is
       * reference start of end1 minus reference start of end2.
       *
       * <======= head ================>
       *
       * left                                             right
       * one_s              two_s      one_e              two_e
       * [------------end1------[======]------end2------------]
       */
      truncate_cigar_r(merged.cigar, head);
      const uint32_t one_seq_len = cigar_qseq_ops(merged.cigar);
      merged.cigar += two.cigar;
      merge_equal_neighbor_cigar_ops(merged.cigar);
      // ADS: need to take care of soft clipping in between;
      merged.seq.resize(one_seq_len);
      merged.seq += revcomp(two.seq);
      // if (merged.qual != "*") {
      //   merged.qual.resize(one_seq_len);
      //   merged.qual += two.qual;
      // }
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       * left                                       right
       * two_s            one_s    two_e            one_e
       * [--end2----------[============]----------end1--]
       */
      const int overlap = two_e - one_s;
      if (overlap >= 0) {
        truncate_cigar_r(merged.cigar, overlap);
        const uint32_t overlap_qlen = cigar_qseq_ops(merged.cigar);
        merged.seq.resize(overlap_qlen);
        // if (merged.qual != "*")
        //   merged.qual.resize(overlap_qlen);
      }
    }
  }
  merged.rnext = "*";
  merged.pnext = 0;
  merged.tlen = 0;
  return two_e - one_s;
}

/********Above are functions for merging pair-end reads********/

// ADS: there is a bug somewhere when a value of 0 is given for
// suffix_len
static string
remove_suff(const string &x, const size_t suffix_len) {
  return x.size() > suffix_len ? x.substr(0, x.size() - suffix_len) : x;
}

bool
bsmap_get_rc(const string &strand_tag) {
  return strand_tag.size() > 5 && strand_tag[5] == '-';
}

bool
bsmap_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 6 && richness_tag[6] == '-';
}

bool
bismark_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 5 && richness_tag[5] == 'G';
}

static void
standardize_format(const string &input_format, sam_rec &aln) {

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    auto z_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "ZS:") == 0;
                             });
    if (z_tag_itr == end(aln.tags))
      throw runtime_error("record appears to be invalid for bsmap");
    const string z_tag = *z_tag_itr;
    aln.tags.erase(z_tag_itr);

    aln.add_tag(bsmap_get_a_rich(z_tag) ? "CV:A:A" : "CV:A:T");
    if (is_rc(aln)) revcomp_inplace(aln.seq);
  }

  if (input_format == "bismark") {
    // remove everything after _ in read name
    aln.qname = aln.qname.substr(0, aln.qname.find_first_of("_"));
    auto xr_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "XR:") == 0;
                             }),
         nm_tag_itr = find_if(begin(aln.tags), end(aln.tags),
                             [](const string &t) {
                               return t.compare (0, 3, "NM:") == 0;
                             });

    if (xr_tag_itr == end(aln.tags))
      throw runtime_error("record appears to be invalid for bismark");

    const string xr_tag = *xr_tag_itr,
                 nm_tag = *nm_tag_itr;

    aln.tags.clear();
    aln.add_tag(nm_tag);
    aln.add_tag(bismark_get_a_rich(xr_tag) ? "CV:A:A" : "CV:A:T");

    if (is_rc(aln)) revcomp_inplace(aln.seq);
  }

  // doesn't depend on mapper
  aln.qual = "*";
}


static vector<string>
load_read_names(const string &mapped_reads_file, const size_t n_reads) {
  SAMReader sam_reader(mapped_reads_file);
  sam_rec aln;
  vector<string> names;
  size_t count = 0;
  while (sam_reader >> aln && count++ < n_reads)
    names.push_back(std::move(aln.qname));
  return names;
}


static size_t
get_max_repeat_count(const vector<string> &names, const size_t suff_len) {
  // assume "suff_len" is shorter than the shortest entry in "names"
  size_t repeat_count = 0;
  size_t tmp_repeat_count = 0;
  for (size_t i = 1; i < names.size() && repeat_count < 1; ++i) {
    if (names[i-1].size() == names[i].size() &&
        equal(begin(names[i-1]), end(names[i-1]) - suff_len,
              begin(names[i]), end(names[i]) - suff_len)) {
      ++tmp_repeat_count;
    }
    else tmp_repeat_count = 0;
    repeat_count = max(repeat_count, tmp_repeat_count);
  }
  return repeat_count;
}


static bool
check_suffix_length(const string &mapped_reads_file, const size_t suff_len,
                    const size_t n_names_to_check) {
  auto names(load_read_names(mapped_reads_file, n_names_to_check));
  // get the minimum read name length
  size_t min_name_len = std::numeric_limits<size_t>::max();
  for (auto &&i: names)
    min_name_len = min(min_name_len, i.size());
  if (min_name_len <= suff_len)
    throw runtime_error("given suffix length exceeds min read name length" );
  sort(begin(names), end(names));
  return get_max_repeat_count(names, suff_len) < 2;
}


static size_t
guess_suffix_length(const string &mapped_reads_file,
                    const size_t n_names_to_check,
                    size_t &repeat_count) {

  // ADS: assuming copy elision but should test it
  auto names(load_read_names(mapped_reads_file, n_names_to_check));

  // get the minimum read name length
  size_t min_name_len = std::numeric_limits<size_t>::max();
  for (auto &&i: names)
    min_name_len = min(min_name_len, i.size());
  assert(min_name_len > 0);

  sort(begin(names), end(names));

  // these will be returned
  size_t suff_len = 0;
  repeat_count = 0;

  // check the possible read name suffix lengths; if any causes a
  // repeat count of more than 1 (here this means == 2), all greater
  // suffix lengths will also
  const size_t max_suff_len = min_name_len - 1;
  while (suff_len < max_suff_len && repeat_count == 0) {
    // check current suffix length guess
    repeat_count = get_max_repeat_count(names, suff_len);
    // we want to lag by one iteration
    if (repeat_count == 0)
      ++suff_len;
  }
  // repeat_count should be equal to the greatest value of
  // tmp_repeat_count over all inner iterations above. If this value
  // is not 1, it will indicate whether we have exactly hit a max
  // repeat count of 2, indicating mates, or exceeded it, indicating
  // there seems not to be a good suffix length to remove for
  // identifying mates

  return suff_len;
}


static bool
check_mates_consecutive(const string &mapped_reads_file,
                        const size_t suff_len, size_t n_reads) {
  // In order to check if mates are consecutive we need to check if a
  // given end has a mate and that mate is not adjacent. This requires
  // storing previous reads, not simply checking for adjacent pairs.

  auto names(load_read_names(mapped_reads_file, n_reads));
  for (auto &&i : names)
    i = remove_suff(i, suff_len);

  unordered_map<string, size_t> mate_lookup;
  for (size_t i = 0; i < names.size(); ++i) {
    auto the_mate = mate_lookup.find(names[i]);
    if (the_mate == end(mate_lookup)) // 1st time seeing this one
      mate_lookup[names[i]] = i;
    else if (the_mate->second != i - 1)
      return false;
  }

  // made it here: all reads with mates are consecutive
  return true;
}


static bool
check_input_file(const string &input_filename) {
  // because this isn't so convenient with our sam file wrapper
  std::ifstream in(input_filename);
  if (!in)
    return false;
  return true;
}


static bool
check_format_matches_sam_header(string input_format,
                                const string &mapped_reads_file) {
  SAMReader sam_reader(mapped_reads_file);
  string header = sam_reader.get_header();
  std::transform(begin(header), end(header), begin(header),
                 [](unsigned char c){return std::tolower(c);});
  std::transform(begin(input_format), end(input_format), begin(input_format),
                 [](unsigned char c){return std::tolower(c);});
  return header.find(input_format) != std::string::npos;
}


int
main_format_reads(int argc, const char **argv) {

  try {

    size_t n_reads_to_check = 1000000;

    string outfile;
    string input_format;
    int max_frag_len = std::numeric_limits<int>::max();
    size_t suff_len = 0;
    bool single_end = false;
    bool VERBOSE = false;

    const string description = "convert SAM/BAM mapped bs-seq reads "
      "to standard dnmtools format";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<sam/bam-file>", 1);
    opt_parse.add_opt("format", 'f', "input format (abismal, bsmap, bismark)",
                      false, input_format);
    opt_parse.add_opt("output", 'o', "output file name",
                      false, outfile);
    opt_parse.add_opt("suff", 's', "read name suffix length",
                      false, suff_len);
    opt_parse.add_opt("single-end", '\0',
                      "assume single-end [do not use with -suff]",
                      false, single_end);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("check", 'c',
                      "check this many reads to validate read name suffix",
                      false, n_reads_to_check);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc < 3 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_FAILURE;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (suff_len != 0 && single_end) {
      cerr << "incompatible arguments specified" << endl
           << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE) {
      cerr << "[input file: " << mapped_reads_file << "]" << endl;
      if (single_end)
        cerr << "[assuming reads are single-end]" << endl;
      cerr << "[output file: "
           << (outfile.empty() ? "stdout" : outfile) << "]" << endl;
    }

    // prep the output stream and make sure it's valid
    std::ofstream of;
    if (!outfile.empty()) {
      of.open(outfile);
      if (!of)
        throw runtime_error("failed opening output file: " + outfile);
    }
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    // check the input and parameters
    if (!check_input_file(mapped_reads_file))
      throw runtime_error("failed to open input file: " + mapped_reads_file);

    if (VERBOSE)
      if (!check_format_matches_sam_header(input_format, mapped_reads_file))
        cerr << "[warning: specified input format not found in header "
             << "(" << input_format << ", "
             << mapped_reads_file << ")]" << endl;

    if (!single_end) {
      if (suff_len == 0) {
        size_t repeat_count = 0;
        suff_len = guess_suffix_length(mapped_reads_file, n_reads_to_check,
                                       repeat_count);
        if (repeat_count > 1)
          throw runtime_error("failed to identify read name suffix length\n"
                              "verify reads are not single-end\n"
                              "specify read name suffix length directly");
        if (VERBOSE)
          cerr << "[read name suffix length guess: " << suff_len << "]" << endl;
      }
      else if (!check_suffix_length(mapped_reads_file, suff_len,
                                    n_reads_to_check))
        throw runtime_error("incorrect read name suffix length in: " +
                            mapped_reads_file);

      if (!check_mates_consecutive(mapped_reads_file, suff_len,
                                   n_reads_to_check))
        throw runtime_error("mates are not consecutive in: " +
                            mapped_reads_file);
    }

    // deal with the header
    SAMReader sam_reader(mapped_reads_file);
    out << sam_reader.get_header(); // includes newline
    // ADS: need to check why the command is quoted and has an extra
    // space at the end
    write_pg_line(argc, argv, "FORMAT_READS", VERSION, out);

    // now process the reads
    sam_rec aln, prev_aln;
    string read_name, prev_name;
    bool previous_was_merged = false;

    while (sam_reader >> aln) {
      standardize_format(input_format, aln);

      read_name = remove_suff(aln.qname, suff_len);
      if (read_name == prev_name) {
        if (!is_rc(aln)) // essentially check for dovetail
          swap(prev_aln, aln);
        sam_rec merged;
        const int frag_len = merge_mates(max_frag_len, prev_aln, aln, merged);
        if (frag_len > 0 && frag_len < max_frag_len) {
          out << merged << '\n';
        }
        else {
          out << prev_aln << '\n'
              << aln << '\n';
        }
        previous_was_merged = true;
      }
      else {
        if (!prev_name.empty() && !previous_was_merged) {
          if (is_a_rich(prev_aln))
            flip_conversion(prev_aln);
          out << prev_aln << '\n';
        }
        previous_was_merged = false;
      }
      prev_name = std::move(read_name);
      prev_aln = std::move(aln);
    }

    if (!prev_name.empty() && !previous_was_merged) {
      if (is_a_rich(prev_aln))
        flip_conversion(prev_aln);
      out << prev_aln << '\n';
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
