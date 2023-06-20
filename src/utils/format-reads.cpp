/* format: a command within dnmtools to ensure SAM and BAM format
 * reads are conforming to expectations of dnmtools software
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
#include <sstream>
#include <cstring>

#include <config.h>

// from HTSlib
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

// from smithlab
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
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


template <class T> void
revcomp_inplace(T first, T last) {
  std::transform(first, last, first, complement);
  std::reverse(first, last);
}


// ADS: do we need to check that both mates are on the same strand? Or
// that they are on opposite strands?
static uint32_t
abismal_get_nm(const sam_rec &aln) {
  auto the_entry = find_if(begin(aln.tags), end(aln.tags),
                           [](const string &t) {
                             return t.compare (0, 3, "NM:") == 0;
                           });
  if (the_entry != end(aln.tags))
    return std::atoi(the_entry->c_str() + 5);
  return 0;
}


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


inline bool
format_is_bam_or_sam(htsFile *hts) {
  const htsFormat *fmt = hts_get_format(hts);
  return fmt->category == sequence_data &&
    (fmt->format == bam || fmt->format == sam);
}


static bam1_t *
get_read(samFile *hts, sam_hdr_t *hdr) {
  bam1_t *b = bam_init1();
  const int result = sam_read1(hts, hdr, b);
  if (result >= 0) return b;

  if (result < -1)
    throw runtime_error("error reading file: " + string(hts->fn));
  else // -1 should mean EOF, so we free this read
    bam_destroy1(b);
  return 0;
}


static int
get_sam_record(sam_hdr_t *hdr, bam1_t *b, sam_rec &sr) {
  int fmt_ret = 0;
  kstring_t line = { 0, 0, NULL };
  if ((fmt_ret = sam_format1(hdr, b, &line)) <= 0)
    throw runtime_error("failed parsing sam record");
  sr = sam_rec(line.s);
  return 0;
}


static int
get_n_cigar(const sam_rec &sr) {
  return count_if(begin(sr.cigar), end(sr.cigar),
                  [](unsigned char c) { return std::isalpha(c);});
}


static uint32_t *
get_cigar(const sam_rec &sr, size_t n_cigar) {
  uint32_t *a_cigar = (uint32_t *)malloc(n_cigar*sizeof(uint32_t));
  ssize_t res = sam_parse_cigar(sr.cigar.c_str(), NULL, &a_cigar, &n_cigar);
  return a_cigar;
}


static int
sam_rec_to_bam(sam_hdr_t *hdr, const sam_rec &sr, bam1_t **b) {
  const size_t n_cigar = get_n_cigar(sr);
  uint32_t *a_cigar = get_cigar(sr, n_cigar);
  int32_t tid = sam_hdr_name2tid(hdr, sr.rname.c_str());
  int32_t mtid = sam_hdr_name2tid(hdr, sr.rnext.c_str());

  *b = bam_init1();
  bam_set1(*b, sr.qname.size(), sr.qname.c_str(), sr.flags, tid,
           sr.pos - 1, sr.mapq, n_cigar, a_cigar, mtid, 0, sr.tlen,
           sr.seq.size(), sr.seq.c_str(), 0, 0ul);

  uint32_t nm = abismal_get_nm(sr);
  int res = bam_aux_append(*b, "NM", 'i', sizeof(uint32_t), (uint8_t*)&nm);

  uint8_t *cv_data = (uint8_t *)calloc(1, sizeof(char));
  *cv_data = abismal_is_a_rich(sr) ? 'A' : 'T';
  res = bam_aux_append(*b, "CV", 'A', sizeof(char), cv_data);

  free(a_cigar);
  return 0;
}


static void
write_bam(samFile *out, sam_hdr_t *hdr, const sam_rec &aln) {
  bam1_t *b = 0;
  int r = sam_rec_to_bam(hdr, aln, &b);
  int res = sam_write1(out, hdr, b);
  if (res < 0)
    throw runtime_error("failed writing sam/bam record");
  bam_destroy1(b);
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


// ADS: will move to using this function once it is written
static void
standardize_format(const string &input_format, bam1_t *aln) {

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    int res = 0; // result of htslib calls
    kstring_t s = KS_INITIALIZE;
    res = bam_aux_get_str(aln, "ZS", &s);
    if (res != 1)
      throw runtime_error("record appears to be invalid for bsmap");
    const bool a_rich = s.l > 6 && s.s[6] == '-';

    uint8_t *zs_tag = bam_aux_get(aln, "ZS");
    if (!zs_tag)
      throw runtime_error("error interpreting bam record");
    res = bam_aux_del(aln, zs_tag);
    if (res < 0)
      throw runtime_error("error updating tag");

    const uint8_t cv_val = a_rich ? 'A' : 'T';
    res = bam_aux_append(aln, "CV", 'A', sizeof(char), &cv_val);
    if (res == -1)
      throw runtime_error("error updating tag");

    if (bam_is_rev(aln)) {
      uint8_t *begin_seq = bam_get_seq(aln);
      uint8_t *end_seq = begin_seq + aln->core.l_qseq;
      revcomp_inplace(begin_seq, end_seq);
    }
  }
  // doesn't depend on mapper
  // aln.qual = "*";
}


static vector<string>
load_read_names(const string &inputfile, const size_t n_reads) {
  samFile *hts = hts_open(inputfile.c_str(), "r");
  if (!hts)
    throw runtime_error("error opening file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr)
    throw runtime_error("failed to read header: " + inputfile);

  bam1_t *aln = 0;
  vector<string> names;
  size_t count = 0;
  while ((aln = get_read(hts, hdr)) && count++ < n_reads)
    names.push_back(string(bam_get_qname(aln)));
  return names;
}


static size_t
get_max_repeat_count(const vector<string> &names, const size_t suff_len) {
  // assume "suff_len" is shorter than the shortest entry in "names"
  size_t repeat_count = 0;
  size_t tmp_repeat_count = 0;
  // allow the repeat_count to go to 2, which might not be the "max"
  // but still would indicate that this suffix length is too long and
  // would result in more that two reads identified mutually as mates.
  for (size_t i = 1; i < names.size() && repeat_count < 2; ++i) {
    if (names[i-1].size() == names[i].size() &&
        equal(begin(names[i-1]), end(names[i-1]) - suff_len, begin(names[i])))
      ++tmp_repeat_count;
    else tmp_repeat_count = 0;
    repeat_count = max(repeat_count, tmp_repeat_count);
  }
  return repeat_count;
}


static bool
check_suff_len(const string &inputfile, const size_t suff_len,
                    const size_t n_names_to_check) {
  auto names(load_read_names(inputfile, n_names_to_check));
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
guess_suff_len(const string &inputfile,
                    const size_t n_names_to_check,
                    size_t &repeat_count) {

  // ADS: assuming copy elision but should test it
  auto names(load_read_names(inputfile, n_names_to_check));

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
check_sorted(const string &inputfile,
                        const size_t suff_len, size_t n_reads) {
  // In order to check if mates are consecutive we need to check if a
  // given end has a mate and that mate is not adjacent. This requires
  // storing previous reads, not simply checking for adjacent pairs.

  auto names(load_read_names(inputfile, n_reads));
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
  samFile* hts = hts_open(input_filename.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("error opening file: " + input_filename);
  const htsFormat *fmt = hts_get_format(hts);
  if (fmt->category != sequence_data)
    throw runtime_error("not sequence data: " + input_filename);
  if (fmt->format != bam && fmt->format != sam)
    throw runtime_error("not SAM/BAM file: " + input_filename);
  return true;
}


static bool
check_format_in_header(string input_format,
                                const string &inputfile) {
  samFile* hts = hts_open(inputfile.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("error opening file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr)
    throw runtime_error("failed to read header: " + inputfile);

  auto begin_hdr = sam_hdr_str(hdr);
  auto end_hdr = begin_hdr + std::strlen(begin_hdr);

  auto it = std::search(begin_hdr, end_hdr,
                        begin(input_format), end(input_format),
                        [](const unsigned char a, const unsigned char b) {
                          return std::toupper(a) == std::toupper(b);
                        });
  return it != end_hdr;
}


static void
format(const string &cmd, const size_t n_threads,
       const string &inputfile,
       const string &outfile,
       const bool bam_format, const string &input_format,
       const size_t suff_len, const size_t max_frag_len) {

  samFile* hts = hts_open(inputfile.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("bad htslib file: " + inputfile);

  htsThreadPool the_thread_pool{hts_tpool_init(n_threads), 0};
  if (hts_set_thread_pool(hts, &the_thread_pool) < 0)
    throw runtime_error("error setting threads");

    if (!format_is_bam_or_sam(hts))
      throw runtime_error("bad file format: " + inputfile);

    sam_hdr_t *hdr = sam_hdr_read(hts);
    if (!hdr)
      throw runtime_error("failed to read header: " + inputfile);

    // open the output file
    samFile *out = hts_open(outfile.c_str(), bam_format ? "wb" : "w");

    if (hts_set_thread_pool(out, &the_thread_pool) < 0)
      throw runtime_error("error setting threads");

    // take care of the output file's header
    sam_hdr_t *hdr_out = bam_hdr_dup(hdr);
    if (sam_hdr_add_line(hdr_out, "PG", "ID",
                         "DNMTOOLS", "VN", VERSION, "CL", cmd.c_str(), NULL))
      throw runtime_error("failed to format header");
    if (sam_hdr_write(out, hdr_out))
      throw runtime_error("failed to output header");
    bam_hdr_destroy(hdr_out);
    hdr_out = 0; // don't use it again...

    // now process the reads
    sam_rec aln, prev_aln;
    string read_name, prev_name;
    bool previous_was_merged = false;

    bam1_t *b;
    while (b = get_read(hts, hdr)) {
      get_sam_record(hdr, b, aln);
      standardize_format(input_format, aln);
      read_name = remove_suff(aln.qname, suff_len);
      if (read_name == prev_name) {
        if (!is_rc(aln)) // essentially check for dovetail
          swap(prev_aln, aln);
        sam_rec merged;
        const int frag_len = merge_mates(max_frag_len, prev_aln, aln, merged);
        if (frag_len > 0 && frag_len < max_frag_len) {
          write_bam(out, hdr, merged);
        }
        else {
          write_bam(out, hdr, prev_aln);
          write_bam(out, hdr, aln);
        }
        previous_was_merged = true;
      }
      else {
        if (!prev_name.empty() && !previous_was_merged) {
          if (is_a_rich(prev_aln))
            flip_conversion(prev_aln);
          write_bam(out, hdr, prev_aln);
        }
        previous_was_merged = false;
      }
      prev_name = std::move(read_name);
      prev_aln = std::move(aln);
    }

    if (!prev_name.empty() && !previous_was_merged) {
      if (is_a_rich(prev_aln))
        flip_conversion(prev_aln);
      write_bam(out, hdr, prev_aln);
    }

    // turn off the lights
    bam_hdr_destroy(hdr);
    hts_close(out);
    hts_close(hts);
    hts_tpool_destroy(the_thread_pool.pool);
}


int
main_format(int argc, const char **argv) {

  try {

    size_t n_reads_to_check = 1000000;

    bool bam_format = false;
    bool use_stdout = false;

    string input_format;
    string outfile;
    int max_frag_len = std::numeric_limits<int>::max();
    size_t suff_len = 0;
    bool single_end = false;
    bool VERBOSE = false;
    size_t n_threads = 1;

    const string description = "convert SAM/BAM mapped bs-seq reads "
      "to standard dnmtools format";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<sam/bam-file> [out-file]", 2);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("bam", 'B', "output in BAM format", false, bam_format);
    opt_parse.add_opt("stdout", '\0',
                      "write to standard output", false, use_stdout);
    opt_parse.add_opt("format", 'f', "input format {abismal, bsmap, bismark}",
                      false, input_format);
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
    if (opt_parse.about_requested() || opt_parse.help_requested() ||
        leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_FAILURE;
    }
    if (suff_len != 0 && single_end) {
      cerr << "incompatible arguments specified" << endl
           << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    if ((leftover_args.size() == 1 && !use_stdout) ||
        (leftover_args.size() == 2 && use_stdout)) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_FAILURE;
    }
    const string inputfile(leftover_args.front());
    if (leftover_args.size() == 2 && !use_stdout)
      outfile = leftover_args.back();
    else
      outfile = string("-"); // so htslib can write to stdout
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char*>(cmd, " "));

    check_input_file(inputfile);

    if (VERBOSE)
      if (!check_format_in_header(input_format, inputfile))
        cerr << "[warning: input format not found in header "
             << "(" << input_format << ", " << inputfile << ")]" << endl;

    if (!single_end) {
      if (suff_len == 0) {
        size_t repeat_count = 0;
        suff_len = guess_suff_len(inputfile, n_reads_to_check, repeat_count);
        if (repeat_count > 1)
          throw runtime_error("failed to identify read name suffix length\n"
                              "verify reads are not single-end\n"
                              "specify read name suffix length directly");
        if (VERBOSE)
          cerr << "[read name suffix length guess: " << suff_len << "]" << endl;
      }
      else if (!check_suff_len(inputfile, suff_len, n_reads_to_check))
        throw runtime_error("wrong read name suffix length [" +
                            to_string(suff_len) + "] in: " + inputfile);
      if (!check_sorted(inputfile, suff_len, n_reads_to_check))
        throw runtime_error("mates not consecutive in: " + inputfile);
    }

    if (VERBOSE)
      cerr << "[input file: " << inputfile << "]" << endl
           << "[mapper: " << input_format << "]" << endl
           << "[configuration: " << (single_end ? "SE" : "PE") << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: " << (bam_format ? "B" : "S") << "AM]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl
           << "[readname suffix length: " << suff_len << "]" << endl;

    format(cmd.str(), n_threads, inputfile, outfile,
           bam_format, input_format, suff_len, max_frag_len);

  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
