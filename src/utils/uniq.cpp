/* uniq: remove duplicate reads from a file of mapped reads in the
 * dnmtools format (as output from format_reads), based on identical
 * mapping location and alignment to the reference.
 *
 * Copyright (C) 2013-2023 University of Southern California and
 *                         Andrew D. Smith
 *
 * Author: Andrew D. Smith
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

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

#include <config.h>

#include <htslib/sam.h>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::runtime_error;
using std::to_string;


inline string
read_name(const bam1_t *b) {return string(bam_get_qname(b));}


inline int32_t
get_tid(bam1_t *b) {return b->core.tid;}


inline size_t
qlen(const bam1_t *r) {return bam_cigar2qlen(r->core.n_cigar, bam_get_cigar(r));}


inline bool
precedes_by_start(const bam1_t *a, const bam1_t *b) {
  return a->core.tid == b->core.tid && a->core.pos < b->core.pos;
}


inline bool
precedes_by_end_and_strand(const bam1_t *a, const bam1_t *b) {
  const size_t end_a = bam_endpos(a), end_b = bam_endpos(b);
  return (end_a < end_b || (end_a == end_b && bam_is_rev(a) < bam_is_rev(b)));
}


inline bool
equivalent_chrom_and_start(const bam1_t *a, const bam1_t *b) {
  return a->core.pos == b->core.pos && a->core.tid == b->core.tid;
}


inline bool
equivalent_end_and_strand(const bam1_t *a, const bam1_t *b) {
  return bam_endpos(a) == bam_endpos(b) && bam_is_rev(a) == bam_is_rev(b);
}


static void
write_stats_output(const size_t reads_in, const size_t reads_out,
                   const size_t good_bases_in, const size_t good_bases_out,
                   const size_t reads_with_dups, const string &statfile) {
  if (!statfile.empty()) {
    const size_t reads_removed = reads_in - reads_out;
    const double non_dup_frac =
      (reads_out - reads_with_dups)/static_cast<double>(reads_in);
    const double dup_rate =
      (reads_removed + reads_with_dups)/static_cast<double>(reads_with_dups);
    ofstream out_stat(statfile);
    if (!out_stat) throw runtime_error("bad stats output file");
    out_stat << "total_reads: " << reads_in << endl
             << "total_bases: " << good_bases_in << endl
             << "unique_reads: " << reads_out << endl
             << "unique_read_bases: " << good_bases_out << endl
             << "non_duplicate_fraction: " << non_dup_frac << endl
             << "duplicate_reads: " << reads_with_dups << endl
             << "reads_removed: " << reads_removed << endl
             << "duplication_rate: " << dup_rate << endl;
  }
}


static void
write_hist_output(const vector<size_t> &hist, const string &histfile) {
  if (!histfile.empty()) {
    ofstream out_hist(histfile);
    if (!out_hist) throw runtime_error("bad hist output file");
    for (size_t i = 0; i < hist.size(); ++i)
      if (hist[i] > 0)
        out_hist << i << '\t' << hist[i] << '\n';
  }
}


/* The "inner" buffer corresponds to all reads sharing chrom, start,
 * end and strand, and is a contiguous subset of the "outer" buffer
 * that shares the same end and strand.
 */
static void
process_inner_buffer(const vector<bam1_t*>::const_iterator it,
                     const vector<bam1_t*>::const_iterator jt,
                     sam_hdr_t *hdr, samFile *out,
                     size_t &reads_out,
                     size_t &good_bases_out,
                     size_t &reads_with_dups,
                     vector<size_t> &hist) {
  const size_t n_reads = std::distance(it, jt);
  const size_t selected = rand() % n_reads;
  if (sam_write1(out, hdr, *(it + selected)) < 0)
    throw runtime_error("failed writing bam record");
  if (hist.size() <= n_reads)
    hist.resize(n_reads + 1);
  hist[n_reads]++;
  good_bases_out += qlen(*(it + selected));
  ++reads_out;
  reads_with_dups += (n_reads > 1);
}


/* The buffer corresponds to reads sharing the same mapping chromosome
   and start position. These are gathered and then processed together. */
static void
process_buffer(size_t &reads_out, size_t &good_bases_out,
               size_t &reads_with_dups, vector<size_t> &hist,
               vector<bam1_t*> &buffer, sam_hdr_t *hdr, samFile *out) {
  sort(begin(buffer), end(buffer), precedes_by_end_and_strand);
  auto it(begin(buffer));
  auto jt = it + 1;
  for (; jt != end(buffer); ++jt)
    if (!equivalent_end_and_strand(*it, *jt)) {
      process_inner_buffer(it, jt, hdr, out, reads_out, good_bases_out,
                           reads_with_dups, hist);
      it = jt;
    }
  process_inner_buffer(it, jt, hdr, out, reads_out, good_bases_out,
                       reads_with_dups, hist);

  // free the bam1_t pointers before clearing the buffer
  for (size_t i = 0; i < buffer.size(); ++i)
    if (buffer[i] != 0) {
      bam_destroy1(buffer[i]);
      buffer[i] = 0;
    }
  buffer.clear();
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


static void
uniq(const bool VERBOSE, const size_t n_threads,
     const string &cmd, const string &infile,
     const string &statfile, const string &histfile,
     const bool bam_format, const string &outfile) {

  samFile* hts = hts_open(infile.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("bad htslib file: " + infile);

  if (n_threads > 1 && hts_set_threads(hts, n_threads/2) < 0)
    throw runtime_error("error setting threads");

  if (hts_get_format(hts)->category != sequence_data)
    throw runtime_error("bad file format: " + infile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr)
    throw runtime_error("failed to read header: " + infile);

  // open the output file
  samFile *out = bam_format ? hts_open(outfile.c_str(), "wb") :
    hts_open(outfile.c_str(), "w");

  if (n_threads > 1 && hts_set_threads(out, (n_threads + 1)/2) < 0)
    throw runtime_error("error setting threads");

  // take care of the output file's header
  sam_hdr_t *hdr_out = bam_hdr_dup(hdr);
  if (sam_hdr_add_line(hdr_out, "PG", "ID",
                       "DNMTOOLS", "VN", VERSION, "CL", cmd.c_str(), NULL))
    throw runtime_error("failed to format header");
  if (sam_hdr_write(out, hdr_out))
    throw runtime_error("failed to output header");
  bam_hdr_destroy(hdr_out);

  // try to load the first read
  bam1_t *aln = get_read(hts, hdr);
  if (!aln)
    throw runtime_error("failed parsing read from input file");

  // values to tabulate stats; these cost almost nothing
  vector<size_t> hist;
  size_t reads_in = 1; // count the one we just got
  size_t good_bases_in = qlen(aln); // count its good bases
  size_t reads_out = 0;
  size_t good_bases_out = 0;
  size_t reads_with_dups = 0;

  vector<bam1_t*> buffer(1, aln); // select output from this buffer

  // to check that reads are sorted properly
  vector<bool> chroms_seen(hdr->n_targets, false);
  int32_t cur_chrom = get_tid(aln);

  while (aln = get_read(hts, hdr)) {
    ++reads_in;
    good_bases_in += qlen(aln);

    // below works because buffer is reset every chrom
    if (precedes_by_start(aln, buffer.front()))
      throw runtime_error("input not properly sorted:\n" +
                          read_name(buffer[0]) + "\n" + read_name(aln));
    const int32_t chrom = get_tid(aln);
    if (chrom != cur_chrom) {
      if (chroms_seen[chrom]) throw runtime_error("input not sorted");
      chroms_seen[chrom] = true;
      cur_chrom = chrom;
    }

    if (!equivalent_chrom_and_start(buffer.front(), aln))
      process_buffer(reads_out, good_bases_out, reads_with_dups,
                     hist, buffer, hdr, out);
    buffer.push_back(aln);
  }
  process_buffer(reads_out, good_bases_out, reads_with_dups,
                 hist, buffer, hdr, out);

  // remember to close these
  bam_hdr_destroy(hdr);
  hts_close(out);
  hts_close(hts);

  write_stats_output(reads_in, reads_out, good_bases_in, good_bases_out,
                     reads_with_dups, statfile);

  write_hist_output(hist, histfile);
}


int
main_uniq(int argc, const char **argv) {

  try {

    bool VERBOSE = false;

    bool bam_format = false;
    bool use_stdout = false;

    // ADS: Not recommended to change this seed. It shouldn't matter
    // at all, and we want results to behave as deterministic.
    size_t the_seed = 408;
    string outfile;
    string statfile;
    string histfile;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "program to remove "
                           "duplicate reads from sorted mapped reads",
                           "<in-file> [out-file]", 2);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("stats", 'S', "statistics output file", false, statfile);
    opt_parse.add_opt("hist", '\0', "histogram output file for library"
                      " complexity analysis", false, histfile);
    opt_parse.add_opt("bam", 'B', "output in BAM format", false, bam_format);
    opt_parse.add_opt("stdout", '\0',
                      "write to standard output", false, use_stdout);
    opt_parse.add_opt("seed", 's', "random seed", false, the_seed);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.help_requested()) {
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
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() == 1 && !use_stdout) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() == 2 && use_stdout) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const string infile(leftover_args.front());
    if (leftover_args.size() == 2 && !use_stdout)
      outfile = leftover_args.back();
    else
      outfile = string("-"); // so htslib can write to stdout
    /****************** END COMMAND LINE OPTIONS *****************/

    // ADS: Random here is because we choose randomly when keeping one
    // among duplicate reads.
    srand(the_seed);

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char*>(cmd, " "));

    if (VERBOSE)
      cerr << "[output file: " << outfile << "]" << endl
           << "[output format: " << (bam_format ? "B" : "S") << "AM]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl
           << "[random number seed: " << the_seed << "]" << endl;

    uniq(VERBOSE, n_threads, cmd.str(),
         infile, statfile, histfile, bam_format, outfile);
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
