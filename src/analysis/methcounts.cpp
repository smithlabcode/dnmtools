/* methcounts: a program for counting the methylated and unmethylated
 * reads mapping over each CpG or C
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Song Qiang and Guilherme Sena
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
#include <fstream>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <unordered_set>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MappedRead.hpp"
#include "MSite.hpp"
#include "bsutils.hpp"
#include "dnmt_error.hpp"
#include "cigar_utils.hpp"

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::end;
using std::to_string;


/* The three functions below here should probably be moved into
   bsutils.hpp. I am not sure if the DDG function is needed, but it
   seems like if one considers strand, and the CHH is not symmetric,
   then one needs this. Also, Qiang should be consulted on this
   because he spent much time thinking about it in the context of
   plants. */
static bool
is_chh(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    is_cytosine(s[i]) &&
    !is_guanine(s[i + 1]) &&
    !is_guanine(s[i + 2]);
}


static bool
is_ddg(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    !is_cytosine(s[i]) &&
    !is_cytosine(s[i + 1]) &&
    is_guanine(s[i + 2]);
}


static bool
is_c_at_g(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) &&
    is_cytosine(s[i]) &&
    !is_cytosine(s[i + 1]) &&
    !is_guanine(s[i + 1]) &&
    is_guanine(s[i + 2]);
}


/* Right now the CountSet objects below are much larger than they need
   to be, for the things we are computing. However, it's not clear
   that the minimum information would really put the memory
   requirement of the program into a more reasonable range, so keeping
   all the information seems reasonable. */

template <class count_type>
struct CountSet {

  string tostring() const {
    std::ostringstream oss;
    oss << pA << '\t' << pC << '\t' << pG << '\t' << pT << '\t'
        << nA << '\t' << nC << '\t' << nG << '\t' << nT << '\t' << N;
    return oss.str();
  }
  void add_count_pos(const char x) {
    if (x == 'T') ++pT; // conditions ordered for efficiency
    else if (x == 'C') ++pC;
    else if (x == 'G') ++pG;
    else if (x == 'A') ++pA;
    else ++N;
  }
  void add_count_neg(const char x) {
    if (x == 'T') ++nT; // conditions ordered for efficiency
    else if (x == 'C') ++nC;
    else if (x == 'G') ++nG;
    else if (x == 'A') ++nA;
    else ++N;
  }
  count_type pos_total() const {return pA + pC + pG + pT;}
  count_type neg_total() const {return nA + nC + nG + nT;}

  count_type unconverted_cytosine() const {return pC;}
  count_type converted_cytosine() const {return pT;}
  count_type unconverted_guanine() const {return nC;}
  count_type converted_guanine() const {return nT;}

  count_type pA{0}, pC{0}, pG{0}, pT{0};
  count_type nA{0}, nC{0}, nG{0}, nT{0};
  count_type N{0};
};


/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * size_t.
 */
static string
get_methylation_context_tag_from_genome(const string &s, const size_t pos) {
  if (is_cytosine(s[pos])) {
    if (is_cpg(s, pos)) return "CpG";
    else if (is_chh(s, pos)) return "CHH";
    else if (is_c_at_g(s, pos)) return "CXG";
    else return "CCG";
  }
  if (is_guanine(s[pos])) {
    if (is_cpg(s, pos - 1)) return "CpG";
    else if (is_ddg(s, pos - 2)) return "CHH";
    else if (is_c_at_g(s, pos - 2)) return "CXG";
    else return "CCG";
  }
  return "N";
}


/* This "has_mutated" function looks on the opposite strand to see
 * if the apparent conversion from C->T was actually already in the
 * DNA because of a mutation or SNP.
 */
template <class count_type>
static bool
has_mutated(const char base, const CountSet<count_type> &cs) {
  static const double MUTATION_DEFINING_FRACTION = 0.5;
  return is_cytosine(base) ?
    (cs.nG < MUTATION_DEFINING_FRACTION*(cs.neg_total())) :
    (cs.pG < MUTATION_DEFINING_FRACTION*(cs.pos_total()));
}

inline static bool
is_cpg_site(const string &s, const size_t pos) {
  return (is_cytosine(s[pos]) ? is_guanine(s[pos+1]) :
          (is_guanine(s[pos]) ?
           (pos > 0 && is_cytosine(s[pos - 1])) : false));
}

template <class count_type, class output_type>
static void
write_output(output_type &out,
             const string &chrom_name, const string &chrom,
             const vector<CountSet<count_type> > &counts,
             bool CPG_ONLY) {

  for (size_t i = 0; i < counts.size(); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {
      MSite the_site;
      the_site.chrom = chrom_name;
      the_site.pos = i;
      the_site.strand = is_cytosine(base) ? '+' : '-';
      const double unconverted = is_cytosine(base) ?
        counts[i].unconverted_cytosine() : counts[i].unconverted_guanine();
      const double converted = is_cytosine(base) ?
        counts[i].converted_cytosine() : counts[i].converted_guanine();
      the_site.n_reads = unconverted + converted;
      the_site.meth = the_site.n_reads > 0 ?
        unconverted/(converted + unconverted) : 0.0;
      the_site.context = get_methylation_context_tag_from_genome(chrom, i) +
        (has_mutated(base, counts[i]) ? "x" : "");
      if (!CPG_ONLY || is_cpg_site(chrom, i))
        out << the_site << "\n";
    }
  }
}


inline static size_t
get_chrom_id(const string &chrom_name,
             const unordered_map<string, size_t> &cl) {

  auto the_chrom(cl.find(chrom_name));
  if (the_chrom == end(cl))
    throw dnmt_error("could not find chrom: " + chrom_name);

  return the_chrom->second;
}



template <class count_type>
static void
write_output(sam_hdr_t *hdr, BGZF *out,
             const int32_t tid, const string &chrom,
             const vector<CountSet<count_type> > &counts,
             bool CPG_ONLY) {

  for (size_t i = 0; i < counts.size(); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {
      MSite the_site;
      the_site.chrom = string(sam_hdr_tid2name(hdr, tid));
      the_site.pos = i;
      the_site.strand = is_cytosine(base) ? '+' : '-';
      const double unconverted = is_cytosine(base) ?
        counts[i].unconverted_cytosine() : counts[i].unconverted_guanine();
      const double converted = is_cytosine(base) ?
        counts[i].converted_cytosine() : counts[i].converted_guanine();
      the_site.n_reads = unconverted + converted;
      the_site.meth = the_site.n_reads > 0 ?
        unconverted/(converted + unconverted) : 0.0;
      the_site.context = get_methylation_context_tag_from_genome(chrom, i) +
        (has_mutated(base, counts[i]) ? "x" : "");
      if (!CPG_ONLY || is_cpg_site(chrom, i)) {
        std::ostringstream oss;
        oss << the_site << '\n';
        const ssize_t err_code =
          bgzf_write(out, oss.str().c_str(), oss.str().size());
        if (err_code < 0) throw dnmt_error(err_code, "error writing output");
      }
    }
  }
}


#define get_tid(b) ((b)->core.tid)
#define get_rlen(b) (bam_cigar2rlen((b)->core.n_cigar, bam_get_cigar(b)))
#define get_pos(b) ((b)->core.pos)
#define get_qlen(b) ((b)->core.l_qseq)


static string
bam_seq_inflate(const bam1_t *b) {
  const size_t qlen = get_qlen(b);
  const auto bam_seq = bam_get_seq(b);
  string seq(qlen, '\0');
  for (size_t i = 0; i < qlen; ++i)
    seq[i] = seq_nt16_str[bam_seqi(bam_seq, i)];
  return seq;
}


static string
bam_cigar_inflate(const bam1_t *b) {
  string cigar;
  cigar.reserve(256);
  auto c_itr = bam_get_cigar(b);
  const auto c_end = c_itr + b->core.n_cigar;
  for (; c_itr != c_end; ++c_itr) {
    cigar += to_string(bam_cigar_oplen(*c_itr));
    cigar += bam_cigar_opchr(*c_itr);
  }
  return cigar;
}


template <class count_type> static void
count_states_pos(const bam1_t *aln, vector<CountSet<count_type> > &counts) {

  const size_t width = get_rlen(aln);
  size_t position = get_pos(aln);
  string seq(bam_seq_inflate(aln));
  const string cigar(bam_cigar_inflate(aln));

  apply_cigar(cigar, seq);
  assert(seq.size() == width);

  const size_t chrom_len = counts.size(); // the counts should have
                                          // one entry per position
  if (chrom_len < position) {
    //throw dnmt_error("mapped past chrom end: " + string(bam_get_qname(aln)));
    position = chrom_len;
  }

  for (size_t i = 0; i < width && position < chrom_len; ++i, ++position)
    counts[position].add_count_pos(seq[i]);
}


template <class count_type> static void
count_states_neg(const bam1_t *aln, vector<CountSet<count_type> > &counts) {

  const size_t width = get_rlen(aln);
  size_t position = get_pos(aln) + width;
  string seq(bam_seq_inflate(aln));
  const string cigar(bam_cigar_inflate(aln));

  revcomp_inplace(seq);
  apply_cigar(cigar, seq);
  revcomp_inplace(seq);
  assert(seq.size() == width);

  const size_t chrom_len = counts.size(); // the counts should have
                                          // one entry per position

  if (chrom_len < position) {
    //throw dnmt_error("mapped past chrom end: " + string(bam_get_qname(aln)));
    position = chrom_len;
  }
  // skip past any part of read not overlapping chrom; condition above
  // ensures value for i remains valid as it is incremented in this
  // first loop below
  size_t i = 0;
  for (; position > chrom_len; ++i, --position);

  for (; i < width && position > 0; ++i)
    counts[--position].add_count_neg(seq[i]);
}


static inline bool
has_gz_ext(const string &fn) {
  static const char *ext = ".gz";
  static const size_t ext_sz = 3;
  return fn.size() > ext_sz &&
    (fn.compare(fn.size() - ext_sz, ext_sz, ext) == 0);
}


static void
process_reads(const bool VERBOSE, const size_t n_threads,
              const string &infile, const string &outfile,
              const string &chroms_file, const bool CPG_ONLY) {

  int err_code = 0;

  vector<string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  if (VERBOSE)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

  unordered_map<string, size_t> name_to_idx;
  vector<size_t> chrom_sizes(chroms.size(), 0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = chroms[i].size();
  }

  // open the hts input file
  samFile *hts = hts_open(infile.c_str(), "r");
  if (!hts) throw dnmt_error("failed to open input file");
  // load the input file's header
  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr) throw dnmt_error("failed to read header");

  unordered_map<int32_t, size_t> tid_to_idx;
  for (int32_t i = 0; i < hdr->n_targets; ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to get
    // "tid_to_idx" using "name_to_idx"
    const string curr_name(hdr->target_name[i]);
    const auto name_itr(name_to_idx.find(curr_name));
    if (name_itr == end(name_to_idx))
      throw dnmt_error("failed to find chrom: " + curr_name);
    tid_to_idx[i] = name_itr->second;
  }
  //// ADS: cross-check the chromosome sizes
  // copy(begin(chrom_sizes), end(chrom_sizes),
  // std::ostream_iterator<size_t>(cout, "\n"));

  /* set the threads */
  htsThreadPool tp{hts_tpool_init(n_threads), 0};
  err_code = hts_set_thread_pool(hts, &tp);
  if (err_code < 0) throw dnmt_error(err_code, "error setting threads");

  const string output_mode = has_gz_ext(outfile) ? "w" : "wu";
  BGZF* out = bgzf_open(outfile.c_str(), output_mode.c_str());
  if (!out) throw dnmt_error("error opening output file: " + outfile);

  err_code = bgzf_thread_pool(out, tp.pool, tp.qsize);
  if (err_code) throw dnmt_error(err_code, "error setting bgzf threads");

  // now process the reads
  bam1_t *aln = bam_init1();
  int32_t prev_tid = -1;

  // this is where all the counts are accumulated
  vector<CountSet<uint32_t> > counts;

  unordered_set<int32_t> chroms_seen;
  vector<string>::const_iterator chrom_itr;
  while ((err_code = sam_read1(hts, hdr, aln)) >= 0) {

    // if chrom changes, output results, get the next one
    const int32_t tid = get_tid(aln);
    if (tid != prev_tid) {

      // write output if there is any; should fail only once
      if (!counts.empty())
        write_output(hdr, out, prev_tid, *chrom_itr, counts, CPG_ONLY);

      // make sure all reads from same chrom are consecutive
      if (chroms_seen.find(tid) != end(chroms_seen))
        throw dnmt_error("reads in SAM file not sorted");
      chroms_seen.insert(tid);

      // get the next chrom to process
      auto chrom_idx(tid_to_idx.find(tid));
      if (chrom_idx == end(tid_to_idx))
        throw dnmt_error("chromosome not found: " +
                         string(sam_hdr_tid2name(hdr, tid)));

      if (VERBOSE)
        cerr << "processing " << sam_hdr_tid2name(hdr, tid) << endl;

      prev_tid = tid;
      chrom_itr = begin(chroms) + chrom_idx->second;

      // reset the counts
      counts.clear();
      counts.resize(chrom_sizes[chrom_idx->second]);
    }

    // do the work for this mapped read, depending on strand
    if (bam_is_rev(aln))
      count_states_neg(aln, counts);
    else
      count_states_pos(aln, counts);
  }
  if (!counts.empty())
    write_output(hdr, out, prev_tid, *chrom_itr, counts, CPG_ONLY);

  // turn off the lights
  bam_destroy1(aln);
  bam_hdr_destroy(hdr);
  err_code = hts_close(hts);
  if (err_code < 0) throw dnmt_error(err_code, "process_reads:hts_close");
  err_code = bgzf_close(out);
  if (err_code) throw dnmt_error(err_code, "process_reads:bgzf_close");
  // do this after the files have been closed
  hts_tpool_destroy(tp.pool);
}


int
main_methcounts(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool CPG_ONLY = false;

    string chroms_file;
    string outfile;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "get methylation levels from "
                           "mapped WGBS reads", "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true , chroms_file);
    opt_parse.add_opt("cpg-only", 'n', "print only CpG context cytosines",
                      false, CPG_ONLY);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
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
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!outfile.empty() && !is_valid_output_file(outfile))
      throw dnmt_error("bad output file: " + outfile);

    process_reads(VERBOSE, n_threads, mapped_reads_file,
                  outfile, chroms_file, CPG_ONLY);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
