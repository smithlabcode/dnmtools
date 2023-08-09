/* counts: a program for counting the methylated and unmethylated
 * reads mapping over each CpG or C
 *
 * Copyright (C) 2011-2023 University of Southern California and
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
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <cstdint> // for [u]int[0-9]+_t

#include "OptionParser.hpp"
// #include "GenomicRegion.hpp"
/* ADS: This code writes MSite objects to files, but does not use
   MSite to do it. Possiby dangerous, but currently much faster. If
   MSite has a way to serialize into a char[] directly, then we should
   use it. */
// #include "MSite.hpp"
#include "bsutils.hpp"
#include "dnmt_error.hpp"
#include "bam_record_utils.hpp"

/* HTSlib */
#include <htslib/sam.h>

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::unordered_set;
using std::unordered_map;

using bamxx::bam_rec;

struct quick_buf : public std::ostringstream,
                   public std::basic_stringbuf<char> {
  // ADS: By user ecatmur on SO; very fast. Seems to work...
  quick_buf() {
    // ...but this seems to depend on data layout
    static_cast<std::basic_ios<char>&>(*this).rdbuf(this);
  }
  void clear() {
    // reset buffer pointers (member functions)
    setp(pbase(), pbase());
  }
  char const* c_str() {
    /* between c_str and insertion make sure to clear() */
    *pptr() = '\0';
    return pbase();
  }
};


// ADS: we should never have to worry about coverage over > 32767 in
// any downstream analysis, so using "int16_t" here would allow to
// detect wrap around and report it as some kind of weird thing, maybe
// zeroing it and flagging the output. As it is now, if we really need
// >65535-fold coverage, we can make the change here.
typedef uint16_t count_type;


static inline bool
eats_ref(const uint32_t c) {return bam_cigar_type(bam_cigar_op(c)) & 2;}


static inline bool
eats_query(const uint32_t c) {return bam_cigar_type(bam_cigar_op(c)) & 1;}


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
struct CountSet {

  string tostring() const {
    std::ostringstream oss;
    oss << pA << '\t' << pC << '\t' << pG << '\t' << pT << '\t'
        << nA << '\t' << nC << '\t' << nG << '\t' << nT;
    // << '\t' << N; /* not used */
    return oss.str();
  }
  void add_count_pos(const char x) {
    if (x == 'T') ++pT; // conditions ordered for efficiency
    else if (x == 'C') ++pC;
    else if (x == 'G') ++pG;
    else if (x == 'A') ++pA;
    // else ++N; /* not used */
  }
  void add_count_neg(const char x) {
    if (x == 'T') ++nT; // conditions ordered for efficiency(??)
    else if (x == 'C') ++nC;
    else if (x == 'G') ++nG;
    else if (x == 'A') ++nA;
    // else ++N; /* not used */
  }
  count_type pos_total() const {return pA + pC + pG + pT;}
  count_type neg_total() const {return nA + nC + nG + nT;}

  count_type unconverted_cytosine() const {return pC;}
  count_type converted_cytosine() const {return pT;}
  count_type unconverted_guanine() const {return nC;}
  count_type converted_guanine() const {return nT;}

  count_type pA{0}, pC{0}, pG{0}, pT{0};
  count_type nA{0}, nC{0}, nG{0}, nT{0};
  // count_type N; /* this wasn't used and breaks alignment */
};


/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * size_t.
 */
static uint32_t
get_tag_from_genome(const string &s, const size_t pos) {
  if (is_cytosine(s[pos])) {
    if (is_cpg(s, pos)) return 0;
    else if (is_chh(s, pos)) return 1;
    else if (is_c_at_g(s, pos)) return 2;
    else return 3;
  }
  if (is_guanine(s[pos])) {
    if (is_cpg(s, pos - 1)) return 0;
    else if (is_ddg(s, pos - 2)) return 1;
    else if (is_c_at_g(s, pos - 2)) return 2;
    else return 3;
  }
  return 4; // shouldn't be used for anything
}


/* This "has_mutated" function looks on the opposite strand to see
 * if the apparent conversion from C->T was actually already in the
 * DNA because of a mutation or SNP.
 */
static bool
has_mutated(const char base, const CountSet &cs) {
  static const double MUTATION_DEFINING_FRACTION = 0.5;
  return is_cytosine(base) ?
    (cs.nG < MUTATION_DEFINING_FRACTION*(cs.neg_total())) :
    (cs.pG < MUTATION_DEFINING_FRACTION*(cs.pos_total()));
}


static inline bool
is_cpg_site(const string &s, const size_t pos) {
  return (is_cytosine(s[pos]) ? is_guanine(s[pos + 1]) :
          (is_guanine(s[pos]) ?
           (pos > 0 && is_cytosine(s[pos - 1])) : false));
}


static inline size_t
get_chrom_id(const string &chrom_name,
             const unordered_map<string, size_t> &cl) {

  auto the_chrom(cl.find(chrom_name));
  if (the_chrom == end(cl))
    throw dnmt_error("could not find chrom: " + chrom_name);

  return the_chrom->second;
}


static const char *tag_values[] = {
  "CpG", // 0
  "CHH", // 1
  "CXG", // 2
  "CCG", // 3
  "N",   // 4
  "CpGx",// 5 <---- MUT_OFFSET
  "CHHx",// 6
  "CXGx",// 7
  "CCGx",// 8
  "Nx"   // 9
};
static const uint32_t MUT_OFFSET = 5;


static inline uint32_t
tag_with_mut(const uint32_t tag, const bool mut) {
  return tag + (mut ? MUT_OFFSET : 0);
}


static void
write_output(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
             const int32_t tid, const string &chrom,
             const vector<CountSet> &counts, bool CPG_ONLY) {

  quick_buf buf; // keep underlying buffer space?

  for (size_t i = 0; i < counts.size(); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {

      const uint32_t the_tag = get_tag_from_genome(chrom, i);
      if (CPG_ONLY && the_tag != 0) continue;

      const bool is_c = is_cytosine(base);
      const double unconverted = is_c ?
        counts[i].unconverted_cytosine() : counts[i].unconverted_guanine();
      const double converted = is_c ?
        counts[i].converted_cytosine() : counts[i].converted_guanine();
      const bool mut = has_mutated(base, counts[i]);
      const size_t n_reads = unconverted + converted;
      buf.clear();
      // ADS: here is where we make an MSite, but not using MSite
      buf << sam_hdr_tid2name(hdr, tid) << '\t'
          << i << '\t'
          << (is_c ? '+' : '-') << '\t'
          << tag_values[tag_with_mut(the_tag, mut)] << '\t'
          << (n_reads > 0 ? unconverted/n_reads : 0.0) << '\t'
          << n_reads << '\n';
      if (!out.write(buf.c_str(), buf.tellp()))
        throw dnmt_error("error writing output");
    }
  }
}


/* ADS: was tempted to use these until I passed in a "b++"... */
// #define get_tid(b) ((b)->core.tid)
// #define get_rlen(b) (bam_cigar2rlen((b)->core.n_cigar, bam_get_cigar(b)))
// #define get_pos(b) ((b)->core.pos)
// #define get_qlen(b) ((b)->core.l_qseq)







static void
count_states_pos(const bam_rec &aln, vector<CountSet> &counts) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = 0; // to match type with b->core.l_qseq
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos) {
        // ADS: beware!!! bam_seqi is a macro, so no "qpos++" inside
        // its arguments! Why macros?!?! Just make sure the compiler
        // inliens it properly ffs!
        counts[rpos++].add_count_pos(seq_nt16_str[bam_seqi(seq, qpos)]);
      }
    }
    else if (eats_query(op)) {
      qpos += n;
    }
    else if (eats_ref(op)) {
      rpos += n;
    }
  }
  // ADS: somehow previous code included a correction for rpos going
  // past the end of the chromosome; this should result at least in a
  // soft-clip by any mapper. I'm not checking it here as even if it
  // happens I don't want to terminate.
  assert(qpos == get_l_qseq(aln));
}


static void
count_states_neg(const bam_rec &aln, vector<CountSet> &counts) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  size_t rpos = get_pos(aln);
  size_t qpos = get_l_qseq(aln); // to match type with b->core.l_qseq
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const size_t end_qpos = qpos - n; // to match type with qpos
      for (; qpos > end_qpos; --qpos) // beware ++ in macro below!!!
        counts[rpos++].add_count_neg(seq_nt16_str[bam_seqi(seq, qpos-1)]);
    }
    else if (eats_query(op)) {
      qpos -= n;
    }
    else if (eats_ref(op)) {
      rpos += n;
    }
  }
  /* qpos is unsigned; would wrap around if < 0 */
  // ADS: Same as count_states_pos; see comment there
  assert(qpos == 0);
}


static void
process_reads(const bool VERBOSE,
              const bool compress_output, const size_t n_threads,
              const string &infile, const string &outfile,
              const string &chroms_file, const bool CPG_ONLY) {

  vector<string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &&i: chroms)
    transform(begin(i), end(i), begin(i),
              [](const char c){return std::toupper(c);});
  if (VERBOSE)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;

  unordered_map<string, size_t> name_to_idx;
  vector<size_t> chrom_sizes(chroms.size(), 0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = chroms[i].size();
  }

  bamxx::bam_tpool tp(n_threads); // Must be destroyed after hts

  // open the hts SAM/BAM input file and get the header
  bamxx::bam_in hts(infile);
  if (!hts) throw dnmt_error("failed to open input file");
  // load the input file's header
  bamxx::bam_header hdr(hts);
  if (!hdr) throw dnmt_error("failed to read header");

  unordered_map<int32_t, size_t> tid_to_idx;
  for (int32_t i = 0; i < hdr.h->n_targets; ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to jump
    // through "name_to_idx" and get "tid_to_idx"
    const string curr_name(hdr.h->target_name[i]);
    const auto name_itr(name_to_idx.find(curr_name));
    if (name_itr == end(name_to_idx))
      throw dnmt_error("failed to find chrom: " + curr_name);
    tid_to_idx[i] = name_itr->second;
  }
  //// ADS: really should cross-check the chromosome sizes

  // open the output file
  const string output_mode = compress_output ? "w" : "wu";
  bamxx::bgzf_file out(outfile, output_mode);
  if (!out) throw dnmt_error("error opening output file: " + outfile);

  /* set the threads for the input file decompression */
  if (n_threads > 1) {
    tp.set_io(hts);
    tp.set_io(out);
  }

  /* now iterate over the reads, switching chromosomes and writing
     output as needed */
  bam_rec aln;
  int32_t prev_tid = -1;

  // this is where all the counts are accumulated
  vector<CountSet> counts;

  unordered_set<int32_t> chroms_seen;
  vector<string>::const_iterator chrom_itr;
  while (hts.read(hdr, aln)) {

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

}


int
main_counts(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool CPG_ONLY = false;
    bool compress_output = false;

    string chroms_file;
    string outfile;
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "get methylation levels from "
                           "mapped bisulfite sequencing reads",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "threads to use (few needed)",
                      false, n_threads);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true , chroms_file);
    opt_parse.add_opt("cpg-only", 'n', "print only CpG context cytosines",
                      false, CPG_ONLY);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
      return EXIT_SUCCESS;
    }
    const string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (n_threads < 0)
      throw dnmt_error("thread count cannot be negative");

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char*>(cmd, " "));

    // file types from HTSlib use "-" for the filename to go to stdout
    if (outfile.empty())
      outfile = "-";

    if (VERBOSE)
      cerr << "[input BAM/SAM file: " << mapped_reads_file << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: "
           << (compress_output ? "bgzf" : "text") << "]" << endl
           << "[genome file: " << chroms_file << "]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[CpG only mode: " << (CPG_ONLY ? "yes" : "no") << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl;

    process_reads(VERBOSE, compress_output, n_threads, mapped_reads_file,
                  outfile, chroms_file, CPG_ONLY);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
