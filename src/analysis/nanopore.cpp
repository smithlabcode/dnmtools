/* nanopore: counts for nanopore data
 *
 * Copyright (C) 2011-2025 University of Southern California and
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

#include <cstdint>  // for [u]int[0-9]+_t
#include <iostream>
#include <numeric>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"
#include "dnmt_error.hpp"

/* HTSlib */
#include <htslib/sam.h>

using std::cerr;
using std::endl;
using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

using bamxx::bam_header;
using bamxx::bam_rec;
using bamxx::bgzf_file;

struct quick_buf : public std::ostringstream,
                   public std::basic_stringbuf<char> {
  // ADS: By user ecatmur on SO; very fast. Seems to work...
  quick_buf() {
    // ...but this seems to depend on data layout
    static_cast<std::basic_ios<char> &>(*this).rdbuf(this);
  }
  void
  clear() {
    // reset buffer pointers (member functions)
    setp(pbase(), pbase());
  }
  char const *
  c_str() {
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
eats_ref(const uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

static inline bool
eats_query(const uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 1;
}

/* The three functions below here should probably be moved into
   bsutils.hpp. I am not sure if the DDG function is needed, but it
   seems like if one considers strand, and the CHH is not symmetric,
   then one needs this. Also, Qiang should be consulted on this
   because he spent much time thinking about it in the context of
   plants. */
static bool
is_chh(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) && is_cytosine(s[i]) && !is_guanine(s[i + 1]) &&
         !is_guanine(s[i + 2]);
}

static bool
is_ddg(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) && !is_cytosine(s[i]) &&
         !is_cytosine(s[i + 1]) && is_guanine(s[i + 2]);
}

static bool
is_c_at_g(const std::string &s, size_t i) {
  return (i < (s.length() - 2)) && is_cytosine(s[i]) &&
         !is_cytosine(s[i + 1]) && !is_guanine(s[i + 1]) &&
         is_guanine(s[i + 2]);
}

/* Right now the CountSet objects below are much larger than they need
   to be, for the things we are computing. However, it's not clear
   that the minimum information would really put the memory
   requirement of the program into a more reasonable range, so keeping
   all the information seems reasonable. */
struct CountSet {

  string
  tostring() const {
    std::ostringstream oss;
    oss << pos_prob << '\t' << neg_prob << '\t' << pos_n << '\t' << neg_n;
    return oss.str();
  }
  void
  add_count_pos(const std::uint8_t x) {
    pos_prob += x;
    ++pos_n;
  }
  void
  add_count_neg(const std::uint8_t x) {
    neg_prob += x;
    ++neg_n;
  }
  count_type
  pos_total() const {
    return pos_n;
  }
  count_type
  neg_total() const {
    return neg_n;
  }

  count_type pos_prob{0}, neg_prob{0};
  count_type pos_n{0}, neg_n{0};
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
    if (is_cpg(s, pos))
      return 0;
    else if (is_chh(s, pos))
      return 1;
    else if (is_c_at_g(s, pos))
      return 2;
    else
      return 3;
  }
  if (is_guanine(s[pos])) {
    if (is_cpg(s, pos - 1))
      return 0;
    else if (is_ddg(s, pos - 2))
      return 1;
    else if (is_c_at_g(s, pos - 2))
      return 2;
    else
      return 3;
  }
  return 4;  // shouldn't be used for anything
}

static const char *tag_values[] = {
  "CpG",   // 0
  "CHH",   // 1
  "CXG",   // 2
  "CCG",   // 3
  "N",     // 4
  "CpGx",  // 5 <---- MUT_OFFSET
  "CHHx",  // 6
  "CXGx",  // 7
  "CCGx",  // 8
  "Nx"     // 9
};
static const uint32_t MUT_OFFSET = 5;

static inline uint32_t
tag_with_mut(const uint32_t tag, const bool mut) {
  return tag + (mut ? MUT_OFFSET : 0);
}

template <const bool require_covered = false>
static void
write_output(const bam_header &hdr, bgzf_file &out, const int32_t tid,
             const string &chrom, const vector<CountSet> &counts,
             bool CPG_ONLY) {

  quick_buf buf;  // keep underlying buffer space?
  for (size_t i = 0; i < std::size(chrom); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {

      const uint32_t the_tag = get_tag_from_genome(chrom, i);
      if (CPG_ONLY && the_tag != 0)
        continue;

      const bool is_c = is_cytosine(base);
      const double n_reads = is_c ? counts[i].pos_n : counts[i].neg_n;
      double n_meth = is_c ? counts[i].pos_prob : counts[i].neg_prob;
      n_meth = n_meth > 0.0 ? (n_meth / 256.0) : 0.0;

      if (require_covered && n_reads == 0)
        continue;

      const bool mut = false;
      buf.clear();
      // ADS: here is where we make an MSite, but not using MSite
      buf << sam_hdr_tid2name(hdr, tid) << '\t' << i << '\t'
          << (is_c ? '+' : '-') << '\t'
          << tag_values[tag_with_mut(the_tag, mut)] << '\t'
          << (n_reads > 0.0 ? n_meth / n_reads : 0.0) << '\t' << n_reads
          << '\n';
      if (!out.write(buf.c_str(), buf.tellp()))
        throw dnmt_error("error writing output");
    }
  }
}

static inline std::vector<std::uint8_t>
get_ml_values(const bamxx::bam_rec &aln) {
  static constexpr auto min_ml_size = 7;  // ML:B:C
  // Get the ML values
  kstring_t s = {0, 0, nullptr};
  if (bam_aux_get_str(aln.b, "ML", &s) <= 0) {
    ks_free(&s);
    return {};
  }
  std::string ml(ks_str(&s));
  ks_free(&s);
  if (std::size(ml) < min_ml_size)
    return {};
  ml = ml.substr(min_ml_size);
  std::vector<std::string> parts = smithlab::split(ml, ",", false);
  std::vector<uint8_t> probs;
  for (const auto &p : parts)
    probs.push_back(atoi(p.data()));

  auto prob_itr = std::cbegin(probs) + std::size(probs) / 2;
  std::vector<std::uint8_t> vals(get_l_qseq(aln), 0);
  const auto seq = bam_get_seq(aln);
  for (std::int32_t i = 0; i + 1 < get_l_qseq(aln); ++i)
    if (seq_nt16_str[bam_seqi(seq, i)] == 'C' &&
        seq_nt16_str[bam_seqi(seq, i + 1)] == 'G')
      vals[i] += *prob_itr++;
  assert(prob_itr == std::cend(probs));

  return vals;
}

static inline std::vector<std::uint8_t>
get_ml_values_rev(const bamxx::bam_rec &aln) {
  static constexpr auto min_ml_size = 7;  // ML:B:C
  // Get the ML values
  kstring_t s = {0, 0, nullptr};
  if (bam_aux_get_str(aln.b, "ML", &s) <= 0) {
    ks_free(&s);
    return {};
  }
  std::string ml(ks_str(&s));
  ks_free(&s);
  if (std::size(ml) < min_ml_size)
    return {};
  ml = ml.substr(min_ml_size);
  std::vector<std::string> parts = smithlab::split(ml, ",", false);
  std::vector<uint8_t> probs;
  for (const auto &p : parts)
    probs.push_back(atoi(p.data()));

  auto prob_itr = std::cbegin(probs) + std::size(probs) / 2;
  const std::uint32_t n = get_l_qseq(aln);
  std::vector<std::uint8_t> vals(n, 0);
  const auto seq = bam_get_seq(aln);
  for (std::int32_t i = 0; i + 1 < get_l_qseq(aln); ++i)
    if (seq_nt16_str[bam_seqi(seq, i)] == 'C' &&
        seq_nt16_str[bam_seqi(seq, i + 1)] == 'G')
      vals[n - i - 2] += *prob_itr++;
  assert(prob_itr == std::cend(probs));

  return vals;
}

static void
count_states_pos(const bam_rec &aln, vector<CountSet> &counts) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);

  auto rpos = get_pos(aln);
  auto qpos = 0;  // to match type with b->core.l_qseq

  const auto probs = get_ml_values(aln);
  if (probs.empty())
    return;

  auto prob_itr = std::cbegin(probs);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos)
        counts[rpos++].add_count_pos(*prob_itr++);
    }
    else if (eats_query(op)) {
      qpos += n;
      prob_itr += n;
    }
    else if (eats_ref(op)) {
      rpos += n;
    }
  }
  // ADS: somehow previous code included a correction for rpos going
  // past the end of the chromosome; this should result at least in a
  // soft-clip by any mapper. I'm not checking it here as even if it
  // happens I don't want to terminate.
  // assert(qpos == get_l_qseq(aln));
  if (qpos != get_l_qseq(aln)) {
    std::cerr << qpos << '\t' << get_l_qseq(aln) << std::endl;
  }
}

[[maybe_unused]] static void
count_states_neg(const bam_rec &aln, vector<CountSet> &counts) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  size_t rpos = get_pos(aln);
  size_t qpos = get_l_qseq(aln);  // to match type with b->core.l_qseq

  const auto probs = get_ml_values_rev(aln);
  if (probs.empty())
    return;
  assert(std::size(probs) == qpos);

  auto prob_itr = std::cend(probs);

  if (qpos == 0)
    return;
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      assert(qpos >= n);
      const size_t end_qpos = qpos - n;  // to match type with qpos
      for (; qpos > end_qpos; --qpos)
        counts[rpos++].add_count_neg(*--prob_itr);
    }
    else if (eats_query(op)) {
      qpos -= n;
      prob_itr -= n;
    }
    else if (eats_ref(op)) {
      rpos += n;
    }
  }
  /* qpos is unsigned; would wrap around if < 0 */
  // ADS: Same as count_states_pos; see comment there
  assert(qpos == 0);
}

static unordered_map<int32_t, size_t>
get_tid_to_idx(const bam_header &hdr,
               const unordered_map<string, size_t> name_to_idx) {
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
  return tid_to_idx;
}

template <class CS>
static void
output_skipped_chromosome(const bool CPG_ONLY, const int32_t tid,
                          const unordered_map<int32_t, size_t> &tid_to_idx,
                          const bam_header &hdr,
                          const vector<string>::const_iterator chroms_beg,
                          const vector<size_t> &chrom_sizes, vector<CS> &counts,
                          bgzf_file &out) {

  // get the index of the next chrom sequence
  const auto chrom_idx = tid_to_idx.find(tid);
  if (chrom_idx == cend(tid_to_idx))
    throw dnmt_error("chrom not found: " + sam_hdr_tid2name(hdr, tid));

  const auto chrom_itr = chroms_beg + chrom_idx->second;

  // reset the counts
  counts.clear();
  counts.resize(chrom_sizes[chrom_idx->second]);

  write_output<false>(hdr, out, tid, *chrom_itr, counts, CPG_ONLY);
}

static bool
consistent_targets(const bam_header &hdr,
                   const unordered_map<int32_t, size_t> &tid_to_idx,
                   const vector<string> &names, const vector<size_t> &sizes) {
  const size_t n_targets = hdr.h->n_targets;
  if (n_targets != names.size())
    return false;

  for (size_t tid = 0; tid < n_targets; ++tid) {
    const string tid_name_sam = sam_hdr_tid2name(hdr, tid);
    const size_t tid_size_sam = sam_hdr_tid2len(hdr, tid);
    const auto idx_itr = tid_to_idx.find(tid);
    if (idx_itr == cend(tid_to_idx))
      return false;
    const auto idx = idx_itr->second;
    if (tid_name_sam != names[idx] || tid_size_sam != sizes[idx])
      return false;
  }
  return true;
}

template <const bool require_covered = false>
static void
process_reads(const bool VERBOSE, const bool show_progress,
              const bool compress_output, const bool include_header,
              const size_t n_threads, const string &infile,
              const string &outfile, const string &chroms_file,
              const bool CPG_ONLY) {
  // first get the chromosome names and sequences from the FASTA file
  vector<string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(cbegin(i), cend(i), begin(i),
              [](const char c) { return std::toupper(c); });
  if (VERBOSE)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;
  const auto chroms_beg = cbegin(chroms);

  unordered_map<string, size_t> name_to_idx;
  vector<size_t> chrom_sizes(chroms.size(), 0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = chroms[i].size();
  }

  bamxx::bam_tpool tp(n_threads);  // Must be destroyed after hts

  // open the hts SAM/BAM input file and get the header
  bamxx::bam_in hts(infile);
  if (!hts)
    throw dnmt_error("failed to open input file");
  // load the input file's header
  bam_header hdr(hts);
  if (!hdr)
    throw dnmt_error("failed to read header");

  unordered_map<int32_t, size_t> tid_to_idx = get_tid_to_idx(hdr, name_to_idx);

  if (!consistent_targets(hdr, tid_to_idx, names, chrom_sizes))
    throw dnmt_error("inconsistent reference genome information");

  // open the output file
  const string output_mode = compress_output ? "w" : "wu";
  bgzf_file out(outfile, output_mode);
  if (!out)
    throw dnmt_error("error opening output file: " + outfile);

  // set the threads for the input file decompression
  if (n_threads > 1) {
    tp.set_io(hts);
    tp.set_io(out);
  }

  if (include_header)
    write_counts_header_from_bam_header(hdr, out);

  // now iterate over the reads, switching chromosomes and writing
  // output as needed
  bam_rec aln;
  int32_t prev_tid = -1;

  // this is where all the counts are accumulated
  vector<CountSet> counts;

  vector<string>::const_iterator chrom_itr{};

  while (hts.read(hdr, aln)) {
    const int32_t tid = get_tid(aln);
    if (get_l_qseq(aln) == 0)
      continue;

    if (tid == -1)  // ADS: skip reads that have no tid -- they are not mapped
      continue;
    if (tid == prev_tid) {
      if (bam_is_rev(aln))
        count_states_neg(aln, counts);
      else
        count_states_pos(aln, counts);
    }
    else {  // chrom has changed, so output results and get the next chrom

      // write output if there is any; counts is empty only for first chrom
      if (!counts.empty())
        write_output<require_covered>(hdr, out, prev_tid, *chrom_itr, counts,
                                      CPG_ONLY);
      // make sure reads are sorted chrom tid number in header
      if (tid < prev_tid) {
        const std::string message = "SAM file is not sorted "
                                    "previous tid: " +
                                    std::to_string(prev_tid) +
                                    " current tid: " + std::to_string(tid);
        throw dnmt_error(message);
      }

      if (!require_covered)
        for (auto i = prev_tid + 1; i < tid; ++i)
          output_skipped_chromosome(CPG_ONLY, i, tid_to_idx, hdr, chroms_beg,
                                    chrom_sizes, counts, out);

      // get the next chrom to process
      auto chrom_idx(tid_to_idx.find(tid));
      if (chrom_idx == end(tid_to_idx))
        throw dnmt_error("chromosome not found: " +
                         string(sam_hdr_tid2name(hdr, tid)));
      if (show_progress)
        cerr << "processing " << sam_hdr_tid2name(hdr, tid) << endl;

      prev_tid = tid;
      chrom_itr = chroms_beg + chrom_idx->second;

      // reset the counts
      counts.clear();
      counts.resize(chrom_sizes[chrom_idx->second]);
    }
  }
  if (!counts.empty())
    write_output<require_covered>(hdr, out, prev_tid, *chrom_itr, counts,
                                  CPG_ONLY);

  // ADS: if some chroms might not be covered by reads, we have to
  // iterate over what remains
  if (!require_covered)
    for (auto i = prev_tid + 1; i < hdr.h->n_targets; ++i)
      output_skipped_chromosome(CPG_ONLY, i, tid_to_idx, hdr, chroms_beg,
                                chrom_sizes, counts, out);
}

int
main_nanopore(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool show_progress = false;
    bool CPG_ONLY = false;
    bool require_covered = false;
    bool compress_output = false;
    bool include_header = false;

    string chroms_file;
    string outfile;
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "get methylation levels from mapped nanopore reads",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "threads to use (few needed)", false,
                      n_threads);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true, chroms_file);
    opt_parse.add_opt("cpg-only", 'n', "print only CpG context cytosines",
                      false, CPG_ONLY);
    opt_parse.add_opt("require-covered", 'r', "only output covered sites",
                      false, require_covered);
    opt_parse.add_opt("header", 'H', "add a header to identify the reference",
                      false, include_header);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
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
    copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    // file types from HTSlib use "-" for the filename to go to stdout
    if (outfile.empty())
      outfile = "-";

    if (VERBOSE)
      cerr << "[input BAM/SAM file: " << mapped_reads_file << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: " << (compress_output ? "bgzf" : "text") << "]"
           << endl
           << "[genome file: " << chroms_file << "]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[CpG only mode: " << (CPG_ONLY ? "yes" : "no") << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl;

    if (require_covered)
      process_reads<true>(VERBOSE, show_progress, compress_output,
                          include_header, n_threads, mapped_reads_file, outfile,
                          chroms_file, CPG_ONLY);
    else
      process_reads(VERBOSE, show_progress, compress_output, include_header,
                    n_threads, mapped_reads_file, outfile, chroms_file,
                    CPG_ONLY);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
