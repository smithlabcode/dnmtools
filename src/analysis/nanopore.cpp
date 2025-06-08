/* nanocount: methylation counts for nanopore data
 *
 * Copyright (C) 2025 Andrew D. Smith
 *
 * Author: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include <array>
#include <cstdint>  // for [u]int[0-9]+_t
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
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

// ADS: here the uint16_t allows for up to 256 reads, each contributing up to
// 256 "counts" in the probability encoding.
typedef uint16_t count_type;

static inline bool
eats_ref(const std::uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

static inline bool
eats_query(const std::uint32_t c) {
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
  std::string
  tostring() const {
    // clang-format off
    std::ostringstream oss;
    oss << hydroxy_pos << '\t' << hydroxy_neg << '\t'
        << methyl_pos << '\t' << methyl_neg << '\t'
        << n_reads_pos << '\t' << n_reads_neg << '\n';
    // clang-format on
    return oss.str();
  }
  void
  add_count_pos(const std::uint8_t h, const std::uint8_t m) {
    hydroxy_pos += h;
    methyl_pos += m;
    ++n_reads_pos;
  }
  void
  add_count_neg(const std::uint8_t h, const std::uint8_t m) {
    hydroxy_neg += h;
    methyl_neg += m;
    ++n_reads_neg;
  }
  double
  get_hydroxy(const bool is_c) const {
    return is_c ? hydroxy_pos : hydroxy_neg;
  }
  double
  get_methyl(const bool is_c) const {
    return is_c ? methyl_pos : methyl_neg;
  }
  double
  get_mods(const bool is_c) const {
    return get_hydroxy(is_c) + get_methyl(is_c);
  }
  double
  get_n_reads(const bool is_c) const {
    return is_c ? n_reads_pos : n_reads_neg;
  }
  count_type hydroxy_pos{0};
  count_type hydroxy_neg{0};
  count_type methyl_pos{0};
  count_type methyl_neg{0};
  count_type n_reads_pos{0};
  count_type n_reads_neg{0};
};

/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * size_t.
 */
static std::uint32_t
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
  "CpG",  // 0
  "CHH",  // 1
  "CXG",  // 2
  "CCG",  // 3
  "N",    // 4
};

template <const bool require_covered = false>
static void
write_output(const bam_header &hdr, bgzf_file &out, const std::int32_t tid,
             const string &chrom, const vector<CountSet> &counts,
             bool CPG_ONLY) {

  quick_buf buf;  // keep underlying buffer space?
  for (size_t i = 0; i < std::size(chrom); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {

      const std::uint32_t the_tag = get_tag_from_genome(chrom, i);
      if (CPG_ONLY && the_tag != 0)
        continue;

      const bool is_c = is_cytosine(base);

      const double n_reads = counts[i].get_n_reads(is_c);
      if (require_covered && n_reads < 1.0)
        continue;
      const double denom = std::max(n_reads, 1.0);
      const double hydroxy = counts[i].get_hydroxy(is_c) / 256.0 / denom;
      const double methyl = counts[i].get_methyl(is_c) / 256.0 / denom;
      const double mods = counts[i].get_mods(is_c) / 256.0 / denom;

      buf.clear();
      // ADS: here is where we make an MSite, but not using MSite
      // clang-format off
      buf << sam_hdr_tid2name(hdr, tid) << '\t' << i << '\t'
          << (is_c ? '+' : '-') << '\t' << tag_values[the_tag] << '\t'
          << mods << '\t'
          << n_reads << '\t'
          << hydroxy << '\t'
          << methyl
          << '\n';
      // clang-format on
      if (!out.write(buf.c_str(), buf.tellp()))
        throw std::runtime_error("error writing output");
    }
  }
}

template <typename T>
static std::tuple<T, T>
get_hydroxy_sites(T mod_pos_beg, T mod_pos_end) {
  const char hydroxy_tag[] = "C+h?";
  const auto hydroxy_tag_size = 4;
  auto hydroxy_beg = strstr(mod_pos_beg, hydroxy_tag);
  if (hydroxy_beg == mod_pos_end)
    return std::tuple<T, T>{};
  auto hydroxy_end = std::find(hydroxy_beg, mod_pos_end, ';');
  if (hydroxy_end == mod_pos_end)
    return std::tuple<T, T>{};
  hydroxy_beg += hydroxy_tag_size;
  return std::tuple<T, T>(hydroxy_beg, hydroxy_end);
}

template <typename T>
static std::tuple<T, T>
get_methyl_sites(T mod_pos_beg, T mod_pos_end) {
  const char methyl_tag[] = "C+h?";
  const auto methyl_tag_size = 4;
  auto methyl_beg = strstr(mod_pos_beg, methyl_tag);
  if (methyl_beg == mod_pos_end)
    return std::tuple<T, T>{};
  auto methyl_end = std::find(methyl_beg, mod_pos_end, ';');
  if (methyl_end == mod_pos_end)
    return std::tuple<T, T>{};
  methyl_beg += methyl_tag_size;
  return std::tuple<T, T>(methyl_beg, methyl_end);
}

static std::tuple<char *, char *>
get_modification_positions(const bamxx::bam_rec &aln) {
  const auto mm_aux = bam_aux_get(aln.b, "MM");
  if (mm_aux == nullptr)
    return {};
  auto mod_pos_beg = bam_aux2Z(mm_aux);
  auto mod_pos_end = mod_pos_beg + std::strlen(mod_pos_beg);
  return {mod_pos_beg, mod_pos_end};
}

struct mod_prob_buffer {
  static constexpr auto init_capacity{128 * 1024};
  std::vector<std::uint8_t> hydroxy_probs;
  std::vector<std::uint8_t> methyl_probs;
  mod_prob_buffer() {
    methyl_probs.reserve(init_capacity);
    hydroxy_probs.reserve(init_capacity);
  }
  bool
  set_probs(const bamxx::bam_rec &aln) {
    const auto get_next_mod_pos = [](auto &b, const auto e) -> std::int32_t {
      const auto isdig = [](const auto x) {
        return std::isdigit(static_cast<unsigned char>(x));
      };
      b = std::find_if(b, e, isdig);
      if (b == e)
        return -1;
      auto r = atoi(b);
      b = std::find_if_not(b, e, isdig);
      return r;
    };

    auto [mod_pos_beg, mod_pos_end] = get_modification_positions(aln);
    auto [hydroxy_beg, hydroxy_end] =
      get_hydroxy_sites(mod_pos_beg, mod_pos_end);
    auto [methyl_beg, methyl_end] = get_methyl_sites(mod_pos_beg, mod_pos_end);

    // check if any are false
    if (mod_pos_beg == nullptr || hydroxy_beg == nullptr ||
        methyl_beg == nullptr)
      return false;

    // assume that hydroxy and methyl both point to CpG sites
    if (!std::equal(hydroxy_beg, hydroxy_end, methyl_beg, methyl_end))
      return false;

    const auto mod_prob = bam_aux_get(aln.b, "ML");
    if (mod_prob == nullptr)
      return false;

    // number of commas is number of hydroxy substrates = CpG sites
    const auto n_cpgs = std::count(hydroxy_beg, hydroxy_end, ',');

    const auto qlen = get_l_qseq(aln);
    const auto seq = bam_get_seq(aln);

    methyl_probs.clear();
    methyl_probs.resize(qlen, 0);

    hydroxy_probs.clear();
    hydroxy_probs.resize(qlen, 0);

    std::int32_t delta = get_next_mod_pos(hydroxy_beg, hydroxy_end);

    auto hydroxy_prob_idx = 0;      // start of modifications
    auto methyl_prob_idx = n_cpgs;  // start methyl after hydroxy

    if (bam_is_rev(aln)) {
      for (auto i = 0; i < qlen; ++i) {
        const auto nuc = seq_nt16_str[bam_seqi(seq, qlen - i - 1)];
        if (nuc == 'G') {
          if (seq_nt16_str[bam_seqi(seq, qlen - i - 2)] == 'C') {
            // assume that when delta hits 0 we have a CpG site
            if (delta != 0)
              return false;
            methyl_probs[i] = bam_auxB2i(mod_prob, methyl_prob_idx++);
            hydroxy_probs[i] = bam_auxB2i(mod_prob, hydroxy_prob_idx++);
          }
          --delta;
          if (delta < 0)
            delta = get_next_mod_pos(hydroxy_beg, hydroxy_end);
        }
      }
    }
    else {
      for (auto i = 0; i + 1 < qlen; ++i) {
        const auto nuc = seq_nt16_str[bam_seqi(seq, i)];
        if (nuc == 'C') {
          if (seq_nt16_str[bam_seqi(seq, i + 1)] == 'G') {
            // assume that when delta hits 0 we have a CpG site
            if (delta != 0)
              return false;
            methyl_probs[i] = bam_auxB2i(mod_prob, methyl_prob_idx++);
            hydroxy_probs[i] = bam_auxB2i(mod_prob, hydroxy_prob_idx++);
          }
          --delta;
          if (delta < 0)
            delta = get_next_mod_pos(hydroxy_beg, hydroxy_end);
        }
      }
    }
    return true;
  }
};

// clang-format off
static const std::uint8_t enc[] = {
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 16
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 32
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 48
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 64
  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,  // 80
  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 96
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 112
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 128
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 144
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 160
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 176
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 192
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 208
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 224
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,  // 240
  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4   // 256
};
// clang-format on

static constexpr auto n_dinucs = 16u;
// clang-format off
static const auto dinucs = std::vector{
  "AA", "AC", "AG", "AT",
  "CA", "CC", "CG", "CT",
  "GA", "GC", "GG", "GT",
  "TA", "TC", "TG", "TT",
};
// clang-format on

struct match_counter {
  std::array<std::uint64_t, 256> c{};

  void
  add_pos(const std::uint8_t r1, const std::uint8_t r2, const std::uint8_t q1,
          const std::uint8_t q2) {
    const auto r1e = enc[r1];
    const auto r2e = enc[r2];
    const auto q1e = enc[q1];
    const auto q2e = enc[q2];
    if ((r1e | r2e | q1e | q2e) & 4u)
      return;
    c[(r1e * 64) + (r2e * 16) + (q1e * 4) + (q2e * 1)]++;
  }

  std::string
  tostring() const {
    std::ostringstream oss;
    for (auto i = 0u; i < n_dinucs; ++i) {
      const auto beg = std::cbegin(c) + (i * n_dinucs);
      oss << dinucs[i];
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << '\t' << *itr;
      oss << '\n';
    }
    return oss.str();
  }

  std::string
  tostring_frac() const {
    static constexpr auto width = 8;
    static constexpr auto prec = 3;
    std::ostringstream oss;
    oss.precision(prec);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    for (auto i = 0u; i < n_dinucs; ++i) {
      oss << dinucs[i];
      const auto beg = std::cbegin(c) + (i * n_dinucs);
      const auto tot = std::max(std::accumulate(beg, beg + n_dinucs, 0.0), 1.0);
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << std::setw(width) << *itr / tot;
      oss << '\n';
    }
    return oss.str();
  }
};

static void
count_states_pos(const bam_rec &aln, vector<CountSet> &counts,
                 mod_prob_buffer &mod_buf, const std::string &chrom,
                 match_counter &mc) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);

  auto rpos = get_pos(aln);
  auto r_itr = std::cbegin(chrom) + rpos;
  const auto r_itr_lim = std::cend(chrom);

  if (!mod_buf.set_probs(aln))
    return;

  auto qpos = 0;                     // to match type with b->core.l_qseq
  auto q_lim = get_l_qseq(aln) - 1;  // if here, this is >= 0

  auto hydroxy_prob_itr = std::cbegin(mod_buf.hydroxy_probs);
  auto methyl_prob_itr = std::cbegin(mod_buf.methyl_probs);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos) {
        const auto r_nuc = *r_itr++;
        mc.add_pos(r_nuc, r_itr == r_itr_lim ? 'N' : *r_itr,
                   seq_nt16_str[bam_seqi(seq, qpos)],
                   qpos == q_lim ? 'N' : seq_nt16_str[bam_seqi(seq, qpos + 1)]);
        counts[rpos++].add_count_pos(*hydroxy_prob_itr, *methyl_prob_itr);
        ++methyl_prob_itr;
        ++hydroxy_prob_itr;
      }
    }
    else if (eats_query(op)) {
      qpos += n;
      methyl_prob_itr += n;
      hydroxy_prob_itr += n;
    }
    else if (eats_ref(op)) {
      rpos += n;
      r_itr += n;
    }
  }
  // ADS: somehow previous code included a correction for rpos going
  // past the end of the chromosome; this should result at least in a
  // soft-clip by any mapper. I'm not checking it here as even if it
  // happens I don't want to terminate.
  assert(qpos == get_l_qseq(aln));
}

[[maybe_unused]] static void
count_states_neg(const bam_rec &aln, vector<CountSet> &counts,
                 mod_prob_buffer &mod_buf, const std::string &chrom,
                 match_counter &mc) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  size_t rpos = get_pos(aln);
  auto r_itr = std::cbegin(chrom) + rpos;
  const auto r_lim = std::cend(chrom);

  size_t qpos = get_l_qseq(aln);  // to match type with b->core.l_qseq
  size_t q_idx = 0;
  size_t l_qseq = get_l_qseq(aln);

  if (!mod_buf.set_probs(aln))
    return;

  auto hydroxy_prob_itr = std::cend(mod_buf.hydroxy_probs);
  auto methyl_prob_itr = std::cend(mod_buf.methyl_probs);

  if (qpos == 0)
    return;
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      assert(qpos >= n);
      const size_t end_qpos = qpos - n;  // to match type with qpos
      for (; qpos > end_qpos; --qpos) {
        const auto r_nuc = *r_itr++;
        const auto q_nuc = seq_nt16_str[bam_seqi(seq, q_idx)];
        ++q_idx;
        mc.add_pos(r_nuc, r_itr == r_lim ? 'N' : *r_itr, q_nuc,
                   q_idx == l_qseq ? 'N' : seq_nt16_str[bam_seqi(seq, q_idx)]);
        --methyl_prob_itr;
        --hydroxy_prob_itr;
        counts[rpos++].add_count_neg(*hydroxy_prob_itr, *methyl_prob_itr);
      }
    }
    else if (eats_query(op)) {
      qpos -= n;
      q_idx += n;
      methyl_prob_itr -= n;
      hydroxy_prob_itr -= n;
    }
    else if (eats_ref(op)) {
      rpos += n;
      r_itr += n;
    }
  }

  /* qpos is unsigned; would wrap around if < 0 */
  // ADS: Same as count_states_pos; see comment there
  assert(qpos == 0);
}

static std::unordered_map<std::int32_t, size_t>
get_tid_to_idx(const bam_header &hdr,
               const std::unordered_map<string, size_t> name_to_idx) {
  std::unordered_map<std::int32_t, size_t> tid_to_idx;
  for (std::int32_t i = 0; i < hdr.h->n_targets; ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to jump
    // through "name_to_idx" and get "tid_to_idx"
    const string curr_name(hdr.h->target_name[i]);
    const auto name_itr(name_to_idx.find(curr_name));
    if (name_itr == end(name_to_idx))
      throw std::runtime_error("failed to find chrom: " + curr_name);
    tid_to_idx[i] = name_itr->second;
  }
  return tid_to_idx;
}

template <class CS>
static void
output_skipped_chromosome(
  const bool CPG_ONLY, const std::int32_t tid,
  const std::unordered_map<std::int32_t, size_t> &tid_to_idx,
  const bam_header &hdr, const vector<string>::const_iterator chroms_beg,
  const vector<size_t> &chrom_sizes, vector<CS> &counts, bgzf_file &out) {

  // get the index of the next chrom sequence
  const auto chrom_idx = tid_to_idx.find(tid);
  if (chrom_idx == std::cend(tid_to_idx))
    throw std::runtime_error("chrom not found: " + sam_hdr_tid2name(hdr, tid));

  const auto chrom_itr = chroms_beg + chrom_idx->second;

  // reset the counts
  counts.clear();
  counts.resize(chrom_sizes[chrom_idx->second]);

  write_output<false>(hdr, out, tid, *chrom_itr, counts, CPG_ONLY);
}

static bool
consistent_targets(const bam_header &hdr,
                   const std::unordered_map<std::int32_t, size_t> &tid_to_idx,
                   const vector<string> &names, const vector<size_t> &sizes) {
  const size_t n_targets = hdr.h->n_targets;
  if (n_targets != names.size())
    return false;

  for (size_t tid = 0; tid < n_targets; ++tid) {
    const string tid_name_sam = sam_hdr_tid2name(hdr, tid);
    const size_t tid_size_sam = sam_hdr_tid2len(hdr, tid);
    const auto idx_itr = tid_to_idx.find(tid);
    if (idx_itr == std::cend(tid_to_idx))
      return false;
    const auto idx = idx_itr->second;
    if (tid_name_sam != names[idx] || tid_size_sam != sizes[idx])
      return false;
  }
  return true;
}

template <const bool require_covered = false>
static std::tuple<match_counter, match_counter>
process_reads(const bool VERBOSE, const bool show_progress,
              const bool compress_output, const bool include_header,
              const size_t n_threads, const string &infile,
              const string &outfile, const string &chroms_file,
              const bool CPG_ONLY, const int strand) {
  // first get the chromosome names and sequences from the FASTA file
  vector<string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(std::cbegin(i), std::cend(i), std::begin(i),
              [](const char c) { return std::toupper(c); });
  if (VERBOSE)
    cerr << "[n chroms in reference: " << chroms.size() << "]" << endl;
  const auto chroms_beg = std::cbegin(chroms);

  std::unordered_map<string, size_t> name_to_idx;
  vector<size_t> chrom_sizes(chroms.size(), 0);
  for (size_t i = 0; i < chroms.size(); ++i) {
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = chroms[i].size();
  }

  bamxx::bam_tpool tp(n_threads);  // Must be destroyed after hts

  // open the hts SAM/BAM input file and get the header
  bamxx::bam_in hts(infile);
  if (!hts)
    throw std::runtime_error("failed to open input file");
  // load the input file's header
  bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::unordered_map<std::int32_t, size_t> tid_to_idx =
    get_tid_to_idx(hdr, name_to_idx);

  if (!consistent_targets(hdr, tid_to_idx, names, chrom_sizes))
    throw std::runtime_error("inconsistent reference genome information");

  // open the output file
  const string output_mode = compress_output ? "w" : "wu";
  bgzf_file out(outfile, output_mode);
  if (!out)
    throw std::runtime_error("error opening output file: " + outfile);

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
  std::int32_t prev_tid = -1;

  // this is where all the counts are accumulated
  vector<CountSet> counts;

  vector<string>::const_iterator chrom_itr{};

  match_counter mc_pos;
  match_counter mc_neg;
  mod_prob_buffer mod_buf;

  while (hts.read(hdr, aln)) {
    const std::int32_t tid = get_tid(aln);
    if (get_l_qseq(aln) == 0)
      continue;

    if (tid == -1)  // ADS: skip reads that have no tid -- they are not mapped
      continue;
    if (tid == prev_tid) {
      const bool is_rev = bam_is_rev(aln);
      if (is_rev && strand != 1)
        count_states_neg(aln, counts, mod_buf, *chrom_itr, mc_neg);
      if (!is_rev && strand != 2)
        count_states_pos(aln, counts, mod_buf, *chrom_itr, mc_pos);
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
        throw std::runtime_error(message);
      }

      if (!require_covered)
        for (auto i = prev_tid + 1; i < tid; ++i)
          output_skipped_chromosome(CPG_ONLY, i, tid_to_idx, hdr, chroms_beg,
                                    chrom_sizes, counts, out);

      // get the next chrom to process
      auto chrom_idx(tid_to_idx.find(tid));
      if (chrom_idx == end(tid_to_idx))
        throw std::runtime_error("chromosome not found: " +
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
  return {mc_pos, mc_neg};
}

int
main_nanocount(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool show_progress = false;
    bool CPG_ONLY = false;
    bool require_covered = false;
    bool compress_output = false;
    bool include_header = false;
    int strand = 0;

    string chroms_file;
    string outfile;
    string stats_file;
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
    opt_parse.add_opt("strand", '\0',
                      "use strand (1=positive, 2=negative, 0=default)", false,
                      strand);
    opt_parse.add_opt("stats", 's', "output match/mismatch stats to this file",
                      false, stats_file);
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
      throw std::runtime_error("thread count cannot be negative");

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

    const auto [mc_pos, mc_neg] = [&] {
      if (require_covered)
        return process_reads<true>(VERBOSE, show_progress, compress_output,
                                   include_header, n_threads, mapped_reads_file,
                                   outfile, chroms_file, CPG_ONLY, strand);
      return process_reads(VERBOSE, show_progress, compress_output,
                           include_header, n_threads, mapped_reads_file,
                           outfile, chroms_file, CPG_ONLY, strand);
    }();

    if (!stats_file.empty()) {
      std::ofstream stats_out(stats_file);
      if (!stats_out)
        std::cerr << "Error opening stats file" << std::endl;
      else {
        stats_out << mc_pos.tostring_frac() << std::endl;
        stats_out << mc_neg.tostring_frac() << std::endl;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
