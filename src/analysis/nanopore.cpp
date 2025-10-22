/* nanocount: methylation counts for nanopore data; see docs for 'counts'
 * command for details.
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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "counts_header.hpp"
#include "dnmt_error.hpp"

/* HTSlib */
#include <htslib/sam.h>

[[nodiscard]] inline bool
is_cytosine(const char c) {
  return c == 'c' || c == 'C';
}

[[nodiscard]] inline bool
is_guanine(const char c) {
  return c == 'g' || c == 'G';
}

[[nodiscard]] inline bool
is_thymine(const char c) {
  return c == 't' || c == 'T';
}

[[nodiscard]] inline bool
is_adenine(const char c) {
  return c == 'a' || c == 'A';
}

[[nodiscard]] inline bool
is_cpg(const std::string &s, const std::size_t i) {
  return i + 1 < std::size(s) && is_cytosine(s[i]) && is_guanine(s[i + 1]);
}

static void
read_fasta_file(const std::string &filename, std::vector<std::string> &names,
                std::vector<std::string> &sequences) {
  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("error opening genome file: " + filename);
  names.clear();
  sequences.clear();

  std::string line;
  while (getline(in, line)) {
    if (line[0] == '>') {
      const auto first_space = line.find_first_of(" \t", 1);
      if (first_space == std::string::npos)
        names.push_back(line.substr(1));
      else
        names.emplace_back(std::cbegin(line) + 1,
                           std::cbegin(line) + line.find_first_of(" \t", 1));
      sequences.emplace_back();
    }
    else
      sequences.back() += line;
  }
}

[[nodiscard]] static std::string
get_basecall_model(const bamxx::bam_header &hdr) {
  kstring_t ks{};

  ks = {0, 0, nullptr};
  const auto rg_ret = sam_hdr_find_line_pos(hdr.h, "RG", 0, &ks);
  if (rg_ret == -2)
    throw std::runtime_error("failure of sam_hdr_find_line_pos");
  if (rg_ret == -1)
    return {};
  ks_free(&ks);

  ks = {0, 0, nullptr};
  const auto ds_ret = sam_hdr_find_tag_pos(hdr.h, "RG", 0, "DS", &ks);
  if (ds_ret == -2)
    throw std::runtime_error("failure of sam_hdr_find_tag_pos");
  if (ds_ret == -1)
    return {};
  std::istringstream buffer(ks.s);
  ks_free(&ks);

  const std::vector<std::string> parts(
    (std::istream_iterator<std::string>(buffer)), {});

  std::string basecall_model;
  for (const auto &key_val : parts) {
    const auto equal_pos = key_val.find('=');
    if (equal_pos == std::string::npos)
      throw std::runtime_error("malformed DS key=val pair: " + key_val);
    if (key_val.substr(0, equal_pos) == "basecall_model")
      return key_val.substr(equal_pos + 1);
  }

  return {};
}

// ADS: here the std::uint16_t allows for up to 256 reads, each contributing up
// to 256 "counts" in the probability encoding.
typedef std::uint16_t count_type;

[[nodiscard]] static inline bool
eats_ref(const std::uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 2;
}

[[nodiscard]] static inline bool
eats_query(const std::uint32_t c) {
  return bam_cigar_type(bam_cigar_op(c)) & 1;
}

/* The three functions below here should probably be moved into
   bsutils.hpp. I am not sure if the DDG function is needed, but it
   seems like if one considers strand, and the CHH is not symmetric,
   then one needs this. Also, Qiang should be consulted on this
   because he spent much time thinking about it in the context of
   plants. */
[[nodiscard]] static bool
is_chh(const std::string &s, const std::size_t i) {
  return i + 2 < std::size(s) && is_cytosine(s[i]) && !is_guanine(s[i + 1]) &&
         !is_guanine(s[i + 2]);
}

[[nodiscard]] static bool
is_ddg(const std::string &s, const std::size_t i) {
  return i + 2 < std::size(s) && !is_cytosine(s[i]) && !is_cytosine(s[i + 1]) &&
         is_guanine(s[i + 2]);
}

[[nodiscard]] static bool
is_c_at_g(const std::string &s, const std::size_t i) {
  return i + 2 < std::size(s) && is_cytosine(s[i]) && !is_cytosine(s[i + 1]) &&
         !is_guanine(s[i + 1]) && is_guanine(s[i + 2]);
}

/* Right now the CountSet objects below are much larger than they need
   to be, for the things we are computing. However, it's not clear
   that the minimum information would really put the memory
   requirement of the program into a more reasonable range, so keeping
   all the information seems reasonable. */
struct CountSet {
  static constexpr auto max_prob_repr = 256.0;
  [[nodiscard]] std::string
  tostring() const {
    // clang-format off
    std::ostringstream oss;
    oss << hydroxy_fwd << '\t' << hydroxy_rev << '\t'
        << methyl_fwd << '\t' << methyl_rev << '\t'
        << n_reads_fwd << '\t' << n_reads_rev << '\n';
    // clang-format on
    return oss.str();
  }
  void
  add_count_fwd(const std::uint8_t h, const std::uint8_t m) {
    hydroxy_fwd += h;
    methyl_fwd += m;
    ++n_reads_fwd;
  }
  void
  add_count_rev(const std::uint8_t h, const std::uint8_t m) {
    hydroxy_rev += h;
    methyl_rev += m;
    ++n_reads_rev;
  }
  [[nodiscard]] double
  get_hydroxy(const bool is_c) const {
    return (is_c ? hydroxy_fwd : hydroxy_rev) / max_prob_repr;
  }
  [[nodiscard]] double
  get_methyl(const bool is_c) const {
    return (is_c ? methyl_fwd : methyl_rev) / max_prob_repr;
  }
  [[nodiscard]] double
  get_mods(const bool is_c) const {
    return get_hydroxy(is_c) + get_methyl(is_c);
  }
  [[nodiscard]] double
  get_n_reads(const bool is_c) const {
    return is_c ? n_reads_fwd : n_reads_rev;
  }
  count_type hydroxy_fwd{0};
  count_type hydroxy_rev{0};
  count_type methyl_fwd{0};
  count_type methyl_rev{0};
  count_type n_reads_fwd{0};
  count_type n_reads_rev{0};
};

/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * std::size_t.
 */
[[nodiscard]] static std::uint32_t
get_tag_from_genome(const std::string &s, const std::size_t pos) {
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

template <typename T>
[[nodiscard]] static std::tuple<T, T>
get_hydroxy_sites(const T mod_pos_beg, const T mod_pos_end) {
  const char hydroxy_tag[] = "C+h?";
  const auto hydroxy_tag_size = 4;
  if (mod_pos_beg == mod_pos_end)
    return {};
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
[[nodiscard]] static std::tuple<T, T>
get_methyl_sites(const T mod_pos_beg, const T mod_pos_end) {
  const char methyl_tag[] = "C+h?";
  const auto methyl_tag_size = 4;
  if (mod_pos_beg == mod_pos_end)
    return {};
  auto methyl_beg = strstr(mod_pos_beg, methyl_tag);
  if (methyl_beg == mod_pos_end)
    return std::tuple<T, T>{};
  auto methyl_end = std::find(methyl_beg, mod_pos_end, ';');
  if (methyl_end == mod_pos_end)
    return std::tuple<T, T>{};
  methyl_beg += methyl_tag_size;
  return std::tuple<T, T>(methyl_beg, methyl_end);
}

[[nodiscard]] static std::tuple<char *, char *>
get_modification_positions(const bamxx::bam_rec &aln) {
  const auto mm_aux = bam_aux_get(aln.b, "MM");
  if (mm_aux == nullptr)
    return {};
  auto mod_pos_beg = bam_aux2Z(mm_aux);
  auto mod_pos_end = mod_pos_beg + std::strlen(mod_pos_beg);
  return {mod_pos_beg, mod_pos_end};
}

struct prob_counter {
  std::array<std::uint64_t, 256> meth_hist{};
  std::array<std::uint64_t, 256> hydro_hist{};
  std::string
  json() const {
    std::ostringstream oss;
    oss << R"({"methyl_hist":[)";
    for (auto i = 0; i < 256; ++i) {
      if (i > 0)
        oss << ',';
      oss << R"(")" << meth_hist[i] << R"(")";
    }
    oss << R"(],"hydroxy_hist":[)";
    for (auto i = 0; i < 256; ++i) {
      if (i > 0)
        oss << ',';
      oss << R"(")" << hydro_hist[i] << R"(")";
    }
    oss << R"(]})";
    return oss.str();
  }
};

struct mod_prob_buffer {
  static constexpr auto init_capacity{128 * 1024};

  std::vector<std::uint8_t> hydroxy_probs;
  std::vector<std::uint8_t> methyl_probs;
  prob_counter pc;

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
            const std::uint8_t m_val = bam_auxB2i(mod_prob, methyl_prob_idx++);
            methyl_probs[i] = m_val;
            pc.meth_hist[m_val]++;

            const std::uint8_t h_val = bam_auxB2i(mod_prob, hydroxy_prob_idx++);
            hydroxy_probs[i] = h_val;
            pc.hydro_hist[h_val]++;
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
            const std::uint8_t m_val = bam_auxB2i(mod_prob, methyl_prob_idx++);
            methyl_probs[i] = m_val;
            pc.meth_hist[m_val]++;

            const std::uint8_t h_val = bam_auxB2i(mod_prob, hydroxy_prob_idx++);
            hydroxy_probs[i] = h_val;
            pc.hydro_hist[h_val]++;
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
static const std::uint8_t encoding[] = {
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
static constexpr auto n_dinucs = 16u;
static const auto dinucs = std::vector{
  "AA", "AC", "AG", "AT",
  "CA", "CC", "CG", "CT",
  "GA", "GC", "GG", "GT",
  "TA", "TC", "TG", "TT",
};
// clang-format on

struct match_counter {
  static constexpr auto dinuc_combos = 256;
  typedef std::array<std::uint64_t, dinuc_combos> container;
  container pos{};
  container neg{};

  void
  add_fwd(const std::uint8_t r1, const std::uint8_t r2, const std::uint8_t q1,
          const std::uint8_t q2) {
    const auto r1e = encoding[r1];
    const auto r2e = encoding[r2];
    const auto q1e = encoding[q1];
    const auto q2e = encoding[q2];
    if ((r1e | r2e | q1e | q2e) & 4u)
      return;
    pos[(r1e * 64) + (r2e * 16) + (q1e * 4) + (q2e * 1)]++;
  }

  void
  add_rev(const std::uint8_t r1, const std::uint8_t r2, const std::uint8_t q1,
          const std::uint8_t q2) {
    const auto r1e = encoding[r1];
    const auto r2e = encoding[r2];
    const auto q1e = encoding[q1];
    const auto q2e = encoding[q2];
    if ((r1e | r2e | q1e | q2e) & 4u)
      return;
    neg[(r1e * 64) + (r2e * 16) + (q1e * 4) + (q2e * 1)]++;
  }

  [[nodiscard]] std::string
  tostring() const {
    std::ostringstream oss;
    oss << "target_query_fwd" << '\n';
    for (auto i = 0u; i < n_dinucs; ++i) {
      const auto beg = std::cbegin(pos) + (i * n_dinucs);
      oss << dinucs[i];
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << '\t' << *itr;
      oss << '\n';
    }
    oss << "target_query_rev" << '\n';
    for (auto i = 0u; i < n_dinucs; ++i) {
      const auto beg = std::cbegin(neg) + (i * n_dinucs);
      oss << dinucs[i];
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << '\t' << *itr;
      oss << '\n';
    }
    return oss.str();
  }

  [[nodiscard]] std::string
  json_row(const container::const_iterator beg) const {
    std::ostringstream oss;
    oss << "{";
    auto j = 0;
    for (auto itr = beg; itr != beg + n_dinucs; ++itr) {
      if (j > 0)
        oss << ",";
      oss << R"(")" << dinucs[j++] << R"(":")" << *itr << R"(")";
    }
    oss << "}";
    return oss.str();
  }

  [[nodiscard]] std::string
  json_table(const container &a) const {
    std::ostringstream oss;
    const auto beg = std::cbegin(a);
    oss << "{";
    for (auto i = 0u; i < n_dinucs; ++i) {
      if (i > 0)
        oss << ",";
      oss << R"(")" << dinucs[i] << R"(":)" << json_row(beg + (i * n_dinucs));
    }
    oss << "}";
    return oss.str();
  }

  [[nodiscard]] std::string
  json() const {
    // clang-format off
    return (std::ostringstream() << "{"
            << R"("map_fwd":)" << json_table(pos) << ","
            << R"("map_rev":)" << json_table(neg)
            << "}").str();
    // clang-format on
  }

  [[nodiscard]] std::string
  tostring_frac() const {
    static constexpr auto width = 8;
    static constexpr auto prec = 3;
    std::ostringstream oss;
    oss.precision(prec);
    oss.setf(std::ios::fixed, std::ios::floatfield);
    oss << "map_fwd\n";
    for (auto i = 0u; i < n_dinucs; ++i) {
      oss << dinucs[i];
      const auto beg = std::cbegin(pos) + (i * n_dinucs);
      const auto tot = std::max(std::accumulate(beg, beg + n_dinucs, 0.0), 1.0);
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << std::setw(width) << *itr / tot;
      oss << '\n';
    }
    oss << "map_ref\n";
    for (auto i = 0u; i < n_dinucs; ++i) {
      oss << dinucs[i];
      const auto beg = std::cbegin(neg) + (i * n_dinucs);
      const auto tot = std::max(std::accumulate(beg, beg + n_dinucs, 0.0), 1.0);
      for (auto itr = beg; itr != beg + n_dinucs; ++itr)
        oss << std::setw(width) << *itr / tot;
      oss << '\n';
    }
    return oss.str();
  }
};

static void
count_states_fwd(const bamxx::bam_rec &aln, std::vector<CountSet> &counts,
                 mod_prob_buffer &mod_buf, const std::string &chrom,
                 match_counter &mc) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  const auto end_ref = std::cend(chrom);

  auto rpos = get_pos(aln);
  auto ref_itr = std::cbegin(chrom) + rpos;

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
        const auto r_nuc = *ref_itr++;
        if (qpos + 1 < end_qpos) {
          mc.add_fwd(r_nuc, ref_itr == end_ref ? 'N' : *ref_itr,
                     seq_nt16_str[bam_seqi(seq, qpos)],
                     qpos == q_lim ? 'N'
                                   : seq_nt16_str[bam_seqi(seq, qpos + 1)]);
        }
        counts[rpos++].add_count_fwd(*hydroxy_prob_itr, *methyl_prob_itr);
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
      ref_itr += n;
    }
  }
  // ADS: somehow previous code included a correction for rpos going
  // past the end of the chromosome; this should result at least in a
  // soft-clip by any mapper. I'm not checking it here as even if it
  // happens I don't want to terminate.
  assert(qpos == get_l_qseq(aln));
}

static void
count_states_rev(const bamxx::bam_rec &aln, std::vector<CountSet> &counts,
                 mod_prob_buffer &mod_buf, const std::string &chrom,
                 match_counter &mc) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  const auto end_ref = std::cend(chrom);
  auto rpos = get_pos(aln);
  auto ref_itr = std::cbegin(chrom) + rpos;

  std::size_t qpos = get_l_qseq(aln);  // to match type with b->core.l_qseq
  std::size_t q_idx = 0;
  std::size_t l_qseq = get_l_qseq(aln);

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
      const std::size_t end_qpos = qpos - n;  // to match type with qpos
      for (; qpos > end_qpos; --qpos) {
        const auto r_nuc = *ref_itr++;
        const auto q_nuc = seq_nt16_str[bam_seqi(seq, q_idx)];
        ++q_idx;
        if (end_qpos + 1 < qpos)
          mc.add_rev(r_nuc, ref_itr == end_ref ? 'N' : *ref_itr, q_nuc,
                     q_idx == l_qseq ? 'N'
                                     : seq_nt16_str[bam_seqi(seq, q_idx)]);
        --methyl_prob_itr;
        --hydroxy_prob_itr;
        counts[rpos++].add_count_rev(*hydroxy_prob_itr, *methyl_prob_itr);
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
      ref_itr += n;
    }
  }

  /* qpos is unsigned; would wrap around if < 0 */
  // ADS: Same as count_states_pos; see comment there
  assert(qpos == 0);
}

[[nodiscard]] static std::tuple<std::map<std::int32_t, std::size_t>,
                                std::set<std::int32_t>>
get_tid_to_idx(
  const bamxx::bam_header &hdr,
  const std::unordered_map<std::string, std::size_t> &name_to_idx) {
  std::set<std::int32_t> missing_tids;
  std::map<std::int32_t, std::size_t> tid_to_idx;
  for (std::int32_t i = 0; i < hdr.h->n_targets; ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to jump
    // through "name_to_idx" and get "tid_to_idx"
    const std::string curr_name(hdr.h->target_name[i]);
    const auto name_itr(name_to_idx.find(curr_name));
    if (name_itr == std::cend(name_to_idx))
      missing_tids.insert(i);
    else
      tid_to_idx[i] = name_itr->second;
  }
  return std::tuple<std::map<std::int32_t, std::size_t>,
                    std::set<std::int32_t>>{tid_to_idx, missing_tids};
}

[[nodiscard]] static bool
consistent_targets(const bamxx::bam_header &hdr,
                   const std::map<std::int32_t, std::size_t> &tid_to_idx,
                   const std::vector<std::string> &names,
                   const std::vector<std::size_t> &sizes) {
  const std::size_t n_targets = hdr.h->n_targets;
  if (n_targets != std::size(names))
    return false;

  for (std::size_t tid = 0; tid < n_targets; ++tid) {
    const std::string tid_name_sam = sam_hdr_tid2name(hdr, tid);
    const std::size_t tid_size_sam = sam_hdr_tid2len(hdr, tid);
    const auto idx_itr = tid_to_idx.find(tid);
    if (idx_itr == std::cend(tid_to_idx))
      return false;
    const auto idx = idx_itr->second;
    if (tid_name_sam != names[idx] || tid_size_sam != sizes[idx])
      return false;
  }
  return true;
}

[[nodiscard]] static bool
consistent_existing_targets(
  const bamxx::bam_header &hdr,
  const std::map<std::int32_t, std::size_t> &tid_to_idx,
  const std::vector<std::string> &names,
  const std::vector<std::size_t> &sizes) {
  const std::size_t n_targets = hdr.h->n_targets;
  for (std::size_t tid = 0; tid < n_targets; ++tid) {
    const auto idx_itr = tid_to_idx.find(tid);
    if (idx_itr == std::cend(tid_to_idx))
      continue;
    const std::string tid_name_sam = sam_hdr_tid2name(hdr, tid);
    const std::size_t tid_size_sam = sam_hdr_tid2len(hdr, tid);
    const auto idx = idx_itr->second;
    if (tid_name_sam != names[idx] || tid_size_sam != sizes[idx])
      return false;
  }
  return true;
}

struct read_processor {
  static constexpr auto default_expected_basecall_model =
    "dna_r10.4.1_e8.2_400bps_hac@v4.3.0";

  bool verbose{};
  bool show_progress{};
  bool compress_output{};
  bool include_header{};
  bool require_covered{};
  bool symmetric{};
  bool cpg_only{};
  bool force{};
  std::uint32_t n_threads{1};
  int strand{0};
  std::string expected_basecall_model{};

  [[nodiscard]] std::string
  tostring() const {
    std::ostringstream oss;
    oss << std::boolalpha;
    oss << "[verbose: " << verbose << "]\n"
        << "[show_progress: " << show_progress << "]\n"
        << "[compress_output: " << compress_output << "]\n"
        << "[include_header: " << include_header << "]\n"
        << "[require_covered: " << require_covered << "]\n"
        << "[symmetric: " << symmetric << "]\n"
        << "[cpg_only: " << cpg_only << "]\n"
        << "[force: " << force << "]\n"
        << "[n_threads: " << n_threads << "]\n"
        << "[strand: " << strand << "]\n"
        << "[expected_basecall_model: " << expected_basecall_model_str()
        << "]\n";
    return oss.str();
  }

  read_processor() : expected_basecall_model{default_expected_basecall_model} {}

  [[nodiscard]] std::string
  expected_basecall_model_str() const {
    return expected_basecall_model.empty() ? "NA" : expected_basecall_model;
  }

  void
  output_skipped_chromosome(
    const std::int32_t tid,
    const std::map<std::int32_t, std::size_t> &tid_to_idx,
    const bamxx::bam_header &hdr,
    const std::vector<std::string>::const_iterator chroms_beg,
    const std::vector<std::size_t> &chrom_sizes, std::vector<CountSet> &counts,
    bamxx::bgzf_file &out) const {

    // get the index of the next chrom sequence
    const auto chrom_idx = tid_to_idx.find(tid);
    if (chrom_idx == std::cend(tid_to_idx)) {
      if (force)
        return;
      else
        throw std::runtime_error("chrom not found: " +
                                 sam_hdr_tid2name(hdr, tid));
    }

    const auto chrom_itr = chroms_beg + chrom_idx->second;

    // reset the counts
    counts.clear();
    counts.resize(chrom_sizes[chrom_idx->second]);

    // ADS: 'false' below for require_covered b/c these sites can't be covered.
    write_output(hdr, out, tid, *chrom_itr, counts);
  }

  void
  write_output_all(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
                   const std::int32_t tid, const std::string &chrom,
                   const std::vector<CountSet> &counts) const {
    static constexpr auto out_fmt = "%ld\t%c\t%s\t%.6g\t%d\t%.6g\t%.6g\n";
    static constexpr auto buf_size = 1024;
    char buffer[buf_size];

    // Put chrom name in buffer and then skip that part for each site because
    // it doesn't change.
    const auto chrom_name = sam_hdr_tid2name(hdr.h, tid);
    if (chrom_name == nullptr)
      throw std::runtime_error("failed to identify chrom for tid: " +
                               std::to_string(tid));
    const auto chrom_name_offset =
      std::snprintf(buffer, buf_size, "%s\t", chrom_name);
    if (chrom_name_offset < 0)
      throw std::runtime_error("failed to write to output buffer");
    auto buffer_after_chrom = buffer + chrom_name_offset;

    for (auto chrom_posn = 0ul; chrom_posn < std::size(chrom); ++chrom_posn) {
      const char base = chrom[chrom_posn];
      if (is_cytosine(base) || is_guanine(base)) {

        const std::uint32_t the_tag = get_tag_from_genome(chrom, chrom_posn);
        if (cpg_only && the_tag != 0)
          continue;

        const bool is_c = is_cytosine(base);

        const std::uint32_t n_reads = counts[chrom_posn].get_n_reads(is_c);
        if (require_covered && n_reads == 0)
          continue;
        const double denom = std::max(n_reads, 1u);
        const double hydroxy = counts[chrom_posn].get_hydroxy(is_c) / denom;
        const double methyl = counts[chrom_posn].get_methyl(is_c) / denom;
        const double mods = counts[chrom_posn].get_mods(is_c) / denom;

        // clang-format off
        const int r = std::snprintf(buffer_after_chrom,
                                    buf_size - chrom_name_offset,
                                    out_fmt,
                                    chrom_posn,
                                    is_c ? '+' : '-',
                                    tag_values[the_tag],
                                    mods,
                                    n_reads,
                                    hydroxy,
                                    methyl);
        // clang-format on

        if (r < 0)
          throw std::runtime_error("failed to write to output buffer");
        out.write(buffer);
      }
    }
  }

  void
  write_output_sym(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
                   const std::int32_t tid, const std::string &chrom,
                   const std::vector<CountSet> &counts) const {
    static constexpr auto out_fmt = "%ld\t+\tCpG\t%.6g\t%d\t%.6g\t%.6g\n";
    static constexpr auto buf_size = 1024;
    char buffer[buf_size];

    // Put chrom name in buffer and then skip that part for each site because
    // it doesn't change.
    const auto chrom_name = sam_hdr_tid2name(hdr.h, tid);
    if (chrom_name == nullptr)
      throw std::runtime_error("failed to identify chrom for tid: " +
                               std::to_string(tid));
    const auto chrom_name_offset =
      std::snprintf(buffer, buf_size, "%s\t", chrom_name);
    if (chrom_name_offset < 0)
      throw std::runtime_error("failed to write to output buffer");
    auto buffer_after_chrom = buffer + chrom_name_offset;

    bool prev_was_c{false};
    std::uint32_t n_reads_pos{};
    double hydroxy_pos{};
    double methyl_pos{};

    for (auto chrom_posn = 0ul; chrom_posn < std::size(chrom); ++chrom_posn) {
      const char base = chrom[chrom_posn];
      if (is_cytosine(base)) {
        prev_was_c = true;
        n_reads_pos = counts[chrom_posn].get_n_reads(true);
        hydroxy_pos = counts[chrom_posn].get_hydroxy(true);
        methyl_pos = counts[chrom_posn].get_methyl(true);
      }
      else {
        if (is_guanine(base) && prev_was_c) {
          const std::uint32_t n_reads_neg =
            counts[chrom_posn].get_n_reads(false);
          const double hydroxy_neg = counts[chrom_posn].get_hydroxy(false);
          const double methyl_neg = counts[chrom_posn].get_methyl(false);

          const auto n_reads = n_reads_pos + n_reads_neg;
          if (require_covered && n_reads == 0)
            continue;

          const double denom = std::max(n_reads, 1u);
          const auto hydroxy = (hydroxy_neg + hydroxy_pos) / denom;
          const auto methyl = (methyl_neg + methyl_pos) / denom;
          const auto mods = hydroxy + methyl;

          // clang-format off
          const int r = std::snprintf(buffer_after_chrom,
                                      buf_size - chrom_name_offset,
                                      out_fmt,
                                      chrom_posn - 1,  // for previous position
                                      mods,
                                      n_reads,
                                      hydroxy,
                                      methyl);
          // clang-format on

          if (r < 0)
            throw std::runtime_error("failed to write to output buffer");
          out.write(buffer);
        }
        prev_was_c = false;
      }
    }
  }

  void
  write_output(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
               const std::int32_t tid, const std::string &chrom,
               const std::vector<CountSet> &counts) const {
    if (symmetric)
      write_output_sym(hdr, out, tid, chrom, counts);
    else
      write_output_all(hdr, out, tid, chrom, counts);
  }

  [[nodiscard]] std::tuple<match_counter, prob_counter>
  operator()(const std::string &infile, const std::string &outfile,
             const std::string &chroms_file) const {
    // first get the chromosome names and sequences from the FASTA file
    std::vector<std::string> chroms, names;
    read_fasta_file(chroms_file, names, chroms);
    for (auto &i : chroms)
      std::transform(std::cbegin(i), std::cend(i), std::begin(i),
                     [](const char c) { return std::toupper(c); });
    if (verbose)
      std::cerr << "[n chroms in reference: " << std::size(chroms) << "]"
                << std::endl;
    const auto chroms_beg = std::cbegin(chroms);

    std::unordered_map<std::string, std::size_t> name_to_idx;
    std::vector<std::size_t> chrom_sizes(std::size(chroms), 0);
    for (std::size_t i = 0; i < std::size(chroms); ++i) {
      name_to_idx[names[i]] = i;
      chrom_sizes[i] = std::size(chroms[i]);
    }

    bamxx::bam_tpool tp(n_threads);  // Must be destroyed after hts

    // open the hts SAM/BAM input file and get the header
    bamxx::bam_in hts(infile);
    if (!hts)
      throw std::runtime_error("failed to open input file");
    // load the input file's header
    bamxx::bam_header hdr(hts);
    if (!hdr)
      throw std::runtime_error("failed to read header");

    auto [tid_to_idx, missing_tids] = get_tid_to_idx(hdr, name_to_idx);
    if (verbose)
      std::cerr << "[targets found: " << std::size(tid_to_idx) << "]\n"
                << "[missing targets: " << std::size(missing_tids) << "]\n";

    if (!force && !missing_tids.empty())
      throw std::runtime_error(
        "missing targets in reference and force not specified");

    if (!force && !consistent_targets(hdr, tid_to_idx, names, chrom_sizes))
      throw std::runtime_error("inconsistent reference genome information");

    if (force &&
        !consistent_existing_targets(hdr, tid_to_idx, names, chrom_sizes))
      throw std::runtime_error("inconsistent reference genome information");

    // open the output file
    const std::string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw std::runtime_error("error opening output file: " + outfile);

    // set the threads for the input file decompression
    if (n_threads > 1) {
      tp.set_io(hts);
      tp.set_io(out);
    }

    // validate the basecall model
    const auto basecall_model = get_basecall_model(hdr);
    if (verbose)
      std::cerr << "[observed basecall model: "
                << (basecall_model.empty() ? "NA" : basecall_model) << "]"
                << std::endl;
    if (!expected_basecall_model.empty() &&
        basecall_model != expected_basecall_model) {
      std::cerr << "failed to match basecall model:" << "\n"
                << "observed="
                << (basecall_model.empty() ? "NA" : basecall_model) << "\n"
                << "expected=" << expected_basecall_model_str() << std::endl;
      return {{}, {}};
    }

    if (include_header)
      write_counts_header_from_bam_header(hdr, out);

    // this is where all the counts are accumulated
    std::vector<CountSet> counts;

    // now iterate over the reads, switching chromosomes and writing
    // output as needed
    bamxx::bam_rec aln;
    std::int32_t prev_tid = -1;

    std::vector<std::string>::const_iterator chrom_itr{};

    match_counter mc;
    mod_prob_buffer mod_buf;

    bool current_target_present = false;

    while (hts.read(hdr, aln)) {
      const std::int32_t tid = get_tid(aln);
      if (get_l_qseq(aln) == 0)
        continue;
      if (tid == -1)  // ADS: skip reads that have no tid -- they are not mapped
        continue;
      if (tid != prev_tid) {  // chrom changed, output results, get next chrom
        // write output if any; counts is empty on first and missing chroms
        if (!counts.empty())
          write_output(hdr, out, prev_tid, *chrom_itr, counts);
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
            output_skipped_chromosome(i, tid_to_idx, hdr, chroms_beg,
                                      chrom_sizes, counts, out);

        // get the next chrom to process
        const auto chrom_idx = tid_to_idx.find(tid);
        current_target_present = (chrom_idx != std::cend(tid_to_idx));
        if (!force && !current_target_present)
          throw std::runtime_error("chromosome not found: " +
                                   std::string(sam_hdr_tid2name(hdr, tid)));

        if (show_progress && current_target_present)
          std::cerr << "processing " << sam_hdr_tid2name(hdr, tid) << std::endl;

        // reset the counts
        counts.clear();
        if (current_target_present) {
          // set size of counts if current target is present
          chrom_itr = chroms_beg + chrom_idx->second;
          counts.resize(chrom_sizes[chrom_idx->second]);

          const bool is_rev = bam_is_rev(aln);
          if (is_rev && strand != 1)
            count_states_rev(aln, counts, mod_buf, *chrom_itr, mc);
          if (!is_rev && strand != 2)
            count_states_fwd(aln, counts, mod_buf, *chrom_itr, mc);
        }
        prev_tid = tid;
      }
      if (current_target_present) {
        const bool is_rev = bam_is_rev(aln);
        if (is_rev && strand != 1)
          count_states_rev(aln, counts, mod_buf, *chrom_itr, mc);
        if (!is_rev && strand != 2)
          count_states_fwd(aln, counts, mod_buf, *chrom_itr, mc);
      }
    }
    if (!counts.empty())
      write_output(hdr, out, prev_tid, *chrom_itr, counts);

    // ADS: if some chroms might not be covered by reads, we have to iterate
    // over what remains
    if (!require_covered)
      for (auto i = prev_tid + 1; i < hdr.h->n_targets; ++i)
        output_skipped_chromosome(i, tid_to_idx, hdr, chroms_beg, chrom_sizes,
                                  counts, out);
    return std::tuple<match_counter, prob_counter>{mc, mod_buf.pc};
  }
};

int
main_nanocount(int argc, char *argv[]) {

  try {

    read_processor rp;

    std::string chroms_file;
    std::string outfile;
    std::string stats_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(std::filesystem::path{argv[0]}.filename(),
                           "get methylation levels from mapped nanopore reads",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "threads to use (few needed)", false,
                      rp.n_threads);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true, chroms_file);
    opt_parse.add_opt("cpg-only", 'n', "print only CpG context cytosines",
                      false, rp.cpg_only);
    opt_parse.add_opt("sym", '\0', "collapse symmetric CpGs (implies -n)",
                      false, rp.symmetric);
    opt_parse.add_opt("require-covered", 'r', "only output covered sites",
                      false, rp.require_covered);
    opt_parse.add_opt("strand", '\0',
                      "strand-specific (1=forward, 2=reverse, 0=both)", false,
                      rp.strand);
    opt_parse.add_opt("stats", 's', "output summary stats in json format",
                      false, stats_file);
    opt_parse.add_opt("force", '\0', "skip consistency checks", false,
                      rp.force);
    opt_parse.add_opt("header", 'H', "add a header to identify the reference",
                      false, rp.include_header);
    opt_parse.add_opt("zip", 'z', "output gzip format", false,
                      rp.compress_output);
    opt_parse.add_opt("progress", '\0', "show progress", false,
                      rp.show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, rp.verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.about_requested() || opt_parse.help_requested() ||
        leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << std::endl
                << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (rp.force)
      rp.expected_basecall_model.clear();

    if (rp.symmetric)
      rp.cpg_only = true;

    if (rp.n_threads <= 0)
      throw std::runtime_error("thread count cannot be negative");

    std::ostringstream cmd;
    std::copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    // file types from HTSlib use "-" for the filename to go to stdout
    if (outfile.empty())
      outfile = "-";

    if (rp.verbose)
      std::cerr << "[input BAM/SAM file: " << mapped_reads_file << "]\n"
                << "[output file: " << outfile << "]\n"
                << "[genome file: " << chroms_file << "]\n"
                << rp.tostring();

    const auto [mc, pc] = rp(mapped_reads_file, outfile, chroms_file);
    if (!stats_file.empty()) {
      std::ofstream stats_out(stats_file);
      if (!stats_out)
        throw std::runtime_error("Error opening stats file: " + stats_file);
      stats_out << R"({"matches":)" << mc.json() << ","
                << R"("probabilities":)" << pc.json() << "}\n";
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
