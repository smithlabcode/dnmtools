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

#ifdef BUILD_NANOPORE

#include "counts_header.hpp"

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"

#include "bamxx.hpp"

#include "nlohmann/json.hpp"

#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

// clang-format off
static constexpr std::array<std::uint8_t, 256> encoding = {
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

static constexpr auto n_nucs = 4u;

[[nodiscard]] inline bool
is_cytosine(const char c) {
  return c == 'c' || c == 'C';
}

[[nodiscard]] inline bool
is_guanine(const char c) {
  return c == 'g' || c == 'G';
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

  const auto parts =
    std::vector<std::string>(std::istream_iterator<std::string>(buffer), {});

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
[[nodiscard]] static inline bool
is_chh(const std::string &s, const std::size_t i) {
  return i + 2 < std::size(s) && is_cytosine(s[i]) && !is_guanine(s[i + 1]) &&
         !is_guanine(s[i + 2]);
}

[[nodiscard]] static inline bool
is_ddg(const std::string &s, const std::size_t i) {
  return i + 2 < std::size(s) && !is_cytosine(s[i]) && !is_cytosine(s[i + 1]) &&
         is_guanine(s[i + 2]);
}

[[nodiscard]] static inline bool
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
  static constexpr auto max_prob_repr = 255.0;
  // ADS: accepting int16_t because of using -1 for unknown prob vs. 0 prob.
  void
  add_count_fwd(const std::int16_t h, const std::int16_t m) {
    hydroxy_fwd += static_cast<std::uint8_t>(h);
    methyl_fwd += static_cast<std::uint8_t>(m);
    ++n_reads_fwd;
  }
  void
  add_count_rev(const std::int16_t h, const std::int16_t m) {
    hydroxy_rev += static_cast<std::uint8_t>(h);
    methyl_rev += static_cast<std::uint8_t>(m);
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

enum class missing_code : std::uint8_t {
  unmodified = 0,
  unknown = 1,
};

// NOLINTBEGIN(*-avoid-c-arrays)
static const char *const tag_values[] = {
  "CpG",  // 0
  "CHH",  // 1
  "CXG",  // 2
  "CCG",  // 3
  "N",    // 4
};
// NOLINTEND(*-avoid-c-arrays)

struct mod_prob_buffer {
  static constexpr auto init_capacity{128 * 1024};
  static constexpr auto max_mods = 10;
  // scratch
  std::array<hts_base_mod, max_mods> mods{};
  std::unique_ptr<hts_base_mod_state, void (*)(hts_base_mod_state *)> m;

  std::vector<std::int16_t> hydroxy_probs;
  std::vector<std::int16_t> methyl_probs;

  mod_prob_buffer() : m{hts_base_mod_state_alloc(), &hts_base_mod_state_free} {
    methyl_probs.reserve(init_capacity);
    hydroxy_probs.reserve(init_capacity);
  }

  [[nodiscard]] bool
  set_probs(const bamxx::bam_rec &aln) {
    static constexpr auto h_idx = 0;
    static constexpr auto m_idx = 1;
    const auto qlen = get_l_qseq(aln);
    const auto d = mods.data();
    const auto is_rev = bam_is_rev(aln);

    if (!bam_aux_get(aln.b, "MM") && !bam_aux_get(aln.b, "Mm"))
      return false;

    methyl_probs.clear();
    methyl_probs.resize(qlen, -1);

    hydroxy_probs.clear();
    hydroxy_probs.resize(qlen, -1);

    bam_parse_basemod(aln.b, m.get());
    // ADS: or bam_parse_basemod2(aln.b, m, HTS_MOD_REPORT_UNCHECKED)

    int n_types{};
    const auto types = bam_mods_recorded(m.get(), &n_types);
    // NOLINTNEXTLINE(*-pointer-arithmetic)
    if (n_types < 2 || (types[h_idx] != 'h' && types[m_idx] != 'm'))
      return false;

    int pos{};
    int n{};
    while ((n = bam_next_basemod(aln.b, m.get(), d, max_mods, &pos)) > 0) {
      if (n < m_idx)
        continue;
      pos = is_rev ? qlen - 1 - pos : pos;
      methyl_probs[pos] = mods[m_idx].qual;
      hydroxy_probs[pos] = mods[h_idx].qual;
    }
    return true;
  }
};

static void
count_states_fwd(const bamxx::bam_rec &aln, std::vector<CountSet> &counts,
                 mod_prob_buffer &mod_buf, const std::string &chrom) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig =
    beg_cig + get_n_cigar(aln);  // NOLINT(*-pointer-arithmetic)

  auto rpos = get_pos(aln);
  auto ref_itr = std::cbegin(chrom) + rpos;

  if (!mod_buf.set_probs(aln))
    return;

  // position in the query sequence
  auto qpos = 0;

  auto hydroxy_prob_itr = std::cbegin(mod_buf.hydroxy_probs);
  auto methyl_prob_itr = std::cbegin(mod_buf.methyl_probs);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = qpos + n;
      for (; qpos < end_qpos; ++qpos, ++rpos) {
        if (*hydroxy_prob_itr >= 0 && *methyl_prob_itr >= 0)
          counts[rpos].add_count_fwd(*hydroxy_prob_itr, *methyl_prob_itr);
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
                 mod_prob_buffer &mod_buf, const std::string &chrom) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig =
    beg_cig + get_n_cigar(aln);  // NOLINT(*-pointer-arithmetic)

  auto rpos = get_pos(aln);
  auto ref_itr = std::cbegin(chrom) + rpos;

  if (!mod_buf.set_probs(aln))
    return;

  // position in the query sequence
  auto qpos = get_l_qseq(aln);
  if (qpos == 0)
    return;

  auto hydroxy_prob_itr = std::cend(mod_buf.hydroxy_probs);
  auto methyl_prob_itr = std::cend(mod_buf.methyl_probs);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = qpos - n;
      for (; qpos > end_qpos; --qpos, ++rpos) {
        --methyl_prob_itr;
        --hydroxy_prob_itr;
        if (*hydroxy_prob_itr >= 0 && *methyl_prob_itr >= 0)
          counts[rpos].add_count_rev(*hydroxy_prob_itr, *methyl_prob_itr);
      }
    }
    else if (eats_query(op)) {
      qpos -= n;
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
  for (std::int32_t i = 0; i < get_n_targets(hdr); ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to jump
    // through "name_to_idx" and get "tid_to_idx"
    // NOLINTNEXTLINE(*-pointer-arithmetic)
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
  const std::size_t n_targets = get_n_targets(hdr);
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
  const std::size_t n_targets = get_n_targets(hdr);
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

struct mod_prob_stats {
  static constexpr auto max_mods = 10;
  static constexpr auto n_values = 256;
  // scratch
  std::array<hts_base_mod, max_mods> mods{};
  std::unique_ptr<hts_base_mod_state, void (*)(hts_base_mod_state *)> m;

  std::array<std::array<std::uint64_t, n_values>, n_nucs> methyl_fwd{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> methyl_rev{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> hydroxy_fwd{};
  std::array<std::array<std::uint64_t, n_values>, n_nucs> hydroxy_rev{};

  mod_prob_stats() : m{hts_base_mod_state_alloc(), &hts_base_mod_state_free} {};

  auto
  operator()(const bamxx::bam_rec &aln) {
    static constexpr auto h_idx = 0;
    static constexpr auto m_idx = 1;
    const auto qlen = get_l_qseq(aln);
    const auto seq = bam_get_seq(aln);
    const auto d = mods.data();
    const auto is_rev = bam_is_rev(aln);

    bam_parse_basemod(aln.b, m.get());
    // ADS: or bam_parse_basemod2(aln.b, m, HTS_MOD_REPORT_UNCHECKED)

    int n_types{};
    const auto types = bam_mods_recorded(m.get(), &n_types);
    // NOLINTNEXTLINE(*-pointer-arithmetic)
    if (n_types < 2 || (types[h_idx] != 'h' && types[m_idx] != 'm'))
      return;

    int pos{};
    int n{};
    while ((n = bam_next_basemod(aln.b, m.get(), d, max_mods, &pos)) > 0) {
      if (n < m_idx)
        continue;
      const auto other_nuc =
        is_rev ? (pos > 0 ? seq_nt16_str[bam_seqi(seq, pos - 1)] : '\0')
               : (pos + 1 < qlen ? seq_nt16_str[bam_seqi(seq, pos + 1)] : '\0');
      // NOLINTBEGIN(*-constant-array-index)
      const auto other_enc = encoding[static_cast<std::uint8_t>(other_nuc)];
      if (other_enc == n_nucs)
        continue;
      if (is_rev) {
        hydroxy_rev[other_enc][mods[h_idx].qual]++;
        methyl_rev[other_enc][mods[m_idx].qual]++;
      }
      else {
        hydroxy_fwd[other_enc][mods[h_idx].qual]++;
        methyl_fwd[other_enc][mods[m_idx].qual]++;
      }
      // NOLINTEND(*-constant-array-index)
    }
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(mod_prob_stats, methyl_fwd, methyl_rev,
                                 hydroxy_fwd, hydroxy_rev)
};

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
  std::int32_t n_threads{1};
  int strand{};
  std::string expected_basecall_model{};

  [[nodiscard]] std::string
  tostring() const {
    const auto strand_str = strand == 0 ? "both" : strand == 1 ? "fwd" : "rev";
    std::ostringstream oss;
    oss << std::boolalpha << "[verbose: " << verbose << "]\n"
        << "[show_progress: " << show_progress << "]\n"
        << "[compress_output: " << compress_output << "]\n"
        << "[include_header: " << include_header << "]\n"
        << "[require_covered: " << require_covered << "]\n"
        << "[symmetric: " << symmetric << "]\n"
        << "[cpg_only: " << cpg_only << "]\n"
        << "[force: " << force << "]\n"
        << "[n_threads: " << n_threads << "]\n"
        << "[strand: " << strand_str << "]\n"
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
    std::array<char, buf_size> buffer{};

    // Put chrom name in buffer and then skip that part for each site because
    // it doesn't change.
    const auto chrom_name = sam_hdr_tid2name(hdr.h, tid);
    if (chrom_name == nullptr)
      throw std::runtime_error("failed to identify chrom for tid: " +
                               std::to_string(tid));
    const auto chrom_name_offset =
      std::snprintf(buffer.data(), buf_size, "%s\t", chrom_name);
    if (chrom_name_offset < 0)
      throw std::runtime_error("failed to write to output buffer");
    auto buffer_after_chrom = buffer.data() + chrom_name_offset;

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
                                    tag_values[the_tag],  // NOLINT
                                    mods,
                                    n_reads,
                                    hydroxy,
                                    methyl);
        // clang-format on

        if (r < 0)
          throw std::runtime_error("failed to write to output buffer");
        out.write(buffer.data());
      }
    }
  }

  void
  write_output_sym(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
                   const std::int32_t tid, const std::string &chrom,
                   const std::vector<CountSet> &counts) const {
    static constexpr auto out_fmt = "%ld\t+\tCpG\t%.6g\t%d\t%.6g\t%.6g\n";
    static constexpr auto buf_size = 1024;
    std::array<char, buf_size> buffer{};

    // Put chrom name in buffer and then skip that part for each site because
    // it doesn't change.
    const auto chrom_name = sam_hdr_tid2name(hdr.h, tid);
    if (chrom_name == nullptr)
      throw std::runtime_error("failed to identify chrom for tid: " +
                               std::to_string(tid));
    const auto chrom_name_offset =
      std::snprintf(buffer.data(), buf_size, "%s\t", chrom_name);
    if (chrom_name_offset < 0)
      throw std::runtime_error("failed to write to output buffer");
    auto buffer_after_chrom = buffer.data() + chrom_name_offset;

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
          out.write(buffer.data());
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

  [[nodiscard]] mod_prob_stats
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
                << '\n';
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
                << (basecall_model.empty() ? "NA" : basecall_model) << "]\n";
    if (!expected_basecall_model.empty() &&
        basecall_model != expected_basecall_model) {
      std::cerr << "failed to match basecall model:\n"
                << "observed="
                << (basecall_model.empty() ? "NA" : basecall_model) << "\n"
                << "expected=" << expected_basecall_model_str() << '\n';
      return {};
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

    mod_prob_stats mps;
    mod_prob_buffer mod_buf;

    bool current_target_present = false;

    while (hts.read(hdr, aln)) {
      const std::int32_t tid = get_tid(aln);
      if (tid == -1)  // ADS: skip reads that have no tid -- they are not mapped
        continue;

      mps(aln);

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
          std::cerr << "processing " << sam_hdr_tid2name(hdr, tid) << '\n';

        // reset the counts
        counts.clear();
        if (current_target_present) {
          // set size of counts if current target is present
          chrom_itr = chroms_beg + chrom_idx->second;
          counts.resize(chrom_sizes[chrom_idx->second]);
        }
        prev_tid = tid;
      }
      if (current_target_present) {
        const bool is_rev = bam_is_rev(aln);
        if (is_rev && strand != 1)
          count_states_rev(aln, counts, mod_buf, *chrom_itr);
        if (!is_rev && strand != 2)
          count_states_fwd(aln, counts, mod_buf, *chrom_itr);
      }
    }
    if (!counts.empty())
      write_output(hdr, out, prev_tid, *chrom_itr, counts);

    // ADS: if some chroms might not be covered by reads, we have to iterate
    // over what remains
    if (!require_covered)
      for (auto i = prev_tid + 1; i < get_n_targets(hdr); ++i)
        output_skipped_chromosome(i, tid_to_idx, hdr, chroms_beg, chrom_sizes,
                                  counts, out);
    return mps;
  }
};

[[nodiscard]] static auto
valid_modification_types(const std::string &infile,
                         const std::uint32_t n_reads_to_check)
  -> std::pair<bool, std::string> {
  using mstate = hts_base_mod_state;
  std::unique_ptr<mstate, void (*)(mstate *)> m(hts_base_mod_state_alloc(),
                                                &hts_base_mod_state_free);
  bamxx::bam_in hts(infile);
  if (!hts)
    throw std::runtime_error("failed to open input file");
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::string message;

  std::uint32_t read_count{};
  bool valid_types{true};
  bamxx::bam_rec aln;
  for (; valid_types && hts.read(hdr, aln) && read_count < n_reads_to_check;
       ++read_count) {
    if (!bam_aux_get(aln.b, "MM") && !bam_aux_get(aln.b, "Mm"))
      continue;
    bam_parse_basemod(aln.b, m.get());
    // ADS: or bam_parse_basemod2(aln.b, m, HTS_MOD_REPORT_UNCHECKED)
    int n_types{};
    const auto types = bam_mods_recorded(m.get(), &n_types);
    // clang-format off
    // NOLINTBEGIN(*-pointer-arithmetic)
    valid_types = ((n_types == 0) ||
                   (n_types == 1 && types[0] == 'C') ||
                   (n_types >= 2 && types[0] == 'h' && types[1] == 'm'));
    // clang-format on
    if (!valid_types) {
      message = "n_types: " + std::to_string(n_types) + "\n";
      for (auto i = 0; i < n_types; ++i)
        message += "type[" + std::to_string(i) +
                   "]=" + std::to_string(static_cast<char>(types[i])) + "\n";
    }
    // NOLINTEND(*-pointer-arithmetic)
  }
  return std::make_pair(valid_types, message);
}

[[nodiscard]] static auto
check_modification_sites(const std::string &infile,
                         const std::uint32_t n_reads_to_check) -> bool {
  static constexpr auto max_mods = 10;
  std::array<hts_base_mod, max_mods> mods{};

  using mstate = hts_base_mod_state;
  std::unique_ptr<mstate, void (*)(mstate *)> m(hts_base_mod_state_alloc(),
                                                &hts_base_mod_state_free);

  bamxx::bam_in hts(infile);
  if (!hts)
    throw std::runtime_error("failed to open input file");
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::uint32_t read_count{};
  std::uint32_t reads_processed{};
  std::uint32_t only_cpgs_counter{};
  const auto d = mods.data();

  bamxx::bam_rec aln;
  for (; hts.read(hdr, aln) && read_count < n_reads_to_check; ++read_count) {
    if (!bam_aux_get(aln.b, "MM") && !bam_aux_get(aln.b, "Mm"))
      continue;

    const auto qlen = get_l_qseq(aln);
    const auto seq = bam_get_seq(aln);
    const auto is_rev = bam_is_rev(aln);

    bam_parse_basemod(aln.b, m.get());

    bool only_cpgs{true};
    int pos{};
    while (bam_next_basemod(aln.b, m.get(), d, max_mods, &pos) > 0) {
      const auto nuc = seq_nt16_str[bam_seqi(seq, pos)];
      const auto other_nuc =
        is_rev ? (pos > 0 ? seq_nt16_str[bam_seqi(seq, pos - 1)] : '\0')
               : (pos + 1 < qlen ? seq_nt16_str[bam_seqi(seq, pos + 1)] : '\0');
      if (!other_nuc)
        continue;
      only_cpgs = is_rev ? other_nuc == 'C' && nuc == 'G'
                         : nuc == 'C' && other_nuc == 'G';
    }
    ++reads_processed;
    only_cpgs_counter += only_cpgs;
  }
  return only_cpgs_counter == reads_processed;
}

int
main_nanocount(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  static constexpr auto n_reads_to_check = 1000;
  try {
    read_processor rp;

    std::string chroms_file;
    std::string outfile;
    std::string stats_file;

    // ADS: add a message to indicate that we assume either all modifications
    // are either called or no modifications are called, for each site.
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "get methylation levels from mapped nanopore reads",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "threads to use (few needed)", false,
                      rp.n_threads);
    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
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
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (rp.force)
      rp.expected_basecall_model.clear();

    if (rp.symmetric)
      rp.cpg_only = true;

    if (rp.n_threads <= 0)
      throw std::runtime_error("thread count cannot be zero");

    std::ostringstream cmd;
    std::copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    const auto [is_valid, message] =
      valid_modification_types(mapped_reads_file, n_reads_to_check);
    if (!is_valid) {
      std::cerr << "modification types are not valid. Violation:\n"
                << message << '\n';
      return EXIT_FAILURE;
    }

    const auto mods_at_cpgs =
      check_modification_sites(mapped_reads_file, n_reads_to_check);

    if (rp.verbose)
      std::cerr << std::boolalpha
                << "[input BAM/SAM file: " << mapped_reads_file << "]\n"
                << "[output file: " << outfile << "]\n"
                << "[genome file: " << chroms_file << "]\n"
                << "[mods only at CpGs: " << mods_at_cpgs << "]\n"
                << rp.tostring();

    // const auto [mc, mps] = rp(mapped_reads_file, outfile, chroms_file);
    const auto mps = rp(mapped_reads_file, outfile, chroms_file);
    if (!stats_file.empty()) {
      std::ofstream stats_out(stats_file);
      if (!stats_out)
        throw std::runtime_error("Error opening stats file: " + stats_file);
      stats_out << nlohmann::json(mps).dump(4) << "\n";
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)

#endif
