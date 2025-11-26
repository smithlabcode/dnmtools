/* summary: does bsrates, levels, sym, states and xsym all in one
 *
 * Copyright (C) 2025 Andrew D. Smith
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

/* ADS: This code writes MSite objects to files, but does not use
   MSite to do it. Possiby dangerous, but currently much faster. If
   MSite has a way to serialize into a char[] directly, then we should
   use it. */
// #include "MSite.hpp"
#include "LevelsCounter.hpp"
#include "bam_record_utils.hpp"
#include "bsutils.hpp"
#include "counts_header.hpp"

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include <bamxx.hpp>

#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

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
  data() {
    /* between c_str and insertion make sure to clear() */
    *pptr() = '\0';
    return pbase();
  }
};

static const char b2c[] = "TNGNNNCNNNNNNNNNNNNA";  // NOLINT

static void
collect_cpgs(const std::string &s, std::vector<std::uint64_t> &cpgs) {
  cpgs.clear();
  const std::uint64_t lim = std::size(s) - 1;
  for (auto i = 0u; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}

template <class BidirIt, class OutputIt>
// constexpr // since C++20
OutputIt
revcomp_copy(BidirIt first, BidirIt last, OutputIt d_first) {
  for (; first != last; ++d_first)
    *d_first = b2c[*(--last) - 'A'];  // NOLINT
  return d_first;
}

static bool
convert_meth_states_pos(const std::vector<std::uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bamxx::bam_rec &aln,
                        std::uint64_t &first_cpg_index, std::string &states) {
  states.clear();

  const std::uint64_t seq_start = get_pos(aln);
  const std::uint64_t width = rlen_from_cigar(aln);
  const std::uint64_t seq_end = seq_start + width;

  std::string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');

  if (std::size(seq_str) != width)
    throw std::runtime_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start
  auto cpg_itr = lower_bound(std::cbegin(cpgs), std::cend(cpgs), seq_start);
  auto first_cpg_itr = std::cend(cpgs);

  if (cpg_itr == std::cend(cpgs))
    return false;

  for (; cpg_itr != std::cend(cpgs) && *cpg_itr < seq_end; cpg_itr++) {
    const char x = seq_str[*cpg_itr - seq_start];
    states += (x == 'T') ? 'T' : ((x == 'C') ? 'C' : 'N');
    if (first_cpg_itr == std::cend(cpgs))
      first_cpg_itr = cpg_itr;
  }

  if (first_cpg_itr != std::cend(cpgs))
    first_cpg_index = distance(std::cbegin(cpgs), first_cpg_itr);

  return states.find_first_of("CT") != std::string::npos;
}

static bool
convert_meth_states_neg(const std::vector<std::uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bamxx::bam_rec &aln,
                        std::uint64_t &first_cpg_index, std::string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const std::uint64_t seq_start = get_pos(aln);
  const std::uint64_t width = rlen_from_cigar(aln);
  const std::uint64_t seq_end = seq_start + width;

  std::string orig_seq;
  get_seq_str(aln, orig_seq);

  std::string seq_str;
  seq_str.resize(orig_seq.size());
  revcomp_copy(std::cbegin(orig_seq), std::cend(orig_seq), std::begin(seq_str));
  apply_cigar(aln, seq_str, 'N');

  if (seq_str.size() != width)
    throw std::runtime_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start - 1
  // the -1 is because we look for G in the read corresponding to a
  // CpG in chromosome, which are indexed in cpgs based on the position of C
  auto cpg_itr = lower_bound(std::cbegin(cpgs), std::cend(cpgs),
                             seq_start > 0 ? seq_start - 1 : 0);
  auto first_cpg_itr = std::cend(cpgs);

  if (cpg_itr == std::cend(cpgs)) {
    return false;
  }
  else {
    for (; cpg_itr != std::cend(cpgs) && *cpg_itr < seq_end - 1; cpg_itr++) {
      const char x = seq_str[*cpg_itr - seq_start + 1];
      states += (x == 'G') ? 'C' : ((x == 'A') ? 'T' : 'N');
      if (first_cpg_itr == std::cend(cpgs))
        first_cpg_itr = cpg_itr;
    }
  }

  if (first_cpg_itr != std::cend(cpgs)) {
    first_cpg_index = distance(std::cbegin(cpgs), first_cpg_itr);
  }

  return states.find_first_of("CT") != std::string::npos;
}

struct bsrate_summary {
  // converted_count_positive is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the positive strand.
  std::uint64_t converted_count_positive{};
  // total_count_positive is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the positive
  // strand.
  std::uint64_t total_count_positive{};
  // bisulfite_conversion_rate_positive is equal to
  // converted_count_positive divided by total_count_positive, a value
  // that is always between 0 and 1. When total_count_positive is 0,
  // then bisulfite_conversion_rate_positive is given a value of 0.
  double bisulfite_conversion_rate_positive{};

  // converted_count_negative is the number of nucleotides covering a
  // cytosine in the reference that show a thymine in the read, and
  // for reads mapping to the negative strand.
  std::uint64_t converted_count_negative{};
  // total_count_negative is the number of nucleotides covering a
  // cytosine in the reference that show either a cytosine or a
  // thymine in the read, and for reads mapping to the negative
  // strand.
  std::uint64_t total_count_negative{};
  // bisulfite_conversion_rate_negative is equal to
  // converted_count_negative divided by total_count_negative, a value
  // that is always between 0 and 1. When total_count_negative is 0,
  // then bisulfite_conversion_rate_negative is given a value of 0.
  double bisulfite_conversion_rate_negative{};

  // converted_count is equal to the sum of converted_count_positive
  // and converted_count_negative.
  std::uint64_t converted_count{};
  // total_count is equal to the sum of total_count_positive and
  // total_count_negative
  std::uint64_t total_count{};
  // bisulfite_conversion_rate is equal to converted_count divided by
  // total_count, a value that is always between 0 and 1. When
  // total_count is 0, then bisulfite_conversion_rate is given a value
  // of 0.
  double bisulfite_conversion_rate{};

  // error_count is the number of nucleotides covering a cytosine in
  // the reference shows either an A or a G in the read.
  std::uint64_t error_count{};
  // valid_count is the number of nucleotides covering a cytosine in
  // the reference shows any nucleotide that is not an N in the read.
  std::uint64_t valid_count{};
  // error_rate is equal to error_count divided by valid_count, and is
  // a value that is always between 0 and 1. When valid_count is 0,
  // then error_rate is given a value of 0.
  double error_rate{};

  void
  update_pos(const char nt) {
    if (nt == 'C' || nt == 'T') {
      ++total_count_positive;
      converted_count_positive += (nt == 'T');
    }
    else if (nt != 'N')
      ++error_count;
  }

  void
  update_neg(const char nt) {
    if (nt == 'C' || nt == 'T') {
      ++total_count_negative;
      converted_count_negative += (nt == 'T');
    }
    else if (nt != 'N')
      ++error_count;
  }

  bsrate_summary &
  operator+=(const bsrate_summary &rhs) {
    // ADS: the "rates" are set to 0.0 here to ensure that they are
    // computed properly using the full data after any accumulation of
    // integer values has completed.
    converted_count_positive += rhs.converted_count_positive;
    total_count_positive += rhs.total_count_positive;
    bisulfite_conversion_rate_positive = 0.0;

    converted_count_negative += rhs.converted_count_negative;
    total_count_negative += rhs.total_count_negative;
    bisulfite_conversion_rate_negative = 0.0;

    converted_count += rhs.converted_count;
    total_count += rhs.total_count;
    bisulfite_conversion_rate = 0.0;

    error_count += rhs.error_count;
    valid_count += rhs.valid_count;
    error_rate = 0.0;
    return *this;
  }

  void
  set_values() {
    bisulfite_conversion_rate_positive =
      static_cast<double>(converted_count_positive) /
      static_cast<double>(
        std::max(total_count_positive, static_cast<std::uint64_t>(1)));

    bisulfite_conversion_rate_negative =
      static_cast<double>(converted_count_negative) /
      static_cast<double>(
        std::max(total_count_negative, static_cast<std::uint64_t>(1)));

    converted_count = converted_count_positive + converted_count_negative;
    total_count = total_count_positive + total_count_negative;

    bisulfite_conversion_rate =
      static_cast<double>(converted_count) /
      static_cast<double>(std::max(total_count, static_cast<std::uint64_t>(1)));

    valid_count = total_count + error_count;

    error_rate =
      static_cast<double>(error_count) /
      static_cast<double>(std::max(valid_count, static_cast<std::uint64_t>(1)));
  }

  std::string
  tostring_as_row() const {
    static constexpr auto precision_val = 5u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    // clang-format off
    oss << total_count_positive << '\t'
        << converted_count_positive << '\t'
        << bisulfite_conversion_rate_positive << '\t';
    oss << total_count_negative << '\t'
        << converted_count_negative << '\t'
        << bisulfite_conversion_rate_negative << '\t';
    oss << total_count << '\t'
        << converted_count << '\t'
        << bisulfite_conversion_rate << '\t';
    oss << error_count << '\t'
        << valid_count << '\t'
        << error_rate;
    // clang-format on
    return oss.str();
  }

  std::string
  tostring_as_yaml_list(const std::uint32_t position) const {
    static constexpr auto precision_val = 5u;
    std::ostringstream oss;
    oss.precision(precision_val);
    oss.setf(std::ios_base::fixed, std::ios_base::floatfield);
    // clang-format off
    oss << "  - base: " << position << '\n'
        << "    ptot: " << total_count_positive << '\n'
        << "    pconv: " << converted_count_positive << '\n'
        << "    prate: " << bisulfite_conversion_rate_positive << '\n'
        << "    ntot: " << total_count_negative << '\n'
        << "    nconv: " << converted_count_negative << '\n'
        << "    nrate: " << bisulfite_conversion_rate_negative << '\n'
        << "    bthtot: " << total_count << '\n'
        << "    bthconv: " << converted_count << '\n'
        << "    bthrate: " << bisulfite_conversion_rate << '\n'
        << "    err: " << error_count << '\n'
        << "    all: " << valid_count << '\n'
        << "    errrate: " << error_rate << '\n';
    // clang-format on
    return oss.str();
  }

  std::string
  tostring() const {
    static constexpr auto sep = ": ";
    std::ostringstream oss;
    oss << "converted_count_positive" << sep << converted_count_positive
        << '\n';

    oss << "total_count_positive" << sep << total_count_positive << '\n';

    oss << "bisulfite_conversion_rate_positive" << sep
        << bisulfite_conversion_rate_positive << '\n';

    oss << "converted_count_negative" << sep << converted_count_negative
        << '\n';

    oss << "total_count_negative" << sep << total_count_negative << '\n';

    oss << "bisulfite_conversion_rate_negative" << sep
        << bisulfite_conversion_rate_negative << '\n';

    oss << "converted_count" << sep << converted_count << '\n';
    oss << "total_count" << sep << total_count << '\n';

    oss << "bisulfite_conversion_rate" << sep << bisulfite_conversion_rate
        << '\n';

    oss << "error_count" << sep << error_count << '\n';
    oss << "valid_count" << sep << valid_count << '\n';
    oss << "error_rate" << sep << error_rate;

    return oss.str();
  }
};

inline bsrate_summary
operator+(bsrate_summary lhs, const bsrate_summary &rhs) {
  lhs += rhs;
  return lhs;
}

static void
write_summary(const std::string &summary_file,
              std::vector<bsrate_summary> &summaries) {
  auto s = reduce(std::cbegin(summaries), std::cend(summaries));
  s.set_values();

  std::ofstream out(summary_file);
  if (!out)
    throw std::runtime_error("failed to open file: " + summary_file);
  out << s.tostring() << '\n';
}

static std::pair<std::uint32_t, std::uint32_t>
count_states_pos(const bool INCLUDE_CPGS, const std::string &chrom,
                 const bamxx::bam_rec &aln,
                 std::vector<bsrate_summary> &summaries, std::size_t &hanging) {
  // NOLINTBEGIN(*-pointer-arithmetic)
  std::uint32_t n_conv = 0, n_uconv = 0;

  /* iterate through reference, query/read and fragment */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig =
    beg_cig + get_n_cigar(aln);  // NOLINT(*-pointer-arithmetic)
  auto rpos = get_pos(aln);
  auto qpos = 0;
  auto fpos = 0;

  const decltype(rpos) chrom_lim =
    static_cast<std::int64_t>(std::size(chrom) - 1);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);  // NOLINT(*-narrowing-conversions)
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = static_cast<std::int32_t>(qpos + n);
      for (; qpos < end_qpos; ++qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim)
          ++hanging;
        if (is_cytosine(chrom[rpos]) &&
            (rpos >= chrom_lim || !is_guanine(chrom[rpos + 1]) ||
             INCLUDE_CPGS)) {
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos)];
          summaries[fpos].update_pos(nt);
          n_conv += (nt == 'T');
          n_uconv += (nt == 'C');
        }
      }
    }
    else {
      if (cigar_eats_query(op))
        qpos += n;  // NOLINT
      if (cigar_eats_ref(op))
        rpos += n;  // NOLINT
      if (cigar_eats_frag(op))
        fpos += n;  // NOLINT
    }
  }
  assert(qpos == get_l_qseq(aln));
  return {n_conv, n_conv + n_uconv};
  // NOLINTEND(*-pointer-arithmetic)
}

static std::pair<std::uint32_t, std::uint32_t>
count_states_neg(const bool INCLUDE_CPGS, const std::string &chrom,
                 const bamxx::bam_rec &aln,
                 std::vector<bsrate_summary> &summaries, std::size_t &hanging) {
  // NOLINTBEGIN(*-pointer-arithmetic)
  std::uint32_t n_conv = 0, n_uconv = 0;

  /* iterate backward over query/read positions but forward over
     reference and fragment positions */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig = beg_cig + get_n_cigar(aln);
  auto rpos = get_pos(aln);
  auto qpos = get_l_qseq(aln);
  auto fpos = 0;

  const decltype(rpos) chrom_lim =
    static_cast<std::int64_t>(std::size(chrom) - 1);

  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const auto op = bam_cigar_op(*c_itr);
    const auto n = bam_cigar_oplen(*c_itr);
    if (cigar_eats_ref(op) && cigar_eats_query(op)) {
      const decltype(qpos) end_qpos = static_cast<std::int32_t>(qpos - n);
      for (; qpos > end_qpos; --qpos, ++rpos, ++fpos) {
        if (rpos > chrom_lim)
          ++hanging;
        if (is_guanine(chrom[rpos]) &&
            (rpos == 0 || !is_cytosine(chrom[rpos - 1]) || INCLUDE_CPGS)) {
          const auto nt = seq_nt16_str[bam_seqi(seq, qpos - 1)];
          summaries[fpos].update_neg(nt);
          n_conv += (nt == 'T');
          n_uconv += (nt == 'C');
        }
      }
    }
    else {
      if (cigar_eats_query(op))
        qpos -= static_cast<std::int32_t>(n);
      if (cigar_eats_ref(op))
        rpos += static_cast<std::int32_t>(n);
      if (cigar_eats_frag(op))
        fpos += static_cast<std::int32_t>(n);
    }
  }
  assert(qpos == 0);
  return {n_conv, n_conv + n_uconv};
  // NOLINTEND(*-pointer-arithmetic)
}

// ADS: we should never have to worry about coverage over > 32767 in
// any downstream analysis, so using "int16_t" here would allow to
// detect wrap around and report it as some kind of weird thing, maybe
// zeroing it and flagging the output. As it is now, if we really need
// >65535-fold coverage, we can make the change here.
typedef std::uint16_t count_type;

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
is_chh(const std::string &s, std::size_t i) {
  return (i < (s.length() - 2)) && is_cytosine(s[i]) && !is_guanine(s[i + 1]) &&
         !is_guanine(s[i + 2]);
}

static bool
is_ddg(const std::string &s, std::size_t i) {
  return (i < (s.length() - 2)) && !is_cytosine(s[i]) &&
         !is_cytosine(s[i + 1]) && is_guanine(s[i + 2]);
}

static bool
is_c_at_g(const std::string &s, std::size_t i) {
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
    std::ostringstream oss;
    oss << pA << '\t' << pC << '\t' << pG << '\t' << pT << '\t' << nA << '\t'
        << nC << '\t' << nG << '\t' << nT;
    // << '\t' << N; /* not used */
    return oss.str();
  }
  void
  add_count_pos(const char x) {
    if (x == 'T')
      ++pT;  // conditions ordered for efficiency
    else if (x == 'C')
      ++pC;
    else if (x == 'G')
      ++pG;
    else if (x == 'A')
      ++pA;
    // else ++N; /* not used */
  }
  void
  add_count_neg(const char x) {
    if (x == 'T')
      ++nT;  // conditions ordered for efficiency(??)
    else if (x == 'C')
      ++nC;
    else if (x == 'G')
      ++nG;
    else if (x == 'A')
      ++nA;
    // else ++N; /* not used */
  }
  count_type
  pos_total() const {
    return pA + pC + pG + pT;
  }
  count_type
  neg_total() const {
    return nA + nC + nG + nT;
  }

  count_type
  unconverted_cytosine() const {
    return pC;
  }
  count_type
  converted_cytosine() const {
    return pT;
  }
  count_type
  unconverted_guanine() const {
    return nC;
  }
  count_type
  converted_guanine() const {
    return nT;
  }

  count_type pA{0}, pC{0}, pG{0}, pT{0};
  count_type nA{0}, nC{0}, nG{0}, nT{0};
  // count_type N; /* this wasn't used and breaks alignment */
};

/* The "tag" returned by this function should be exclusive, so that
 * the order of checking conditions doesn't matter. There is also a
 * bit of a hack in that the unsigned "pos" could wrap, but this still
 * works as long as the chromosome size is not the maximum size of a
 * std::size_t.
 */
static std::uint32_t
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

static inline bool
is_cpg(const std::uint32_t x) {
  return x == 0;
}

static inline bool
is_chh(const std::uint32_t x) {
  return x == 1;
}

static inline bool
is_cxg(const std::uint32_t x) {
  return x == 2;
}

static inline bool
is_ccg(const std::uint32_t x) {
  return x == 3;
}

/* This "has_mutated" function looks on the opposite strand to see
 * if the apparent conversion from C->T was actually already in the
 * DNA because of a mutation or SNP.
 */
static bool
has_mutated(const char base, const CountSet &cs) {
  static const double MUTATION_DEFINING_FRACTION = 0.5;
  return is_cytosine(base)
           ? (cs.nG < MUTATION_DEFINING_FRACTION * (cs.neg_total()))
           : (cs.pG < MUTATION_DEFINING_FRACTION * (cs.pos_total()));
}

// NOLINTBEGIN(*-avoid-c-arrays)
static const char *const tag_values[] = {
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
// NOLINTEND(*-avoid-c-arrays)
static constexpr std::uint32_t MUT_OFFSET = 5;

static inline std::uint32_t
tag_with_mut(const std::uint32_t tag, const bool mut) {
  return tag + (mut ? MUT_OFFSET : 0);
}

// NOLINTBEGIN(cert-err58-cpp,*-avoid-non-const-global-variables)
LevelsCounter cpg("cpg");
LevelsCounter cpg_symmetric("cpg_symmetric");
LevelsCounter chh("chh");
LevelsCounter cxg("cxg");
LevelsCounter ccg("ccg");
LevelsCounter cytosines("cytosines");
// NOLINTEND(cert-err58-cpp,*-avoid-non-const-global-variables)

static void
update_lc(const bool mut, const std::uint32_t n_reads,
          const std::uint32_t n_meth, LevelsCounter &lc) {
  static constexpr auto alpha = 0.05;
  if (mut) {
    ++lc.mutations;
  }
  if (n_reads > 0) {
    ++lc.sites_covered;
    lc.max_depth = std::max(lc.max_depth, static_cast<std::uint64_t>(n_reads));
    lc.total_c += n_meth;
    lc.total_t += n_reads - n_meth;
    const double meth = n_meth / static_cast<double>(n_reads);
    lc.total_meth += meth;
    double lower{};
    double upper{};
    wilson_ci_for_binomial(alpha, n_reads, meth, lower, upper);
    lc.called_meth += (lower > 0.5);    // NOLINT(*-avoid-magic-numbers)
    lc.called_unmeth += (upper < 0.5);  // NOLINT(*-avoid-magic-numbers)
  }
  ++lc.total_sites;
}

static void
write_output(const bamxx::bam_header &hdr, bamxx::bgzf_file &out,
             bamxx::bgzf_file &out_sym, const int32_t tid,
             const std::string &chrom, const std::vector<CountSet> &counts) {
  quick_buf buf;  // keep underlying buffer space?
  buf.clear();

  buf << sam_hdr_tid2name(hdr, tid) << '\n';
  if (!out.write(buf.data(), buf.tellp()))
    throw std::runtime_error("error writing output");

  std::uint32_t prev_unconverted{};
  std::uint32_t prev_n_reads{};
  std::size_t prev_idx{};
  bool prev_is_c{};
  bool prev_mut{};
  std::uint32_t offset{};

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const char base = chrom[i];
    if (is_cytosine(base) || is_guanine(base)) {
      const std::uint32_t the_tag = get_tag_from_genome(chrom, i);
      const bool is_c = is_cytosine(base);
      const auto unconverted = is_c ? counts[i].unconverted_cytosine()
                                    : counts[i].unconverted_guanine();
      const auto converted =
        is_c ? counts[i].converted_cytosine() : counts[i].converted_guanine();
      const auto n_reads = static_cast<std::uint32_t>(unconverted + converted);

      const bool mut = has_mutated(base, counts[i]);

      update_lc(mut, n_reads, unconverted, cytosines);
      // cytosines.update(site);

      if (is_cpg(the_tag)) {
        update_lc(mut, n_reads, unconverted, cpg);
        if (prev_idx + 1 == i && !is_c && prev_is_c)
          update_lc(mut || prev_mut, n_reads + prev_n_reads,
                    unconverted + prev_unconverted, cpg_symmetric);
      }
      else if (is_chh(the_tag))
        update_lc(mut, n_reads, unconverted, chh);
      else if (is_ccg(the_tag))
        update_lc(mut, n_reads, unconverted, ccg);
      else if (is_cxg(the_tag))
        update_lc(mut, n_reads, unconverted, cxg);

      if (n_reads > 0) {
        buf.clear();
        // ADS: here is where we make an MSite, but not using MSite
        buf << i - offset << '\t' << unconverted << '\t' << converted << '\n';
        if (!out.write(buf.data(), buf.tellp()))
          throw std::runtime_error("error writing output");
        offset = i;
      }

      if (is_cpg(the_tag) && !is_c && prev_is_c) {
        buf.clear();
        // ADS: here is where we make an MSite, but not using MSite

        // NOLINTBEGIN(*-constant-array-index)
        buf << sam_hdr_tid2name(hdr, tid) << '\t' << prev_idx << '\t' << '+'
            << '\t' << tag_values[tag_with_mut(the_tag, mut || prev_mut)]
            << '\t'
            << ((n_reads + prev_n_reads) > 0
                  ? (unconverted + prev_unconverted) / (n_reads + prev_n_reads)
                  : 0.0)
            << '\t' << n_reads + prev_n_reads << '\n';
        // NOLINTEND(*-constant-array-index)
        if (!out_sym.write(buf.data(), buf.tellp()))
          throw std::runtime_error("error writing output");
      }
      prev_n_reads = n_reads;
      prev_unconverted = unconverted;
      prev_idx = i;
      prev_is_c = is_c;
      prev_mut = mut;
    }
  }
}

static void
count_states_pos(const bamxx::bam_rec &aln, std::vector<CountSet> &counts) {
  /* Move through cigar, reference and read positions without
     inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig =
    beg_cig + get_n_cigar(aln);  // NOLINT(*-pointer-arithmetic)
  auto rpos = get_pos(aln);
  auto qpos = 0;  // to match type with b->core.l_qseq
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);  // NOLINT(*-narrowing-conversions)
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const decltype(qpos) end_qpos = static_cast<std::int32_t>(qpos + n);
      for (; qpos < end_qpos; ++qpos) {
        // ADS: beware!!! bam_seqi is a macro, so no "qpos++" inside
        // its arguments! Why macros?!?! Just make sure the compiler
        // inliens it properly ffs!
        counts[rpos++].add_count_pos(seq_nt16_str[bam_seqi(seq, qpos)]);
      }
    }
    else if (eats_query(op)) {
      qpos += n;  // NOLINT
    }
    else if (eats_ref(op)) {
      rpos += n;  // NOLINT
    }
  }
  // ADS: somehow previous code included a correction for rpos going
  // past the end of the chromosome; this should result at least in a
  // soft-clip by any mapper. I'm not checking it here as even if it
  // happens I don't want to terminate.
  assert(qpos == get_l_qseq(aln));
}

static void
count_states_neg(const bamxx::bam_rec &aln, std::vector<CountSet> &counts) {
  /* Move through cigar, reference and (*backward*) through read
     positions without inflating cigar or read sequence */
  const auto seq = bam_get_seq(aln);
  const auto beg_cig = bam_get_cigar(aln);
  const auto end_cig =
    beg_cig + get_n_cigar(aln);  // NOLINT(*-pointer-arithmetic)
  std::size_t rpos = get_pos(aln);
  std::size_t qpos = get_l_qseq(aln);  // to match type with b->core.l_qseq
  for (auto c_itr = beg_cig; c_itr != end_cig; ++c_itr) {
    const char op = bam_cigar_op(*c_itr);  // NOLINT(*-narrowing-conversions)
    const std::uint32_t n = bam_cigar_oplen(*c_itr);
    if (eats_ref(op) && eats_query(op)) {
      const std::size_t end_qpos = qpos - n;  // to match type with qpos
      for (; qpos > end_qpos; --qpos)         // beware ++ in macro below!!!
        counts[rpos++].add_count_neg(seq_nt16_str[bam_seqi(seq, qpos - 1)]);
    }
    else if (eats_query(op)) {
      qpos -= n;  // NOLINT
    }
    else if (eats_ref(op)) {
      rpos += n;  // NOLINT
    }
  }
  /* qpos is unsigned; would wrap around if < 0 */
  // ADS: Same as count_states_pos; see comment there
  assert(qpos == 0);
}

static std::unordered_map<int32_t, std::size_t>
get_tid_to_idx(
  const bamxx::bam_header &hdr,
  const std::unordered_map<std::string, std::size_t> &name_to_idx) {
  std::unordered_map<int32_t, std::size_t> tid_to_idx;
  for (int32_t i = 0; i < hdr.h->n_targets; ++i) {
    // "curr_name" gives a "tid_to_name" mapping allowing to jump
    // through "name_to_idx" and get "tid_to_idx"
    const std::string curr_name(
      hdr.h->target_name[i]);  // NOLINT(*-pointer-arithmetic)
    const auto name_itr(name_to_idx.find(curr_name));
    if (name_itr == std::cend(name_to_idx))
      throw std::runtime_error("failed to find chrom: " + curr_name);
    tid_to_idx[i] = name_itr->second;
  }
  return tid_to_idx;
}

template <class CS>
static void
output_skipped_chromosome(
  const int32_t tid, const std::unordered_map<int32_t, std::size_t> &tid_to_idx,
  const bamxx::bam_header &hdr,
  const std::vector<std::string>::const_iterator chroms_beg,
  const std::vector<std::size_t> &chrom_sizes, std::vector<CS> &counts,
  bamxx::bgzf_file &out, bamxx::bgzf_file &out_sym) {
  // get the index of the next chrom sequence
  const auto chrom_idx = tid_to_idx.find(tid);
  if (chrom_idx == std::cend(tid_to_idx))
    throw std::runtime_error("chrom not found: " + sam_hdr_tid2name(hdr, tid));

  const auto chrom_itr =
    chroms_beg + chrom_idx->second;  // NOLINT(*-narrowing-conversions)

  // reset the counts
  counts.clear();
  counts.resize(chrom_sizes[chrom_idx->second]);

  write_output(hdr, out, out_sym, tid, *chrom_itr, counts);
}

static bool
consistent_targets(const bamxx::bam_header &hdr,
                   const std::unordered_map<int32_t, std::size_t> &tid_to_idx,
                   const std::vector<std::string> &names,
                   const std::vector<std::size_t> &sizes) {
  const std::int32_t n_targets = hdr.h->n_targets;
  if (n_targets != std::size(names))
    return false;

  for (std::int32_t tid = 0; tid < n_targets; ++tid) {
    const std::string tid_name_sam = sam_hdr_tid2name(hdr, tid);  // NOLINT
    const std::size_t tid_size_sam = sam_hdr_tid2len(hdr, tid);   // NOLINT
    const auto idx_itr = tid_to_idx.find(tid);
    if (idx_itr == std::cend(tid_to_idx))
      return false;
    const auto idx = idx_itr->second;
    if (tid_name_sam != names[idx] || tid_size_sam != sizes[idx])
      return false;
  }
  return true;
}

static void
process_reads(const bool VERBOSE, const bool show_progress,
              const bool compress_output, const bool include_header,
              const std::int32_t n_threads, const std::string &infile,
              const std::string &chroms_file, const std::string &counts_outfile,
              const std::string &sym_outfile, const std::string &levels_outfile,
              const std::string &bsrate_outfile,
              const std::string &states_outfile) {
  // assumed maximum length of a fragment for bsrate
  static constexpr const std::size_t output_size = 10000;
  std::size_t hanging = 0;

  // first get the chromosome names and sequences from the FASTA file
  std::vector<std::string> chroms, names;
  read_fasta_file_short_names(chroms_file, names, chroms);
  for (auto &i : chroms)
    transform(std::cbegin(i), std::cend(i), std::begin(i),
              [](const char c) { return std::toupper(c); });
  if (VERBOSE)
    std::cerr << "[n chroms in reference: " << chroms.size() << "]\n";
  const auto chroms_beg = cbegin(chroms);

  std::unordered_map<std::string, std::size_t> name_to_idx;
  std::vector<std::size_t> chrom_sizes(chroms.size(), 0);
  for (std::size_t i = 0; i < chroms.size(); ++i) {
    name_to_idx[names[i]] = i;
    chrom_sizes[i] = chroms[i].size();
  }

  std::vector<bsrate_summary> summaries(output_size);

  bamxx::bam_tpool tp(n_threads);  // Must be destroyed after hts

  // open the hts SAM/BAM input file and get the header
  bamxx::bam_in hts(infile);
  if (!hts)
    throw std::runtime_error("failed to open input file");
  // load the input file's header
  bamxx::bam_header hdr(hts);
  if (!hdr)
    throw std::runtime_error("failed to read header");

  std::unordered_map<int32_t, std::size_t> tid_to_idx =
    get_tid_to_idx(hdr, name_to_idx);

  if (!consistent_targets(hdr, tid_to_idx, names, chrom_sizes))
    throw std::runtime_error("inconsistent reference genome information");

  // open the output file
  const std::string output_mode = compress_output ? "w" : "wu";
  bamxx::bgzf_file out(counts_outfile, output_mode);
  if (!out)
    throw std::runtime_error("error opening output file: " + counts_outfile);

  // open the sym output file
  bamxx::bgzf_file out_sym(sym_outfile, output_mode);
  if (!out_sym)
    throw std::runtime_error("error opening output file: " + sym_outfile);

  quick_buf buf;  // keep underlying buffer space?
  buf.clear();

  // open the states output file
  bamxx::bgzf_file out_states(states_outfile, output_mode);
  if (!out_states)
    throw std::runtime_error("error opening output file: " + states_outfile);

  // set the threads for the input file decompression
  if (n_threads > 1) {
    tp.set_io(hts);
    tp.set_io(out);
    tp.set_io(out_sym);
  }

  if (include_header)
    write_counts_header_from_bam_header(hdr, out);

  // now iterate over the reads, switching chromosomes and writing
  // output as needed
  bamxx::bam_rec aln;
  int32_t prev_tid = -1;

  // this is where all the counts are accumulated
  std::vector<CountSet> counts;
  std::vector<std::uint64_t> cpgs;

  std::vector<std::string>::const_iterator chrom_itr;

  while (hts.read(hdr, aln)) {
    if (get_tid(aln) == prev_tid) {
      // do the work for this mapped read, depending on strand
      if (bam_is_rev(aln))
        count_states_neg(aln, counts);
      else
        count_states_pos(aln, counts);
    }
    else {  // chrom changes, output results, get the next one
      // write output if there is any; should fail only once
      if (!counts.empty())
        write_output(hdr, out, out_sym, prev_tid, *chrom_itr, counts);

      const int32_t tid = get_tid(aln);
      if (tid < prev_tid)
        throw std::runtime_error("SAM file is not sorted");

      for (auto i = prev_tid + 1; i < tid; ++i)
        output_skipped_chromosome(i, tid_to_idx, hdr, chroms_beg, chrom_sizes,
                                  counts, out, out_sym);

      // get the next chrom to process
      auto chrom_idx(tid_to_idx.find(tid));
      if (chrom_idx == std::cend(tid_to_idx))
        throw std::runtime_error("chromosome not found: " +
                                 std::string(sam_hdr_tid2name(hdr, tid)));
      if (show_progress)
        std::cerr << "processing " << sam_hdr_tid2name(hdr, tid) << '\n';

      prev_tid = tid;
      chrom_itr =
        chroms_beg + chrom_idx->second;  // NOLINT(*-narrowing-conversions)

      // reset the counts
      counts.clear();
      counts.resize(chrom_sizes[chrom_idx->second]);
      // use_this_chrom = seq_to_use.empty() || chrom_idx == chrom_idx_to_use;

      collect_cpgs(*chrom_itr, cpgs);
    }

    static constexpr auto use_this_chrom = true;

    if (use_this_chrom) {
      static constexpr auto reads_are_a_rich = false;
      // do the work for the current mapped read
      if (reads_are_a_rich)
        flip_conversion(aln);
      static constexpr auto INCLUDE_CPGS = false;
      if (bam_is_rev(aln))
        count_states_neg(INCLUDE_CPGS, *chrom_itr, aln, summaries, hanging);
      else
        count_states_pos(INCLUDE_CPGS, *chrom_itr, aln, summaries, hanging);
    }

    std::uint64_t first_cpg_index = std::numeric_limits<std::uint64_t>::max();
    std::string seq;

    const bool has_cpgs =
      bam_is_rev(aln)
        ? convert_meth_states_neg(cpgs, hdr, aln, first_cpg_index, seq)
        : convert_meth_states_pos(cpgs, hdr, aln, first_cpg_index, seq);

    if (has_cpgs) {
      buf.clear();
      buf << sam_hdr_tid2name(hdr, aln) << '\t' << first_cpg_index << '\t'
          << seq << '\n';
      if (!out_states.write(buf.data(), buf.tellp())) {
        throw std::runtime_error{"failure writing states output"};
      }
    }
  }
  if (!counts.empty())
    write_output(hdr, out, out_sym, prev_tid, *chrom_itr, counts);

  // ADS: if some chroms might not be covered by reads, we have to
  // iterate over what remains
  for (auto i = prev_tid + 1; i < hdr.h->n_targets; ++i)
    output_skipped_chromosome(i, tid_to_idx, hdr, chroms_beg, chrom_sizes,
                              counts, out, out_sym);

  // before writing the output we need to ensure the derived
  // quantities in bsrate_summary have been calculated
  std::for_each(std::begin(summaries), std::end(summaries),
                [](bsrate_summary &s) { s.set_values(); });
  write_summary(bsrate_outfile, summaries);

  std::ofstream out_levels(levels_outfile);
  if (!out_levels)
    throw std::runtime_error{"bad output file: " + levels_outfile};
  out_levels << cytosines << '\n'
             << cpg << '\n'
             << cpg_symmetric << '\n'
             << chh << '\n'
             << ccg << '\n'
             << cxg << '\n';
}

int
main_summary(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    bool VERBOSE = false;
    bool show_progress = false;
    bool require_covered = false;
    bool compress_output = false;
    bool include_header = false;

    std::string chroms_file;
    std::string levels_outfile;
    std::string bsrate_outfile;
    std::string counts_outfile;
    std::string sym_outfile;
    std::string states_outfile;
    std::int32_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "do bsrate, levels, sym, xsym, xcouns and states",
                           "-c <chroms> <mapped-reads>");
    opt_parse.add_opt("threads", 't', "threads to use (few needed)", false,
                      n_threads);
    opt_parse.add_opt("levels-out", '\0', "levels output file", true,
                      levels_outfile);
    opt_parse.add_opt("bsrate-out", '\0', "bsrate output file", true,
                      bsrate_outfile);
    opt_parse.add_opt("counts-out", '\0', "counts output file", true,
                      counts_outfile);
    opt_parse.add_opt("sym-out", '\0', "sym output file", true, sym_outfile);
    opt_parse.add_opt("states-out", '\0', "states output file", true,
                      states_outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome file (FASTA format)",
                      true, chroms_file);
    opt_parse.add_opt("require-covered", 'r', "only output covered sites",
                      false, require_covered);
    opt_parse.add_opt("header", 'H', "add a header to identify the reference",
                      false, include_header);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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

    if (n_threads < 0)
      throw std::runtime_error("thread count cannot be negative");

    std::ostringstream cmd;
    std::copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    if (VERBOSE)
      std::cerr << "[input BAM/SAM file: " << mapped_reads_file << "]\n"
                << "[counts file: " << counts_outfile << "]\n"
                << "[sym file:    " << sym_outfile << "]\n"
                << "[bsrate file: " << bsrate_outfile << "]\n"
                << "[levels file: " << levels_outfile << "]\n"
                << "[states file: " << states_outfile << "]\n"
                << "[genome file: " << chroms_file << "]\n"
                << "[threads requested: " << n_threads << "]\n"
                << "[command line: \"" << cmd.str() << "\"]\n";

    process_reads(VERBOSE, show_progress, compress_output, include_header,
                  n_threads, mapped_reads_file, chroms_file, counts_outfile,
                  sym_outfile, levels_outfile, bsrate_outfile, states_outfile);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
