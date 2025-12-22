/* Copyright (C) 2012-2025 Andrew D. Smith
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

[[maybe_unused]] static constexpr auto about = R"(
dmr: compute DMRs based on HMRs and probability of differences at each site
)";

#include "Interval6.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <numeric>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

static bool
parse_methdiff_line(const char *c, const std::size_t c_sz, std::string &chrom,
                    std::uint32_t &pos, char &strand, std::string &context,
                    double &diffscore, std::uint32_t &n_meth_a,
                    std::uint32_t &n_unmeth_a, std::uint32_t &n_meth_b,
                    std::uint32_t &n_unmeth_b) {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  // NOLINTBEGIN(*-pointer-arithmetic)
  const auto c_end = c + c_sz;
  auto field_s = c;
  auto field_e = std::find_if(field_s + 1, c_end, is_sep);
  bool failed = field_e == c_end;

  // chromosome name
  {
    const std::uint32_t d = std::distance(field_s, field_e);
    chrom = std::string{field_s, d};
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || field_e == c_end;

  // position
  {
    const auto [ptr, ec] = std::from_chars(field_s, field_e, pos);
    failed = failed || ec != std::errc();
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  // below because strand is 1 base wide
  failed = failed || field_e != field_s + 1 || field_e == c_end;

  // strand
  strand = *field_s;
  failed = failed || (strand != '-' && strand != '+');

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || field_e == c_end;

  // context
  {
    const std::uint32_t d = std::distance(field_s, field_e);
    context = std::string{field_s, d};
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || field_e == c_end;

  // score for difference in methylation (contingency table p-value)
  {
#ifdef __APPLE__
    const int ret = std::sscanf(field_s, "%lf", &diffscore);
    failed = failed || ret < 1;
#else
    const auto [ptr, ec] = std::from_chars(field_s, field_e, diffscore);
    failed = failed || ec != std::errc();
#endif
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || (field_e == c_end);

  // counts methylated in methylome "a"
  {
    const auto [ptr, ec] = std::from_chars(field_s, c_end, n_meth_a);
    failed = failed || ec != std::errc();
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || field_e == c_end;

  // counts unmethylated in methylome "a"
  {
    const auto [ptr, ec] = std::from_chars(field_s, c_end, n_unmeth_a);
    failed = failed || ec != std::errc();
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);
  field_e = std::find_if(field_s + 1, c_end, is_sep);
  failed = failed || field_e == c_end;

  // counts methylated in methylome "b"
  {
    const auto [ptr, ec] = std::from_chars(field_s, c_end, n_meth_b);
    failed = failed || ec != std::errc();
  }

  field_s = std::find_if(field_e + 1, c_end, not_sep);

  // counts unmethylated in methylome "a"
  {
    const auto [ptr, ec] = std::from_chars(field_s, c_end, n_unmeth_b);
    // final field needs to fail if we haven't reached the end
    failed = failed || ec != std::errc() || ptr != c_end;
  }

  // NOLINTEND(*-pointer-arithmetic)

  return !failed;
}

static std::vector<MSite>
read_diffs_file(const std::string &diffs_file) {
  bamxx::bgzf_file in(diffs_file, "r");
  if (!in)
    throw std::runtime_error("could not open file: " + diffs_file);

  std::string chrom, name;
  char strand{};
  double diffscore{};
  std::uint32_t pos{}, meth_a{}, unmeth_a{}, meth_b{}, unmeth_b{};

  std::vector<MSite> cpgs;
  std::string line;
  while (getline(in, line)) {
    // NOLINTBEGIN(*-pointer-arithmetic)
    if (!parse_methdiff_line(line.data(), std::size(line), chrom, pos, strand,
                             name, diffscore, meth_a, unmeth_a, meth_b,
                             unmeth_b))
      throw std::runtime_error("bad methdiff line: " + line);
    // NOLINTEND(*-pointer-arithmetic)

    cpgs.emplace_back(chrom, pos, strand, name, diffscore, 1);
  }
  return cpgs;
}

static void
complement(const std::size_t max_end, const std::vector<Interval6> &a,
           const std::size_t start, const std::size_t end,
           std::vector<Interval6> &cmpl) {
  cmpl.push_back(Interval6(a[start]));
  cmpl.back().start = 0;
  for (std::size_t i = start; i < end; ++i) {
    cmpl.back().stop = a[i].start;
    cmpl.push_back(Interval6(a[i]));
    cmpl.back().start = a[i].stop;
  }
  cmpl.back().stop = max_end;
}

static std::vector<std::size_t>
get_chrom_ends(const std::vector<Interval6> &r) {
  std::vector<std::size_t> ends;
  for (std::size_t i = 0; i + 1 < std::size(r); ++i)
    if (r[i].chrom != r[i + 1].chrom)
      ends.push_back(i + 1);
  ends.push_back(std::size(r));
  return ends;
}

static std::vector<Interval6>
complement(const std::size_t max_end, const std::vector<Interval6> &r) {
  std::vector<std::size_t> r_chroms = get_chrom_ends(r);
  std::vector<Interval6> r_cmpl;
  std::size_t t{};
  for (std::size_t i = 0; i < std::size(r_chroms); ++i) {
    complement(max_end, r, t, r_chroms[i], r_cmpl);
    t = r_chroms[i];
  }
  return r_cmpl;
}

static bool
check_no_overlap(const std::vector<Interval6> &regions) {
  for (std::size_t i = 1; i < std::size(regions); ++i)
    if (regions[i].chrom == regions[i - 1].chrom &&
        regions[i].start < regions[i - 1].stop)
      return false;
  return true;
}

static inline MSite
get_left_msite(const Interval6 &r) {
  return {r.chrom, r.start, r.strand, r.name, 0.0, 1u};
}

static inline MSite
get_right_msite(const Interval6 &r) {
  return {r.chrom, r.stop, r.strand, r.name, 0.0, 1u};
}

static std::vector<std::pair<std::size_t, std::size_t>>
separate_sites(const std::vector<Interval6> &dmrs,
               const std::vector<MSite> &sites) {
  std::vector<std::pair<std::size_t, std::size_t>> sep_sites;
  for (const auto &dmr : dmrs) {
    const auto a = get_left_msite(dmr);
    const auto b = get_right_msite(dmr);
    const auto a_insert = lower_bound(std::cbegin(sites), std::cend(sites), a);
    const auto b_insert = lower_bound(std::cbegin(sites), std::cend(sites), b);
    sep_sites.emplace_back(distance(std::cbegin(sites), a_insert),
                           distance(std::cbegin(sites), b_insert));
  }
  return sep_sites;
}

static inline double
pval_from_msite(const MSite &s) {
  return s.meth;  // abused as a p-value here
}

static void
get_cpg_stats(const bool LOW_CUTOFF, const double sig_cutoff,
              const std::vector<MSite> &cpgs, const std::size_t start_idx,
              const std::size_t end_idx, std::size_t &total_cpgs,
              std::size_t &total_sig) {
  total_cpgs = end_idx - start_idx;
  for (std::size_t i = start_idx; i < end_idx; ++i) {
    const auto pval = pval_from_msite(cpgs[i]);
    if ((LOW_CUTOFF && (pval < sig_cutoff)) ||
        (!LOW_CUTOFF && (pval > 1.0 - sig_cutoff)))
      ++total_sig;
  }
}

template <class T>
[[nodiscard]] auto
intersection_by_base(const std::vector<T> &a,
                     const std::vector<T> &b) -> std::vector<T> {
  const auto end_precedes = [](const T &x, const T &y) {
    const auto cmp = x.chrom.compare(y.chrom);
    return cmp < 0 || (cmp == 0 && x.stop < y.stop);
  };
  const auto overlaps = [](const T &x, const T &y) {
    const auto cmp = x.chrom.compare(y.chrom);
    return cmp == 0 && std::max(x.start, x.start) < std::min(x.stop, y.stop);
  };

  auto a_itr = std::cbegin(a);
  auto a_end = std::cend(a);
  auto b_itr = std::cbegin(b);
  auto b_end = std::cend(b);

  std::vector<T> c;
  while (a_itr != a_end && b_itr != b_end) {
    if (overlaps(*a_itr, *b_itr))
      c.emplace_back(a_itr->chrom, std::max(a_itr->start, b_itr->start),
                     std::min(a_itr->stop, b_itr->stop), std::string{}, 0.0,
                     '+');
    if (end_precedes(*a_itr, *b_itr))
      ++a_itr;
    else
      ++b_itr;
  }
  return c;
}

int
main_dmr(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const std::string description =
      R"(
computes DMRs based on HMRs and probability of differences at
individual CpGs";
)";

    bool VERBOSE = false;
    double sig_cutoff = 0.05;  // NOLINT(*-avoid-magic-numbers)

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description,
                           "<methdiffs_1_gt_2> <hmr_1> <hmr_2> "
                           "<dmr_1_lt_2> <dmr_2_lt_1>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("cutoff", 'c', "Significance cutoff", false, sig_cutoff);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (std::size(leftover_args) != 5) {  // NOLINT(*-avoid-magic-numbers)
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string diffs_file = leftover_args[0];
    const std::string hmr1_file = leftover_args[1];
    const std::string hmr2_file = leftover_args[2];
    const std::string outfile_a = leftover_args[3];
    const std::string outfile_b = leftover_args[4];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      std::cerr << "[LOADING HMRS] " << hmr1_file << '\n';

    const auto regions_a = read_intervals6(hmr1_file);
    if (!std::is_sorted(std::cbegin(regions_a), std::cend(regions_a)))
      throw std::runtime_error("regions not sorted in file: " + hmr1_file);
    if (!check_no_overlap(regions_a))
      throw std::runtime_error("regions overlap in file: " + hmr1_file);

    if (VERBOSE)
      std::cerr << "[LOADING HMRS] " << hmr2_file << '\n';

    const auto regions_b = read_intervals6(hmr2_file);
    if (!std::is_sorted(std::cbegin(regions_b), std::cend(regions_b)))
      throw std::runtime_error("regions not sorted in file: " + hmr2_file);
    if (!check_no_overlap(regions_b))
      throw std::runtime_error("regions overlap in file: " + hmr2_file);

    if (VERBOSE)
      std::cerr << "[COMPUTING SYMMETRIC DIFFERENCE]" << '\n';

    const auto get_max_stop = [](const auto &x) {
      const auto acc = [](const auto s, const auto &r) {
        return std::max(s, r.stop);
      };
      return std::accumulate(std::cbegin(x), std::cend(x), 0u, acc);
    };
    const auto max_stop =
      std::max(get_max_stop(regions_a), get_max_stop(regions_b));

    const auto a_cmpl = complement(max_stop, regions_a);
    const auto b_cmpl = complement(max_stop, regions_b);

    auto dmrs_a = intersection_by_base(regions_a, b_cmpl);
    auto dmrs_b = intersection_by_base(regions_b, a_cmpl);

    // separate the regions by chrom and by desert
    if (VERBOSE)
      std::cerr << "[READING CPG METH DIFFS]" << '\n';
    const auto cpgs = read_diffs_file(diffs_file);
    if (VERBOSE)
      std::cerr << "[read " << std::size(cpgs) << " sites from " + diffs_file
                << "]\n";

    if (!std::is_sorted(std::cbegin(cpgs), std::cend(cpgs)))
      throw std::runtime_error("CpGs not sorted in: " + diffs_file);
    if (VERBOSE)
      std::cerr << "[TOTAL CPGS]: " << std::size(cpgs) << '\n';

    auto sep_sites = separate_sites(dmrs_a, cpgs);

    for (std::size_t i = 0; i < std::size(dmrs_a); ++i) {
      std::size_t total_cpgs{}, total_sig{};
      get_cpg_stats(true, sig_cutoff, cpgs, sep_sites[i].first,
                    sep_sites[i].second, total_cpgs, total_sig);
      dmrs_a[i].name += ":" + std::to_string(total_cpgs);
      dmrs_a[i].score = static_cast<double>(total_sig);
    }

    sep_sites = separate_sites(dmrs_b, cpgs);

    for (std::size_t i = 0; i < std::size(dmrs_b); ++i) {
      std::size_t total_cpgs{}, total_sig{};
      get_cpg_stats(false, sig_cutoff, cpgs, sep_sites[i].first,
                    sep_sites[i].second, total_cpgs, total_sig);
      dmrs_b[i].name += ":" + std::to_string(total_cpgs);
      dmrs_b[i].score = static_cast<double>(total_sig);
    }

    std::ofstream out_a(outfile_a);
    if (!out_a)
      throw std::runtime_error("failed to open output file: " + outfile_a);
    std::copy(std::cbegin(dmrs_a), std::cend(dmrs_a),
              std::ostream_iterator<Interval6>(out_a, "\n"));

    std::ofstream out_b(outfile_b);
    if (!out_b)
      throw std::runtime_error("failed to open output file: " + outfile_b);
    std::copy(std::cbegin(dmrs_b), std::cend(dmrs_b),
              std::ostream_iterator<Interval6>(out_b, "\n"));

    if (VERBOSE)
      std::cerr << "[OUTPUT FORMAT] COL4=NAME:N_COVERED_CPGS COL5=N_SIG_CPGS\n";
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
