/* dmr: computes DMRs based on HMRs and probability of differences at
 * individual CpGs
 *
 * Copyright (C) 2012-2023 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include "GenomicRegion.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cassert>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::find_if;
using std::from_chars;
using std::ifstream;
using std::max;
using std::pair;
using std::runtime_error;
using std::string;
using std::vector;

static bool
parse_methdiff_line(const char *c, const char *c_end, string &chrom,
                    uint32_t &pos, char &strand, string &context,
                    double &diffscore, uint32_t &n_meth_a, uint32_t &n_unmeth_a,
                    uint32_t &n_meth_b, uint32_t &n_unmeth_b) {
  constexpr auto is_sep = [](const char x) { return x == ' ' || x == '\t'; };
  constexpr auto not_sep = [](const char x) { return x != ' ' && x != '\t'; };

  // NOLINTBEGIN(*-pointer-arithmetic)
  auto field_s = c;
  auto field_e = std::find_if(field_s + 1, c_end, is_sep);
  bool failed = field_e == c_end;

  // chromosome name
  {
    const uint32_t d = std::distance(field_s, field_e);
    chrom = string{field_s, d};
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
    const uint32_t d = std::distance(field_s, field_e);
    context = string{field_s, d};
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

static vector<MSite>
read_diffs_file(const string &diffs_file) {
  bamxx::bgzf_file in(diffs_file, "r");
  if (!in)
    throw runtime_error("could not open file: " + diffs_file);

  string chrom, name;
  char strand{};
  double diffscore{};
  uint32_t pos{}, meth_a{}, unmeth_a{}, meth_b{}, unmeth_b{};

  vector<MSite> cpgs;
  string line;
  while (getline(in, line)) {
    // NOLINTBEGIN(*-pointer-arithmetic)
    if (!parse_methdiff_line(line.data(), line.data() + size(line), chrom, pos,
                             strand, name, diffscore, meth_a, unmeth_a, meth_b,
                             unmeth_b))
      throw runtime_error("bad methdiff line: " + line);
    // NOLINTEND(*-pointer-arithmetic)

    cpgs.emplace_back(chrom, pos, strand, name, diffscore, 1);
  }
  return cpgs;
}

static void
complement_regions(const size_t max_end, const vector<GenomicRegion> &a,
                   const size_t start, const size_t end,
                   vector<GenomicRegion> &cmpl) {
  cmpl.push_back(GenomicRegion(a[start]));
  cmpl.back().set_start(0);
  for (size_t i = start; i < end; ++i) {
    cmpl.back().set_end(a[i].get_start());
    cmpl.push_back(GenomicRegion(a[i]));
    cmpl.back().set_start(a[i].get_end());
  }
  cmpl.back().set_end(max_end);
}

static vector<size_t>
get_chrom_ends(const vector<GenomicRegion> &r) {
  vector<size_t> ends;
  for (size_t i = 0; i < r.size() - 1; ++i)
    if (!r[i].same_chrom(r[i + 1]))
      ends.push_back(i + 1);
  ends.push_back(r.size());
  return ends;
}

static vector<GenomicRegion>
complement_regions(const size_t max_end, const vector<GenomicRegion> &r) {
  vector<size_t> r_chroms = get_chrom_ends(r);
  vector<GenomicRegion> r_cmpl;
  size_t t = 0;
  for (size_t i = 0; i < size(r_chroms); ++i) {
    complement_regions(max_end, r, t, r_chroms[i], r_cmpl);
    t = r_chroms[i];
  }
  return r_cmpl;
}

static bool
check_no_overlap(const vector<GenomicRegion> &regions) {
  for (size_t i = 1; i < regions.size(); ++i)
    if (regions[i].same_chrom(regions[i - 1]) &&
        regions[i].get_start() < regions[i - 1].get_end())
      return false;
  return true;
}

static inline MSite
get_left_msite(const GenomicRegion &r) {
  return {r.get_chrom(), r.get_start(), r.get_strand(), r.get_name(), 0.0, 1u};
}

static inline MSite
get_right_msite(const GenomicRegion &r) {
  return {r.get_chrom(), r.get_end(), r.get_strand(), r.get_name(), 0.0, 1u};
}

static vector<pair<size_t, size_t>>
separate_sites(const vector<GenomicRegion> &dmrs, const vector<MSite> &sites) {
  vector<pair<size_t, size_t>> sep_sites;
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
              const vector<MSite> &cpgs, const size_t start_idx,
              const size_t end_idx, size_t &total_cpgs, size_t &total_sig) {
  total_cpgs = end_idx - start_idx;
  for (size_t i = start_idx; i < end_idx; ++i) {
    const auto pval = pval_from_msite(cpgs[i]);
    if ((LOW_CUTOFF && (pval < sig_cutoff)) ||
        (!LOW_CUTOFF && (pval > 1.0 - sig_cutoff)))
      ++total_sig;
  }
}

int
main_dmr(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const string description =
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
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << '\n'
           << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 5) {  // NOLINT(*-avoid-magic-numbers)
      cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const string diffs_file = leftover_args[0];
    const string hmr1_file = leftover_args[1];
    const string hmr2_file = leftover_args[2];
    const string outfile_a = leftover_args[3];
    const string outfile_b = leftover_args[4];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[LOADING HMRS] " << hmr1_file << '\n';

    vector<GenomicRegion> regions_a;
    ReadBEDFile(hmr1_file, regions_a);
    assert(check_sorted(regions_a));
    if (!check_sorted(regions_a))
      throw runtime_error("regions not sorted in file: " + hmr1_file);
    if (!check_no_overlap(regions_a))
      throw runtime_error("regions overlap in file: " + hmr1_file);

    if (VERBOSE)
      cerr << "[LOADING HMRS] " << hmr2_file << '\n';

    vector<GenomicRegion> regions_b;
    ReadBEDFile(hmr2_file, regions_b);
    assert(check_sorted(regions_b));
    if (!check_sorted(regions_b))
      throw runtime_error("regions not sorted in file: " + hmr2_file);
    if (!check_no_overlap(regions_b))
      throw runtime_error("regions overlap in file: " + hmr2_file);

    if (VERBOSE)
      cerr << "[COMPUTING SYMMETRIC DIFFERENCE]" << '\n';

    size_t max_end = 0;
    for (const auto &r : regions_a)
      max_end = max(max_end, r.get_end());  // cppcheck-suppress useStlAlgorithm
    for (const auto &r : regions_b)
      max_end = max(max_end, r.get_end());  // cppcheck-suppress useStlAlgorithm

    const auto a_cmpl = complement_regions(max_end, regions_a);
    const auto b_cmpl = complement_regions(max_end, regions_b);

    vector<GenomicRegion> dmrs_a, dmrs_b;
    genomic_region_intersection_by_base(regions_a, b_cmpl, dmrs_a);
    genomic_region_intersection_by_base(regions_b, a_cmpl, dmrs_b);

    // separate the regions by chrom and by desert
    if (VERBOSE)
      cerr << "[READING CPG METH DIFFS]" << '\n';
    const auto cpgs = read_diffs_file(diffs_file);
    if (VERBOSE)
      cerr << "[read " << size(cpgs) << " sites from " + diffs_file << "]"
           << '\n';

    if (!check_sorted(cpgs))
      throw runtime_error("CpGs not sorted in: " + diffs_file);
    if (VERBOSE)
      cerr << "[TOTAL CPGS]: " << cpgs.size() << '\n';

    auto sep_sites = separate_sites(dmrs_a, cpgs);

    for (size_t i = 0; i < dmrs_a.size(); ++i) {
      size_t total_cpgs = 0, total_sig = 0;
      get_cpg_stats(true, sig_cutoff, cpgs, sep_sites[i].first,
                    sep_sites[i].second, total_cpgs, total_sig);
      dmrs_a[i].set_name(dmrs_a[i].get_name() + ":" + toa(total_cpgs));
      dmrs_a[i].set_score(static_cast<float>(total_sig));
    }

    sep_sites = separate_sites(dmrs_b, cpgs);

    for (size_t i = 0; i < dmrs_b.size(); ++i) {
      size_t total_cpgs = 0, total_sig = 0;
      get_cpg_stats(false, sig_cutoff, cpgs, sep_sites[i].first,
                    sep_sites[i].second, total_cpgs, total_sig);
      dmrs_b[i].set_name(dmrs_b[i].get_name() + ":" + toa(total_cpgs));
      dmrs_b[i].set_score(static_cast<float>(total_sig));
    }

    std::ofstream out_a(outfile_a);
    std::copy(std::cbegin(dmrs_a), std::cend(dmrs_a),
              std::ostream_iterator<GenomicRegion>(out_a, "\n"));

    std::ofstream out_b(outfile_b);
    std::copy(std::cbegin(dmrs_b), std::cend(dmrs_b),
              std::ostream_iterator<GenomicRegion>(out_b, "\n"));

    if (VERBOSE)
      cerr << "[OUTPUT FORMAT] COL4=NAME:N_COVERED_CPGS COL5=N_SIG_CPGS\n";
  }
  catch (const std::exception &e) {
    cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
