/* roimethstat: average methylation in each of a set of regions
 *
 * Copyright (C) 2014-2025 Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Masaru Nakajima
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

#include "Interval6.hpp"
#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "xcounts_utils.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

static void
update(LevelsCounter &lc, const xcounts_entry &xse) {
  static constexpr auto one_half = 0.5;
  const std::uint64_t n_reads = xse.n_meth + xse.n_unmeth;
  if (n_reads > 0) {
    ++lc.sites_covered;
    lc.max_depth = std::max(lc.max_depth, n_reads);
    lc.total_c += xse.n_meth;
    lc.total_t += xse.n_unmeth;
    const auto meth = static_cast<double>(xse.n_unmeth) / n_reads;
    lc.total_meth += meth;
    double lower{};
    double upper{};
    wilson_ci_for_binomial(lc.alpha, n_reads, meth, lower, upper);
    lc.called_meth += (lower > one_half);
    lc.called_unmeth += (upper < one_half);
  }
  ++lc.total_sites;
}

static void
process_chrom(const bool report_more_info, const char level_code,
              const std::vector<Interval6> &intervals,
              const std::vector<xcounts_entry> &sites, std::ostream &out) {
  std::uint64_t j = 0;
  for (auto i = 0ul; i < std::size(intervals); ++i) {
    while (j < std::size(sites) && sites[j].pos < intervals[i].start)
      ++j;
    LevelsCounter lc;
    while (j < std::size(sites) && sites[j].pos < intervals[i].stop)
      update(lc, sites[j++]);

    Interval6 interval(intervals[i]);
    interval.score =
      level_code == 'w'
        ? lc.mean_meth_weighted()
        : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());
    interval.name +=
      "_" + std::to_string(
              level_code == 'w'
                ? lc.coverage()
                : (level_code == 'u' ? lc.sites_covered : lc.total_called()));
    out << to_string(interval);
    if (report_more_info)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

static void
process_chrom(const bool report_more_info,
              const std::vector<Interval6> &intervals, std::ostream &out) {
  LevelsCounter lc;
  const std::string lc_formatted = format_levels_counter(lc);
  for (const auto &r : intervals) {
    out << to_string(r);
    if (report_more_info)
      out << '\t' << lc_formatted;
    out << '\n';
  }
}

static void
process_from_xcounts(const std::uint32_t n_threads, const bool report_more_info,
                     const char level_code, const std::string &xsym_file,
                     const std::vector<Interval6> &intervals_in,
                     std::ostream &out) {
  const auto sites_by_chrom = read_xcounts_by_chrom(n_threads, xsym_file);
  // const auto intervals = get_Interval6s(intervals_file);

  std::vector<std::vector<Interval6>> intervals_by_chrom;
  std::string prev_chrom;
  for (auto i = 0u; i < std::size(intervals_in); ++i) {
    if (intervals_in[i].chrom != prev_chrom) {
      intervals_by_chrom.push_back(std::vector<Interval6>());
      prev_chrom = intervals_in[i].chrom;
    }
    intervals_by_chrom.back().push_back(intervals_in[i]);
  }

  for (const auto &intervals : intervals_by_chrom) {
    const auto chrom_name = intervals.front().chrom;
    const auto sites = sites_by_chrom.find(chrom_name);
    if (sites != std::cend(sites_by_chrom))
      process_chrom(report_more_info, level_code, intervals, sites->second,
                    out);
    else
      process_chrom(report_more_info, intervals, out);
  }
}

[[nodiscard]] bool
cmp_within_chrom(const Interval6 &r1, const Interval6 &r2) {
  return r1.start < r2.start ||
         (r1.start == r2.start &&
          (r1.stop < r2.stop ||
           (r1.stop == r2.stop && (r1.strand < r2.strand))));
}

[[nodiscard]] bool
cmp_within_chrom(const MSite &s1, const MSite &s2) {
  return s1.pos < s2.pos;
}

template <class T>
[[nodiscard]] typename T::const_iterator
get_chrom_order(const T &order, const MSite &s) {
  return order.find(s.chrom);
}

template <class T>
[[nodiscard]] typename T::const_iterator
get_chrom_order(const T &order, const Interval6 &r) {
  return order.find(r.chrom);
}

struct cmp_chrom_order {
  explicit cmp_chrom_order(
    const std::unordered_map<std::string, std::uint32_t> &m) : order{m} {}

  template <class T>
  [[nodiscard]] bool
  operator()(const T &a, const T &b) const {
    const auto ac = get_chrom_order(order, a);
    if (ac == std::cend(order))
      return false;
    const auto bc = get_chrom_order(order, b);
    if (bc == std::cend(order))
      return true;
    const auto c = static_cast<int>(ac->second) - static_cast<int>(bc->second);
    return c < 0 || (c == 0 && cmp_within_chrom(a, b));
  }

  const std::unordered_map<std::string, std::uint32_t> &order;  // NOLINT
};

[[nodiscard]] static std::vector<std::string>
get_chroms(const std::string &filename) {
  // this function might look slow but it seems not to bottleneck
  constexpr auto not_in_chrom = [](const std::uint8_t c) {
    return c == ' ' || c == '\t' || c == '\n';
  };

  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("cannot open file: " + filename);

  std::string prev_chrom;
  std::string line;
  std::vector<std::string> chroms_seen;
  while (getline(in, line)) {
    const auto chrom_end =
      std::find_if(std::cbegin(line), std::cend(line), not_in_chrom);
    if (chrom_end == std::cend(line))
      throw std::runtime_error("bad counts line: " + line);
    if (!std::equal(std::cbegin(prev_chrom), std::cend(prev_chrom),
                    std::cbegin(line), chrom_end)) {
      const std::string chrom{std::cbegin(line), chrom_end};
      auto x = find(std::cbegin(chroms_seen), std::cend(chroms_seen), chrom);
      if (x != std::cend(chroms_seen))
        throw std::runtime_error("chroms not sorted");
      chroms_seen.push_back(chrom);
      prev_chrom = chrom;
    }
  }
  return chroms_seen;
}

template <class ForwardIt, class T>
[[nodiscard]] static inline std::pair<ForwardIt, ForwardIt>
region_bounds(const std::unordered_map<std::string, std::uint32_t> &chrom_order,
              ForwardIt first, ForwardIt last, const T &region) {
  // ADS: in the sites below, a and b, the strand is set to '+'
  // always. Otherwise the ordering on MSite would force strand to be
  // used. The consequence is that sites on the positive strand would not be
  // included if they begin at the same nucleotide where a region starts if
  // that region is on the negative strand const char
  // strand(region.get_strand());

  const std::string chrom(region.chrom);
  const MSite a(chrom, region.start, '+', "", 0, 0);
  const MSite b(chrom, region.stop, '+', "", 0, 0);

  // ADS: This function seems like std::equal_range but both elements
  // below use "lower_bound" because "b" is strictly following valid
  // elements. Otherwise we would include in our range all elements
  // comparing equal to "b", which are invalid.
  const cmp_chrom_order cmp{chrom_order};
  return {lower_bound(first, last, a, cmp), lower_bound(first, last, b, cmp)};
}

[[nodiscard]] static bool
is_sorted_within_chrom(const std::vector<MSite> &sites) {
  constexpr auto pos_cmp = [](const MSite &s1, const MSite &s2) {
    return s1.pos < s2.pos;
  };
  auto a = std::cbegin(sites);
  while (a != std::cend(sites)) {
    const auto a_chrom = a->chrom;
    const auto b = find_if(a, std::cend(sites), [&a_chrom](const MSite &s) {
      return s.chrom != a_chrom;
    });
    if (!std::is_sorted(a, b, pos_cmp))
      return false;
    a = b;
  }
  return true;
}

[[nodiscard]] static std::vector<MSite>
read_sites(const std::string &filename) {
  std::vector<MSite> sites;
  bamxx::bgzf_file in(filename, "r");
  if (in) {
    MSite s;
    while (read_site(in, s))
      sites.push_back(s);
  }
  return sites;
}

static void
process_preloaded(
  const bool verbose, const bool report_more_info, const char level_code,
  const std::string &sites_file,
  const std::unordered_map<std::string, std::uint32_t> &chrom_order,
  const std::vector<Interval6> &regions, std::ostream &out) {
  const auto sites = read_sites(sites_file);
  if (sites.empty())
    throw std::runtime_error("failed to read sites: " + sites_file);

  if (!is_sorted_within_chrom(sites))
    throw std::runtime_error("sites not sorted: " + sites_file);

  if (verbose)
    std::cerr << "[n_sites=" << std::size(sites) << "]\n";

  for (const auto &r : regions) {
    LevelsCounter lc;
    const auto b =
      region_bounds(chrom_order, std::cbegin(sites), std::cend(sites), r);
    for (auto c = b.first; c != b.second; ++c)
      lc.update(*c);

    const double score =
      level_code == 'w'
        ? lc.mean_meth_weighted()
        : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());
    Interval6 r_scored{r};
    r_scored.score = score;
    out << to_string(r_scored);
    if (report_more_info)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

template <class T>
[[nodiscard]] static inline bool
not_after(const T &chrom_order, const MSite &site,
          const std::uint32_t &chrom_idx, const std::size_t end_pos) {
  const auto sc_itr = chrom_order.find(site.chrom);
  if (sc_itr == std::cend(chrom_order))
    throw std::runtime_error("lookup failure: " + site.chrom);
  const auto sc = sc_itr->second;
  return (sc == chrom_idx && site.pos < end_pos) || sc < chrom_idx;
}

[[nodiscard]] static inline bool
at_or_after(const MSite &site, const std::string &chrom,
            const std::size_t start_pos) {
  // not using dictionary in this function because equality of chrom
  // can be tested using the strings regardless of index
  return (start_pos <= site.pos && site.chrom == chrom);
}

[[nodiscard]] static LevelsCounter
calc_site_stats(
  std::ifstream &sites_in, const Interval6 &region,
  const std::unordered_map<std::string, std::uint32_t> &chrom_order) {
  const std::string chrom{region.chrom};
  const auto chrom_itr = chrom_order.find(chrom);
  if (chrom_itr == std::cend(chrom_order))
    throw std::runtime_error("lookup failure: " + chrom);
  const auto chrom_idx = chrom_itr->second;
  const auto start_pos = region.start;
  const auto end_pos = region.stop;
  find_offset_for_msite(chrom_order, chrom, start_pos, sites_in);

  LevelsCounter lc;

  MSite site;
  while (sites_in >> site && not_after(chrom_order, site, chrom_idx, end_pos))
    if (at_or_after(site, chrom, start_pos))
      lc.update(site);

  sites_in.clear();  // clear here otherwise an eof means we skip later
                     // intervals, which could be valid
  return lc;
}

static void
process_on_disk(
  const bool report_more_info, const char level_code,
  const std::string &sites_file,
  const std::unordered_map<std::string, std::uint32_t> &chrom_order,
  const std::vector<Interval6> &regions, std::ostream &out) {
  std::ifstream in(sites_file);
  if (!in)
    throw std::runtime_error("failed to open file: " + sites_file);

  for (const auto &region : regions) {
    const auto lc = calc_site_stats(in, region, chrom_order);
    const double score =
      level_code == 'w'
        ? lc.mean_meth_weighted()
        : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());
    Interval6 r{region};
    r.score = score;
    out << to_string(r);
    if (report_more_info)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

[[nodiscard]] static std::size_t
get_bed_columns(const std::string &regions_file) {
  std::ifstream in(regions_file);
  if (!in)
    throw std::runtime_error("cannot open file: " + regions_file);
  std::string line;
  getline(in, line);
  std::istringstream iss(line);
  std::string token;
  std::size_t n_columns = 0;
  while (iss >> token)
    ++n_columns;
  return n_columns;
}

int
main_roimethstat(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const std::string description = R"""(
Compute average site methylation levels in each interval from a given
set of genomic intervals. The 5th column (the "score" column in BED
format) is determined by the '-l' or '-level' argument.

Columns (beyond the first 6) in the BED format output:
(7) weighted mean methylation
(8) unweighted mean methylation
(9) fractional methylation
(10) number of sites in the region
(11) number of sites covered at least once
(12) number of observations in reads indicating methylation
(13) total number of observations from reads in the region
)""";
    constexpr auto valid_n_cols = [](const auto n_cols) {
      return n_cols != 3 && n_cols < 6;  // NOLINT(*-avoid-magic-numbers)
    };

    static const std::string default_name_prefix = "X";

    bool verbose = false;
    bool print_numeric_only = false;
    bool preload = false;
    bool report_more_info = false;
    bool sort_data_if_needed = false;
    bool allow_extra_fields = false;
    std::uint32_t n_threads = 1;

    std::string level_code = "w";

    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<intervals-bed> <methylation-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "output file", true, outfile);
    opt_parse.add_opt("numeric", 'N', "print numeric values only (not NAs)",
                      false, print_numeric_only);
    opt_parse.add_opt("preload", 'L', "load all site sites", false, preload);
    opt_parse.add_opt("sort", 's', "sort data if needed", false,
                      sort_data_if_needed);
    opt_parse.add_opt("level", 'l',
                      "the level to report as score column "
                      "in bed format output (w, u or f)",
                      false, level_code);
    opt_parse.add_opt("more-levels", 'M', "report more methylation information",
                      false, report_more_info);
    opt_parse.add_opt("threads", 't', "threads to use (if input compressed)",
                      false, n_threads);
    opt_parse.add_opt("relaxed", '\0',
                      "input has extra fields (used for nanopore)", false,
                      allow_extra_fields);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << opt_parse.about_message_raw()
                << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message_raw() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_FAILURE;
    }
    if (std::size(leftover_args) != 2) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_FAILURE;
    }
    if (level_code != "w" && level_code != "u" && level_code != "f") {
      std::cerr << "selected level must be in {w, u, f}\n";
      return EXIT_FAILURE;
    }
    const std::string regions_file = leftover_args.front();
    const std::string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    MSite::no_extra_fields = (allow_extra_fields == false);

    const bool is_xcounts = get_is_xcounts_file(sites_file);
    if (!is_msite_file(sites_file) && !is_xcounts)
      throw std::runtime_error("dnmtools counts or xcounts format required: " +
                               sites_file);

    // make a map that specifies their order; otherwise we can't
    // ensure regions are sorted in the same way
    std::unordered_map<std::string, std::uint32_t> chrom_order;
    if (!is_xcounts)
      for (const auto &i : get_chroms(sites_file))
        chrom_order.emplace(i, std::size(chrom_order));

    if (verbose)
      std::cerr << "loading regions" << '\n';

    if (!std::filesystem::is_regular_file(regions_file))
      // otherwise we could not read the file twice
      throw std::runtime_error("regions file must be regular file");

    const auto n_columns = get_bed_columns(regions_file);
    if (valid_n_cols(n_columns))
      throw std::runtime_error("format must be 3 or 6+ column bed: " +
                               regions_file);
    if (is_msite_file(regions_file))
      throw std::runtime_error("seems to be a counts file: " + regions_file +
                               "\ncheck order of input arguments");

    auto regions = read_intervals6(regions_file);
    if (!std::is_sorted(std::begin(regions), std::end(regions))) {
      if (sort_data_if_needed) {
        if (verbose)
          std::cerr << "sorting regions\n";
        std::sort(std::begin(regions), std::end(regions),
                  cmp_chrom_order(chrom_order));
      }
      else
        throw std::runtime_error("not sorted: " + regions_file +
                                 " (consider -s)");
    }

    // if columns are 3 we name the regions otherwise that column will
    // be empty and the output format will be wrong
    if (n_columns == 3)
      for (auto i = 0u; i < std::size(regions); ++i)
        regions[i].name = default_name_prefix + std::to_string(i);

    if (verbose)
      std::cerr << "[n_regions=" << std::size(regions) << "]\n";

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open outfile: " + outfile);

    if (is_xcounts)
      process_from_xcounts(n_threads, report_more_info, level_code[0],
                           sites_file, regions, out);
    else if (preload)
      process_preloaded(verbose, report_more_info, level_code[0], sites_file,
                        chrom_order, regions, out);
    else
      process_on_disk(report_more_info, level_code[0], sites_file, chrom_order,
                      regions, out);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
