/* roimethstat: average methylation in each of a set of regions
 *
 * Copyright (C) 2014-2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Masaru Nakajima
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

#include <bamxx.hpp>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <filesystem>

#include "GenomicRegion.hpp"
#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "bsutils.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::distance;
using std::endl;
using std::ifstream;
using std::is_sorted;
using std::pair;
using std::runtime_error;
using std::string;
using std::to_string;
using std::unordered_map;
using std::vector;

using bamxx::bgzf_file;

namespace fs = std::filesystem;

bool
cmp_within_chrom(const GenomicRegion &r1, const GenomicRegion &r2) {
  return (r1.get_start() < r2.get_start() ||
          (r1.get_start() == r2.get_start() &&
           (r1.get_end() < r2.get_end() ||
            (r1.get_end() == r2.get_end() &&
             (r1.get_strand() < r2.get_strand())))));
}

bool
cmp_within_chrom(const MSite &s1, const MSite &s2) {
  return s1.pos < s2.pos;
}

template<class T> typename T::const_iterator
get_chrom_order(const T &order, const MSite &s) {
  return order.find(s.chrom);
}

template<class T> typename T::const_iterator
get_chrom_order(const T &order, const GenomicRegion &r) {
  return order.find(r.get_chrom());
}

struct cmp_chrom_order {
  cmp_chrom_order(const unordered_map<string, uint32_t> &m): order{m} {}

  template<class T> bool operator()(const T &a, const T &b) const {
    const auto ac = get_chrom_order(order, a);
    if (ac == cend(order)) return false;
    const auto bc = get_chrom_order(order, b);
    if (bc == cend(order)) return true;
    const auto c = static_cast<int>(ac->second) - static_cast<int>(bc->second);
    return c < 0 || (c == 0 && cmp_within_chrom(a, b));
  }

  const unordered_map<string, uint32_t> &order;
};

static vector<string>
get_chroms(const string &filename) {
  // this function might look slow but it seems not to bottleneck
  constexpr auto not_in_chrom = [](const uint8_t c) {
    return c == ' ' || c == '\t' || c == '\n';
  };

  bgzf_file in(filename, "r");
  if (!in) throw runtime_error("cannot open file: " + filename);

  string prev_chrom;
  string line;
  vector<string> chroms_seen;
  while (getline(in, line)) {
    const auto chrom_end = find_if(cbegin(line), cend(line), not_in_chrom);
    if (chrom_end == cend(line))
      throw runtime_error("bad counts line: " + line);
    if (!equal(cbegin(prev_chrom), cend(prev_chrom), cbegin(line), chrom_end)) {
      const string chrom{cbegin(line), chrom_end};
      auto x = find(cbegin(chroms_seen), cend(chroms_seen), chrom);
      if (x != cend(chroms_seen)) throw runtime_error("chroms not sorted");
      chroms_seen.push_back(chrom);
      prev_chrom = std::move(chrom);
    }
  }
  return chroms_seen;
}

template <class ForwardIt, class T>
static inline pair<ForwardIt, ForwardIt>
region_bounds(const unordered_map<string, uint32_t> &chrom_order,
              ForwardIt first, ForwardIt last, const T &region) {
  // ADS: in the sites below, a and b, the strand is set to '+'
  // always. Otherwise the ordering on MSite would force strand to be
  // used. The consequence is that sites on the positive strand would
  // not be included if they begin at the same nucleotide where a
  // region starts if that region is on the negative strand
  // const char strand(region.get_strand());

  const string chrom(region.get_chrom());
  const MSite a(chrom, region.get_start(), '+', "", 0, 0);
  const MSite b(chrom, region.get_end(), '+', "", 0, 0);

  // ADS: This function seems like std::equal_range but both elements
  // below use "lower_bound" because "b" is strictly following valid
  // elements. Otherwise we would include in our range all elements
  // comparing equal to "b", which are invalid.
  const cmp_chrom_order cmp{chrom_order};
  return {lower_bound(first, last, a, cmp), lower_bound(first, last, b, cmp)};
}

static string
format_levels_counter(const LevelsCounter &lc) {
  // ...
  // (7) weighted mean methylation
  // (8) unweighted mean methylation
  // (9) fractional methylation
  // (10) number of sites in the region
  // (11) number of sites covered at least once
  // (12) number of observations in reads indicating methylation
  // (13) total number of observations from reads in the region
  std::ostringstream oss;
  // clang-format off
  oss << lc.mean_meth_weighted() << '\t'
      << lc.mean_meth() << '\t'
      << lc.fractional_meth() << '\t'
      << lc.total_sites << '\t'
      << lc.sites_covered << '\t'
      << lc.total_c << '\t'
      << (lc.total_c + lc.total_t);
  // clang-format on
  return oss.str();
}

static bool
is_sorted_within_chrom(const vector<MSite> &sites) {
  constexpr auto pos_cmp = [](const MSite &s1, const MSite &s2) {
    return s1.pos < s2.pos;
  };
  auto a = cbegin(sites);
  while (a != cend(sites)) {
    const auto b = find_if(a, cend(sites),
                           [a](const MSite &s) { return s.chrom != a->chrom; });
    if (!is_sorted(a, b, pos_cmp)) return false;
    a = b;
  }
  return true;
}

static vector<MSite>
read_sites(const string &filename) {
  vector<MSite> sites;
  bgzf_file in(filename, "r");
  if (in) {
    MSite s;
    while (read_site(in, s)) sites.push_back(s);
  }
  return sites;
}

static void
process_preloaded(const bool VERBOSE, const bool report_more_information,
                  const char level_code, const string &sites_file,
                  const unordered_map<string, uint32_t> &chrom_order,
                  const vector<GenomicRegion> &regions, std::ostream &out) {

  const auto sites = read_sites(sites_file);
  if (sites.empty()) throw runtime_error("failed to read sites: " + sites_file);

  if (!is_sorted_within_chrom(sites))
    throw runtime_error("sites not sorted: " + sites_file);

  if (VERBOSE) cerr << "[n_sites=" << sites.size() << "]" << endl;

  for (auto &r : regions) {
    LevelsCounter lc;
    const auto b = region_bounds(chrom_order, cbegin(sites), cend(sites), r);
    for (auto c = b.first; c != b.second; ++c) lc.update(*c);

    const double score =
      level_code == 'w'
        ? lc.mean_meth_weighted()
        : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());
    GenomicRegion r_scored{r};
    r_scored.set_score(score);
    out << r_scored;
    if (report_more_information)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

template<class T> static inline bool
not_after(const T &chrom_order, const MSite &site, const uint32_t &chrom_idx,
          const size_t end_pos) {
  const auto sc_itr = chrom_order.find(site.chrom);
  if (sc_itr == cend(chrom_order))
    throw runtime_error("lookup failure: " + site.chrom);
  const auto sc = sc_itr->second;
  return (sc == chrom_idx && site.pos < end_pos) || sc < chrom_idx;
}

static inline bool
at_or_after(const MSite &site, const string &chrom, const size_t start_pos) {
  // not using dictionary in this function because equality of chrom
  // can be tested using the strings regardless of index
  return (start_pos <= site.pos && site.chrom == chrom);
}

static LevelsCounter
calc_site_stats(ifstream &sites_in, const GenomicRegion &region,
                const unordered_map<string, uint32_t> &chrom_order) {
  const string chrom{region.get_chrom()};
  const auto chrom_itr = chrom_order.find(chrom);
  if (chrom_itr == cend(chrom_order))
    throw runtime_error("lookup failure: " + chrom);
  const auto chrom_idx = chrom_itr->second;
  const auto start_pos = region.get_start();
  const auto end_pos = region.get_end();
  find_offset_for_msite(chrom_order, chrom, start_pos, sites_in);

  LevelsCounter lc;

  MSite site;
  while (sites_in >> site && not_after(chrom_order, site, chrom_idx, end_pos))
    if (at_or_after(site, chrom, start_pos)) lc.update(site);

  sites_in.clear();  // clear here otherwise an eof means we skip later
                     // intervals, which could be valid
  return lc;
}

static void
process_on_disk(const bool report_more_information, const char level_code,
                const string &sites_file,
                const unordered_map<string, uint32_t> &chrom_order,
                const vector<GenomicRegion> &regions, std::ostream &out) {
  ifstream in(sites_file);
  if (!in) throw runtime_error("failed to open file: " + sites_file);

  for (auto &region : regions) {
    const auto lc = calc_site_stats(in, region, chrom_order);
    const double score =
      level_code == 'w'
        ? lc.mean_meth_weighted()
        : (level_code == 'u' ? lc.mean_meth() : lc.fractional_meth());
    GenomicRegion r{region};
    r.set_score(score);
    out << r;
    if (report_more_information)
      out << '\t' << format_levels_counter(lc);
    out << '\n';
  }
}

static size_t
get_bed_columns(const string &regions_file) {
  ifstream in(regions_file);
  if (!in) throw runtime_error("failed to open file: " + regions_file);

  string line;
  getline(in, line);

  std::istringstream iss(line);
  string token;
  size_t n_columns = 0;
  while (iss >> token) ++n_columns;

  return n_columns;
}


int
main_roimethstat(int argc, const char **argv) {
  try {
    // ADS: this information should be somehow presented to the user
    // when the tool is run without arguments.
    static const string description = R"""(

Compute average site methylation levels in each interval from a given
set of genomic intervals

Columns (beyond the first 6) in the BED format output:
(7) weighted mean methylation
(8) unweighted mean methylation
(9) fractional methylation
(10) number of sites in the region
(11) number of sites covered at least once
(12) number of observations in reads indicating methylation
(13) total number of observations from reads in the region

)""";

    static const string default_name_prefix = "X";

    bool VERBOSE = false;
    bool print_numeric_only = false;
    bool preload = false;
    bool report_more_information = false;
    bool sort_data_if_needed = false;

    string level_code = "w";

    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "Compute average site "
                           "methylation levels in each interval from "
                           "a given set of genomic intervals",
                           "<intervals-bed> <methylation-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
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
                      false, report_more_information);
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
      return EXIT_FAILURE;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    if (level_code != "w" && level_code != "u" && level_code != "f") {
      cerr << "selected level must be in {w, u, f}" << endl;
      return EXIT_FAILURE;
    }
    const string regions_file = leftover_args.front();
    const string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!is_msite_file(sites_file))
      throw runtime_error("dnmtools counts format required: " + sites_file);

    // make a map that specifies their order; otherwise we can't
    // ensure regions are sorted in the same way
    unordered_map<string, uint32_t> chrom_order;
    for (auto &i : get_chroms(sites_file))
      chrom_order.emplace(i, chrom_order.size());

    if (VERBOSE) cerr << "loading regions" << endl;

    if (!fs::is_regular_file(regions_file))
      // otherwise we could not read the file twice
      throw runtime_error("regions file must be regular file");

    // MAGIC: below allow for ==3 or >=6 columns in the bed format
    const auto n_columns = get_bed_columns(regions_file);
    if (n_columns != 3 && n_columns < 6)
      throw runtime_error("format must be 3 or 6+ column bed: " + regions_file);
    if (is_msite_file(regions_file))
      throw runtime_error("seems to be a counts file: " + regions_file + "\n" +
                          "check order of input arguments");

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!is_sorted(begin(regions), end(regions))) {
      if (sort_data_if_needed) {
        if (VERBOSE) cerr << "sorting regions" << endl;
        sort(begin(regions), end(regions), cmp_chrom_order(chrom_order));
      }
      else
        throw runtime_error("not sorted: " + regions_file + " (consider -s)");
    }

    // if columns are 3 we name the regions otherwise that column will
    // be empty and the output format will be wrong
    if (n_columns == 3)
      for (auto i = 0u; i < regions.size(); ++i)
        regions[i].set_name(default_name_prefix + to_string(i));

    if (VERBOSE) cerr << "[n_regions=" << regions.size() << "]" << endl;

    std::ofstream of;
    if (!outfile.empty()) {
      of.open(outfile);
      if (!of) throw runtime_error("failed to open outfile: " + outfile);
    }
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    if (preload)
      process_preloaded(VERBOSE, report_more_information, level_code[0],
                        sites_file, chrom_order, regions, out);
    else
      process_on_disk(report_more_information, level_code[0], sites_file,
                      chrom_order, regions, out);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
