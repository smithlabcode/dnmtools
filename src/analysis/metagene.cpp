/* Copyright (C) 2010-2025 Andrew D. Smith
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
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

[[maybe_unused]] static constexpr auto about = R"(
metagene (tsscpgplot): get data to plot methylation level around a TSS
)";

#include "Interval6.hpp"
#include "LevelsCounter.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

static MSite
tss_from_gene(const Interval6 &r) {
  MSite s;
  s.chrom = r.chrom;
  s.pos = r.strand == '+' ? r.start : r.stop;
  s.strand = r.strand;
  return s;
}

static void
process_chrom(const std::uint32_t region_size,
              const std::vector<Interval6> &genes,
              const std::pair<std::uint32_t, std::uint32_t> &bounds,
              const std::vector<MSite> &sites,
              std::vector<LevelsCounter> &levels) {
  const auto cmp = [](const MSite &a, const MSite &b) { return a.pos < b.pos; };
  const std::uint32_t twice_rs = 2 * region_size;

  for (auto i = bounds.first; i < bounds.second; ++i) {
    const MSite tss = tss_from_gene(genes[i]);
    MSite left = tss;
    left.pos = left.pos > region_size ? left.pos - region_size : 0;

    MSite right = tss;
    right.pos += region_size;

    const auto lim =
      std::upper_bound(std::cbegin(sites), std::cend(sites), right, cmp);
    auto itr =
      std::lower_bound(std::cbegin(sites), std::cend(sites), left, cmp);
    for (; itr != lim; ++itr) {
      const auto dist = itr->pos - left.pos;
      const auto base_idx = tss.strand == '+' ? dist : twice_rs - dist;
      levels[base_idx].update(*itr);
    }
  }
}

template <typename T>
static void
collapse_bins(const std::uint32_t bin_size, std::vector<T> &v) {
  const std::uint32_t n_bins =
    std::ceil(static_cast<double>(std::size(v)) / bin_size);
  std::vector<T> vv(n_bins);
  for (auto i = 0u; i < std::size(v); ++i)
    vv[i / bin_size] += v[i];
  v.swap(vv);
}

int
metagene(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  constexpr auto description =
    R"(
Compute the information needed for metagene plots of DNA methylation
levels. The columns in the output correspond to the fields calculated
globally by the `levels` and per-region by the `roi` command. Input
for features is in BED format, and when present the 6th column is used
to indicate strand. For features of non-zero width (where the 2nd and
3rd columns are not identical) the negative strand will indicate that
3rd column should be used. This means, for example, if the features are
genes, and the promoters are of interest, the strand will be used
correctly.
)";
  try {
    std::string outfile;
    std::uint32_t region_size = 5000;  // NOLINT(*-avoid-magic-numbers)
    bool verbose = false;
    bool show_progress = false;
    std::uint32_t bin_size = 50;  // NOLINT(*-avoid-magic-numbers)

    /****************** GET COMMAND LINE ARGUMENTS ***************************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<features-bed> <counts>");
    opt_parse.add_opt("output", 'o', "output file (default: terminal)", false,
                      outfile);
    opt_parse.add_opt("size", 's', "analyze this size in both directions",
                      false, region_size);
    opt_parse.add_opt("bin", 'b', "bin size", false, bin_size);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
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
    if (std::size(leftover_args) != 2) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string features_file_name = leftover_args.front();
    const std::string cpg_file_name = leftover_args.back();
    /**********************************************************************/

    if (verbose)
      std::cerr << "[loading feature annotations data]\n";
    auto features = read_intervals6(features_file_name);
    std::sort(std::begin(features), std::end(features));
    if (verbose)
      std::cerr << "[number of features: " << std::size(features) << "]\n";

    // identify the start and end of ranges for each chromosome
    using pos_pair = std::pair<std::uint32_t, std::uint32_t>;
    std::unordered_map<std::string, pos_pair> lookup;

    std::string chrom_name;
    auto prev_idx = 0u;
    for (auto i = 0u; i < std::size(features); ++i)
      if (features[i].chrom != chrom_name) {
        if (!chrom_name.empty())
          lookup.insert({chrom_name, {prev_idx, i}});
        prev_idx = i;
        chrom_name = features[i].chrom;
      }
    lookup.insert({chrom_name, {prev_idx, std::size(features)}});
    if (verbose)
      std::cerr << "[number of chroms with features: " << std::size(lookup)
                << "]\n";

    const auto pair_diff = [&lookup](const auto x) {
      return (x != std::cend(lookup)) ? x->second.second - x->second.first : 0u;
    };

    std::vector<LevelsCounter> levels(2 * region_size);

    bamxx::bgzf_file cpgin(cpg_file_name, "r");
    if (!cpgin)
      throw std::runtime_error("failed to open file: " + cpg_file_name);

    std::vector<MSite> sites;
    std::string line;
    chrom_name.clear();
    while (getline(cpgin, line)) {
      const auto the_site = MSite(line);
      if (the_site.chrom != chrom_name) {
        if (!sites.empty()) {
          const auto bounds = lookup.find(chrom_name);
          if (bounds != std::cend(lookup))
            process_chrom(region_size, features, bounds->second, sites, levels);
          const auto n_features = pair_diff(bounds);
          if (show_progress)
            std::cerr << "[sites=" << std::size(sites)
                      << " features=" << n_features << "]" << '\n';
          sites.clear();
        }
        if (show_progress)
          std::cerr << "[processing: " << the_site.chrom << "]";
        chrom_name = the_site.chrom;
      }
      sites.push_back(the_site);
    }

    if (!sites.empty()) {
      const auto bounds = lookup.find(chrom_name);
      if (bounds != std::cend(lookup))
        process_chrom(region_size, features, bounds->second, sites, levels);
      const auto n_features = pair_diff(bounds);
      if (show_progress)
        std::cerr << "[sites=" << std::size(sites) << " features=" << n_features
                  << "]\n";
    }

    collapse_bins(bin_size, levels);

    if (verbose)
      std::cerr << "output columns:\n"
                << LevelsCounter::format_header() << '\n';

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open file: " + outfile);
    for (auto i = 0u; i < std::size(levels); ++i)
      out << i * bin_size << '\t' << levels[i].format_row() << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
