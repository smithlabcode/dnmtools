/* Copyright (C) 2025 Andrew D Smith
 *
 * Author: Andrew D Smith
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "MSite.hpp"

#include "OptionParser.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

constexpr auto description =
  R"(This program computes statistics on the autocorrelation of methylation
levels. The input file must be in "counts" format from dnmtools. The output
file has one line for each distance between sites. Sites are assumed to be CpG
sites, but need not be. Each line of output has the follwing values:

    distance, correlation, N, sdX, sdY, covXY

The value of N is the number of observations contributing.)";

struct genomic_interval {
  std::string chrom;
  std::uint32_t start_pos{};
  std::uint32_t end_pos{};

  auto
  operator<(const genomic_interval &other) const -> bool {
    const auto x = chrom.compare(other.chrom);
    return (x < 0 || (x == 0 && (start_pos < other.start_pos ||
                                 (start_pos == other.start_pos &&
                                  (end_pos < other.end_pos)))));
  }
};

auto
operator>>(std::istream &in, genomic_interval &gi) -> std::istream & {
  if (std::string line; getline(in, line)) {
    std::istringstream iss(line);
    if (!(iss >> gi.chrom >> gi.start_pos >> gi.end_pos))
      in.setstate(std::ios_base::failbit);
  }
  return in;
}

static auto
region_precedes_site(const genomic_interval &region,
                     const MSite &site) -> bool {
  // check if the region precedes a the site; [a, b) doesn't contain x
  const int x = region.chrom.compare(site.chrom);
  return x < 0 || (x == 0 && region.end_pos <= site.pos);
}

static auto
region_contains_site(const genomic_interval &region,
                     const MSite &site) -> bool {
  // check if a given site is contained in a given region location
  // Containment is for half open intervals [a, b)
  return (region.chrom == site.chrom && region.start_pos <= site.pos &&
          site.pos < region.end_pos);
}

// Find the position of the region boundary (inclusive counting of
// bases) to the right of the current site; **assume the region
// contains the site**. Index of the region is updated.
static auto
boundary_position(const std::vector<genomic_interval> &regions,
                  std::uint32_t &idx, const MSite &site) -> std::uint32_t {
  // move index of regions so the region doesn't entirely precede it
  while (idx < std::size(regions) && region_precedes_site(regions[idx], site))
    ++idx;

  // check if the region contains the site
  if (idx < std::size(regions) && region_contains_site(regions[idx], site))
    // the right side (closed interval) of containing region; need to return
    // the site that is within the region, since it will be used to know when
    // to stop, and that must be inclusive
    return regions[idx].end_pos - 1;
  // by default return a very far boundary position
  return std::numeric_limits<std::uint32_t>::max();
}

static auto
strands_are_good(const int strand, const MSite &a, const MSite &b) -> bool {
  return strand == 0 || (strand == 1 && a.strand == b.strand) ||
         (strand == -1 && a.strand != b.strand);
}

struct sum_stats {
  std::uint32_t N{};  // counts
  double X{}, Y{};    // sums
  double XY{};        // cross product
  double XX{}, YY{};  // sums of squares

  auto
  update(const double x, const double y) -> void {
    ++N;
    X += x;
    Y += y;
    XY += x * y;
    XX += x * x;
    YY += y * y;
  }

  auto
  get_values(double &sdX, double &sdY, double &covXY) const -> double {
    // Sum XY - N.mu(X).mu(Y) = Sum XY - Sum(X)Sum(Y)/N
    covXY = XY - (X * Y) / N;
    // sqrt[SSX - N.mu(X).mu(X)]
    sdX = std::sqrt(XX - (X * X) / N);
    // sqrt[SSY - N.mu(Y).mu(Y)]
    sdY = std::sqrt(YY - (Y * Y) / N);
    // Pearson correlation
    const double rXY = covXY / (sdX * sdY);

    const double size_factor = 1.0 / std::sqrt(N - 1.0);
    sdX *= size_factor;
    sdY *= size_factor;

    return rXY;
  }
};

std::ostream &
operator<<(std::ostream &out, const sum_stats &s) {
  double sdX = 0.0, sdY = 0.0, covXY = 0.0;
  const auto rXY = s.get_values(sdX, sdY, covXY);
  out << rXY << '\t' << s.N << '\t' << sdX << '\t' << sdY << '\t' << covXY;
  return out;
}

static auto
site_allowed(const std::vector<genomic_interval> &regions, const MSite &site,
             std::uint32_t &idx) -> bool {
  // check if a site is allowed to be used for correlation calculation
  // depending on whether you want to exlude or include certain regions
  while (idx < size(regions) && region_precedes_site(regions[idx], site))
    ++idx;
  return idx < size(regions) && region_contains_site(regions[idx], site);
}

static auto
process_chrom(const bool require_same_region,
              const std::vector<genomic_interval> &regions,
              std::uint32_t &region_idx,  // region_idx changes along chrom
              const std::uint32_t min_reads, const std::uint32_t max_dist,
              const std::uint32_t window_size, const std::vector<MSite> &sites,
              const int strand, std::vector<sum_stats> &the_stats) -> void {
  // early exit if there are not enough sites
  if (std::size(sites) <= 1)
    return;

  // assign CpGs in the given std::vector to the distance tables a & b;
  // each row corresponds to a gap distance in bp
  const auto j_lim = std::end(sites);
  const auto i_lim = j_lim - 1;
  for (auto i = std::begin(sites); i != i_lim; ++i) {
    // check if current site is covered by enough reads
    if (i->n_reads >= min_reads) {
      // determine the limit of sites to consider
      std::uint32_t position_limit = i->pos + max_dist;
      if (require_same_region)
        position_limit =
          std::min(position_limit, boundary_position(regions, region_idx, *i));

      // start j at the 1st neighbor
      for (auto j = i + 1; j != j_lim && j->pos <= position_limit; ++j) {
        // check other site has enough reads and right orientation
        if (j->n_reads >= min_reads && strands_are_good(strand, *i, *j)) {
          // get range of posns current sites should contribute to
          const std::uint32_t the_dist = j->pos - i->pos;
          std::uint32_t d =
            (the_dist > window_size) ? the_dist - window_size : 1;
          const std::uint32_t d_lim =
            std::min(the_dist + window_size, max_dist);

          // include data for current sites in each relevant bin
          while (d <= d_lim)
            the_stats[d++].update(i->meth, j->meth);
        }
      }
    }
  }
}

static auto
process_sites_all_neighbors(const bool verbose, std::ifstream &in,
                            const bool require_same_region,
                            const std::vector<genomic_interval> &regions,
                            const std::uint32_t min_reads,
                            const std::uint32_t max_dist,
                            const std::uint32_t window_size, const int strand,
                            std::vector<sum_stats> &the_stats) -> bool {
  std::vector<MSite> sites;
  MSite the_site;
  std::string prev_chrom;

  std::uint32_t region_idx_add{};
  std::uint32_t region_idx_process{};

  while (in >> the_site) {
    if (the_site.chrom != prev_chrom) {
      if (!sites.empty())
        process_chrom(require_same_region, regions, region_idx_process,
                      min_reads, max_dist, window_size, sites, strand,
                      the_stats);
      if (verbose)
        std::cerr << "processing: " << the_site.chrom << '\n';
      sites.clear();
    }
    if (regions.empty() || site_allowed(regions, the_site, region_idx_add))
      sites.push_back(the_site);
    prev_chrom = std::move(the_site.chrom);
  }
  if (!in.eof())
    return false;
  if (!sites.empty())
    process_chrom(require_same_region, regions, region_idx_process, min_reads,
                  max_dist, window_size, sites, strand, the_stats);
  return true;
}

static auto
process_chrom_kth(const bool require_same_region,
                  const std::vector<genomic_interval> &regions,
                  std::uint32_t &region_idx,  // region_idx changes along chrom
                  const std::uint32_t min_reads, const std::uint32_t max_dist,
                  const std::uint32_t window_size,
                  const std::vector<MSite> &sites, const int strand,
                  const std::uint32_t the_neighbor,
                  std::vector<sum_stats> &the_stats) -> void {
  // early exit if there are not enough sites for specified neighbor
  if (std::size(sites) <= the_neighbor)
    return;

  // limit of outer iteration to allow for k-th neighbors
  const auto j_lim = std::end(sites);
  const auto i_lim = j_lim - the_neighbor - 1;

  // iterate over "left" site to get all pairs
  for (auto i = std::begin(sites); i != i_lim; ++i) {
    // check if current site is covered by enough reads
    if (i->n_reads >= min_reads) {
      // determine the limit of "right" sites to consider
      std::uint32_t position_limit = i->pos + max_dist;
      if (require_same_region)
        position_limit =
          std::min(position_limit, boundary_position(regions, region_idx, *i));

      const auto j = i + the_neighbor;  // get "right" neighbor site
      // make sure right neighbor is in range
      if (j != j_lim && j->pos <= position_limit) {
        // check if other site has enough reads
        if (j->n_reads >= min_reads && strands_are_good(strand, *i, *j)) {
          // get range of posns current sites should contribute to
          const std::uint32_t the_dist = j->pos - i->pos;
          std::uint32_t d =
            (the_dist > window_size) ? the_dist - window_size : 1;
          const std::uint32_t d_lim =
            std::min(the_dist + window_size, max_dist);
          while (d <= d_lim)
            the_stats[d++].update(i->meth, j->meth);
        }
      }
    }
  }
}

static auto
process_sites_kth_neighbor(const bool verbose, std::ifstream &in,
                           const bool require_same_region,
                           const std::vector<genomic_interval> &regions,
                           const std::uint32_t min_reads,
                           const std::uint32_t max_dist,
                           const std::uint32_t window_size, const int strand,
                           const std::uint32_t the_neighbor,
                           std::vector<sum_stats> &the_stats) -> bool {
  std::vector<MSite> sites;
  MSite the_site;
  std::string prev_chrom;

  std::uint32_t region_idx_add = 0;
  std::uint32_t region_idx_process = 0;

  while (in >> the_site) {
    if (the_site.chrom != prev_chrom) {
      if (!sites.empty())
        process_chrom_kth(require_same_region, regions, region_idx_process,
                          min_reads, max_dist, window_size, sites, strand,
                          the_neighbor, the_stats);
      if (verbose)
        std::cerr << "processing: " << the_site.chrom << '\n';
      sites.clear();
    }
    if (regions.empty() || site_allowed(regions, the_site, region_idx_add))
      sites.push_back(the_site);
    prev_chrom = std::move(the_site.chrom);
  }
  if (!in.eof())
    return false;
  if (!sites.empty())
    process_chrom_kth(require_same_region, regions, region_idx_process,
                      min_reads, max_dist, window_size, sites, strand,
                      the_neighbor, the_stats);
  return true;
}

static auto
load_genomic_intervals(std::ifstream &in) -> std::vector<genomic_interval> {
  std::vector<genomic_interval> regions;
  for (genomic_interval r; in >> r;)
    regions.emplace_back(r);

  return regions;
}

static auto
newline_terminated_file(const std::string &filename) -> bool {
  std::ifstream in(filename);
  if (!in)
    return false;
  in.seekg(-1, std::ios_base::end);  // go to one position before EOF
  auto c = '\0';
  in.get(c);
  return c == '\n';
}

auto
main_autocorr(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  static constexpr auto default_max_dist = 4000;
  static constexpr auto default_min_reads = 10;
  static constexpr auto default_min_sites = 500;

  static constexpr auto header = "distance\tcorrelation\tN\tsdX\tsdY\tcovXY\n";

  const auto megabytes = [](const auto x) {
    constexpr auto mb = 1024 * 1024;
    return std::to_string(static_cast<double>(x) / mb) + "MB";
  };

  std::string input_filename;
  std::string regions_file;
  std::string outfile;

  std::uint32_t max_dist = default_max_dist;
  std::uint32_t min_reads = default_min_reads;
  std::uint32_t min_sites = default_min_sites;
  std::uint32_t the_neighbor = 0;  // 0=nearest
  std::uint32_t window_size = 0;
  int strand = 0;  // code: same=1, different=-1 and any=0

  bool verbose = false;
  bool require_same_region = false;
  bool include_header = false;

  /****************** COMMAND LINE OPTIONS ********************/
  OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                         description, "<methylation-file>");
  opt_parse.set_prog_descr_raw();
  opt_parse.add_opt("input", 'i', "input file name", true, input_filename);
  opt_parse.add_opt("output", 'o', "output file name", true, outfile);
  opt_parse.add_opt("max-dist", 'd', "maximum distance for pairs of sites",
                    false, max_dist);
  opt_parse.add_opt("min-sites", 's', "minimum sites needed for correlation",
                    false, min_sites);
  opt_parse.add_opt("min-coverage", 'c', "minimum coverage needed at a site",
                    false, min_reads);
  opt_parse.add_opt("nth-neighbor", 'n',
                    "use only nth neighbor (0=any, 1=nearest, ...)", false,
                    the_neighbor);
  opt_parse.add_opt("window-size", 'w', "smoothing window for distances", false,
                    window_size);
  opt_parse.add_opt("strand", 't',
                    "use neighbors by strand: same=1, different=-1, any=0",
                    false, strand);
  opt_parse.add_opt("regions", 'R', "file of regions to process (BED format)",
                    false, regions_file);
  opt_parse.add_opt("same-region", 'S', "require both sites in same region",
                    false, require_same_region);
  opt_parse.add_opt("header", 'H', "put header in output (R data frame format)",
                    false, include_header);
  opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
  opt_parse.set_show_defaults();
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
  /****************** END COMMAND LINE OPTIONS *****************/

  const auto input_filesize = std::filesystem::file_size(input_filename);

  if (verbose)
    std::cerr << "min_reads: " << min_reads << '\n'
              << "max_dist: " << max_dist << '\n'
              << "min_sites: " << min_sites << '\n'
              << "the_neighbor: " << the_neighbor << '\n'
              << "input_filename: " << input_filename << '\n'
              << "input_filesize: " << megabytes(input_filesize) << '\n';

  std::vector<genomic_interval> regions;
  if (!regions_file.empty()) {
    if (verbose)
      std::cerr << "loading regions: " << regions_file << '\n';

    std::ifstream regions_in(regions_file);
    if (!regions_in) {
      std::cerr << "bad regions file: " << regions_file << '\n';
      return EXIT_FAILURE;
    }

    regions = load_genomic_intervals(regions_in);
    if (!is_sorted(std::begin(regions), std::end(regions))) {
      std::cerr << "regions not sorted: " << regions_file << '\n';
      return EXIT_FAILURE;
    }
    if (verbose)
      std::cerr << "total regions: " << std::size(regions) << '\n';
  }

  std::ifstream sites_in(input_filename);
  if (!sites_in) {
    std::cout << "bad input file: " << input_filename << '\n';
    return EXIT_FAILURE;
  }

  if (!newline_terminated_file(input_filename)) {
    std::cout << "truncated input file: " << input_filename << '\n';
    return EXIT_FAILURE;
  }

  std::vector<sum_stats> the_stats(max_dist + 1);
  const auto success = [&] {
    if (the_neighbor == 0)
      return process_sites_all_neighbors(verbose, sites_in, require_same_region,
                                         regions, min_reads, max_dist,
                                         window_size, strand, the_stats);
    return process_sites_kth_neighbor(verbose, sites_in, require_same_region,
                                      regions, min_reads, max_dist, window_size,
                                      strand, the_neighbor, the_stats);
  }();

  if (!success) {
    std::cerr << "failure parsing file: " << input_filename << '\n';
    return EXIT_FAILURE;
  }

  std::ofstream out(outfile);
  if (!out) {
    std::cerr << "failed to open output: " << outfile << '\n';
    return EXIT_FAILURE;
  }

  if (include_header)
    out << header;

  auto i = 0u;
  for (const auto &s : the_stats) {
    if (s.N >= min_sites)
      out << i << '\t' << s << '\n';
    ++i;
  }

  return EXIT_SUCCESS;
}
