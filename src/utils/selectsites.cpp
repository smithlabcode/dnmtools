/* Copyright (C) 2019-2025 Andrew D. Smith
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

[[maybe_unused]] static constexpr auto about = R"(
selectsites: Select sites, specified in a methcounts format file, that are
contained in given (bed format) intervals.
)";

#include "Interval6.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <htslib/sam.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

struct selectsites_summary {
  // command_line is the command used to produce this summary file and the
  // corresponding results
  std::string command_line{};

  // n_target_regions is the number of target regions provided as
  // input to the command
  std::uint64_t n_target_regions{};

  // target_region_size is the sum of the sizes of each target region
  std::uint64_t target_region_size{};

  // n_target_regions_collapsed is the number of target regions after having
  // collapsed the input target regions merging those that overlap
  std::uint64_t n_target_regions_collapsed{};

  // target_region_collapsed_size is the sum of the sizes of target regions
  // after collapsing
  std::uint64_t target_region_collapsed_size{};

  // n_sites_total is the total number of sites available in the input counts
  // file. This value is displayed as zero if the command specified to process
  // the sites on disk.
  std::uint64_t n_sites_total{};

  // n_sites_selected is the total number of sites in the output file
  std::uint64_t n_sites_selected{};

  template <typename T>
  [[nodiscard]] static auto
  measure_target_regions(const T &t)
    -> std::tuple<std::uint64_t, std::uint64_t> {
    return {
      std::size(t),
      std::accumulate(
        std::cbegin(t), std::cend(t), 0ul,
        [](const std::uint64_t a, const typename T::value_type &v) {
          return a + size(v);
        }),
    };
  }
};

[[nodiscard]] auto
to_string(const selectsites_summary &s) -> std::string {
  std::ostringstream oss;
  oss << "command_line: " << s.command_line << '\n'
      << "n_target_regions: " << s.n_target_regions << '\n'
      << "target_region_size: " << s.target_region_size << '\n'
      << "n_target_regions_collapsed: " << s.n_target_regions_collapsed << '\n'
      << "target_region_collapsed_size: " << s.target_region_collapsed_size
      << '\n'
      << "n_sites_total: " << s.n_sites_total << '\n'
      << "n_sites_selected: " << s.n_sites_selected;
  return oss.str();
}

static auto
write_stats_output(const selectsites_summary &summary,
                   const std::string &summary_file) -> void {
  if (!summary_file.empty()) {
    std::ofstream out_summary(summary_file);
    if (!out_summary)
      throw std::runtime_error("bad summary output file");
    out_summary << to_string(summary) << '\n';
  }
}

static void
collapsebed(std::vector<Interval6> &regions) {
  std::ptrdiff_t j = 0;
  for (const auto &r : regions) {
    if (regions[j].chrom == r.chrom && r.start <= regions[j].stop)
      regions[j].stop = std::max(regions[j].stop, r.stop);
    else
      regions[++j] = r;
  }
  regions.erase(std::begin(regions) + j + 1, std::end(regions));
}

static inline bool
precedes(const Interval6 &r, const MSite &s) {
  return (r.chrom < s.chrom || (r.chrom == s.chrom && r.stop <= s.pos));
}

static inline bool
contains(const Interval6 &r, const MSite &s) {
  return (r.chrom == s.chrom && (r.start <= s.pos && s.pos < r.stop));
}

static auto
process_all_sites(
  const bool verbose, const std::string &sites_file,
  const std::unordered_map<std::string, std::vector<Interval6>> &regions,
  bamxx::bgzf_file &out) -> std::tuple<std::uint64_t, std::uint64_t> {
  bamxx::bgzf_file in(sites_file, "r");
  if (!in)
    throw std::runtime_error("cannot open file: " + sites_file);

  std::uint64_t n_sites_total = 0u;
  std::uint64_t n_sites_selected = 0u;

  MSite the_site, prev_site;
  std::vector<Interval6>::const_iterator i, i_lim;
  bool chrom_is_relevant = false;
  while (read_site(in, the_site)) {
    ++n_sites_total;
    if (the_site.chrom != prev_site.chrom) {
      if (verbose)
        std::cerr << "processing " << the_site.chrom << '\n';
      const auto r = regions.find(the_site.chrom);
      chrom_is_relevant = (r != cend(regions));
      if (chrom_is_relevant) {
        i = cbegin(r->second);
        i_lim = cend(r->second);
      }
    }
    if (chrom_is_relevant) {
      while (i != i_lim && precedes(*i, the_site))
        ++i;
      if (i != i_lim && contains(*i, the_site)) {
        ++n_sites_selected;
        write_site(out, the_site);
      }
    }
    std::swap(prev_site, the_site);
  }
  return {n_sites_selected, n_sites_total};
}

[[nodiscard]] static auto
get_sites_in_region(std::ifstream &site_in, const Interval6 &region,
                    bamxx::bgzf_file &out) -> std::uint64_t {
  const std::string chrom{region.chrom};
  const std::size_t start_pos = region.start;
  const std::size_t end_pos = region.stop;
  find_offset_for_msite(chrom, start_pos, site_in);

  std::uint64_t n_sites_selected = 0u;
  // ADS: should only happen once that "the_site.chrom < chrom" and
  // this is only needed because of end state of binary search on disk
  std::string line;
  while (getline(site_in, line)) {
    const auto the_site = MSite(line);
    if (chrom <= the_site.chrom &&
        (the_site.chrom != chrom || end_pos <= the_site.pos))
      break;
    if (start_pos <= the_site.pos) {
      if (!out.write(line.data(), std::size(line)))
        throw std::runtime_error("error writing output");
      ++n_sites_selected;
    }
  }
  return n_sites_selected;
}

static auto
process_with_sites_on_disk(const std::string &sites_file,
                           const std::vector<Interval6> &regions,
                           bamxx::bgzf_file &out) -> std::uint64_t {
  std::ifstream in(sites_file);
  if (!in)
    throw std::runtime_error("cannot open file: " + sites_file);
  std::uint64_t n_sites_selected = 0ul;
  for (auto i = 0u; i < std::size(regions) && in; ++i)
    n_sites_selected += get_sites_in_region(in, regions[i], out);
  return n_sites_selected;
}

static void
regions_by_chrom(
  std::vector<Interval6> &regions,
  std::unordered_map<std::string, std::vector<Interval6>> &lookup) {
  for (auto &&r : regions)
    lookup[r.chrom].push_back(r);
  regions.clear();
  regions.shrink_to_fit();
}

static bool
is_compressed_file(const std::string &filename) {
  const bamxx::bgzf_file f(filename.data(), "r");
  htsFormat fmt;
  const int ret = hts_detect_format(f.f->fp, &fmt);
  if (ret != 0)
    throw std::runtime_error("failed to detect format: " + filename);
  return fmt.compression != no_compression;
}

static auto
get_command_line(const int argc,
                 const char *const argv[]  // NOLINT(*-avoid-c-arrays)
                 ) -> std::string {
  if (argc == 0)
    return std::string{};
  std::ostringstream cmd;
  cmd << '"';
  std::copy(argv, argv + (argc - 1),
            std::ostream_iterator<const char *>(cmd, " "));
  cmd << argv[argc - 1] << '"';  // NOLINT(*-pointer-arithmetic)
  return cmd.str();
}

int
main_selectsites(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    bool verbose = false;
    bool keep_file_on_disk = false;
    bool compress_output = false;
    bool allow_extra_fields = false;

    std::string outfile("-");
    std::string summary_file;

    const std::string description =
      "Select sites inside a set of genomic intervals. "
      "Sites must be specified in methcounts format. "
      "Intervals must be specified in bed format.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<intervals-bed> <sites-file>", 2);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("disk", 'd',
                      "process sites on disk "
                      "(fast if target intervals are few)",
                      false, keep_file_on_disk);
    opt_parse.add_opt("summary", 'S', "write summary to this file", false,
                      summary_file);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("relaxed", '\0', "input has extra fields", false,
                      allow_extra_fields);
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
    if (std::size(leftover_args) != 2) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string regions_file = leftover_args.front();
    const std::string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    MSite::no_extra_fields = (allow_extra_fields == false);

    selectsites_summary summary;
    summary.command_line = get_command_line(argc, argv);

    if (!std::filesystem::is_regular_file(sites_file))
      throw std::runtime_error("bad input sites file: " + sites_file);

    if (is_compressed_file(sites_file)) {
      keep_file_on_disk = false;
      if (verbose)
        std::cerr << "input file is so must be loaded\n";
    }

    auto regions = read_intervals6(regions_file);
    if (!std::is_sorted(std::cbegin(regions), std::cend(regions)))
      throw std::runtime_error("regions not sorted in file: " + regions_file);
    std::tie(summary.n_target_regions, summary.target_region_size) =
      selectsites_summary::measure_target_regions(regions);

    const std::size_t n_orig_regions = std::size(regions);
    collapsebed(regions);
    if (verbose && n_orig_regions != std::size(regions))
      std::cerr << "[number of regions merged due to overlap: "
                << n_orig_regions - std::size(regions) << "]\n";

    std::tie(summary.n_target_regions_collapsed,
             summary.target_region_collapsed_size) =
      selectsites_summary::measure_target_regions(regions);

    std::unordered_map<std::string, std::vector<Interval6>> regions_lookup;
    if (!keep_file_on_disk)
      regions_by_chrom(regions, regions_lookup);

    // open the output file
    const std::string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw std::runtime_error("error opening output file: " + outfile);

    if (keep_file_on_disk)
      summary.n_sites_selected =
        process_with_sites_on_disk(sites_file, regions, out);
    else
      std::tie(summary.n_sites_selected, summary.n_sites_total) =
        process_all_sites(verbose, sites_file, regions_lookup, out);

    write_stats_output(summary, summary_file);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
