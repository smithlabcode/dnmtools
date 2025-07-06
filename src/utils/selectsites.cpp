/* selectsites: program to select sites, specified in a methcounts
 * format file, that are contained in given (bed format) intervals
 *
 * Copyright (C) 2019-2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include <algorithm>
#include <bamxx.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "GenomicRegion.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::ios_base;
using std::runtime_error;
using std::string;
using std::tie;
using std::tuple;
using std::unordered_map;

using bamxx::bgzf_file;

namespace fs = std::filesystem;

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
  c_str() {
    /* between c_str and insertion make sure to clear() */
    *pptr() = '\0';
    return pbase();
  }
};

struct selectsites_summary {

  // command_line is the command used to produce this summary file and
  // the corresponding results
  std::string command_line{};

  // n_target_regions is the number of target regions provided as
  // input to the command
  std::uint64_t n_target_regions{};

  // target_region_size is the sum of the sizes of each target region
  std::uint64_t target_region_size{};

  // n_target_regions_collapsed is the number of target regions after
  // having collapsed the input target regions merging those that
  // overlap
  std::uint64_t n_target_regions_collapsed{};

  // target_region_collapsed_size is the sum of the sizes of target
  // regions after collapsing
  std::uint64_t target_region_collapsed_size{};

  // n_sites_total is the total number of sites available in the input
  // counts file. This value is displayed as zero if the command
  // specified to process the sites on disk.
  std::uint64_t n_sites_total{};

  // n_sites_selected is the total number of sites in the output file
  std::uint64_t n_sites_selected{};

  template <typename T>
  static auto
  measure_target_regions(const T &t) -> tuple<std::uint64_t, std::uint64_t> {
    return {
      std::size(t),
      accumulate(cbegin(t), cend(t), 0ul,
                 [](const std::uint64_t a, const typename T::value_type &v) {
                   return a + v.get_width();
                 }),
    };
  }

  auto
  tostring() const -> std::string {
    std::ostringstream oss;
    oss << "command_line: " << command_line << '\n'
        << "n_target_regions: " << n_target_regions << '\n'
        << "target_region_size: " << target_region_size << '\n'
        << "n_target_regions_collapsed: " << n_target_regions_collapsed << '\n'
        << "target_region_collapsed_size: " << target_region_collapsed_size
        << '\n'
        << "n_sites_total: " << n_sites_total << '\n'
        << "n_sites_selected: " << n_sites_selected;
    return oss.str();
  }
};

static auto
write_stats_output(const selectsites_summary &summary,
                   const std::string &summary_file) -> void {
  if (!summary_file.empty()) {
    std::ofstream out_summary(summary_file);
    if (!out_summary)
      throw runtime_error("bad summary output file");
    out_summary << summary.tostring() << endl;
  }
}

static void
collapsebed(std::vector<GenomicRegion> &regions) {
  std::size_t j = 0;
  for (std::size_t i = 1; i < std::(regions); ++i) {
    if (regions[j].same_chrom(regions[i]) &&
        regions[i].get_start() <= regions[j].get_end()) {
      regions[j].set_end(std::max(regions[j].get_end(), regions[i].get_end()));
    }
    else {
      regions[++j] = regions[i];
    }
  }
  regions.erase(begin(regions) + j + 1, end(regions));
}

static inline bool
precedes(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() < s.chrom ||
          (r.get_chrom() == s.chrom && r.get_end() <= s.pos));
}

static inline bool
contains(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() == s.chrom &&
          (r.get_start() <= s.pos && s.pos < r.get_end()));
}

static auto
process_all_sites(
  const bool VERBOSE, const std::string &sites_file,
  const unordered_map<std::string, std::vector<GenomicRegion>> &regions,
  bgzf_file &out) -> tuple<std::uint64_t, std::uint64_t> {
  bgzf_file in(sites_file, "r");
  if (!in)
    throw runtime_error("cannot open file: " + sites_file);

  std::uint64_t n_sites_total = 0u;
  std::uint64_t n_sites_selected = 0u;

  MSite the_site, prev_site;
  std::vector<GenomicRegion>::const_iterator i, i_lim;
  bool chrom_is_relevant = false;
  while (read_site(in, the_site)) {
    ++n_sites_total;
    if (the_site.chrom != prev_site.chrom) {
      if (VERBOSE)
        cerr << "processing " << the_site.chrom << endl;
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

static auto
get_sites_in_region(std::ifstream &site_in, const GenomicRegion &region,
                    bgzf_file &out) -> std::uint64_t {
  quick_buf buf;

  const std::string chrom{region.get_chrom()};
  const std::size_t start_pos = region.get_start();
  const std::size_t end_pos = region.get_end();
  find_offset_for_msite(chrom, start_pos, site_in);

  MSite the_site;
  std::uint64_t n_sites_selected = 0u;
  // ADS: should only happen once that "the_site.chrom < chrom" and
  // this is only needed because of end state of binary search on disk
  while (site_in >> the_site &&
         (the_site.chrom < chrom ||
          (the_site.chrom == chrom && the_site.pos < end_pos)))
    if (start_pos <= the_site.pos) {
      // ADS: don't like this -- always needing to "clear()" this
      // struct is bad...
      buf.clear();
      buf << the_site << '\n';
      if (!out.write(buf.data(), buf.tellp()))
        throw runtime_error("error writing output");
      ++n_sites_selected;
    }
  return n_sites_selected;
}

static auto
process_with_sites_on_disk(const std::string &sites_file,
                           std::vector<GenomicRegion> &regions,
                           bgzf_file &out) -> std::uint64_t {
  std::ifstream in(sites_file);
  if (!in)
    throw runtime_error("cannot open file: " + sites_file);

  std::uint64_t n_sites_selected = 0ul;
  for (auto i = 0u; i < std::size(regions) && in; ++i)
    n_sites_selected += get_sites_in_region(in, regions[i], out);
  return n_sites_selected;
}

static void
regions_by_chrom(
  std::vector<GenomicRegion> &regions,
  unordered_map<std::string, std::vector<GenomicRegion>> &lookup) {
  for (auto &&r : regions) {
    const std::string chrom_name(r.get_chrom());
    if (lookup.find(chrom_name) == end(lookup))
      lookup[chrom_name] = std::vector<GenomicRegion>();
    lookup[chrom_name].push_back(r);
  }
  regions.clear();
  regions.shrink_to_fit();
}

inline bool
file_exists(const std::string &filename) {
  return (access(filename.data(), F_OK) == 0);
}

static bool
is_compressed_file(const std::string &filename) {
  const bgzf_file f(filename.data(), "r");
  htsFormat fmt;
  const int ret = hts_detect_format(f.f->fp, &fmt);
  if (ret != 0)
    throw runtime_error("failed to detect format: " + filename);
  return fmt.compression != no_compression;
}

static auto
get_command_line(const int argc, const char *const argv[]) -> std::string {
  if (argc == 0)
    return std::string();
  std::ostringstream cmd;
  cmd << '"';
  copy(argv, argv + (argc - 1), std::ostream_iterator<const char *>(cmd, " "));
  cmd << argv[argc - 1] << '"';
  return cmd.str();
}

int
main_selectsites(int argc, char *argv[]) {

  try {

    bool VERBOSE = false;
    bool keep_file_on_disk = false;
    bool compress_output = false;

    std::string outfile("-");
    std::string summary_file;

    const std::string description =
      "Select sites inside a set of genomic intervals. "
      "Sites must be specified in methcounts format. "
      "Intervals must be specified in bed format.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<intervals-bed> <sites-file>", 2);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("disk", 'd',
                      "process sites on disk "
                      "(fast if target intervals are few)",
                      false, keep_file_on_disk);
    opt_parse.add_opt("summary", 'S', "write summary to this file", false,
                      summary_file);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
    std::vector<std::string> leftover_args;
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
      return EXIT_SUCCESS;
    }
    if (std::size(leftover_args) != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const std::string regions_file = leftover_args.front();
    const std::string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    selectsites_summary summary;
    summary.command_line = get_command_line(argc, argv);

    if (!fs::is_regular_file(sites_file))
      throw runtime_error("bad input sites file: " + sites_file);

    if (is_compressed_file(sites_file)) {
      keep_file_on_disk = false;
      if (VERBOSE)
        cerr << "input file is so must be loaded" << endl;
    }

    std::vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in file: " + regions_file);
    std::tie(summary.n_target_regions, summary.target_region_size) =
      selectsites_summary::measure_target_regions(regions);

    const std::size_t n_orig_regions = std::size(regions);
    collapsebed(regions);
    if (VERBOSE && n_orig_regions != std::size(regions))
      cerr << "[number of regions merged due to overlap: "
           << n_orig_regions - std::size(regions) << "]" << endl;

    std::tie(summary.n_target_regions_collapsed,
             summary.target_region_collapsed_size) =
      selectsites_summary::measure_target_regions(regions);

    unordered_map<std::string, std::vector<GenomicRegion>> regions_lookup;
    if (!keep_file_on_disk)
      regions_by_chrom(regions, regions_lookup);

    // open the output file
    const std::string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw runtime_error("error opening output file: " + outfile);

    if (keep_file_on_disk)
      summary.n_sites_selected =
        process_with_sites_on_disk(sites_file, regions, out);
    else
      std::tie(summary.n_sites_selected, summary.n_sites_total) =
        process_all_sites(VERBOSE, sites_file, regions_lookup, out);

    write_stats_output(summary, summary_file);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
