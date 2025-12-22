/* Copyright (C) 2014-2023 Andrew D. Smith, Benjamin E Decato and Fang Fang
 *
 * Authors: Andrew D. Smith and Benjamin E Decato and Fang Fang
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
amrtester: A program for testing whether a genomic region has allele-specific
methylation
)";

#include "Epiread.hpp"
#include "EpireadStats.hpp"

#include <Interval6.hpp>
#include <OptionParser.hpp>
#include <smithlab_os.hpp>
#include <smithlab_utils.hpp>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

static void
backup_to_start_of_current_record(std::ifstream &in) {
  static constexpr std::size_t assumed_max_valid_line_width = 10000;
  std::size_t count = 0;
  while (in.tellg() > 0 && in.peek() != '\n' && in.peek() != '\r' &&
         count++ < assumed_max_valid_line_width)
    in.seekg(-1, std::ios_base::cur);
  if (count > assumed_max_valid_line_width)
    throw std::runtime_error("file contains a line longer than " +
                             std::to_string(assumed_max_valid_line_width));
}

static std::streampos
find_first_epiread_ending_after_position(const std::string &query_chrom,
                                         const std::size_t query_pos,
                                         std::ifstream &in) {
  in.seekg(0, std::ios_base::end);
  auto high_pos = in.tellg();
  std::size_t eof = in.tellg();
  in.seekg(0, std::ios_base::beg);
  std::streamoff low_pos = 0;

  std::string chrom, seq;
  std::size_t start = 0ul;

  // This is just binary search on disk
  while (high_pos > low_pos + 1) {
    const std::streamoff mid_pos = (low_pos + high_pos) / 2;

    in.seekg(mid_pos);
    backup_to_start_of_current_record(in);

    // we've hit the end of file without finding an epiread
    if (low_pos + 2 == static_cast<std::streamoff>(eof))
      return -1;

    if (!(in >> chrom >> start >> seq))
      throw std::runtime_error("problem loading reads");

    if (chrom < query_chrom ||
        (chrom == query_chrom && start + std::size(seq) <= query_pos))
      low_pos = mid_pos;
    else
      high_pos = mid_pos;
  }
  return low_pos;
}

static void
load_reads(const std::string &reads_file_name, const Interval6 &region,
           std::vector<small_epiread> &the_reads) {
  // open and check the file
  std::ifstream in(reads_file_name.c_str());
  if (!in)
    throw std::runtime_error("cannot open input file " + reads_file_name);

  const std::string query_chrom(region.chrom);
  const std::size_t query_start = region.start;
  const std::size_t query_end = region.stop;
  const std::streampos low_offset =
    find_first_epiread_ending_after_position(query_chrom, query_start, in);

  in.seekg(low_offset, std::ios_base::beg);
  backup_to_start_of_current_record(in);

  std::string chrom, seq;
  std::size_t start = 0ul;
  while ((in >> chrom >> start >> seq) && chrom == query_chrom &&
         start < query_end)
    the_reads.emplace_back(start, seq);
}

static void
convert_coordinates(const std::vector<std::size_t> &cpg_positions,
                    Interval6 &region) {
  const auto lb_pos = [](const auto &v, const auto x) {
    return std::distance(std::cbegin(v),
                         std::lower_bound(std::cbegin(v), std::cend(v), x));
  };
  region.start = lb_pos(cpg_positions, region.start);
  region.stop = lb_pos(cpg_positions, region.stop);
}

inline static bool
is_cpg(const std::string &s, const std::size_t idx) {
  return std::toupper(s[idx]) == 'C' && std::toupper(s[idx + 1]) == 'G';
}

static void
collect_cpgs(const std::string &s, std::vector<std::size_t> &cpgs) {
  const std::size_t lim = s.length() - 1;
  for (std::size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}

static void
clip_reads(const std::size_t start_pos, const std::size_t end_pos,
           std::vector<small_epiread> &r) {
  std::size_t j = 0;
  for (std::size_t i = 0; i < std::size(r); ++i) {
    if (start_pos < r[i].pos + r[i].seq.length() && r[i].pos < end_pos) {
      if (r[i].pos < start_pos) {
        assert(start_pos - r[i].pos < r[i].seq.length());
        r[i].seq = r[i].seq.substr(start_pos - r[i].pos);
        r[i].pos = start_pos;
      }
      if (r[i].end() > end_pos)
        r[i].seq = r[i].seq.substr(0, end_pos - r[i].pos);
      r[j] = r[i];
      ++j;
    }
  }
  r.erase(std::begin(r) + j, std::end(r));
}

// give names to regions if they do not exist
static void
ensure_regions_are_named(std::vector<Interval6> &regions) {
  auto region_name_idx = 0u;
  for (auto &region : regions)
    if (region.name.empty())
      region.name = "region" + std::to_string(++region_name_idx);
}

int
main_amrtester(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr double critical_value = 0.01;

    bool verbose = false;
    bool show_progress = false;
    bool use_bic = false;
    bool correct_for_read_count = true;

    std::string outfile;
    std::string chrom_file;

    // NOLINTBEGIN(*-avoid-magic-numbers)
    std::size_t max_itr{10};
    double high_prob{0.75};
    double low_prob{0.25};
    // NOLINTEND(*-avoid-magic-numbers)

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "resolve epi-alleles",
                           "<bed-regions> <mapped-reads>");
    opt_parse.add_opt("output", 'o', "output file", false, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome fasta file", true,
                      chrom_file);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("nordc", 'r', "turn off read count correction", false,
                      correct_for_read_count);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("progress", 'P', "show progress", false, show_progress);

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
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string regions_file(leftover_args.front());
    const std::string reads_file_name(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    const EpireadStats epistat{low_prob, high_prob, critical_value,
                               max_itr,  use_bic,   correct_for_read_count};

    if (!validate_epiread_file(reads_file_name))
      throw std::runtime_error("invalid states file: " + reads_file_name);

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    std::vector<std::string> chroms;
    std::vector<std::string> chrom_names;
    std::unordered_map<std::string, std::size_t> chrom_lookup;
    read_fasta_file_short_names(chrom_file, chrom_names, chroms);
    for (auto i = 0u; i < std::size(chroms); ++i) {
      if (chrom_lookup.find(chrom_names[i]) != std::cend(chrom_lookup))
        throw std::runtime_error("repeated chromosome name: " + chrom_names[i]);
      chrom_lookup[chrom_names[i]] = i;
    }

    auto regions = read_intervals6(regions_file);
    if (!std::is_sorted(std::cbegin(regions), std::cend(regions)))
      throw std::runtime_error("regions not sorted in: " + regions_file);

    ensure_regions_are_named(regions);

    auto n_regions = std::size(regions);
    if (verbose)
      std::cerr << "number of regions: " << n_regions << '\n';

    std::string chrom_name;
    std::vector<std::size_t> cpg_positions;

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open file: " + outfile);

    bool is_significant = false;

    ProgressBar progress(std::size(regions));
    auto progress_idx = 0u;

    for (auto &region : regions) {
      if (show_progress && progress.time_to_report(progress_idx))
        progress.report(std::cerr, progress_idx);
      ++progress_idx;

      // get the correct chrom if it has changed
      if (region.chrom != chrom_name) {
        chrom_name = region.chrom;
        auto the_chrom = chrom_lookup.find(chrom_name);
        if (the_chrom == end(chrom_lookup))
          throw std::runtime_error("could not find chrom: " + chrom_name);

        cpg_positions.clear();
        collect_cpgs(chroms[the_chrom->second], cpg_positions);
      }

      Interval6 conv_region(region);
      convert_coordinates(cpg_positions, conv_region);

      std::vector<small_epiread> reads;
      load_reads(reads_file_name, conv_region, reads);

      clip_reads(conv_region.start, conv_region.stop, reads);

      const auto score =
        reads.empty() ? 1.0 : epistat.test_asm(reads, is_significant);
      region.score = score;
      region.name += ":" + std::to_string(std::size(reads));
      out << region << '\n';
    }
    if (show_progress)
      std::cerr << "\r100%\n";
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
