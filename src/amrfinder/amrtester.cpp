/* amrtester: A program for testing whether a genomic region has
 * allele-specific methylation
 *
 * Copyright (C) 2014-2023 University of Southern California and
 *                         Benjamin E Decato and Andrew D. Smith and Fang Fang
 *
 * Authors: Andrew D. Smith and Benjamin E Decato and Fang Fang
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
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "Epiread.hpp"
#include "EpireadStats.hpp"

#include <GenomicRegion.hpp>
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

using std::begin;
using std::cerr;
using std::cout;
using std::end;
using std::runtime_error;
using std::streampos;
using std::string;
using std::unordered_map;
using std::vector;

using epi_r = small_epiread;

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

static void
backup_to_start_of_current_record(std::ifstream &in) {
  static const size_t assumed_max_valid_line_width = 10000;
  size_t count = 0;
  while (in.tellg() > 0 && in.peek() != '\n' && in.peek() != '\r' &&
         count++ < assumed_max_valid_line_width)
    in.seekg(-1, std::ios_base::cur);
  if (count > assumed_max_valid_line_width)
    throw runtime_error("file contains a line longer than " +
                        std::to_string(assumed_max_valid_line_width));
}

static streampos
find_first_epiread_ending_after_position(const string &query_chrom,
                                         const size_t query_pos,
                                         std::ifstream &in) {
  in.seekg(0, std::ios_base::end);
  auto high_pos = in.tellg();
  size_t eof = in.tellg();
  in.seekg(0, std::ios_base::beg);
  std::streamoff low_pos = 0;

  string chrom, seq;
  size_t start = 0ul;

  // This is just binary search on disk
  while (high_pos > low_pos + 1) {
    const std::streamoff mid_pos = (low_pos + high_pos) / 2;

    in.seekg(mid_pos);
    backup_to_start_of_current_record(in);

    // we've hit the end of file without finding an epiread
    if (low_pos + 2 == static_cast<std::streamoff>(eof))
      return -1;

    if (!(in >> chrom >> start >> seq)) {
      throw runtime_error("problem loading reads");
    }
    if (chrom < query_chrom ||
        (chrom == query_chrom && start + seq.length() <= query_pos))
      low_pos = mid_pos;
    else
      high_pos = mid_pos;
  }
  return low_pos;
}

static void
load_reads(const string &reads_file_name, const GenomicRegion &region,
           vector<epi_r> &the_reads) {
  // open and check the file
  std::ifstream in(reads_file_name.c_str());
  if (!in)
    throw runtime_error("cannot open input file " + reads_file_name);

  const string query_chrom(region.get_chrom());
  const size_t query_start = region.get_start();
  const size_t query_end = region.get_end();
  const streampos low_offset =
    find_first_epiread_ending_after_position(query_chrom, query_start, in);

  in.seekg(low_offset, std::ios_base::beg);
  backup_to_start_of_current_record(in);

  string chrom, seq;
  size_t start = 0ul;
  while ((in >> chrom >> start >> seq) && chrom == query_chrom &&
         start < query_end)
    the_reads.emplace_back(start, seq);
}

static void
convert_coordinates(const vector<size_t> &cpg_positions,
                    GenomicRegion &region) {
  const size_t start_pos =
    lower_bound(cbegin(cpg_positions), cend(cpg_positions),
                region.get_start()) -
    cbegin(cpg_positions);

  const size_t end_pos =
    lower_bound(cbegin(cpg_positions), cend(cpg_positions), region.get_end()) -
    cbegin(cpg_positions);

  region.set_start(start_pos);
  region.set_end(end_pos);
}

inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}

static void
collect_cpgs(const string &s, vector<size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}

static void
clip_reads(const size_t start_pos, const size_t end_pos, vector<epi_r> &r) {
  size_t j = 0;
  for (size_t i = 0; i < r.size(); ++i) {
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
  r.erase(begin(r) + j, end(r));
}

// give names to regions if they do not exist
static void
ensure_regions_are_named(
  vector<GenomicRegion> &regions  // cppcheck-suppress constParameterReference
) {
  auto region_name_idx = 0u;
  for (auto region : regions)
    if (region.get_name().empty())
      region.set_name("region" + std::to_string(++region_name_idx));
}

int
main_amrtester(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr double critical_value = 0.01;

    bool verbose = false;
    bool show_progress = false;
    bool use_bic = false;
    bool correct_for_read_count = true;

    string outfile;
    string chrom_file;

    size_t max_itr = 10;
    double high_prob = 0.75, low_prob = 0.25;

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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << '\n'
           << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    const string regions_file(leftover_args.front());
    const string reads_file_name(leftover_args.back());
    /****************** END COMMAND LINE OPTIONS *****************/

    const EpireadStats epistat{low_prob, high_prob, critical_value,
                               max_itr,  use_bic,   correct_for_read_count};

    if (!validate_epiread_file(reads_file_name))
      throw runtime_error("invalid states file: " + reads_file_name);

    /* first load in all the chromosome sequences and names, and make
       a map from chromosome name to the location of the chromosome
       itself */
    vector<string> chroms;
    vector<string> chrom_names;
    unordered_map<string, size_t> chrom_lookup;
    read_fasta_file_short_names(chrom_file, chrom_names, chroms);
    for (auto i = 0u; i < size(chroms); ++i) {
      if (chrom_lookup.find(chrom_names[i]) != cend(chrom_lookup))
        throw runtime_error("repeated chromosome name: " + chrom_names[i]);
      chrom_lookup[chrom_names[i]] = i;
    }

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in: " + regions_file);

    ensure_regions_are_named(regions);

    auto n_regions = size(regions);
    if (verbose)
      cerr << "number of regions: " << n_regions << '\n';

    string chrom_name;
    vector<size_t> cpg_positions;

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    bool is_significant = false;

    ProgressBar progress(size(regions));
    auto progress_idx = 0u;

    for (auto &region : regions) {
      if (show_progress && progress.time_to_report(progress_idx))
        progress.report(cerr, progress_idx);
      ++progress_idx;

      // get the correct chrom if it has changed
      if (region.get_chrom() != chrom_name) {
        chrom_name = region.get_chrom();
        auto the_chrom = chrom_lookup.find(chrom_name);
        if (the_chrom == end(chrom_lookup))
          throw runtime_error("could not find chrom: " + chrom_name);

        cpg_positions.clear();
        collect_cpgs(chroms[the_chrom->second], cpg_positions);
      }

      GenomicRegion conv_region(region);
      convert_coordinates(cpg_positions, conv_region);

      vector<epi_r> reads;
      load_reads(reads_file_name, conv_region, reads);

      clip_reads(conv_region.get_start(), conv_region.get_end(), reads);

      const auto score =
        reads.empty() ? 1.0 : epistat.test_asm(reads, is_significant);
      region.set_score(static_cast<float>(score));
      region.set_name(region.get_name() + ":" + toa(reads.size()));
      out << region << '\n';
    }
    if (show_progress)
      cerr << "\r100%\n";
  }
  catch (const std::exception &e) {
    cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-avoid-magic-numbers,*-narrowing-conversions)
