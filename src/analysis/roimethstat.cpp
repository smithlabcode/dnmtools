/* roimethstat: average methylation in each of a set of regions
 *
 * Copyright (C) 2014-2022 Andrew D. Smith
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <stdexcept>
#include <regex>

#include <bamxx.hpp>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"

#include "MSite.hpp"

#include "bsutils.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::make_pair;
using std::ios_base;
using std::runtime_error;
using std::ifstream;
using std::isfinite;
using std::is_sorted;
using std::regex_match;

using bamxx::bgzf_file;

static pair<bool, bool>
meth_unmeth_calls(const size_t n_meth, const size_t n_unmeth) {
  static const double alpha = 0.95;
  // get info for binomial test
  double lower = 0.0, upper = 0.0;
  const size_t total = n_meth + n_unmeth;
  wilson_ci_for_binomial(alpha, total,
                         static_cast<double>(n_meth)/total, lower, upper);
  return make_pair(lower > 0.5, upper < 0.5);
}


static std::pair<size_t, size_t>
region_bounds(const vector<MSite> &sites, const GenomicRegion &region) {
  const string chrom(region.get_chrom());
  const char strand(region.get_strand());

  const MSite a(chrom, region.get_start(), strand, "", 0, 0);
  vector<MSite>::const_iterator a_ins(lower_bound(begin(sites), end(sites), a));

  const MSite b(chrom, region.get_end(), strand, "", 0, 0);
  vector<MSite>::const_iterator b_ins(lower_bound(begin(sites), end(sites), b));

  return make_pair(a_ins - begin(sites), b_ins - begin(sites));
}


/*
   This function is used to make sure the output is done the same way
   in the different places it might be generated. One issue is that
   there is a lot of string and stream manipulation here, which might
   not be desirable if there are lots of regions, as when the genome
   is binned.
*/
static string
format_output_line(const bool report_more_information,
                   const GenomicRegion &r, const double score,
                   const double weighted_mean_meth, const double mean_meth,
                   const double fractional_meth,
                   const size_t total_cpgs, const size_t cpgs_with_reads,
                   const size_t total_meth, const size_t total_reads) {

  std::ostringstream oss;
  oss << r.get_chrom() << '\t'
      << r.get_start() << '\t'
      << r.get_end() << '\t'
      << r.get_name() << '\t';

  if (isfinite(score)) oss << score;
  else oss << "NA";

  oss << '\t' << r.get_strand();

  if (report_more_information) {

    oss << '\t';
    if (isfinite(weighted_mean_meth)) oss << weighted_mean_meth;
    else oss << "NA";

    oss << '\t';
    if (isfinite(mean_meth)) oss << mean_meth;
    else oss << "NA";

    oss << '\t';
    if (isfinite(fractional_meth)) oss << fractional_meth;
    else oss << "NA";

    oss << '\t' << total_cpgs
        << '\t' << cpgs_with_reads
        << '\t' << total_meth
        << '\t' << total_reads;
  }
  return oss.str();
}


static void
process_with_cpgs_loaded(const bool VERBOSE,
                         const bool sort_data_if_needed,
                         const bool PRINT_NUMERIC_ONLY,
                         const bool report_more_information,
                         const char level_code,
                         const string &cpgs_file,
                         vector<GenomicRegion> &regions,
                         std::ostream &out) {

  bgzf_file in(cpgs_file, "r");
  if (!in) throw runtime_error("cannot open file: " + cpgs_file);

  vector<MSite> cpgs;
  MSite the_cpg;
  while (read_site(in, the_cpg))
    cpgs.push_back(the_cpg);

  if (!is_sorted(begin(cpgs), end(cpgs))) {
    if (sort_data_if_needed)
      sort(begin(cpgs), end(cpgs));
    else
      throw runtime_error("data not sorted: " + cpgs_file);
  }

  if (VERBOSE)
    cerr << "[n_cpgs=" << cpgs.size() << "]" << endl;

  for (size_t i = 0; i < regions.size(); ++i) {

    const pair<size_t, size_t> bounds(region_bounds(cpgs, regions[i]));

    size_t total_meth = 0, total_reads = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    double mean_meth = 0.0;

    for (size_t j = bounds.first; j < bounds.second; ++j) {
      if (cpgs[j].n_reads > 0) {
        total_meth += cpgs[j].n_meth();
        total_reads += cpgs[j].n_reads;
        ++cpgs_with_reads;

        auto calls = meth_unmeth_calls(cpgs[j].n_meth(), cpgs[j].n_unmeth());
        called_total += (calls.first || calls.second);
        called_meth += calls.first;

        mean_meth += cpgs[j].meth;
      }
    }

    const size_t total_cpgs = bounds.second - bounds.first;
    const double fractional_meth =
      static_cast<double>(called_meth)/called_total;
    const double weighted_mean_meth =
      static_cast<double>(total_meth)/total_reads;
    const double unweighted_mean_meth = mean_meth/cpgs_with_reads;

    const double score = (level_code == 'w' ? weighted_mean_meth :
                          (level_code == 'u' ? unweighted_mean_meth :
                           fractional_meth));

    if (!PRINT_NUMERIC_ONLY || isfinite(score))
      out << format_output_line(report_more_information, regions[i], score,
                                weighted_mean_meth, unweighted_mean_meth,
                                fractional_meth,
                                total_cpgs, cpgs_with_reads,
                                total_meth, total_reads)
          << endl;
  }
}


/*

/// ADS: using code from MSite.*pp files now which should be in common
/// for all counts files that must be searched on disk. Keeping this
/// here for now.

////////////////////////////////////////////////////////////////////////
///
///  CODE BELOW HERE IS FOR SEARCHING ON DISK
///

static void
move_to_start_of_line(ifstream &in) {
  char next;
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}

static void
find_start_line(const string &chr, const size_t idx, ifstream &cpg_in) {

  cpg_in.seekg(0, ios_base::beg);
  const size_t begin_pos = cpg_in.tellg();
  cpg_in.seekg(0, ios_base::end);
  const size_t end_pos = cpg_in.tellg();

  if (end_pos - begin_pos < 2)
    throw runtime_error("empty meth file");

  size_t step_size = (end_pos - begin_pos)/2;

  cpg_in.seekg(0, ios_base::beg);
  string low_chr;
  size_t low_idx = 0;
  cpg_in >> low_chr >> low_idx;

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  cpg_in.seekg(-2, ios_base::end);
  move_to_start_of_line(cpg_in);
  string high_chr;
  size_t high_idx;
  cpg_in >> high_chr >> high_idx;

  size_t pos = step_size;
  cpg_in.seekg(pos, ios_base::beg);
  move_to_start_of_line(cpg_in);

  while (step_size > 0) {
    string mid_chr;
    size_t mid_idx = 0;
    // ADS: possible here to check that the chroms are in the right
    // order, at least the ones we have seen.
    cpg_in >> mid_chr >> mid_idx;
    step_size /= 2;
    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      std::swap(mid_chr, high_chr);
      std::swap(mid_idx, high_idx);
      pos -= step_size;
    }
    else {
      std::swap(mid_chr, low_chr);
      std::swap(mid_idx, low_idx);
      pos += step_size;
    }
    cpg_in.seekg(pos, ios_base::beg);
    move_to_start_of_line(cpg_in);
  }
}
*/

static bool
cpg_not_past_region(const GenomicRegion &region, const size_t end_pos,
                    const MSite &cpg) {
  return (cpg.chrom == region.get_chrom() && cpg.pos < end_pos) ||
    cpg.chrom < region.get_chrom();
}


static void
get_cpg_stats(ifstream &cpg_in, const GenomicRegion region,
              size_t &total_meth, size_t &total_reads,
              size_t &total_cpgs, size_t &cpgs_with_reads,
              size_t &called_total, size_t &called_meth,
              double &mean_meth) {

  const string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_offset_for_msite(chrom, start_pos, cpg_in);

  MSite cpg;
  while (cpg_in >> cpg && cpg_not_past_region(region, end_pos, cpg)) {
    if (start_pos <= cpg.pos && cpg.chrom == chrom) {
      ++total_cpgs;
      if (cpg.n_reads > 0) {

        total_meth += cpg.n_meth();
        total_reads += cpg.n_reads;
        ++cpgs_with_reads;

        auto calls = meth_unmeth_calls(cpg.n_meth(), cpg.n_unmeth());
        called_total += (calls.first || calls.second);
        called_meth += calls.first;

        mean_meth += cpg.meth;
      }
    }
  }
  cpg_in.clear();
}


static void
process_with_cpgs_on_disk(const bool PRINT_NUMERIC_ONLY,
                          const bool report_more_information,
                          const char level_code,
                          const string &cpgs_file,
                          vector<GenomicRegion> &regions,
                          std::ostream &out) {

  ifstream in(cpgs_file);
  for (size_t i = 0; i < regions.size() && in; ++i) {

    size_t total_meth = 0, total_reads = 0;
    size_t cpgs_with_reads = 0;
    size_t called_total = 0, called_meth = 0;
    size_t total_cpgs = 0;
    double mean_meth = 0.0;

    get_cpg_stats(in, regions[i], total_meth, total_reads, total_cpgs, cpgs_with_reads,
                  called_total, called_meth, mean_meth);

    const double fractional_meth =
      static_cast<double>(called_meth)/called_total;
    const double weighted_mean_meth =
      static_cast<double>(total_meth)/total_reads;
    const double unweighted_mean_meth = mean_meth/cpgs_with_reads;

    const double score = (level_code == 'w' ? weighted_mean_meth :
                          (level_code == 'u' ? unweighted_mean_meth :
                           fractional_meth));

    if (!PRINT_NUMERIC_ONLY || isfinite(score))
      out << format_output_line(report_more_information, regions[i], score,
                                weighted_mean_meth, unweighted_mean_meth,
                                fractional_meth,
                                total_cpgs, cpgs_with_reads,
                                total_meth, total_reads)
          << endl;
  }
}

///
///  END OF CODE FOR SEARCHING ON DISK
///
////////////////////////////////////////////////////////////////////////


static size_t
check_bed_format(const string &regions_file) {

  ifstream in(regions_file);
  if (!in)
    throw runtime_error("cannot open file: " + regions_file);

  string line;
  getline(in, line);

  std::istringstream iss(line);
  string token;
  size_t n_columns = 0;
  while (iss >> token)
    ++n_columns;

  return n_columns;
}


int
main_roimethstat(int argc, const char **argv) {

  try {

    // ADS: this information should be somehow presented to the user
    // when the tool is run without arguments.
    static const string description = R"""(

Compute average CpG methylation levels in each interval from a given
set of genomic intervals

Columns (beyond the first 6) in the BED format output:
(7) weighted mean methylation
(8) unweighted mean methylation
(9) fractional methylation
(10) number of CpGs in the region
(11) number of CpGs covered at least once
(12) number of observations in reads indicating methylation
(13) total number of observations from reads in the region

)""";

    static const string default_name_prefix = "X";

    bool VERBOSE = false;
    bool PRINT_NUMERIC_ONLY = false;
    bool load_entire_file = false;
    bool report_more_information = false;
    bool sort_data_if_needed = false;

    string level_code = "w";

    string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "Compute average CpG "
                           "methylation levels in each interval from "
                           "a given set of genomic intervals",
                           "<intervals-bed> <methylation-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("numeric", 'N', "print numeric values only (not NAs)",
                      false, PRINT_NUMERIC_ONLY);
    opt_parse.add_opt("preload", 'L', "Load all CpG sites",
                      false, load_entire_file);
    opt_parse.add_opt("sort", 's', "sort data if needed",
                      false, sort_data_if_needed);
    opt_parse.add_opt("level", 'l', "the level to report as score column "
                      "in bed format output (w, u or f)", false, level_code);
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
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (sort_data_if_needed && !load_entire_file) {
      cerr << "cannot sort data unless all data is loaded" << endl;
      return EXIT_SUCCESS;
    }
    if (level_code != "w" && level_code != "u" && level_code != "f") {
      cerr << "selected level must be in {w, u, f}" << endl;
      return EXIT_SUCCESS;
    }
    const string regions_file = leftover_args.front();
    const string cpgs_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "loading regions" << endl;

    const size_t n_columns = check_bed_format(regions_file);
    // MAGIC: this allows for exactly 3 or at least 6 columns in the
    // bed format
    if (n_columns != 3 && n_columns < 6)
      throw runtime_error("format must be 3 or 6+ column bed: " + regions_file);
    if (is_msite_file(regions_file)) {
      cerr << opt_parse.help_message() << endl;
      throw runtime_error("The file seems to be a methylation file: " +
          regions_file + "\nCheck the order of the input arguments");
    }
    if (!is_msite_file(cpgs_file)) {
      cerr << opt_parse.help_message() << endl;
      throw runtime_error("The file is not a methylation file: " + cpgs_file);
    }


    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!is_sorted(begin(regions), end(regions))) {
      if (sort_data_if_needed) {
        if (VERBOSE)
          cerr << "sorting regions" << endl;
        sort(begin(regions), end(regions));
      }
      else throw runtime_error("regions not sorted in file: " + regions_file);
    }

    if (n_columns == 3) // then we should name the regions
      for (size_t i = 0; i < regions.size(); ++i)
        regions[i].set_name(default_name_prefix + std::to_string(i));

    if (VERBOSE)
      cerr << "[n_regions=" << regions.size() << "]" << endl;

    std::ofstream of;
    if (!outfile.empty()) {
      of.open(outfile);
      if (!of)
        throw runtime_error("failed to open outfile: " + outfile);
    }
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    if (load_entire_file)
      process_with_cpgs_loaded(VERBOSE, sort_data_if_needed,
                               PRINT_NUMERIC_ONLY,
                               report_more_information,
                               level_code[0],
                               cpgs_file, regions, out);
    else
      process_with_cpgs_on_disk(PRINT_NUMERIC_ONLY,
                                report_more_information,
                                level_code[0],
                                cpgs_file, regions, out);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
