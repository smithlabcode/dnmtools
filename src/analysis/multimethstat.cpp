/* Copyright (C) 2014-2022 Andrew D. Smith
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
multistat: average methylation in each of a set of regions for each methylome
)";

#include "Interval6.hpp"
#include "MSite.hpp"
#include "bsutils.hpp"

#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

static std::pair<bool, bool>
meth_unmeth_calls(const std::size_t n_meth, const std::size_t n_unmeth) {
  static constexpr double alpha = 0.95;
  static constexpr double one_half = 0.5;
  double lower{}, upper{};
  const auto total = n_meth + n_unmeth;
  wilson_ci_for_binomial(alpha, total, static_cast<double>(n_meth) / total,
                         lower, upper);
  return std::make_pair(lower > one_half, upper < one_half);
}

static std::pair<std::size_t, std::size_t>
region_bounds(const std::vector<MSite> &sites, const Interval6 &region) {
  const MSite a(region.chrom, region.start, region.strand, "", 0, 0);
  const auto a_ins = std::lower_bound(std::cbegin(sites), std::cend(sites), a);
  const MSite b(region.chrom, region.stop, region.strand, "", 0, 0);
  const auto b_ins = std::lower_bound(std::cbegin(sites), std::cend(sites), b);
  return std::make_pair(std::distance(std::cbegin(sites), a_ins),
                        std::distance(std::cbegin(sites), b_ins));
}

/*
   This function is used to make sure the output is done the same way in the
   different places it might be generated. One issue is that there is a lot of
   string and stream manipulation here, which might not be desirable if there
   are lots of regions, as when the genome is binned.
*/
static std::string
format_output_line(const bool report_more_information, const Interval6 &r,
                   const std::vector<double> &score,
                   const std::vector<double> &weighted_mean_meth,
                   const std::vector<double> &mean_meth,
                   const std::vector<double> &fractional_meth,
                   const std::size_t total_cpgs,
                   const std::vector<std::size_t> &cpgs_with_reads,
                   const std::vector<std::size_t> &total_meth,
                   const std::vector<std::size_t> &total_reads) {
  std::ostringstream oss;
  oss << r.chrom << '\t' << r.start << '\t' << r.stop << '\t' << r.name << '\t';

  const std::size_t n_samples = std::size(score);

  for (std::size_t i = 0; i < n_samples; ++i) {
    if (std::isfinite(score[i]))
      oss << score[i];
    else
      oss << "NA";

    // GS: commenting this out because we shouldn't be reporting
    // strand several times
    // oss << '\t' << r.strand;

    if (report_more_information) {
      oss << '\t';
      if (std::isfinite(weighted_mean_meth[i]))
        oss << weighted_mean_meth[i];
      else
        oss << "NA";

      oss << '\t';
      if (std::isfinite(mean_meth[i]))
        oss << mean_meth[i];
      else
        oss << "NA";

      oss << '\t';
      if (std::isfinite(fractional_meth[i]))
        oss << fractional_meth[i];
      else
        oss << "NA";

      oss << '\t' << total_cpgs << '\t' << cpgs_with_reads[i] << '\t'
          << total_meth[i] << '\t' << total_reads[i];
    }
    if (i < n_samples - 1)
      oss << '\t';
  }
  return oss.str();
}

static void
get_site_and_values(std::string &line, MSite &site,
                    std::vector<double> &values) {
  // in tabular, chrom/pos/strand is separated by colons
  std::replace(std::begin(line), std::end(line), ':', '\t');
  std::istringstream iss(line);
  iss >> site.chrom >> site.pos >> site.strand >> site.context;
  values.clear();
  double v = 0.0;
  while (iss >> v)
    values.push_back(v);
}

static bool
all_finite(const std::vector<double> &scores) {
  return std::all_of(std::cbegin(scores), std::cend(scores),
                     [](const auto x) { return std::isfinite(x); });
}

static void
process_with_cpgs_loaded(const bool verbose,
                         // const bool sort_data_if_needed,
                         const bool print_numeric_only,
                         const bool report_more_information,
                         const char level_code, const std::string &cpgs_file,
                         const std::vector<Interval6> &regions,
                         std::ostream &out) {
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("cannot open file: " + cpgs_file);

  std::string header;
  getline(in, header);

  std::istringstream iss(header);
  std::string col_name;
  std::vector<std::string> col_names;
  while (iss >> col_name)
    col_names.push_back(col_name);

  const auto end_uniq = std::unique(std::begin(col_names), std::end(col_names));
  const auto n_uniq = std::distance(std::begin(col_names), end_uniq);
  // make sure all columns come in pairs of counts
  if (2 * static_cast<std::size_t>(n_uniq) != std::size(col_names))
    throw std::runtime_error("wrong header format:\n" + header);

  col_names.resize(std::distance(std::begin(col_names), end_uniq));
  const std::size_t n_samples = std::size(col_names);
  if (verbose)
    std::cerr << "[n_samples=" << n_samples << "]\n";

  std::vector<MSite> cpgs;
  std::vector<std::vector<double>> values;
  MSite the_cpg;
  std::string line;
  while (getline(in, line)) {
    values.push_back(std::vector<double>());
    get_site_and_values(line, the_cpg, values.back());
    cpgs.push_back(the_cpg);

    if (std::size(values.back()) != 2 * n_samples)
      throw std::runtime_error(
        "wrong row length. Expected " + std::to_string(2 * n_samples) +
        ", got " + std::to_string(std::size(values.back())) + "\n" + line);
  }

  if (verbose)
    std::cerr << "[n_cpgs=" << std::size(cpgs) << "]\n";

  if (!std::is_sorted(std::cbegin(cpgs), std::cend(cpgs))) {
    /* GS: need a way to sort cpgs and values to preserve the order.
     * Commenting this for now until I know the best way to do this
    if (sort_data_if_needed)
      std::sort(begin(cpgs), end(cpgs));
    else*/
    throw std::runtime_error("data not sorted: " + cpgs_file);
  }

  if (verbose)
    std::cerr << "[n_cpgs=" << std::size(cpgs) << "]\n";

  // write header as first column
  for (std::size_t i = 0; i < n_samples; ++i)
    out << '\t' << col_names[i];
  out << '\n';

  // preallocate entire space for a single row
  std::vector<std::size_t> total_meth(n_samples, 0);
  std::vector<std::size_t> total_reads(n_samples, 0);
  std::vector<std::size_t> cpgs_with_reads(n_samples, 0);
  std::vector<std::size_t> called_total(n_samples, 0);
  std::vector<std::size_t> called_meth(n_samples, 0);
  std::vector<double> weighted_mean_meth(n_samples, 0);
  std::vector<double> fractional_meth(n_samples, 0);
  std::vector<double> unweighted_mean_meth(n_samples, 0);
  std::vector<double> score(n_samples, 0);

  for (std::size_t i = 0; i < std::size(regions); ++i) {
    const std::pair<std::size_t, std::size_t> bounds(
      region_bounds(cpgs, regions[i]));

    for (std::size_t j = 0; j < n_samples; ++j) {
      total_meth[j] = total_reads[j] = cpgs_with_reads[j] = called_total[j] =
        called_meth[j] = 0;
      double mean_meth = 0.0;

      for (std::size_t k = bounds.first; k < bounds.second; ++k) {
        const std::size_t meth = values[k][2 * j + 1];
        const std::size_t tot = values[k][2 * j];
        if (tot > 0) {
          total_reads[j] += tot;
          total_meth[j] += meth;
          ++cpgs_with_reads[j];

          const auto calls = meth_unmeth_calls(meth, tot - meth);
          called_total[j] += (calls.first || calls.second);
          called_meth[j] += calls.first;

          mean_meth += meth / static_cast<double>(tot);
        }
      }

      fractional_meth[j] =
        static_cast<double>(called_meth[j]) / called_total[j];
      weighted_mean_meth[j] =
        static_cast<double>(total_meth[j]) / total_reads[j];
      unweighted_mean_meth[j] = mean_meth / cpgs_with_reads[j];

      score[j] =
        (level_code == 'w' ? weighted_mean_meth[j]
                           : (level_code == 'u' ? unweighted_mean_meth[j]
                                                : fractional_meth[j]));
    }

    const std::size_t total_cpgs = bounds.second - bounds.first;
    if (!print_numeric_only || all_finite(score))
      out << format_output_line(report_more_information, regions[i], score,
                                weighted_mean_meth, unweighted_mean_meth,
                                fractional_meth, total_cpgs, cpgs_with_reads,
                                total_meth, total_reads)
          << '\n';
  }
}

////////////////////////////////////////////////////////////////////////
///
///  CODE BELOW HERE IS FOR SEARCHING ON DISK
///

static void
move_to_start_of_line(std::ifstream &in) {
  char next{};
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}

static void
get_chr_and_idx(std::ifstream &in, std::string &chrom, std::size_t &pos) {
  std::string tmp;
  in >> tmp;
  std::replace(std::begin(tmp), std::end(tmp), ':', ' ');
  std::istringstream iss(tmp);
  iss >> chrom >> pos;
}

static void
find_start_line(const std::string &chr, const std::size_t idx,
                std::ifstream &cpg_in) {
  cpg_in.seekg(0, std::ios_base::beg);
  const std::size_t begin_pos = cpg_in.tellg();
  cpg_in.seekg(0, std::ios_base::end);
  const std::size_t end_pos = cpg_in.tellg();

  if (end_pos - begin_pos < 2)
    throw std::runtime_error("empty meth file");

  std::size_t step_size = (end_pos - begin_pos) / 2;

  cpg_in.seekg(0, std::ios_base::beg);
  std::string low_chr;
  std::size_t low_idx = 0;
  get_chr_and_idx(cpg_in, low_chr, low_idx);

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  cpg_in.seekg(-2, std::ios_base::end);
  move_to_start_of_line(cpg_in);
  std::string high_chr;
  std::size_t high_idx{};
  get_chr_and_idx(cpg_in, high_chr, high_idx);

  std::size_t pos = step_size;
  cpg_in.seekg(pos, std::ios_base::beg);
  move_to_start_of_line(cpg_in);

  while (step_size > 0) {
    std::string mid_chr;
    std::size_t mid_idx = 0;
    // ADS: possible here to check that the chroms are in the right
    // order, at least the ones we have seen.
    get_chr_and_idx(cpg_in, mid_chr, mid_idx);
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
    cpg_in.seekg(pos, std::ios_base::beg);
    move_to_start_of_line(cpg_in);
  }
}

static bool
cpg_not_past_region(const Interval6 &region, const std::size_t end_pos,
                    const MSite &cpg) {
  return (cpg.chrom == region.chrom && cpg.pos < end_pos) ||
         cpg.chrom < region.chrom;
}

static void
get_cpg_stats(std::ifstream &cpg_in, const Interval6 &region,
              std::vector<std::size_t> &total_meth,
              std::vector<std::size_t> &total_reads, std::size_t &total_cpgs,
              std::vector<std::size_t> &cpgs_with_reads,
              std::vector<std::size_t> &called_total,
              std::vector<std::size_t> &called_meth,
              std::vector<double> &mean_meth) {
  const std::string chrom(region.chrom);
  const std::size_t start_pos = region.start;
  const std::size_t end_pos = region.stop;
  find_start_line(chrom, start_pos, cpg_in);

  const std::size_t n_samples = std::size(total_meth);

  MSite cpg;
  std::vector<double> values;
  values.reserve(2 * n_samples);
  std::string line;
  while (getline(cpg_in, line)) {
    values.clear();
    get_site_and_values(line, cpg, values);
    if (std::size(values) != 2 * n_samples)
      throw std::runtime_error("wrong row length. Expected " +
                               std::to_string(2 * n_samples) + ", got " +
                               std::to_string(std::size(values)) + "\n" + line);

    if (cpg_not_past_region(region, end_pos, cpg)) {
      if (start_pos <= cpg.pos && cpg.chrom == chrom) {
        ++total_cpgs;
        for (std::size_t j = 0; j < n_samples; ++j) {
          const std::size_t meth = values[2 * j + 1];
          const std::size_t tot = values[2 * j];
          if (values[2 * j] > 0) {
            total_meth[j] += meth;
            total_reads[j] += tot;
            ++cpgs_with_reads[j];

            const auto calls = meth_unmeth_calls(meth, tot - meth);
            called_total[j] += (calls.first || calls.second);
            called_meth[j] += calls.first;

            mean_meth[j] += meth / static_cast<double>(tot);
          }
        }
      }
    }
  }
  cpg_in.clear();
}

static void
process_with_cpgs_on_disk(const bool print_numeric_only,
                          const bool report_more_information,
                          const char level_code, const std::string &cpgs_file,
                          const std::vector<Interval6> &regions,
                          std::ostream &out) {
  std::ifstream in(cpgs_file);
  std::string header;
  getline(in, header);

  std::istringstream iss(header);
  std::string col_name;
  std::vector<std::string> col_names;
  while (iss >> col_name)
    col_names.push_back(col_name);

  const auto end_uniq = std::unique(std::begin(col_names), std::end(col_names));
  const auto n_uniq = std::distance(std::begin(col_names), end_uniq);
  // make sure all columns come in pairs of counts
  if (2 * static_cast<std::size_t>(n_uniq) != std::size(col_names))
    throw std::runtime_error("wrong header format:\n" + header);

  col_names.resize(std::distance(std::begin(col_names), end_uniq));
  const std::size_t n_samples = std::size(col_names);

  // write header as first column
  for (std::size_t i = 0; i < n_samples; ++i)
    out << '\t' << col_names[i];
  out << '\n';

  // preallocate entire space for a single row
  std::vector<std::size_t> total_meth(n_samples, 0);
  std::vector<std::size_t> total_reads(n_samples, 0);
  std::vector<std::size_t> cpgs_with_reads(n_samples, 0);
  std::vector<std::size_t> called_total(n_samples, 0);
  std::vector<std::size_t> called_meth(n_samples, 0);
  std::vector<double> weighted_mean_meth(n_samples, 0);
  std::vector<double> fractional_meth(n_samples, 0);
  std::vector<double> unweighted_mean_meth(n_samples, 0);
  std::vector<double> score(n_samples, 0.0);

  for (std::size_t i = 0; i < std::size(regions) && in; ++i) {
    std::size_t total_cpgs = 0;
    auto mean_meth = std::vector<double>(n_samples, 0.0);

    for (std::size_t j = 0; j < n_samples; ++j) {
      total_meth[j] = total_reads[j] = cpgs_with_reads[j] = called_total[j] =
        called_meth[j] = 0;
    }

    get_cpg_stats(in, regions[i], total_meth, total_reads, total_cpgs,
                  cpgs_with_reads, called_total, called_meth, mean_meth);

    for (std::size_t j = 0; j < n_samples; ++j) {
      fractional_meth[j] =
        static_cast<double>(called_meth[j]) / called_total[j];
      weighted_mean_meth[j] =
        static_cast<double>(total_meth[j]) / total_reads[j];
      unweighted_mean_meth[j] = mean_meth[j] / cpgs_with_reads[j];
      score[j] =
        (level_code == 'w' ? weighted_mean_meth[j]
                           : (level_code == 'u' ? unweighted_mean_meth[j]
                                                : fractional_meth[j]));
    }

    if (!print_numeric_only || all_finite(score))
      out << format_output_line(report_more_information, regions[i], score,
                                weighted_mean_meth, unweighted_mean_meth,
                                fractional_meth, total_cpgs, cpgs_with_reads,
                                total_meth, total_reads)
          << '\n';
  }
}

// end of code for searching on disk

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
main_multimethstat(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr auto default_name_prefix = "X";
    constexpr auto valid_n_cols = [](const auto n_cols) {
      return n_cols != 3 && n_cols < 6;  // NOLINT(*-avoid-magic-numbers)
    };

    bool verbose{false};
    bool print_numeric_only{false};
    bool load_entire_file{false};
    bool report_more_information{false};
    bool sort_data_if_needed{false};
    std::string level_code = "w";
    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "Compute average CpG "
                           "methylation in each of a set of genomic intervals",
                           "<intervals-bed> <methylation-file>");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("numeric", 'N', "print numeric values only (not NAs)",
                      false, print_numeric_only);
    opt_parse.add_opt("preload", 'L', "Load all CpG sites", false,
                      load_entire_file);
    opt_parse.add_opt("sort", 's', "sort data if needed", false,
                      sort_data_if_needed);
    opt_parse.add_opt("level", 'l',
                      "the level to report as score column "
                      "in bed format output (w, u or f)",
                      false, level_code);
    opt_parse.add_opt("more-levels", 'M', "report more methylation information",
                      false, report_more_information);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
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
    if (sort_data_if_needed && !load_entire_file) {
      std::cerr << "cannot sort data unless all data is loaded\n";
      return EXIT_SUCCESS;
    }
    if (level_code != "w" && level_code != "u" && level_code != "f") {
      std::cerr << "selected level must be in {w, u, f}\n";
      return EXIT_SUCCESS;
    }
    const std::string regions_file = leftover_args.front();
    const std::string cpgs_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (verbose)
      std::cerr << "loading regions\n";

    const auto n_columns = get_bed_columns(regions_file);
    if (valid_n_cols(n_columns))
      throw std::runtime_error("format must be 3 or 6+ column bed: " +
                               regions_file);
    auto regions = read_intervals6(regions_file);
    if (!std::is_sorted(std::cbegin(regions), std::cend(regions))) {
      if (sort_data_if_needed) {
        if (verbose)
          std::cerr << "sorting regions\n";
        std::sort(std::begin(regions), std::end(regions));
      }
      else
        throw std::runtime_error("regions not sorted in file: " + regions_file);
    }

    if (n_columns == 3)  // then we should name the regions
      for (std::size_t i = 0; i < std::size(regions); ++i)
        regions[i].name = default_name_prefix + std::to_string(i);

    if (verbose)
      std::cerr << "[n_regions=" << std::size(regions) << "]\n";

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);
    if (load_entire_file)
      process_with_cpgs_loaded(verbose,  // sort_data_if_needed,
                               print_numeric_only, report_more_information,
                               level_code[0], cpgs_file, regions, out);
    else
      process_with_cpgs_on_disk(print_numeric_only, report_more_information,
                                level_code[0], cpgs_file, regions, out);
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
