/* merge-methcounts: a program for merging methcounts files
 *
 * Copyright (C) 2011-2025 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D Smith
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

#include "MSite.hpp"
#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-owning-memory,*-avoid-magic-numbers,*-narrowing-conversions,*-pointer-arithmetic)

static void
set_invalid(MSite &s) {
  s.pos = std::numeric_limits<std::size_t>::max();
}

static bool
is_valid(const MSite &s) {
  return s.pos != std::numeric_limits<std::size_t>::max();
}

[[nodiscard]] static bool
any_sites_unprocessed(
  const std::vector<std::string> &filenames,
  const std::vector<std::unique_ptr<bamxx::bgzf_file>> &infiles,
  std::vector<bool> &outdated, std::vector<MSite> &sites) {
  const std::size_t n_files = std::size(sites);

  bool sites_remain = false;
  for (std::size_t i = 0; i < n_files; ++i) {
    if (outdated[i]) {
      outdated[i] = false;
      MSite tmp_site;
      if (read_site(*infiles[i], tmp_site)) {
        // ADS: chrom order within a file already tested
        if (tmp_site.pos <= sites[i].pos && tmp_site.chrom == sites[i].chrom)
          throw std::runtime_error("sites not sorted in " + filenames[i]);
        sites_remain = true;
        sites[i] = tmp_site;
      }
      else
        set_invalid(sites[i]);
    }
    else if (is_valid(sites[i]))
      sites_remain = true;
  }
  return sites_remain;
}

[[nodiscard]] static bool
site_precedes(const std::vector<MSite> &sites,
              const std::unordered_map<std::string, std::size_t> &chroms_order,
              const std::size_t a, const std::size_t b) {
  if (chroms_order.empty())
    return sites[a] < sites[b];

  const std::size_t a_chr = chroms_order.find(sites[a].chrom)->second;
  const std::size_t b_chr = chroms_order.find(sites[b].chrom)->second;

  return a_chr < b_chr || (a_chr == b_chr && sites[a].pos < sites[b].pos);
}

[[nodiscard]] static std::size_t
find_minimum_site(
  const std::vector<MSite> &sites,
  const std::unordered_map<std::string, std::size_t> &chroms_order,
  const std::vector<bool> &outdated) {
  const std::size_t n_files = std::size(sites);

  // ms_id is the id of the minimum site, the next one to print
  std::size_t ms_id = std::numeric_limits<std::size_t>::max();

  // make sure there is at least one remaining site to print
  for (std::size_t i = 0;
       i < n_files && ms_id == std::numeric_limits<std::size_t>::max(); ++i)
    if (is_valid(sites[i]) && !outdated[i])
      ms_id = i;

  if (ms_id == std::numeric_limits<std::size_t>::max())
    throw std::runtime_error("failed in a next site to print");

  // now find the earliest site to print among those that could be
  for (std::size_t i = 0; i < n_files; ++i)
    if (!outdated[i] && is_valid(sites[i]) &&
        site_precedes(sites, chroms_order, i, ms_id))
      ms_id = i;

  return ms_id;
}

[[nodiscard]] static bool
same_location(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom && a.pos == b.pos;
}

[[nodiscard]] static std::size_t
collect_sites_to_print(
  const std::vector<MSite> &sites,
  const std::unordered_map<std::string, std::size_t> &chroms_order,
  const std::vector<bool> &outdated, std::vector<bool> &to_print) {
  const std::size_t n_files = std::size(sites);

  const std::size_t min_site_idx =
    find_minimum_site(sites, chroms_order, outdated);

  for (std::size_t i = 0; i < n_files; ++i)
    // condition below covers "is_valid(sites[i])"
    if (same_location(sites[min_site_idx], sites[i]))
      to_print[i] = true;

  return min_site_idx;
}

[[nodiscard]] static bool
any_mutated(const std::vector<bool> &to_print,
            const std::vector<MSite> &sites) {
  std::size_t i = 0;
  while (i < std::size(sites) && !(to_print[i] && sites[i].is_mutated()))
    ++i;
  return (i < std::size(sites));
}

template <std::size_t buffer_size>
static void
write_line_for_tabular(std::array<char, buffer_size> &buffer,
                       const bool write_fractional,
                       const bool report_any_mutated,
                       const std::size_t min_reads, std::ostream &out,
                       const std::vector<bool> &to_print,
                       const std::vector<MSite> &sites, MSite min_site) {
  const std::size_t n_files = std::size(sites);

  min_site.set_unmutated();
  if (report_any_mutated && any_mutated(to_print, sites))
    min_site.set_mutated();

  // ADS: is this the format we want for the row names?
  auto cursor = buffer.data();
  auto bytes_left = buffer_size;
  {
    const auto n_bytes =
      std::snprintf(cursor, bytes_left, "%s:%zu:%c:%s", min_site.chrom.data(),
                    min_site.pos, min_site.strand, min_site.context.data());
    if (n_bytes < 0 || n_bytes == static_cast<std::int32_t>(bytes_left))
      throw std::runtime_error("failed to write output line");
    cursor += n_bytes;
    bytes_left -= n_bytes;
  }

  if (write_fractional) {
    for (std::size_t i = 0; i < n_files; ++i) {
      const std::size_t r = sites[i].n_reads;
      const auto n_bytes = [&] {
        return (to_print[i] && r >= min_reads)
                 ? std::snprintf(cursor, bytes_left, "\t%.6g", sites[i].meth)
                 : std::snprintf(cursor, bytes_left, "\tNA");
      }();
      if (n_bytes < 0 || n_bytes == static_cast<std::int32_t>(bytes_left))
        throw std::runtime_error("failed to write output line");
      cursor += n_bytes;
      bytes_left -= n_bytes;
    }
  }
  else
    for (std::size_t i = 0; i < n_files; ++i) {
      const auto n_bytes = [&] {
        return to_print[i]
                 ? std::snprintf(cursor, bytes_left, "\t%zu\t%.6g",
                                 sites[i].n_reads, sites[i].n_meth_f())
                 : std::snprintf(cursor, bytes_left, "\t0\t0");
      }();
      if (n_bytes < 0 || n_bytes == static_cast<std::int32_t>(bytes_left))
        throw std::runtime_error("failed to write output line");
      cursor += n_bytes;
      bytes_left -= n_bytes;
    }

  if (static_cast<std::ptrdiff_t>(buffer_size) <=
      std::distance(buffer.data(), cursor))
    throw std::runtime_error("failed to write output line");

  *cursor++ = '\n';
  out.write(buffer.data(), std::distance(buffer.data(), cursor));
}

static void
write_line_for_merged_counts(std::ostream &out, const bool report_any_mutated,
                             const std::vector<bool> &to_print,
                             const std::vector<MSite> &sites, MSite min_site) {
  const std::size_t n_files = std::size(sites);

  min_site.set_unmutated();
  if (report_any_mutated && any_mutated(to_print, sites))
    min_site.set_mutated();

  double meth_sum = 0;
  min_site.n_reads = 0;
  for (std::size_t i = 0; i < n_files; ++i)
    if (to_print[i]) {
      meth_sum += sites[i].n_meth();
      min_site.n_reads += sites[i].n_reads;
    }
  min_site.meth = meth_sum / std::max(1ul, min_site.n_reads);

  out << min_site << '\n';
}

[[nodiscard]] static std::string
remove_extension(const std::string &filename) {
  const std::size_t last_dot = filename.find_last_of('.');
  if (last_dot == std::string::npos)
    return filename;
  else
    return filename.substr(0, last_dot);
}

[[nodiscard]] static std::string
remove_suffix(const std::string &suffix, const std::string &filename) {
  if (filename.substr(std::size(filename) - std::size(suffix),
                      std::size(suffix)) == suffix)
    return filename.substr(0, std::size(filename) - std::size(suffix));
  return filename;
}

static void
get_orders_by_file(const std::string &filename,
                   std::vector<std::string> &chroms_order) {
  bamxx::bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("bad file: " + filename);
  chroms_order.clear();

  std::unordered_set<std::string> chroms_seen;
  std::string line;
  std::string prev_chrom;

  while (getline(in, line)) {
    line.resize(line.find_first_of(" \t"));
    if (line != prev_chrom) {
      if (chroms_seen.find(line) != std::cend(chroms_seen))
        throw std::runtime_error("chroms out of order in: " + filename);
      chroms_seen.insert(line);
      chroms_order.push_back(line);
      std::swap(line, prev_chrom);
    }
  }
}

static void
get_chroms_order(const std::vector<std::string> &filenames,
                 std::unordered_map<std::string, std::size_t> &chroms_order) {
  // get order of chroms in each file
  std::vector<std::vector<std::string>> orders(std::size(filenames));
  for (std::size_t i = 0; i < std::size(filenames); ++i)
    get_orders_by_file(filenames[i], orders[i]);

  // get the union of chrom sets
  std::unordered_set<std::string> the_union;
  for (std::size_t i = 0; i < std::size(orders); ++i)
    for (std::size_t j = 0; j < std::size(orders[i]); ++j)
      the_union.insert(orders[i][j]);

  // get an adjacency list and in-degree for each node
  std::unordered_map<std::string, std::vector<std::string>> adj_lists;
  std::unordered_map<std::string, std::size_t> in_degree;
  for (auto &&i : the_union) {
    in_degree[i] = 0;
    adj_lists[i] = std::vector<std::string>();
  }
  for (auto &&i : orders) {
    auto j = std::cbegin(i);
    for (auto k = j + 1; k != std::cend(i); ++j, ++k) {
      adj_lists[*j].push_back(*k);
      ++in_degree[*k];
    }
  }

  std::queue<std::string> q;  // invariant: nodes with no incoming edge
  for (auto &&i : the_union)
    if (in_degree[i] == 0)
      q.push(i);

  while (!q.empty()) {
    const std::string u = q.front();
    q.pop();

    // iterate over the edges (u, v)
    for (auto &&v : adj_lists[u]) {
      --in_degree[v];         // degree should not appear here as 0
      if (in_degree[v] == 0)  // this should only happen once per v
        q.push(v);
    }
    adj_lists[u].clear();  // delete node; already had in_degree 0

    chroms_order.insert(std::make_pair(u, std::size(chroms_order)));
  }

  // finally, make sure we found a consistent order
  for (auto &&i : adj_lists)
    if (!i.second.empty())
      throw std::runtime_error("inconsistent order of chroms between files");
}

/*
  This utility does two things, and they are grouped together here because of
  how they are done, not because the uses are related. (1) merge-methcounts
  can take a set of methcounts output files and combine them into one. There
  are several reasons a user might want to do this. An example is when
  technical replicates are performed, and analyzed separately to understand
  technical variance (e.g., between sequencing runs or library preps). After
  examining the technical variation, subsequent analyses might be best
  conducted on all the data together. So all the methcounts files can be
  combined into one using merge-methcounts. In this case, the coverage at any
  site is the sum of the coverages in the original methcounts files, and the
  methylation level at any site is the weighted mean. (2) merge-methcounts can
  take a set of methcounts output files, and create a table that contains all
  the same information. The table format is helpful if subsequent analyses are
  best done using a data table, for example in R. When producing a tabular
  format, merge-methcounts allows the user to select whether the desired
  output is in counts or fractions.
 */
int
main_merge_methcounts(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  static constexpr auto buffer_size = 65536;
  try {
    static const std::string description = "merge multiple methcounts files";

    std::string outfile;
    bool verbose = false;
    bool report_any_mutated = false;
    bool write_tabular_format = false;
    bool write_fractional = false;
    bool radmeth_format = false;
    bool ignore_chroms_order = false;
    bool add_first_column_header = false;

    std::string header_info;
    std::string column_name_suffix = "RM";
    std::string suffix_to_remove;

    std::size_t min_reads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<methcounts-files>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)", true,
                      outfile);
    opt_parse.add_opt("header", 'H', "header to print (ignored for tabular)",
                      false, header_info);
    opt_parse.add_opt("tabular", 't', "output as table", false,
                      write_tabular_format);
    opt_parse.add_opt("radmeth", '\0',
                      "make radmeth input "
                      "(use with tabular and without fractional)",
                      false, radmeth_format);
    opt_parse.add_opt("remove", '\0',
                      "remove this suffix from filenames when "
                      "making column names; default removes from final dot",
                      false, suffix_to_remove);
    opt_parse.add_opt("suff", 's',
                      "column name suffixes, one for total reads and one for "
                      "methylated reads, to be separated from sample name "
                      "with underscore",
                      false, column_name_suffix);
    opt_parse.add_opt("fractional", 'f', "output fractions (requires tabular)",
                      false, write_fractional);
    opt_parse.add_opt("reads", 'r', "min reads (for fractional)", false,
                      min_reads);
    opt_parse.add_opt("ignore", '\0',
                      "Do not attempt to determine chromosome. "
                      "Lexicographic order will be used.",
                      false, ignore_chroms_order);
    opt_parse.add_opt("mut", 'm',
                      "If any of the sites being merged indicates "
                      "mutated, mark the result has mutated.",
                      false, report_any_mutated);
    opt_parse.add_opt("1st-col-header", '\0', "add name for 1st col in header",
                      false, add_first_column_header);
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
    if (write_fractional && !write_tabular_format) {
      std::cerr << "fractional output only available for tabular format\n";
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (std::size(column_name_suffix) != 2) {
      std::cerr << "column name suffix must be 2 letters\n";
      return EXIT_SUCCESS;
    }
    std::vector<std::string> meth_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    std::unordered_map<std::string, std::size_t> chroms_order;
    if (!ignore_chroms_order) {
      if (verbose)
        std::cerr << "resolving chromosome order\n";
      get_chroms_order(meth_files, chroms_order);
      if (verbose) {
        std::cerr << "chromosome order\n";
        std::vector<std::string> v(std::size(chroms_order));
        for (auto &&i : chroms_order)
          v[i.second] = i.first;
        for (auto &&i : v)
          std::cerr << i << '\n';
      }
    }

    const std::size_t n_files = std::size(meth_files);

    std::vector<std::unique_ptr<bamxx::bgzf_file>> infiles(n_files);
    for (std::size_t i = 0; i < n_files; ++i) {
      infiles[i] = std::make_unique<bamxx::bgzf_file>(meth_files[i], "r");
      if (!(*infiles[i]))
        throw std::runtime_error("cannot open file: " + meth_files[i]);
    }

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);

    // print header if user specifies or if tabular output format
    if (write_tabular_format) {
      std::vector<std::string> colnames(std::size(meth_files));
      std::transform(std::cbegin(meth_files), std::cend(meth_files),
                     std::begin(colnames), [&](const auto &x) {
                       auto fn = std::filesystem::path{x}.filename().string();
                       return suffix_to_remove.empty()
                                ? remove_extension(fn)
                                : remove_suffix(suffix_to_remove, fn);
                     });

      if (!write_fractional && !radmeth_format) {
        std::vector<std::string> tmp;
        for (auto &&i : colnames) {
          tmp.push_back(i + "_" + column_name_suffix[0]);
          tmp.push_back(i + "_" + column_name_suffix[1]);
        }
        std::swap(tmp, colnames);
      }

      if (add_first_column_header)
        out << "site\t";
      std::copy(std::cbegin(colnames), std::cend(colnames),
                std::ostream_iterator<std::string>(out, "\t"));
      out << '\n';
    }
    else if (!header_info.empty())
      out << "#" << header_info << '\n';

    std::vector<MSite> sites(n_files);
    std::vector<bool> outdated(n_files, true);
    std::vector<bool> sites_to_print;  // declared here to keep allocation
    std::vector<std::unordered_set<std::string>> chroms_seen(n_files);

    std::array<char, buffer_size> buffer{};

    while (any_sites_unprocessed(meth_files, infiles, outdated, sites)) {
      sites_to_print.clear();
      sites_to_print.resize(n_files, false);

      // below idx is one index among the sites to print
      const std::size_t idx =
        collect_sites_to_print(sites, chroms_order, outdated, sites_to_print);

      // output the appropriate sites' data
      if (write_tabular_format)
        write_line_for_tabular(buffer, write_fractional, report_any_mutated,
                               min_reads, out, sites_to_print, sites,
                               sites[idx]);
      else
        write_line_for_merged_counts(out, report_any_mutated, sites_to_print,
                                     sites, sites[idx]);

      std::swap(outdated, sites_to_print);
    }
  }
  catch (const std::runtime_error &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-owning-memory,*-avoid-magic-numbers,*-narrowing-conversions,*-pointer-arithmetic)
