/* merge-xsym: a program for merging xcounts format xsym files.
 *
 * Copyright (C) 2024 Andrew D. Smith
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

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <charconv>

#include <bamxx.hpp>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "MSite.hpp"
#include "counts_header.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::copy;
using std::runtime_error;
using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::from_chars;
using std::ostream;
using std::ostream_iterator;

using bamxx::bgzf_file;

template<typename T>
using num_lim = std::numeric_limits<T>;

static void
set_invalid(MSite &s) {s.pos = num_lim<uint64_t>::max();}

static bool
is_valid(const MSite &s) {return s.pos != num_lim<uint64_t>::max();}


inline bamxx::bgzf_file &
read_site(bamxx::bgzf_file &f,
          string &chrom_name,
          uint64_t &offset,
          MSite &s) {
  static constexpr uint32_t buffer_size = 1024;

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, buffer_size);
  if (ret) throw runtime_error("failed to acquire buffer");

  if (!getline(f, line))
    return f;

  assert(line.s[0] != "#");

  if (!std::isdigit(line.s[0])) {
    chrom_name = line.s;
    offset = 0;
    if (!getline(f, line))
      return f;
    assert(line.s[0] != "#");
  }

  uint32_t pos_step = 0, n_meth = 0, n_unmeth = 0;
  const auto end_line = line.s + line.l;
  auto res = from_chars(line.s, end_line, pos_step);
  res = from_chars(res.ptr + 1, end_line, n_meth);
  res = from_chars(res.ptr + 1, end_line, n_unmeth);

  offset += pos_step;

  s = MSite();
  s.chrom = chrom_name;
  s.pos = offset;
  s.n_reads = n_meth + n_unmeth;
  s.meth = static_cast<double>(n_meth)/s.n_reads;
  ks_free(&line);

  return f;
}


static bool
any_sites_unprocessed(const vector<string> &filenames,
                      const vector<bgzf_file*> &infiles,
                      vector<bool> &outdated,
                      vector<string> &chrom_names,
                      vector<uint64_t> &offsets,
                      vector<MSite> &sites) {

  const uint32_t n_files = sites.size();

  bool sites_remain = false;
  for (uint32_t i = 0; i < n_files; ++i) {
    if (outdated[i]) {
      outdated[i] = false;
      MSite tmp_site;
      if (read_site(*infiles[i], chrom_names[i], offsets[i], tmp_site)) {
        // ADS: chrom order within a file already tested
        if (tmp_site.pos <= sites[i].pos && tmp_site.chrom == sites[i].chrom)
          throw runtime_error("sites not sorted in " + filenames[i]);
        assert(tmp_site.chrom != string("#"));
        sites_remain = true;
        sites[i] = tmp_site;
      }
      else set_invalid(sites[i]);
    }
    else if (is_valid(sites[i]))
      sites_remain = true;
  }
  return sites_remain;
}


static bool
site_precedes(const vector<MSite> &sites,
              const unordered_map<string, uint32_t> &chroms_order,
              const uint32_t a, const uint32_t b) {
  if (chroms_order.empty())
    return sites[a] < sites[b];

  const uint32_t a_chr = chroms_order.find(sites[a].chrom)->second;
  const uint32_t b_chr = chroms_order.find(sites[b].chrom)->second;

  return ((a_chr < b_chr) ||
          (a_chr == b_chr && sites[a].pos < sites[b].pos));
}


static uint32_t
find_minimum_site(const vector<MSite> &sites,
                  const unordered_map<string, uint32_t> &chroms_order,
                  const vector<bool> &outdated) {

  const uint32_t n_files = sites.size();

  // ms_id is the id of the minimum site, the next one to print
  uint32_t ms_id = num_lim<uint32_t>::max();

  // make sure there is at least one remaining site to print
  for (uint32_t i = 0; i < n_files && ms_id == num_lim<uint32_t>::max(); ++i)
    if (is_valid(sites[i]) && !outdated[i]) ms_id = i;

  if (ms_id == num_lim<uint32_t>::max())
    throw runtime_error("failed in a next site to print");

  // now find the earliest site to print among those that could be
  for (uint32_t i = 0; i < n_files; ++i)
    if (!outdated[i] && is_valid(sites[i]) &&
        site_precedes(sites, chroms_order, i, ms_id))
      ms_id = i;

  return ms_id;
}


static bool
same_location(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom && a.pos == b.pos;
}


static uint32_t
collect_sites_to_print(const vector<MSite> &sites,
                       const unordered_map<string, uint32_t> &chroms_order,
                       const vector<bool> &outdated,
                       vector<bool> &to_print) {

  const uint32_t n_files = sites.size();

  const uint32_t min_site_idx = find_minimum_site(sites, chroms_order, outdated);

  for (uint32_t i = 0; i < n_files; ++i)
    // condition below covers "is_valid(sites[i])"
    if (same_location(sites[min_site_idx], sites[i]))
      to_print[i] = true;

  return min_site_idx;
}


static bool
any_mutated(const vector<bool> &to_print,
            const vector<MSite> &sites) {
  uint32_t i = 0;
  while (i < sites.size() && !(to_print[i] && sites[i].is_mutated()))
    ++i;
  return (i < sites.size());
}


static void
write_line_for_tabular(const bool write_fractional,
                       const bool report_any_mutated,
                       const uint32_t min_reads,
                       std::ostream &out,
                       const vector<bool> &to_print,
                       const vector<MSite> &sites,
                       MSite min_site) {

  const uint32_t n_files = sites.size();

  min_site.set_unmutated();
  if (report_any_mutated && any_mutated(to_print, sites))
    min_site.set_mutated();

  // ADS: this format assumes the context is CpG and the strand is '+'
  out << min_site.chrom << ':' << min_site.pos;
  if (write_fractional) {
    for (uint32_t i = 0; i < n_files; ++i) {
      const uint32_t r = sites[i].n_reads;
      if (to_print[i] && r >= min_reads)
        out << '\t' << sites[i].meth;
      else
        out << '\t' << "NA";
    }
  }
  else
    for (uint32_t i = 0; i < n_files; ++i) {
      if (to_print[i])
        out << '\t' << sites[i].n_reads << '\t' << sites[i].n_meth();
      else out << '\t' << 0 << '\t'<< 0;
    }
  out << '\n';
}


static void
write_line_for_merged_counts(std::ostream &out,
                             const bool report_any_mutated,
                             const vector<bool> &to_print,
                             const vector<MSite> &sites,
                             MSite min_site) {

  const uint32_t n_files = sites.size();

  min_site.set_unmutated();
  if (report_any_mutated && any_mutated(to_print, sites))
    min_site.set_mutated();

  double meth_sum = 0;
  min_site.n_reads = 0;
  for (uint32_t i = 0; i < n_files; ++i)
    if (to_print[i]) {
      meth_sum += sites[i].n_meth();
      min_site.n_reads += sites[i].n_reads;
    }
  min_site.meth = meth_sum/std::max(1ul, min_site.n_reads);
  min_site.context = "CpG";

  out << min_site << '\n';
}


static string
remove_extension(const std::string &filename) {
  const uint32_t last_dot = filename.find_last_of(".");
  if (last_dot == std::string::npos) return filename;
  else return filename.substr(0, last_dot);
}


static string
remove_suffix(const string &suffix, const std::string &filename) {
  if (filename.substr(filename.size() - suffix.size(), suffix.size()) == suffix)
    return filename.substr(0, filename.size() - suffix.size());
  return filename;
}


static void
read_past_header(bamxx::bgzf_file &in) {
  static constexpr uint32_t buffer_size = 1024;

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, buffer_size);
  if (ret) throw runtime_error("failed to acquire buffer");

  while (getline(in, line) && is_counts_header_line(line.s) &&
         string(line.s) != string("#"))
    ;
  ks_free(&line);
}


static vector<string>
parse_xcounts_header(const string &filename) {
  static constexpr uint32_t buffer_size = 1024;

  bamxx::bgzf_file in(filename, "r");
  if (!in) throw runtime_error("failed to open input file");

  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, buffer_size);
  if (ret) throw runtime_error("failed to acquire buffer");

  vector<string> v;
  while (getline(in, line) && is_counts_header_line(line.s))
    v.emplace_back(line.s);
  ks_free(&line);

  return v;
}


static unordered_map<string, uint32_t>
get_chroms_order(const vector<string> &header) {
  unordered_map<string, uint32_t> chroms_order;
  uint32_t order = 0;
  for (const auto &h: header)
    if (h != "#" && h.find("DNMTOOLS") == string::npos)
      chroms_order.emplace(h.substr(1, h.find_first_of(" ") - 1), order++);
  return chroms_order;
}


int
main_merge_xsym(int argc, const char **argv) {
  /* This utility does two things, and they the same two things
   * explained in the comment just above the main in
   * merge-methcounts.cpp
   */
  try {

    static const string description = "merge multiple xsym files";

    string outfile;
    bool VERBOSE;
    bool report_any_mutated = false;
    bool write_tabular_format = false;
    bool write_fractional = false;
    bool radmeth_format = false;
    bool ignore_chroms_order = false;
    bool add_first_column_header = false;

    string header_info;
    string column_name_suffix = "RM";
    string suffix_to_remove;

    uint32_t min_reads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<xsym-files>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("header", 'H',"header to print (ignored for tabular)",
                      false, header_info);
    opt_parse.add_opt("tabular", 't', "output as table",
                      false, write_tabular_format);
    opt_parse.add_opt("radmeth", '\0', "make radmeth input "
                      "(use with tabular and without fractional)",
                      false, radmeth_format);
    opt_parse.add_opt("remove", '\0', "remove this suffix from filenames when "
                      "making column names; default removes from final dot",
                      false, suffix_to_remove);
    opt_parse.add_opt("suff", 's',
                      "column name suffixes, one for total reads and one for "
                      "methylated reads, to be separated from sample name "
                      "with underscore",
                      false, column_name_suffix);
    opt_parse.add_opt("fractional", 'f', "output fractions (requires tabular)",
                      false, write_fractional);
    opt_parse.add_opt("reads", 'r', "min reads (for fractional)",
                      false, min_reads);
    opt_parse.add_opt("ignore", '\0',
                      "do not check chrom order and assume lexicographic",
                      false, ignore_chroms_order);
    opt_parse.add_opt("mut", 'm',"If any of the sites being merged indicates "
                      "mutated, mark the result has mutated.",
                      false, report_any_mutated);
    opt_parse.add_opt("1st-col-header", '\0',"add name for 1st col in header",
                      false, add_first_column_header);
    opt_parse.add_opt("verbose", 'v',"print more run info", false, VERBOSE);
    opt_parse.set_show_defaults();
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
    if (write_fractional && !write_tabular_format) {
      cerr << "fractional output only available for tabular format" << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (size(column_name_suffix) != 2) {
      cerr << "column name suffix must be 2 letters" << endl;
      return EXIT_SUCCESS;
    }
    vector<string> xsym_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    const auto n_files = size(xsym_files);

    vector<vector<string>> headers;
    for (const auto &fn: xsym_files)
      headers.push_back(parse_xcounts_header(fn));

    for (const auto &h: headers)
      if (h != headers.front())
        throw runtime_error{"inconsistent xsym headers between input files"};

    unordered_map<string, uint32_t> chroms_order;
    if (!ignore_chroms_order)
      chroms_order = get_chroms_order(headers.front());

    if (VERBOSE && !ignore_chroms_order) {
      cerr << "chromosome order" << endl;
      vector<string> v(size(chroms_order));
      for (auto &&i : chroms_order) v[i.second] = i.first;
      copy(cbegin(v), cend(v), ostream_iterator<string>(cerr, "\n"));
    }

    vector<bgzf_file*> infiles(n_files);
    for (uint32_t i = 0; i < n_files; ++i) {
      infiles[i] = new bgzf_file(xsym_files[i], "r");
      if (!(*infiles[i]))
        throw runtime_error("cannot open file: " + xsym_files[i]);
    }

    for (auto &f: infiles)
      read_past_header(*f);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
    if (!out) throw runtime_error("error opening output file: " + outfile);

    // print header if user specifies or if tabular output format
    if (write_tabular_format) {

      vector<string> colnames;
      for (const auto &i : xsym_files)
        colnames.push_back(strip_path(i));

      for (auto &i : colnames)
        i = suffix_to_remove.empty() ?
          remove_extension(i) : remove_suffix(suffix_to_remove, i);

      if (!write_fractional && !radmeth_format) {
        vector<string> tmp;
        for (auto &&i : colnames) {
          tmp.push_back(i + "_" + column_name_suffix[0]);
          tmp.push_back(i + "_" + column_name_suffix[1]);
        }
        swap(tmp, colnames);
      }

      if (add_first_column_header)
        out << "site\t";
      copy(begin(colnames), end(colnames),
           std::ostream_iterator<string>(out, "\t"));
      out << endl;
    }
    else if (!header_info.empty())
      out << "#" << header_info << endl;

    vector<MSite> sites(n_files);
    vector<string> chrom_names(n_files);
    vector<uint64_t> offsets(n_files);
    vector<bool> outdated(n_files, true);
    vector<bool> sites_to_print; // declared here to keep allocation
    vector<unordered_set<string>> chroms_seen(n_files);

    while (any_sites_unprocessed(xsym_files, infiles, outdated, chrom_names,
                                 offsets, sites)) {
      sites_to_print.clear();
      sites_to_print.resize(n_files, false);

      // below idx is one index among the sites to print
      const uint32_t idx =
        collect_sites_to_print(sites, chroms_order, outdated, sites_to_print);

      // output the appropriate sites' data
      if (write_tabular_format)
        write_line_for_tabular(write_fractional, report_any_mutated, min_reads,
                               out, sites_to_print, sites, sites[idx]);
      else
        write_line_for_merged_counts(out, report_any_mutated,
                                     sites_to_print, sites, sites[idx]);

      swap(outdated, sites_to_print);
    }

    for (auto &&f : infiles)
      delete f;
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
