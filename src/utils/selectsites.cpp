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

#include <unistd.h>
#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <unordered_map>

#include <bamxx.hpp>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MSite.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ios_base;
using std::runtime_error;
using std::ifstream;
using std::unordered_map;

using bamxx::bgzf_file;

static void
collapsebed(vector<GenomicRegion> &regions) {
  size_t j = 0;
  for (size_t i = 1; i < regions.size(); ++i) {
    if (regions[j].same_chrom(regions[i]) &&
        regions[i].get_start() <= regions[j].get_end()) {
      regions[j].set_end(std::max(regions[j].get_end(),
                                  regions[i].get_end()));
    }
    else {
      regions[++j] = regions[i];
    }
  }
  regions.erase(begin(regions) + j + 1, end(regions));
}

static bool
precedes(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() < s.chrom ||
          (r.get_chrom() == s.chrom && r.get_end() <= s.pos));
}


static bool
contains(const GenomicRegion &r, const MSite &s) {
  return (r.get_chrom() == s.chrom &&
          (r.get_start() <= s.pos && s.pos < r.get_end()));
}


template <class T>
static void
process_all_sites(const bool VERBOSE,
                  const string &sites_file,
                  const unordered_map<string, vector<GenomicRegion>> &regions,
                  T &out) {

  bgzf_file in(sites_file, "r");
  if (!in) throw runtime_error("cannot open file: " + sites_file);

  MSite the_site, prev_site;
  vector<GenomicRegion>::const_iterator i, i_lim;
  bool chrom_is_relevant = false;
  while (read_site(in, the_site)) {
    if (the_site.chrom != prev_site.chrom) {
      if (VERBOSE)
        cerr << "processing " << the_site.chrom << endl;
      auto r = regions.find(the_site.chrom);
      chrom_is_relevant = (r != end(regions));
      if (chrom_is_relevant) {
        i = begin(r->second);
        i_lim = end(r->second);
      }
    }
    if (chrom_is_relevant) {
      while (i != i_lim && precedes(*i, the_site))
        ++i;

      if (contains(*i, the_site))
        write_site(out, the_site);
    }
    std::swap(prev_site, the_site);
  }
}


static void
get_sites_in_region(ifstream &site_in, const GenomicRegion &region,
                    std::ostream &out) {

  string chrom(region.get_chrom());
  const size_t start_pos = region.get_start();
  const size_t end_pos = region.get_end();
  find_offset_for_msite(chrom, start_pos, site_in);

  MSite the_site;
  // ADS: should only happen once that "the_site.chrom < chrom" and
  // this is only needed because of end state of binary search on disk
  while (site_in >> the_site &&
         (the_site.chrom < chrom ||
          (the_site.chrom == chrom && the_site.pos < end_pos)))
    if (start_pos <= the_site.pos)
      out << the_site << endl;
}


static void
process_with_sites_on_disk(const string &sites_file,
                           vector<GenomicRegion> &regions,
                           std::ostream &out) {

  ifstream in(sites_file);
  if (!in)
    throw runtime_error("cannot open file: " + sites_file);

  for (size_t i = 0; i < regions.size() && in; ++i)
    get_sites_in_region(in, regions[i], out);
}


static void
regions_by_chrom(vector<GenomicRegion> &regions,
                 unordered_map<string, vector<GenomicRegion> > &lookup) {
  for (auto &&r: regions) {
    const string chrom_name(r.get_chrom());
    if (lookup.find(chrom_name) == end(lookup))
      lookup[chrom_name] = vector<GenomicRegion>();
    lookup[chrom_name].push_back(r);
  }
  regions.clear();
  regions.shrink_to_fit();
}

inline bool
file_exists(const string &filename) {
  return (access(filename.c_str(), F_OK) == 0);
}

int
main_selectsites(int argc, const char **argv) {

  try {

    bool VERBOSE = false;
    bool LOAD_ENTIRE_FILE = false;

    string outfile;

    const string description =
      "Select sites inside a set of genomic intervals. "
      "Sites must be specified in methcounts format. "
      "Intervals must be specified in bed format.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<intervals-bed> <sites-file>", 2);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("preload", 'p',
                      "preload sites (use for large target intervals)",
                      false, LOAD_ENTIRE_FILE);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string regions_file = leftover_args.front();
    const string sites_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (isdir(sites_file.c_str()) || !file_exists(sites_file))
      throw runtime_error("bad input sites file: " + sites_file);

    vector<GenomicRegion> regions;
    ReadBEDFile(regions_file, regions);
    if (!check_sorted(regions))
      throw runtime_error("regions not sorted in file: " + regions_file);

    const size_t n_orig_regions = regions.size();
    collapsebed(regions);
    if (VERBOSE && n_orig_regions != regions.size())
      cerr << "[number of regions merged due to overlap: "
           << n_orig_regions - regions.size() << "]" << endl;

    unordered_map<string, vector<GenomicRegion>> regions_lookup;
    if ((outfile.empty() || !has_gz_ext(outfile)) && LOAD_ENTIRE_FILE)
      regions_by_chrom(regions, regions_lookup);

    if (outfile.empty() || !has_gz_ext(outfile)) {
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile);
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
      if (!outfile.empty() && !out)
        throw runtime_error("failed to open output file: " + outfile);

      if (LOAD_ENTIRE_FILE)
        process_all_sites(VERBOSE, sites_file, regions_lookup, out);
      else
        process_with_sites_on_disk(sites_file, regions, out);
    }
    else {
      // not supporting search on disk for gz file
      bgzf_file out(outfile, "w");
      process_all_sites(VERBOSE, sites_file, regions_lookup, out);
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
