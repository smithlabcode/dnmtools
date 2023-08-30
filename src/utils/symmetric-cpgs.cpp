/* symmetric-cpgs: extract the CpG sites from a methcounts output
 * file and produce a new one with the CpGs treated unstranded.
 *
 * Copyright (C) 2014 University of Southern California and
 *                    Andrew D. Smith
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
#include <fstream>
#include <stdexcept>
#include <unordered_set>
#include <filesystem>
#include <bamxx.hpp>

// from smithlab_cpp
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include "MSite.hpp"

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_set;

using bamxx::bgzf_file;

static inline bool
found_symmetric(const MSite &prev_cpg, const MSite &curr_cpg) {
  // assumes check for CpG already done
  return (prev_cpg.strand == '+' && curr_cpg.strand == '-' &&
          prev_cpg.pos + 1 == curr_cpg.pos);
}

template<class T> static bool
process_sites(bgzf_file &in, T &out) {
  MSite prev_site, curr_site;
  bool prev_is_cpg = false;
  if (read_site(in, prev_site))
    if (prev_site.is_cpg()) prev_is_cpg = true;

  unordered_set<string> chroms_seen;
  bool sites_are_sorted = true;

  while (read_site(in, curr_site)) {
    if (curr_site.is_cpg()) {
      if (curr_site.chrom != prev_site.chrom) {
        const auto chrom_itr = chroms_seen.find(curr_site.chrom);
        if (chrom_itr != cend(chroms_seen)) {
          sites_are_sorted = false;
          break;
        }
        else
          chroms_seen.insert(curr_site.chrom);
      }
      else if (curr_site.pos <= prev_site.pos) {
        sites_are_sorted = false;
        break;
      }
      if (prev_is_cpg && found_symmetric(prev_site, curr_site)) {
        prev_site.add(curr_site);
        write_site(out, prev_site);
      }
      prev_is_cpg = true;
    }
    else
      prev_is_cpg = false;
    std::swap(prev_site, curr_site);
  }

  return sites_are_sorted;
}

int
main_symmetric_cpgs(int argc, const char **argv) {
  try {
    string outfile;
    // (not used) bool VERBOSE = false;

    const string description =
      "Get CpG sites and make methylation levels symmetric.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    // opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<string> leftover_args;
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    bgzf_file in(filename, "r");
    if (!in) throw std::runtime_error("could not open file: " + filename);

    bool sites_are_sorted = true;
    if (outfile.empty() || !has_gz_ext(outfile)) {
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      sites_are_sorted = process_sites(in, out);
    }
    else {
      bgzf_file out(outfile, "w");
      sites_are_sorted = process_sites(in, out);
    }

    if (!sites_are_sorted) {
      cerr << "sites are not sorted in: " << filename << endl;
      namespace fs = std::filesystem;
      const fs::path outpath{outfile};
      if (fs::exists(outpath)) fs::remove(outpath);
      return EXIT_FAILURE;
    }
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
