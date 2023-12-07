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

#include <bamxx.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

// from smithlab_cpp
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;
using std::unordered_set;

using bamxx::bgzf_file;

static inline bool
found_symmetric(const MSite &prev, const MSite &curr) {
  // assumes check for CpG already done
  return (prev.pos + 1 == curr.pos && prev.strand == '+' && curr.strand == '-');
}

static inline void
ensure_positive_cpg(MSite &s) {
  s.pos -= (s.strand == '-');
  s.strand = '+';
}

template<class T> static std::tuple<MSite, bool>
get_first_site(T &in, T &out) {
  bool prev_is_cpg = false;
  MSite prev_site;
  string line;
  bool within_header = true;
  while (within_header && getline(in, line)) {
    if (line[0] == '#') {
      line += '\n';
      out.write(line);
    }
    else {
      prev_site.initialize(line.data(), line.data() + size(line));
      if (prev_site.is_cpg()) prev_is_cpg = true;
      within_header = false;
    }
  }
  return {prev_site, prev_is_cpg};
}

template<class T> static bool
process_sites(T &in, T &out) {

  // get the first site while dealing with the header
  auto [prev_site, prev_is_cpg] = get_first_site(in, out);

  // setup to verfiy that chromosomes are together
  unordered_set<string> chroms_seen;
  chroms_seen.insert(prev_site.chrom);
  bool sites_are_sorted = true;

  MSite curr_site;
  while (read_site(in, curr_site)) {
    const bool same_chrom = prev_site.chrom == curr_site.chrom;
    if (same_chrom) {
      if (curr_site.pos <= prev_site.pos) return false;
      if (prev_is_cpg) {
        if (curr_site.is_mate_of(prev_site))
          curr_site.add(prev_site);
        else {
          ensure_positive_cpg(prev_site);
          write_site(out, prev_site);
        }
      }
    }
    else {
      if (chroms_seen.find(curr_site.chrom) != cend(chroms_seen)) return false;
      chroms_seen.insert(curr_site.chrom);

      if (prev_is_cpg) {
        ensure_positive_cpg(prev_site);
        write_site(out, prev_site);
      }
    }
    std::swap(prev_site, curr_site);
    prev_is_cpg = prev_site.is_cpg();
  }
  if (prev_is_cpg) {
    ensure_positive_cpg(prev_site);
    write_site(out, prev_site);
  }

  return sites_are_sorted;
}

int
main_symmetric_cpgs(int argc, const char **argv) {
  try {
    // file types from HTSlib use "-" for the filename to go to stdout
    string outfile{"-"};
    // (not used) bool VERBOSE = false;
    bool compress_output = false;
    int32_t n_threads = 1;

    const string description =
      "Get CpG sites and make methylation levels symmetric.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
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

    if (n_threads <= 0) throw runtime_error("threads must be positive");
    bamxx::bam_tpool tp(n_threads);

    // const bool show_progress = VERBOSE && isatty(fileno(stderr));
    bgzf_file in(filename, "r");
    if (!in) throw runtime_error("could not open file: " + filename);

    // open the output file
    const string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out) throw runtime_error("error opening output file: " + outfile);

    if (n_threads > 1) {
      tp.set_io(in);
      tp.set_io(out);
    }

    const bool sites_are_sorted = process_sites(in, out);

    if (!sites_are_sorted) {
      namespace fs = std::filesystem;
      cerr << "sites are not sorted in: " << filename << endl;
      const fs::path outpath{outfile};
      if (fs::exists(outpath)) fs::remove(outpath);
      return EXIT_FAILURE;
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
