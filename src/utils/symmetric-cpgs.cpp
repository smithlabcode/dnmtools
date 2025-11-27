/* symmetric-cpgs: extract the CpG sites from a methcounts output
 * file and produce a new one with the CpGs treated unstranded.
 *
 * Copyright (C) 2014-2025 Andrew D. Smith
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

#include "MSite.hpp"
#include "OptionParser.hpp"
#include "counts_header.hpp"

#include <bamxx.hpp>

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

static inline void
ensure_positive_cpg(MSite &s) {
  s.pos -= (s.strand == '-');
  s.strand = '+';
}

template <class T>
static std::tuple<MSite, bool>
get_first_site(T &in, T &out) {
  bool prev_is_cpg = false;
  MSite prev_site;
  std::string line;
  bool within_header = true;
  while (within_header && getline(in, line)) {
    if (is_counts_header_line(line)) {
      write_counts_header_line(line, out);
    }
    else {
      // NOLINTNEXTLINE(*-pointer-arithmetic)
      if (!prev_site.initialize(line.data(), line.data() + std::size(line)))
        throw std::runtime_error("failed to parse line:\n" + line);
      if (prev_site.is_cpg())
        prev_is_cpg = true;
      within_header = false;
    }
  }
  return {prev_site, prev_is_cpg};
}

template <class T>
static bool
process_sites(const bool verbose, T &in, T &out) {
  // get the first site while dealing with the header
  auto [prev_site, prev_is_cpg] = get_first_site(in, out);

  // setup to verfiy that chromosomes are together
  std::unordered_set<std::string> chroms_seen;
  chroms_seen.insert(prev_site.chrom);
  bool sites_are_sorted = true;

  MSite curr_site;
  while (read_site(in, curr_site)) {
    const bool same_chrom = prev_site.chrom == curr_site.chrom;
    if (same_chrom) {
      if (curr_site.pos <= prev_site.pos)
        return false;
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
      if (verbose)
        std::cerr << "processing: " << curr_site.chrom << '\n';
      if (chroms_seen.find(curr_site.chrom) != cend(chroms_seen))
        return false;
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
main_symmetric_cpgs(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    // file types from HTSlib use "-" for the filename to go to stdout
    std::string outfile{"-"};
    bool verbose = false;
    bool compress_output = false;
    bool allow_extra_fields = false;
    int32_t n_threads = 1;

    const std::string description =
      "Get CpG sites and make methylation levels symmetric.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("relaxed", '\0', "allow extra fields in input", false,
                      allow_extra_fields);
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
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    MSite::no_extra_fields = (allow_extra_fields == false);

    if (n_threads <= 0)
      throw std::runtime_error("threads must be positive");
    bamxx::bam_tpool tp(n_threads);

    // const bool show_progress = VERBOSE && isatty(fileno(stderr));
    bamxx::bgzf_file in(filename, "r");
    if (!in)
      throw std::runtime_error("could not open file: " + filename);

    // open the output file
    const std::string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw std::runtime_error("error opening output file: " + outfile);

    if (n_threads > 1) {
      if (in.is_bgzf())
        tp.set_io(in);
      tp.set_io(out);
    }

    const bool sites_are_sorted = process_sites(verbose, in, out);

    if (!sites_are_sorted) {
      std::cerr << "sites are not sorted in: " << filename << '\n';
      const std::filesystem::path outpath{outfile};
      if (std::filesystem::exists(outpath) && !std::filesystem::remove(outpath))
        throw std::runtime_error("failed to remove file: " + outpath.string());
      return EXIT_FAILURE;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
