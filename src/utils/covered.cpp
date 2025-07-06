/* covered: filter counts files so they only have covered sites.
 *
 * Copyright (C) 2023 Andrew D. Smith
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

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// from smithlab_cpp
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;

using bamxx::bgzf_file;


static inline bool
get_is_mutated(const kstring_t &line) {
  const auto end_itr = line.s + line.l;
  return std::find(line.s, end_itr, 'x') != end_itr;
}

static inline uint32_t
get_n_reads(const kstring_t &line) {
  const auto end_itr = std::make_reverse_iterator(line.s + line.l);
  const auto beg_itr = std::make_reverse_iterator(line.s);
  auto n_reads_pos = std::find_if(
    end_itr, beg_itr, [](const char c) { return c == ' ' || c == '\t'; });
  ++n_reads_pos;
  return atoi(n_reads_pos.base());
}

int
main_covered(int argc, char *argv[]) {
  try {
    size_t n_threads = 1;

    string outfile{"-"};
    const string description =
      "filter counts files so they only have covered sites";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<counts-file> (\"-\" for standard input)", 1);
    opt_parse.add_opt("output", 'o', "output file (default is standard out)",
                      false, outfile);
    opt_parse.add_opt("threads", 't', "threads for compression (use few)",
                      false, n_threads);
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

    bamxx::bam_tpool tpool(n_threads);

    bgzf_file in(filename, "r");
    if (!in) throw runtime_error("could not open file: " + filename);

    const auto outfile_mode = in.is_compressed() ? "w" : "wu";

    bgzf_file out(outfile, outfile_mode);
    if (!out) throw runtime_error("error opening output file: " + outfile);

    if (n_threads > 1) {
      // ADS: something breaks when we use the thread for the input
      if (in.is_bgzf())
        tpool.set_io(in);
      tpool.set_io(out);
    }

    // use the kstring_t type to more directly use the BGZF file
    kstring_t line{0, 0, nullptr};
    const int ret = ks_resize(&line, 1024);
    if (ret) throw runtime_error("failed to acquire buffer");

    bool write_ok = true;
    while (bamxx::getline(in, line) && write_ok) {
      const bool is_mutated = get_is_mutated(line);
      const uint32_t n_reads = get_n_reads(line);
      if (n_reads > 0u || is_mutated) {
        line.s[line.l++] = '\n';
        write_ok =
          (bgzf_write(out.f, line.s, line.l) == static_cast<int64_t>(line.l));
      }
    }
    if (!write_ok) {
      cerr << "failed writing to: " << outfile << '\n';
      return EXIT_FAILURE;
    }
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
