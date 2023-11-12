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

#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>
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

using bamxx::bgzf_file;

int
main_covered(int argc, const char **argv) {
  try {

    size_t n_threads = 1;

    string outfile;
    const string description =
      "filter counts files so they only have covered sites";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methcounts-file>");
    opt_parse.add_opt("output", 'o', "output file (required)", true, outfile);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
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

    // const bool show_progress = VERBOSE && isatty(fileno(stderr));
    bgzf_file in(filename, "r");
    if (!in) throw std::runtime_error("could not open file: " + filename);

    // open the output file
    bamxx::bgzf_file out(outfile, "w");
    if (!out) throw std::runtime_error("error opening output file: " + outfile);

    if (n_threads > 1) {
      tpool.set_io(in);
      tpool.set_io(out);
    }

    MSite site;
    while (read_site(in, site))
      if (site.n_reads > 0)
        write_site(out, site);
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
