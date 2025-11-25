/*    apxlc: a program for approximately and very quickly counting lines
 *
 *    Copyright (C) 2012-2022 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include "OptionParser.hpp"
#include "smithlab_os.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>  // IWYU pragma: keep
#include <stdexcept>
#include <string>
#include <vector>

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

static std::size_t
get_approx_line_count(const bool VERBOSE, const std::string &filename,
                      const std::size_t n_samples, std::size_t sample_size) {
  static const std::size_t megabyte = (1ul << 20);
  static const std::size_t kilobyte = (1ul << 10);

  const std::size_t filesize = get_filesize(filename);

  if (sample_size == 0)
    sample_size = std::min(megabyte / 10, filesize / n_samples);

  const std::size_t increment =
    std::floor((filesize - sample_size * n_samples) / (n_samples - 1.0)) +
    sample_size;

  assert(filesize > n_samples && filesize > sample_size &&
         filesize > n_samples * sample_size);

  if (VERBOSE) {
    std::cerr << "[PROCESSING FILE: " << filename << "]" << '\n'
              << "[FILESIZE: " << static_cast<double>(filesize) / megabyte
              << "MB]" << '\n'
              << "[CHUNK SIZE: "
              << static_cast<std::size_t>(1.0 * sample_size / kilobyte) << "KB]"
              << '\n'
              << "[NUM CHUNKS: " << n_samples << "]" << '\n'
              << "[TOTAL SAMPLE: " << (1.0 * n_samples * sample_size) / megabyte
              << "MB]" << '\n';
  }
  std::ifstream in(filename, std::ios_base::binary);
  if (!in)
    throw std::runtime_error("cannot open input file " + filename);

  std::vector<char> buffer(sample_size);
  double total_lines = 0.0;
  for (std::size_t i = 0; i < filesize && in.good(); i += increment) {
    in.seekg(i, std::ios_base::beg);
    in.read(&buffer.front(), sample_size);
    if (in.good())
      total_lines += (0.5 + count(buffer.begin(), buffer.end(), '\n'));
  }
  return (filesize * total_lines) / (n_samples * sample_size);
}

int
main_lc_approx(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    std::size_t n_samples = 100;
    std::size_t sample_size = 0;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "approximate line counting in large files",
                           "<file1> <file2> ...");
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("samples", 'n', "number of samples", false, n_samples);
    opt_parse.add_opt("size", 'z', "sample size (bytes)", false, sample_size);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
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
    if (leftover_args.size() < 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_FAILURE;
    }
    std::vector<std::string> filenames(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    for (auto i = 0u; i < std::size(filenames); ++i)
      std::cout << filenames[i] << "\t"
                << get_approx_line_count(VERBOSE, filenames[i], n_samples,
                                         sample_size)
                << '\n';
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-avoid-magic-numbers,*-narrowing-conversions)
