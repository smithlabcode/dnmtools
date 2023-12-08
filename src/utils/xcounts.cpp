/* xcounts: reformat counts so they only give the m and u counts in a
 * wig format
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
#include <charconv>
#include <system_error>

// from smithlab_cpp
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "MSite.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using std::string;
using std::to_chars;
using std::vector;

using bamxx::bgzf_file;

inline auto
getline(bgzf_file &file, kstring_t &line) -> bgzf_file & {
  if (file.f == nullptr) return file;
  const int x = bgzf_getline(file.f, '\n', &line);
  if (x == -1) {
    file.destroy();
    free(line.s);
    line = {0, 0, nullptr};
  }
  if (x < -1) {
    // ADS: this is an error condition and should be handled
    // differently from the EOF above.
    file.destroy();
    free(line.s);
    line = {0, 0, nullptr};
  }
  return file;
}


template<typename T>
static inline uint32_t
fill_output_buffer(const uint32_t offset, const MSite &s, T &buf) {
  auto buf_end = buf.data() + buf.size();
  auto res = to_chars(buf.data(), buf_end, s.pos - offset);
  *res.ptr++ = '\t';
  res = to_chars(res.ptr, buf_end, s.n_meth());
  *res.ptr++ = '\t';
  res = to_chars(res.ptr, buf_end, s.n_unmeth());
  *res.ptr++ = '\n';
  return std::distance(buf.data(), res.ptr);
}


int
main_xcounts(int argc, const char **argv) {
  try {
    size_t n_threads = 1;

    string outfile{"-"};
    const string description =
      "compress counts files by removing context information";

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

    const auto outfile_mode = in.compression() ? "w" : "wu";

    bgzf_file out(outfile, outfile_mode);
    if (!out) throw runtime_error("error opening output file: " + outfile);

    if (n_threads > 1) {
      // ADS: something breaks when we use the thread for the input
      tpool.set_io(in);
      tpool.set_io(out);
    }

    // use the kstring_t type to more directly use the BGZF file
    kstring_t line{0, 0, nullptr};
    const int ret = ks_resize(&line, 1024);
    if (ret) throw runtime_error("failed to acquire buffer");

    vector<char> buf(128);

    uint32_t offset = 0;
    string prev_chrom;
    bool status_ok = true;
    int64_t sz = 0;

    MSite site;
    while (status_ok && getline(in, line)) {
      if (line.s[0] == '#') {
        const auto header_line = string(line.s) + "\n";
        sz = size(header_line);
        status_ok = bgzf_write(out.f, header_line.data(), sz) == sz;
        continue;
      }
      status_ok = site.initialize(line.s, line.s + line.l);
      if (!status_ok) break;

      if (site.chrom != prev_chrom) {
        prev_chrom = site.chrom;
        offset = 0;

        site.chrom += '\n';
        sz = size(site.chrom);
        status_ok = bgzf_write(out.f, site.chrom.data(), sz) == sz;
      }
      if (site.n_reads > 0 || site.is_mutated()) {
        sz = fill_output_buffer(offset, site, buf);
        status_ok = bgzf_write(out.f, buf.data(), sz) == sz;
        offset = site.pos;
      }
    }
    if (!status_ok) {
      cerr << "failed converting "
           << filename << " to " << outfile << endl;
      return EXIT_FAILURE;
    }
  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
