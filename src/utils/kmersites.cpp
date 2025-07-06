/* kmersites: a program to generate a wiggle format file (using the
 * UCSC Genome Browser wiggle format) to indicate the location of
 * sites matching a specific k-mer
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
#include <numeric>
#include <stdexcept>
#include <cstdint> // for [u]int[0-9]+_t
#include <filesystem>
#include <iterator>
#include <algorithm>

#include "OptionParser.hpp"
#include "dnmt_error.hpp"
#include "smithlab_os.hpp"

#include <bamxx.hpp>

namespace fs = std::filesystem;

using bamxx::bgzf_file;

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::to_string;

static inline auto
process_chrom_wig(const string &kmer, const int offset, const string &name,
                  const string &chrom, bgzf_file &out) -> void {

  static const auto variable_step_chrom_header = "variableStep chrom=";

  out.write(variable_step_chrom_header + name + "\n");

  const auto kmer_size = size(kmer);
  const auto chrom_size = size(chrom);
  if (kmer_size > chrom_size)
    throw dnmt_error("kmer size " + to_string(kmer_size) +
                     " larger than chrom size " + to_string(chrom_size));

  const auto beg_kmer = cbegin(kmer);
  const auto end_kmer = cend(kmer);

  const auto end_chrom = cend(chrom);
  auto chrom_itr = cbegin(chrom);
  auto chrom_itr_k = chrom_itr + kmer_size;

  auto pos = 0;
  while (chrom_itr_k != end_chrom) {
    if (std::equal(beg_kmer, end_kmer, chrom_itr++, chrom_itr_k++))
      out.write(to_string(pos + offset) + "\t1\n");
    ++pos;
  }
}

static inline auto
process_chrom_with_named_lines(const string &kmer, const int offset,
                               const string &name, const string &chrom,
                               bgzf_file &out) -> void {

  const auto kmer_size = size(kmer);
  const auto chrom_size = size(chrom);
  if (kmer_size > chrom_size)
    throw dnmt_error("kmer size " + to_string(kmer_size) +
                     " larger than chrom size " + to_string(chrom_size));

  const auto beg_kmer = cbegin(kmer);
  const auto end_kmer = cend(kmer);

  const auto end_chrom = cend(chrom);
  auto chrom_itr = cbegin(chrom);
  auto chrom_itr_k = chrom_itr + kmer_size;

  auto pos = 0;
  while (chrom_itr_k != end_chrom) {
    if (std::equal(beg_kmer, end_kmer, chrom_itr++, chrom_itr_k++))
      out.write(name + "\t" + to_string(pos + offset) + "\t1\n");
    ++pos;
  }
}

auto
kmersites(const int argc, char *argv[]) -> int {
  try {

    bool verbose = false;
    bool show_progress = false;
    bool compress_output = false;
    bool name_each_line = false;

    string kmer = "CG";
    string outfile;
    // int n_threads = 1;
    int offset = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(fs::path(string(*argv)).filename(),
                           "get sites matching kmer",
                           "<fasta-file>");
    // opt_parse.add_opt("threads", 't', "threads to use (few needed)",
    //                   false, n_threads);
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("offset", 'O', "offset within kmer to report",
                      false, offset);
    opt_parse.add_opt("kmer", 'k', "kmer to report", false, kmer);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("name-each-line", '\0', "name each line with chrom",
                      false, name_each_line);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.about_requested() || opt_parse.help_requested() ||
        leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    const string chroms_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (offset < 0)
      throw dnmt_error("offset must be non-negative (specified=" +
                       to_string(offset));

    // if (n_threads < 0)
    //   throw dnmt_error("thread count cannot be negative");

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    // file types from HTSlib use "-" for the filename to go to stdout
    if (outfile.empty()) outfile = "-";

    if (verbose)
      cerr << "[input fastq file: " << chroms_file << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: " << (compress_output ? "bgzf" : "text") << "]"
           << endl
           // << "[threads requested: " << n_threads << "]" << endl
           << "[k-mer to report: " << kmer << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl;

    vector<string> names, chroms;
    read_fasta_file_short_names(chroms_file, names, chroms);
    for (auto &chrom : chroms)
      std::transform(cbegin(chrom), cend(chrom), begin(chrom),
                     [](const char c) { return std::toupper(c); });

    // open the output file
    const auto output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out) throw dnmt_error("error opening output file: " + outfile);

    for (auto i = 0u; i < size(names); ++i) {
      if (show_progress) cerr << "processing: " << names[i] << endl;
      if (name_each_line)
        process_chrom_with_named_lines(kmer, offset, names[i], chroms[i], out);
      else
        process_chrom_wig(kmer, offset, names[i], chroms[i], out);
    }
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
