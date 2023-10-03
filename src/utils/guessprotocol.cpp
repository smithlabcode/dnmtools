/* guessprotocol: a program for guessing whether a wgbs protocol is
 * original, pbat or random pbat
 *
 * Copyright (C) 2019-2023
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
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <filesystem>
#include <sstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::runtime_error;
using std::array;

static string
guess_protocol(const double fraction_t_rich) {
  if (fraction_t_rich >= 0.8) {
    return "original";
  }
  if (fraction_t_rich <= 0.2) {
    return "pbat";
  }
  if (fraction_t_rich >= 0.4 && fraction_t_rich <= 0.6) {
    return "random";
  }
  return "inconclusive";
}

struct guessprotocol_summary {
  string protocol;
  string layout;
  uint64_t n_t_rich_reads{};
  uint64_t n_reads{};
  double fraction_t_rich_reads{};

  void evaluate() {
    fraction_t_rich_reads = static_cast<double>(n_t_rich_reads)/n_reads;
    protocol = guess_protocol(fraction_t_rich_reads);
  }

  string tostring() const {
    std::ostringstream oss;
    oss << "protocol: " << protocol << '\n'
        << "fraction_t_rich_reads: " << fraction_t_rich_reads  << '\n'
        << "t_rich_reads: " << n_t_rich_reads << '\n'
        << "reads_examined: " << n_reads;
    return oss.str();
  }
};

// store each read from one end
struct FASTQRecord {
  string name;
  string seq;
};

// see if two reads from two ends match to each other (they should
// have the same name)
static bool
mates(const size_t to_ignore_at_end, // in case names have #0/1 name ends
      const FASTQRecord &a, const FASTQRecord &b) {
  assert(to_ignore_at_end < std::size(a.name));
  return equal(cbegin(a.name), cend(a.name) - to_ignore_at_end,
               cbegin(b.name));
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
static std::istream&
operator>>(std::istream& s, FASTQRecord &r) {
  constexpr auto n_error_codes = 5u;
  enum err_code { none, bad_name, bad_seq, bad_plus, bad_qual };
  static const array<runtime_error, n_error_codes> error_msg = {
    runtime_error(""),
    runtime_error("failed to parse fastq name line"),
    runtime_error("failed to parse fastq sequence line"),
    runtime_error("failed to parse fastq plus line"),
    runtime_error("failed to parse fastq qual line")
  };

  err_code ec = err_code::none;

  getline(s, r.name);

  if (r.name.empty() || r.name[0] != '@')
    ec = err_code::bad_name;

  r.name.resize(r.name.find_first_of(' '));
  const auto nm_sz = r.name.find_first_of(' ');
  copy(cbegin(r.name) + 1, cbegin(r.name) + nm_sz, begin(r.name));

  if (!getline(s, r.seq))
    ec = err_code::bad_seq;

  string tmp;
  if (!getline(s, tmp))
    ec = err_code::bad_plus;

  if (!getline(s, tmp))
    ec = err_code::bad_qual;

  if (ec != err_code::none)
    throw error_msg[ec];

  return s;
}

int
main_guessprotocol(int argc, const char **argv) {

  try {

    constexpr auto description = "guess bisulfite protocol for a library";

    string outfile;
    size_t reads_to_check = 1000000;
    size_t name_suffix_len = 0;

    namespace fs = std::filesystem;
    const string cmd_name = std::filesystem::path(argv[0]).filename();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(cmd_name, description,
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check",
                      false, reads_to_check);
    opt_parse.add_opt("ignore", 'i', "length of read name suffix "
                      "to ignore when matching", false, name_suffix_len);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested() || leftover_args.size() > 2) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    const vector<string> reads_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    guessprotocol_summary summary;
    if (reads_files.size() == 2) {
      summary.layout = "paired";

      // input: paired-end reads with end1 and end2
      std::ifstream in1(reads_files.front());
      if (!in1)
        throw runtime_error("cannot open input file: " + reads_files.front());

      std::ifstream in2(reads_files.back());
      if (!in2)
        throw runtime_error("cannot open input file: " + reads_files.back());

      FASTQRecord end_one, end_two;
      while (in1 >> end_one && in2 >> end_two &&
             summary.n_reads < reads_to_check) {
        summary.n_reads++;

        // two reads should be in paired-ends
        if (!mates(name_suffix_len, end_one, end_two))
          throw runtime_error("expected mates, got: " +
                              end_one.name + " and " + end_two.name);

        const double end_one_a =
          count(begin(end_one.seq), end(end_one.seq), 'A') +
          count(begin(end_one.seq), end(end_one.seq), 'C');
        const double end_one_t =
          count(begin(end_one.seq), end(end_one.seq), 'T') +
          count(begin(end_one.seq), end(end_one.seq), 'G');

        const double end_two_a =
          count(begin(end_two.seq), end(end_two.seq), 'A') +
          count(begin(end_two.seq), end(end_two.seq), 'C');
        const double end_two_t =
          count(begin(end_two.seq), end(end_two.seq), 'T') +
          count(begin(end_two.seq), end(end_two.seq), 'G');

        const double t_rich_count = (end_one_t + end_two_a);
        const double pbat_count = (end_one_a + end_two_t);

        summary.n_t_rich_reads += (t_rich_count > pbat_count);
      }
    }
    else {
      summary.layout = "single";

      // input: single-end reads
      std::ifstream in(reads_files.front());
      if (!in)
        throw runtime_error("cannot open input file: " + reads_files.front());

      FASTQRecord r;
      while (in >> r && summary.n_reads < reads_to_check) {
        summary.n_reads++;
        const double a = (count(begin(r.seq), end(r.seq), 'A') +
                          count(begin(r.seq), end(r.seq), 'C'));
        const double t = (count(begin(r.seq), end(r.seq), 'T') +
                          count(begin(r.seq), end(r.seq), 'G'));
        summary.n_t_rich_reads += (t > a);
      }
    }

    summary.evaluate();
    if (!outfile.empty()) {
      std::ofstream out(outfile);
      if (!out) throw runtime_error("failed to open output file: " + outfile);
      out << summary.tostring() << endl;
    }
    else cout << summary.tostring() << endl;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
