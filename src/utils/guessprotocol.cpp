/* guessprotocol: a program for guessing whether a wgbs protocol is
 * wgbs, pbat or random pbat
 *
 * Copyright (C) 2019-2023 Andrew D. Smith
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

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <array>

#include "OptionParser.hpp"
#include "numerical_utils.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::array;
using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::runtime_error;
using std::string;
using std::vector;

using bamxx::bgzf_file;

constexpr int nuc_to_idx[] = {
 /*  0*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 16*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 32*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 48*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 64*/  4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 80*/  4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /* 96*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*112*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 /*128*/  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

struct nucleotide_model {
  vector<double> pr{};
  vector<double> lpr{};
  double bisulfite_conversion_rate{};
  bool is_t_rich{};

  nucleotide_model(const vector<double> &bc, const double conv_rate,
                   const bool itr)
      : pr{bc}, bisulfite_conversion_rate{conv_rate}, is_t_rich{itr} {
    auto nuc_from = is_t_rich ? 1 : 2;
    auto nuc_to = is_t_rich ? 3 : 0;
    pr[nuc_to] += bisulfite_conversion_rate * pr[nuc_from];
    pr[nuc_from] *= (1.0 - bisulfite_conversion_rate);
    assert(reduce(cbegin(pr), cend(pr), 0.0) == 1.0);
    lpr.resize(std::size(pr));
    transform(cbegin(pr), cend(pr), begin(lpr),
              [](const double x) { return log(x); });
  }

  double operator()(const string &s) const {
    return accumulate(cbegin(s), cend(s), 0.0,
                      [&](const double x, const char c) {
                        const auto i = nuc_to_idx[static_cast<uint8_t>(c)];
                        return i == 4 ? x : x + lpr[i];
                      });
  };

  string tostring() const {
    std::ostringstream oss;
    oss << "pr:\n";
    for (auto i : pr) oss << i << '\n';
    oss << "log pr:\n";
    for (auto i : lpr) oss << i << '\n';
    oss << bisulfite_conversion_rate << '\n' << is_t_rich;
    return oss.str();
  }
};

struct guessprotocol_summary {
  // protocol is the guessed protocol (wgbs, pbat, rpbat, or inconclusive)
  // based on the content of the reads.
  string protocol;
  // confidence indicates the level of confidence in the guess for the 
  // protocol.
  string confidence;
  // layout indicates whether the reads are paired or single-ended. 
  string layout;
  // n_reads_wgbs is the average number of reads (for single-ended reads) or
  // read pairs (for paired reads) where read1 is T-rich.
  double n_reads_wgbs{};
  // n_reads is the number of evaluated reads or read pairs.
  uint64_t n_reads{};
  // wgbs_fraction is the probability that a read (for single-ended reads) or
  // the read1 of a read pair (for paired reads) is T-rich.
  double wgbs_fraction{};

  void evaluate() {
    const auto frac = n_reads_wgbs / std::max(1ul, n_reads);
    if (frac <= 0.01)      {protocol = "pbat"; confidence = "high";}
    else if (frac <= 0.1)  {protocol = "pbat"; confidence = "low";}
    else if (frac <= 0.2)  {protocol = "rpbat"; confidence = "low";}
    else if (frac <= 0.8)  {protocol = "rpbat"; confidence = "high";}
    else if (frac <= 0.9)  {protocol = "rpbat"; confidence = "low";}
    else if (frac <= 0.99) {protocol = "wgbs"; confidence = "low";}
    else                  {protocol = "wgbs"; confidence = "high";}
    wgbs_fraction = frac;
  }

  string tostring() const {
    std::ostringstream oss;
    oss << "protocol: " << protocol << '\n'
        << "confidence: " << confidence << '\n'
        << "wgbs_fraction: " << wgbs_fraction  << '\n'
        << "n_reads_wgbs: " << n_reads_wgbs << '\n'
        << "n_reads: " << n_reads;
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
static bgzf_file &
operator>>(bgzf_file &s, FASTQRecord &r) {
  constexpr auto n_error_codes = 5u;

  enum err_code { none, bad_name, bad_seq, bad_plus, bad_qual };

  static const array<runtime_error, n_error_codes> error_msg = {
    runtime_error(""), runtime_error("failed to parse fastq name line"),
    runtime_error("failed to parse fastq sequence line"),
    runtime_error("failed to parse fastq plus line"),
    runtime_error("failed to parse fastq qual line")
  };

  err_code ec = err_code::none;

  if (!getline(s, r.name)) return s;

  if (r.name.empty() || r.name[0] != '@') ec = err_code::bad_name;

  const auto nm_end = r.name.find_first_of(" \t");
  const auto nm_sz = (nm_end == string::npos ? r.name.size() : nm_end) - 1;
  r.name.erase(copy_n(cbegin(r.name) + 1, nm_sz, begin(r.name)), cend(r.name));

  if (!getline(s, r.seq)) ec = err_code::bad_seq;

  string tmp;
  if (!getline(s, tmp)) ec = err_code::bad_plus;

  if (!getline(s, tmp)) ec = err_code::bad_qual;

  if (ec != err_code::none) throw error_msg[ec];

  return s;
}

int
main_guessprotocol(int argc, const char **argv) {

  try {

    static const vector<double> human_base_comp = {0.295, 0.205, 0.205, 0.295};
    static const vector<double> flat_base_comp = {0.25, 0.25, 0.25, 0.25};

    constexpr auto description = "guess bisulfite protocol for a library";

    bool verbose;
    bool use_human;
    string outfile;
    size_t reads_to_check = 1000000;
    size_t name_suffix_len = 0;
    double bisulfite_conversion_rate = 0.98;

    namespace fs = std::filesystem;
    const string cmd_name = std::filesystem::path(argv[0]).filename();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(cmd_name, description,
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check",
                      false, reads_to_check);
    opt_parse.add_opt("ignore", 'i', "length of read name suffix "
                      "to ignore when matching", false, name_suffix_len);
    opt_parse.add_opt("bisulfite", 'b', "bisulfite conversion rate",
                      false, bisulfite_conversion_rate);
    opt_parse.add_opt("human", 'H', "assume human genome",
                      false, use_human);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("verbose", 'v',
                      "report available information during the run",
                      false, verbose);
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

    auto base_comp = flat_base_comp;
    if (use_human) base_comp = human_base_comp;

    nucleotide_model t_rich_model(base_comp, bisulfite_conversion_rate, true);
    nucleotide_model a_rich_model(base_comp, bisulfite_conversion_rate, false);

    guessprotocol_summary summary;
    summary.layout = reads_files.size() == 2 ? "paired" : "single";

    if (verbose) {
      if (reads_files.size() == 2)
        cerr << "data layout: "
             << "paired" << '\n'
             << "read1 file: " << reads_files.front() << '\n'
             << "read2 file: " << reads_files.back() << '\n';
      else
        cerr << "data layout: "
             << "single" << '\n'
             << "read file: " << reads_files.front() << '\n';
      cerr << "reads to check: " << reads_to_check << '\n'
           << "read name suffix length: " << name_suffix_len << '\n'
           << "bisulfite conversion: " << bisulfite_conversion_rate << '\n';
    }

    if (reads_files.size() == 2) {

      // input: paired-end reads with end1 and end2
      bgzf_file in1(reads_files.front(), "r");
      if (!in1)
        throw runtime_error("cannot open file: " + reads_files.front());

      bgzf_file in2(reads_files.back(), "r");
      if (!in2)
        throw runtime_error("cannot open file: " + reads_files.back());

      FASTQRecord r1, r2;
      while (in1 >> r1 && in2 >> r2 && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        if (!mates(name_suffix_len, r1, r2))
          throw runtime_error("expected mates: " + r1.name + ", " + r2.name);

        const double ta = t_rich_model(r1.seq) + a_rich_model(r2.seq);
        const double at = a_rich_model(r1.seq) + t_rich_model(r2.seq);

        const auto prob_read1_t_rich = exp(ta - log_sum_log(ta, at));
        summary.n_reads_wgbs += prob_read1_t_rich;
      }
    }
    else {

      // input: single-end reads
      bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw runtime_error("cannot open file: " + reads_files.front());

      FASTQRecord r;
      while (in >> r && summary.n_reads < reads_to_check) {
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;
      }
    }

    summary.evaluate();

    if (!outfile.empty()) {
      std::ofstream out(outfile);
      if (!out) throw runtime_error("failed to open: " + outfile);
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
