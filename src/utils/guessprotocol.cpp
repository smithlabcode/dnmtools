/* guessprotocol: a program for guessing whether a wgbs protocol is
 * wgbs, pbat, rpbat, or rrbs
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
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "OptionParser.hpp"
#include "numerical_utils.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::array;
using std::cbegin;
using std::cend;
using std::cerr;
using std::cout;
using std::endl;
using std::greater;
using std::min;
using std::pair;
using std::priority_queue;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::vector;

using bamxx::bgzf_file;

constexpr int nuc_to_idx[] = {
    // clang-format off
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
    // clang-format on
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
    // adjust probabilities after the conversion
    pr[nuc_to] += bisulfite_conversion_rate * pr[nuc_from];
    pr[nuc_from] *= (1.0 - bisulfite_conversion_rate);
    assert(reduce(cbegin(pr), cend(pr), 0.0) == 1.0);
    lpr.resize(std::size(pr));
    transform(cbegin(pr), cend(pr), begin(lpr),
              [](const double x) { return log(x); });
  }

  // calculate sequence log likelihood
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
    for (auto i : pr)
      oss << i << '\n';
    oss << "log pr:\n";
    for (auto i : lpr)
      oss << i << '\n';
    oss << bisulfite_conversion_rate << '\n' << is_t_rich;
    return oss.str();
  }
};

struct guessprotocol_summary {

  static constexpr auto wgbs_cutoff_confident = 0.99;
  static constexpr auto wgbs_cutoff_unconfident = 0.9;
  static constexpr auto rpbat_cutoff_confident_high = 0.8;
  static constexpr auto rpbat_cutoff_confident_low = 0.2;
  static constexpr auto pbat_cutoff_unconfident = 0.1;
  static constexpr auto pbat_cutoff_confident = 0.01;
  static constexpr auto rrbs_cutoff_confident = 0.9;
  static constexpr auto rrbs_cutoff_unconfident = 0.8;

  // protocol is the guessed protocol (wgbs, pbat, rpbat, rrbs, or inconclusive)
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
  // read_len is the length of each read (assuming equal read length)
  uint64_t read_len{};
  // wgbs_fraction is the probability that a read (for single-ended reads) or
  // the read1 of a read pair (for paired reads) is T-rich.
  double wgbs_fraction{};
  // rrbs_fraction is the fraction of reads containing the most enriched k-mer
  // at the 5' end
  double rrbs_fraction{};
  // bool value indicating whether RRBS or not
  bool is_rrbs = false;
  // rrbs_confidence is the confidence in the guess of rrbs
  string rrbs_confidence = "NA";

  void evaluate() {
    const auto frac = n_reads_wgbs / n_reads;
    protocol = "inconclusive";

    // assigning wgbs (near one)
    if (frac > wgbs_cutoff_confident) {
      protocol = "wgbs";
      confidence = "high";
    } else if (frac > wgbs_cutoff_unconfident) {
      protocol = "wgbs";
      confidence = "low";
    }
    // assigning pbat (near zero)
    else if (frac < pbat_cutoff_confident) {
      protocol = "pbat";
      confidence = "high";
    } else if (frac < pbat_cutoff_unconfident) {
      protocol = "pbat";
      confidence = "low";
    }
    // assigning rpbat (towards middle)
    else if (frac > rpbat_cutoff_confident_low &&
             frac < rpbat_cutoff_confident_high) {
      protocol = "rpbat";
      confidence = "high";
    } else {
      protocol = "rpbat";
      confidence = "low";
    }

    wgbs_fraction = frac;

    if (rrbs_fraction > rrbs_cutoff_confident) {
      is_rrbs = true;
      rrbs_confidence = "high";
    } else if (rrbs_fraction > rrbs_cutoff_unconfident) {
      is_rrbs = true;
      rrbs_confidence = "low";
    }
  }

  string tostring() const {
    std::ostringstream oss;
    oss << "protocol: " << protocol << '\n'
        << "confidence: " << confidence << '\n'
        << "wgbs_fraction: " << wgbs_fraction << '\n'
        << "n_reads_wgbs: " << n_reads_wgbs << '\n'
        << "n_reads: " << n_reads << '\n'
        << "is_rrbs: " << std::boolalpha << is_rrbs << '\n'
        << "rrbs_fraction: " << rrbs_fraction << '\n'
        << "rrbs_confidence: " << rrbs_confidence << '\n';
    return oss.str();
  }
};

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
  return equal(cbegin(a.name), cend(a.name) - to_ignore_at_end, cbegin(b.name));
}

// Read 4 lines one time from fastq and fill in the FASTQRecord structure
static bgzf_file &operator>>(bgzf_file &s, FASTQRecord &r) {
  constexpr auto n_error_codes = 5u;

  enum err_code { none, bad_name, bad_seq, bad_plus, bad_qual };

  static const array<runtime_error, n_error_codes> error_msg = {
      runtime_error(""), runtime_error("failed to parse fastq name line"),
      runtime_error("failed to parse fastq sequence line"),
      runtime_error("failed to parse fastq plus line"),
      runtime_error("failed to parse fastq qual line")};

  err_code ec = err_code::none;

  if (!getline(s, r.name))
    return s;

  if (r.name.empty() || r.name[0] != '@')
    ec = err_code::bad_name;

  const auto nm_end = r.name.find_first_of(" \t");
  const auto nm_sz = (nm_end == string::npos ? r.name.size() : nm_end) - 1;
  r.name.erase(copy_n(cbegin(r.name) + 1, nm_sz, begin(r.name)), cend(r.name));

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

// read FASTA file
static bgzf_file &operator>>(bgzf_file &s, string &r) {
  string line;
  while (getline(s, line)) {
    if (!line.empty() && line[0] != '>' && line[0] != '@' &&
        line.find('N') == string::npos) {
      std::transform(line.begin(), line.end(), line.begin(),
                     [](unsigned char c) { return std::toupper(c); });
      r = line;
      break;
    }
  }
  return s;
}

void count_kmers(unordered_map<string, double> &total_counts,
                 unordered_map<string, double> &head_counts, const bool head,
                 const uint32_t kmer_value, const string &s) {
  if (head) {
    string head_kmer = s.substr(0, kmer_value);
    if (head_kmer.find('N') == string::npos)
      ++head_counts[head_kmer];
  }
  for (auto i = 0u; i < s.size() - kmer_value + 1; ++i) {
    string kmer = s.substr(i, kmer_value);
    if (kmer.find('N') == string::npos)
      ++total_counts[kmer];
  }
}

#include <boost/math/distributions/binomial.hpp>

double binomial_test(uint32_t kmer_count, uint32_t total_positions,
                     double prob) {
  boost::math::binomial_distribution<> binom(total_positions, prob);
  return 1 - boost::math::cdf(binom, kmer_count - 1);
}

// calculate over-represented k-mers at 5' end
unordered_map<string, double>
calculate_pval(const unordered_map<string, double> &total_counts,
               const unordered_map<string, double> &head_counts,
               const uint32_t kmer_value, const uint64_t seq_len,
               const uint64_t n_seqs, const double p_val_cutoff) {
  unordered_map<string, double> p_values;
  for (const auto &k : head_counts) {
    const double total = total_counts.at(k.first);
    const double prob = total / ((seq_len - kmer_value + 1) * n_seqs);
    const double p_value = binomial_test(k.second, n_seqs, prob);
    // cout << "kmer: " << k.first << ", head count = " << k.second
    //      << ", total count = " << total << ", p-value = " << p_value << '\n';
    if (p_value < p_val_cutoff)
      p_values[k.first] = p_value;
  }
  return p_values;
}

// update kmer_counts by converting part of T back to C based on kmer frequency
void ct_convertion(unordered_map<string, double> &kmer_counts,
                   const unordered_map<string, double> &p_values,
                   const unordered_map<string, double> &kmer_freq,
                   const uint64_t n_seqs) {
  for (const auto &k : p_values) {
    string t_kmer = k.first;
    for (char &c : t_kmer) {
      if (c == 'T') {
        uint32_t t_exp = kmer_freq.at(t_kmer) * n_seqs;
        int32_t delta = kmer_counts.at(t_kmer) - t_exp;
        if (delta > 0) {
          kmer_counts.at(t_kmer) -= delta;
          // cout << "kmer = " << t_kmer << ", t_exp = " << t_exp
          //      << ", delta = " << delta
          //      << ", new_t_count = " << kmer_counts.at(t_kmer);
          c = 'C';
          kmer_counts.at(t_kmer) += delta;
          // cout << ", new_c_count = " << kmer_counts.at(t_kmer) << '\n';
        }
      }
    }
  }
}

// find the k-mer with the max fraction at 5' end
pair<double, double>
find_max_frac(unordered_map<string, double> &kmer_counts,
              const unordered_map<string, double> &p_values,
              const uint64_t n_seqs) {
  // pair<frac, p-value>
  pair<double, double> max_enriched(0.0, 0.0);
  for (const auto &k : p_values) {
    double frac = kmer_counts.at(k.first) / static_cast<double>(n_seqs);
    if (frac > max_enriched.first)
      max_enriched.first = frac;
  }
  return max_enriched;
}

pair<double, double>
kmer_enrich(const unordered_map<string, double> &total_counts,
            unordered_map<string, double> &head_counts,
            const unordered_map<string, double> &kmer_freq,
            const uint32_t kmer_value, const uint64_t seq_len,
            const uint64_t n_seqs, const double p_val_cutoff) {
  const auto p_values = calculate_pval(total_counts, head_counts, kmer_value,
                                       seq_len, n_seqs, p_val_cutoff);
  ct_convertion(head_counts, p_values, kmer_freq, n_seqs);
  return find_max_frac(head_counts, p_values, n_seqs);
}

int main_guessprotocol(int argc, const char **argv) {

  try {

    static const vector<double> human_base_comp = {0.295, 0.205, 0.205, 0.295};
    static const vector<double> flat_base_comp = {0.25, 0.25, 0.25, 0.25};

    constexpr auto description = "guess bisulfite protocol for a library";

    bool verbose = false;
    bool use_human = false;
    string outfile;
    size_t reads_to_check = 1000000;
    size_t name_suffix_len = 0;
    double bisulfite_conversion_rate = 0.98;
    uint32_t kmer_value = 3;
    string ref_genome{};
    double p_val_cutoff = 0.01;

    namespace fs = std::filesystem;
    const string cmd_name = std::filesystem::path(argv[0]).filename();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(cmd_name, description,
                           "<end1-fastq> [<end2-fastq>]");
    opt_parse.add_opt("nreads", 'n', "number of reads in initial check", false,
                      reads_to_check);
    opt_parse.add_opt("ignore", 'i',
                      "length of read name suffix "
                      "to ignore when matching",
                      false, name_suffix_len);
    opt_parse.add_opt("bisulfite", 'b', "bisulfite conversion rate", false,
                      bisulfite_conversion_rate);
    opt_parse.add_opt("human", 'H', "assume human genome", false, use_human);
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("verbose", 'v',
                      "report available information during the run", false,
                      verbose);
    opt_parse.add_opt("refgenome", 'r',
                      "reference genome to calculate k-mer frequency", true,
                      ref_genome);
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
    if (use_human)
      base_comp = human_base_comp;

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

    // calculate kmer frequency using reference genome
    bgzf_file ref(ref_genome, "r");
    if (!ref)
      throw runtime_error("cannot open file: " + ref_genome);

    unordered_map<string, double> kmer_freq, tmp;
    uint64_t n_ref = 0;
    uint64_t ref_seq_len = 0;
    string ref_seq;

    while (ref >> ref_seq && n_ref < reads_to_check) {
      if (n_ref == 0)
        ref_seq_len = ref_seq.size();
      n_ref++;
      count_kmers(kmer_freq, tmp, false, kmer_value, ref_seq);
    }
    double total_pos = n_ref * (ref_seq_len - kmer_value + 1);
    std::for_each(
        kmer_freq.begin(), kmer_freq.end(),
        [total_pos](pair<const string, double> &p) { p.second /= total_pos; });

    if (reads_files.size() == 2) {

      // input: paired-end reads with end1 and end2
      bgzf_file in1(reads_files.front(), "r");
      if (!in1)
        throw runtime_error("cannot open file: " + reads_files.front());

      bgzf_file in2(reads_files.back(), "r");
      if (!in2)
        throw runtime_error("cannot open file: " + reads_files.back());

      FASTQRecord r1, r2;
      unordered_map<string, double> total_counts_1{}, total_counts_2{},
          head_counts_1{}, head_counts_2{};

      while (in1 >> r1 && in2 >> r2 && summary.n_reads < reads_to_check) {
        if (summary.n_reads == 0)
          summary.read_len = size(r1.seq);
        summary.n_reads++;

        if (!mates(name_suffix_len, r1, r2))
          throw runtime_error("expected mates: " + r1.name + ", " + r2.name);

        const double ta = t_rich_model(r1.seq) + a_rich_model(r2.seq);
        const double at = a_rich_model(r1.seq) + t_rich_model(r2.seq);

        const auto prob_read1_t_rich = exp(ta - log_sum_log(ta, at));
        summary.n_reads_wgbs += prob_read1_t_rich;

        count_kmers(total_counts_1, head_counts_1, true, kmer_value, r1.seq);
        count_kmers(total_counts_2, head_counts_2, true, kmer_value, r2.seq);
      }
      const auto r1_enrich =
          kmer_enrich(total_counts_1, head_counts_1, kmer_freq, kmer_value,
                      summary.read_len, summary.n_reads, p_val_cutoff);
      const auto r2_enrich =
          kmer_enrich(total_counts_2, head_counts_2, kmer_freq, kmer_value,
                      summary.read_len, summary.n_reads, p_val_cutoff);

      // take the average of the two enrichment
      summary.rrbs_fraction = (r1_enrich.first + r2_enrich.first) / 2;
    } else {
      // input: single-end reads
      bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw runtime_error("cannot open file: " + reads_files.front());

      FASTQRecord r;
      unordered_map<string, double> total_counts{}, head_counts{};

      while (in >> r && summary.n_reads < reads_to_check) {
        if (summary.n_reads == 0)
          summary.read_len = size(r.seq);
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;

        count_kmers(total_counts, head_counts, true, kmer_value, r.seq);
      }
      const auto r_enrich =
          kmer_enrich(total_counts, head_counts, kmer_freq, kmer_value,
                      summary.read_len, summary.n_reads, p_val_cutoff);
      summary.rrbs_fraction = r_enrich.first;
    }

    summary.evaluate();

    if (!outfile.empty()) {
      std::ofstream out(outfile);
      if (!out)
        throw runtime_error("failed to open: " + outfile);
      out << summary.tostring() << endl;
    } else
      cout << summary.tostring() << endl;
  } catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
