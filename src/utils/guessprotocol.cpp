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

constexpr char idx_to_nuc[] = {'A', 'C', 'G', 'T'};

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

template <typename T = string::const_iterator>
uint64_t kmer_to_int(T start, const T stop) {
  uint64_t num = 0;
  uint8_t base = 4;
  for (T it = start; it != stop; ++it) {
    int idx = nuc_to_idx[static_cast<uint8_t>(*it)];
    if (idx == 4) {
      throw std::invalid_argument(
          "Invalid nucleotide character when converting to integer");
    }
    num = num * base + idx;
  }
  return num;
}

string int_to_kmer(uint64_t num, const uint32_t kmer_value) {
  string kmer(kmer_value, 'N');
  uint8_t base = 4;
  for (size_t i = kmer_value; i > 0; --i) {
    kmer[i - 1] = idx_to_nuc[num % base];
    num /= base;
  }
  return kmer;
}

void count_kmers(const bool head, const uint32_t kmer_value, const string &s,
                 int barcode, vector<double> &total_counts,
                 vector<double> &head_counts) {
  try {
    const uint32_t start = std::max(0, barcode);
    if (head) {
      if (s.size() < kmer_value)
        throw std::out_of_range("Read length smaller than k-mer length");

      string head_kmer = s.substr(start, kmer_value);
      if (head_kmer.find('N') == string::npos)
        head_counts[kmer_to_int(head_kmer.begin(), head_kmer.end())]++;
    }
    for (auto i = start; i <= size(s) - start - kmer_value + 1; ++i) {
      string kmer = s.substr(i, kmer_value);
      if (kmer.find('N') == string::npos)
        total_counts[kmer_to_int(kmer.begin(), kmer.end())]++;
    }
  } catch (const std::out_of_range &e) {
    cerr << "Error: " << e.what() << endl;
    cerr << "Sequence: " << s << endl;
  }
}

vector<double> load_kmer_freq(const string &filename,
                              vector<double> &kmer_freq) {
  std::ifstream infile(filename);
  if (!infile)
    throw runtime_error("cannot open file: " + filename);

  string kmer;
  double freq;
  while (infile >> kmer >> freq) {
    uint64_t idx = kmer_to_int(kmer.begin(), kmer.end());
    kmer_freq[idx] = freq;
  }
  return kmer_freq;
}

vector<double> calculate_kmer_freq(const string &ref_genome,
                                   const uint32_t kmer_value,
                                   const uint64_t reads_to_check,
                                   vector<double> &kmer_freq) {
  bgzf_file ref(ref_genome, "r");
  if (!ref)
    throw runtime_error("cannot open file: " + ref_genome);

  uint64_t n_ref = 0;
  uint64_t ref_seq_len = 0;
  string ref_seq;
  vector<double> tmp;

  while (ref >> ref_seq && n_ref < reads_to_check) {
    if (n_ref == 0)
      ref_seq_len = ref_seq.size();
    n_ref++;
    count_kmers(false, kmer_value, ref_seq, 0, kmer_freq, tmp);
  }
  double total_pos = n_ref * (ref_seq_len - kmer_value + 1);
  for (auto &freq : kmer_freq)
    freq /= total_pos;

  return kmer_freq;
}

void save_kmer_freq(const string &filename, const vector<double> &kmer_freq,
                    const uint32_t kmer_value) {
  std::ofstream outfile(filename);
  if (!outfile)
    throw runtime_error("cannot open file for writing: " + filename);

  for (size_t i = 0; i < kmer_freq.size(); ++i) {
    string kmer = int_to_kmer(i, kmer_value);
    outfile << kmer << " " << kmer_freq[i] << "\n";
  }
}

#include <boost/math/distributions/binomial.hpp>

double binomial_test(uint32_t kmer_count, uint32_t total_positions,
                     double prob) {
  boost::math::binomial_distribution<> binom(total_positions, prob);
  return 1 - boost::math::cdf(binom, kmer_count - 1);
}

void ct_conversion(const vector<double> &kmer_freq, const uint64_t n_seqs,
                   const uint32_t kmer_value, vector<double> &kmer_counts) {
  for (size_t i = 0; i < kmer_counts.size(); ++i) {
    string t_kmer = int_to_kmer(i, kmer_value);
    try {
      for (char &c : t_kmer)
        if (c == 'T') {
          uint32_t t_exp = kmer_freq[i] * n_seqs;
          int32_t delta = kmer_counts[i] - t_exp;
          if (delta > 0) {
            kmer_counts[i] -= delta;
            // cout << "kmer = " << t_kmer << ", t_exp = " << t_exp
            //       << ", delta = " << delta
            //       << ", new_t_count = " << kmer_counts[i];
            c = 'C';
            uint64_t c_idx = kmer_to_int(t_kmer.begin(), t_kmer.end());
            kmer_counts[c_idx] += delta;
            // cout << ", new_c_count = " << kmer_counts[c_idx] << '\n';
          }
        }
    } catch (const std::out_of_range &e) {
      cerr << "Error: " << e.what() << endl;
      cerr << "Key not found: " << t_kmer << endl;
    }
  }
}

// calculate over-represented k-mers at 5' end
vector<double> calculate_pval(const vector<double> &total_counts,
                              const vector<double> &head_counts,
                              const uint32_t kmer_value, const uint64_t seq_len,
                              const uint64_t n_seqs) {
  vector<double> p_values(static_cast<size_t>(std::pow(4, kmer_value)),
                          std::numeric_limits<double>::infinity());
  for (size_t i = 0; i < head_counts.size(); ++i) {
    const double total = total_counts[i];
    const double prob = total / ((seq_len - kmer_value + 1) * n_seqs);
    if (head_counts[i] != 0)
      p_values[i] = binomial_test(head_counts[i], n_seqs, prob);
    // cout << "kmer: " << int_to_kmer(i, kmer_value) << ", head count = " <<
    // head_counts[i] << ", total count = " << total << ", p-value = " <<
    // p_value << '\n';
  }
  return p_values;
}

// find the k-mer with the max fraction at 5' end
double check_significance(const vector<double> &kmer_counts,
                          const vector<double> &p_values, const uint64_t n_seqs,
                          const double p_val_cutoff) {
  // pair<frac, p-value>
  double max_enriched = 0.0;
  for (size_t i = 0; i < kmer_counts.size(); ++i) {
    double frac = kmer_counts[i] / n_seqs;
    if (frac > max_enriched && p_values[i] < p_val_cutoff)
      max_enriched = frac;
  }
  return max_enriched;
}

double kmer_enrich(const vector<double> &total_counts,
                   const vector<double> &kmer_freq, const uint32_t kmer_value,
                   const uint64_t seq_len, const uint64_t n_seqs,
                   const double p_val_cutoff, vector<double> &head_counts) {
  ct_conversion(kmer_freq, n_seqs, kmer_value, head_counts);
  const auto p_values =
      calculate_pval(total_counts, head_counts, kmer_value, seq_len, n_seqs);
  return check_significance(head_counts, p_values, n_seqs, p_val_cutoff);
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
    int barcode = 0;

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
    opt_parse.add_opt("barcode", 'x',
                      "barcode length of the 5' end that needs to be omitted "
                      "when determining RRBS",
                      false, barcode);
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

    // look for the kmer frequency file
    string kmer_freq_file = ref_genome + ".gp";
    vector<double> kmer_freq(static_cast<size_t>(std::pow(4, kmer_value)), 0.0);

    if (std::filesystem::exists(kmer_freq_file)) {
      kmer_freq = load_kmer_freq(kmer_freq_file, kmer_freq);
      cout << "Loaded k-mer frequencies from " << kmer_freq_file << '\n';
    } else {
      kmer_freq = calculate_kmer_freq(ref_genome, kmer_value, reads_to_check,
                                      kmer_freq);
      save_kmer_freq(kmer_freq_file, kmer_freq, kmer_value);
      cout << "Calculated and saved k-mer frequencies to " << kmer_freq_file
           << '\n';
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
      vector<double> total_counts_1(
          static_cast<size_t>(std::pow(4, kmer_value)), 0.0);
      vector<double> total_counts_2(
          static_cast<size_t>(std::pow(4, kmer_value)), 0.0);
      vector<double> head_counts_1(static_cast<size_t>(std::pow(4, kmer_value)),
                                   0.0);
      vector<double> head_counts_2(static_cast<size_t>(std::pow(4, kmer_value)),
                                   0.0);

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

        count_kmers(true, kmer_value, r1.seq, barcode, total_counts_1,
                    head_counts_1);
        count_kmers(true, kmer_value, r2.seq, barcode, total_counts_2,
                    head_counts_2);
      }
      const auto r1_enrich =
          kmer_enrich(total_counts_1, kmer_freq, kmer_value, summary.read_len,
                      summary.n_reads, p_val_cutoff, head_counts_1);
      const auto r2_enrich =
          kmer_enrich(total_counts_2, kmer_freq, kmer_value, summary.read_len,
                      summary.n_reads, p_val_cutoff, head_counts_2);

      // take the average of the two enrichment
      summary.rrbs_fraction = (r1_enrich + r2_enrich) / 2;
    } else {
      // input: single-end reads
      bgzf_file in(reads_files.front(), "r");
      if (!in)
        throw runtime_error("cannot open file: " + reads_files.front());

      FASTQRecord r;
      vector<double> total_counts(static_cast<size_t>(std::pow(4, kmer_value)),
                                  0.0);
      vector<double> head_counts(static_cast<size_t>(std::pow(4, kmer_value)),
                                 0.0);

      while (in >> r && summary.n_reads < reads_to_check) {
        if (summary.n_reads == 0)
          summary.read_len = size(r.seq);
        summary.n_reads++;

        const double t = t_rich_model(r.seq);
        const double a = a_rich_model(r.seq);

        const auto prob_t_rich = exp(t - log_sum_log(t, a));
        summary.n_reads_wgbs += prob_t_rich;

        count_kmers(true, kmer_value, r.seq, barcode, total_counts,
                    head_counts);
      }
      const auto r_enrich =
          kmer_enrich(total_counts, kmer_freq, kmer_value, summary.read_len,
                      summary.n_reads, p_val_cutoff, head_counts);
      summary.rrbs_fraction = r_enrich;
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
