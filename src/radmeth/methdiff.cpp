/* methdiff: Computes probability that individual CpGs have higher
 *           methylation in file A than in file B, where files A and B
 *           are specified on the command line.
 *
 * Copyright (C) 2011-2023 Andrew D Smith
 *
 * Author: Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <bamxx.hpp>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "MSite.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::min;
using std::runtime_error;
using std::string;
using std::vector;

using std::ofstream;
using std::ostream_iterator;

using bamxx::bgzf_file;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) { return q; }
  else if (q == 0) { return p; }
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + log(1.0 + exp(smaller - larger));
}

static inline double
lnchoose(const unsigned int n, unsigned int m) {
  if (m == n || m == 0) return 0;
  if (m * 2 > n) m = n - m;
  using std::lgamma;
  return lgamma(n + 1.0) - lgamma(m + 1.0) - lgamma((n - m) + 1.0);
}

static inline double
log_hyper_g_greater(size_t meth_a, size_t unmeth_a, size_t meth_b,
                    size_t unmeth_b, size_t k) {
  return (
    lnchoose(meth_b + unmeth_b - 1, k) +
    lnchoose(meth_a + unmeth_a - 1, meth_a + meth_b - 1 - k) -
    lnchoose(meth_a + unmeth_a + meth_b + unmeth_b - 2, meth_a + meth_b - 1));
}

static double
test_greater_population(const size_t meth_a, const size_t unmeth_a,
                        const size_t meth_b, const size_t unmeth_b) {
  double p = 0;
  for (size_t k = (meth_b > unmeth_a) ? meth_b - unmeth_a : 0; k < meth_b; ++k)
    p = log_sum_log(p,
                    log_hyper_g_greater(meth_a, unmeth_a, meth_b, unmeth_b, k));
  return exp(p);
}

template<class T> T &
write_methdiff_site(T &out, const MSite &a, const MSite &b,
                    const double diffscore) {
  std::ostringstream oss;
  // clang-format off
  oss << a.chrom << '\t'
      << a.pos << '\t'
      << a.strand << '\t'
      << a.context << '\t'
      << diffscore << '\t'
      << a.n_meth() << '\t'
      << a.n_unmeth() << '\t'
      << b.n_meth() << '\t'
      << b.n_unmeth() << '\n';
  // clang-format on
  return out << oss.str();
}

bgzf_file &
write_methdiff_site(bgzf_file &out, const MSite &a, const MSite &b,
                    const double diffscore) {
  std::ostringstream oss;
  // clang-format off
  oss << a.chrom << '\t'
      << a.pos << '\t'
      << a.strand << '\t'
      << a.context << '\t'
      << diffscore << '\t'
      << a.n_meth() << '\t'
      << a.n_unmeth() << '\t'
      << b.n_meth() << '\t'
      << b.n_unmeth() << '\n';
  // clang-format on
  out.write(oss.str());
  return out;
}

static bool
site_precedes(const size_t lhs_chrom_idx, const size_t lhs_pos,
              const size_t rhs_chrom_idx, const size_t rhs_pos) {
  return (lhs_chrom_idx < rhs_chrom_idx ||
          (lhs_chrom_idx == rhs_chrom_idx && (lhs_pos < rhs_pos)));
}

static size_t
get_chrom_id(std::unordered_map<string, size_t> &chrom_order,
             std::unordered_set<string> &chroms_seen, const MSite &s) {
  if (chroms_seen.find(s.chrom) != end(chroms_seen))
    throw runtime_error("unsorted chromosomes found");
  chroms_seen.insert(s.chrom);
  auto x = chrom_order.find(s.chrom);
  if (x == end(chrom_order)) {
    const size_t idx = chrom_order.size();
    chrom_order.insert(std::make_pair(s.chrom, idx));
    return idx;
  }
  else
    return x->second;
}

static string
bad_order(const std::unordered_map<string, size_t> &chrom_order,
          const string prev_chrom, const size_t prev_pos, const string chrom,
          const size_t pos) {
  std::ostringstream oss;
  const size_t chrom_id = chrom_order.find(chrom)->second;
  const size_t prev_chrom_id = chrom_order.find(prev_chrom)->second;
  // clang-format off
  oss << "bad order:\n"
      << "chrom=" << chrom << " [id=" << chrom_id << "] "
      << "pos=" << pos << "\n"
      << "appears after\n"
      << "chrom=" << prev_chrom << " [id=" << prev_chrom_id << "] "
      << "pos=" << prev_pos << "\n";
  // clang-format on
  return oss.str();
}

template<class T> static void
process_sites(const bool VERBOSE, bgzf_file &in_a, bgzf_file &in_b,
              const bool allow_uncovered, const double pseudocount, T &out) {
  // chromosome order in the files
  std::unordered_map<string, size_t> chrom_order;
  std::unordered_set<string> chroms_seen_a, chroms_seen_b;

  MSite a, b;
  string prev_chrom_a, prev_chrom_b;
  size_t chrom_id_a = 0, chrom_id_b = 0;
  size_t prev_chrom_id_a = 0, prev_chrom_id_b = 0;
  size_t prev_pos_a = 0, prev_pos_b = 0;

  string line;
  while (read_site(in_a, a)) {
    if (prev_chrom_a.compare(a.chrom) != 0) {
      prev_chrom_id_a = chrom_id_a;
      chrom_id_a = get_chrom_id(chrom_order, chroms_seen_a, a);
      prev_chrom_a = a.chrom;
      if (VERBOSE) cerr << "processing " << a.chrom << endl;
    }
    if (site_precedes(chrom_id_a, a.pos, prev_chrom_id_a, prev_pos_a))
      throw runtime_error(
        bad_order(chrom_order, prev_chrom_a, prev_pos_a, a.chrom, a.pos));

    bool advance_b = true;
    while (advance_b && read_site(in_b, b)) {
      if (prev_chrom_b.compare(b.chrom) != 0) {
        prev_chrom_id_b = chrom_id_b;
        chrom_id_b = get_chrom_id(chrom_order, chroms_seen_b, b);
        prev_chrom_b = b.chrom;
      }
      if (site_precedes(chrom_id_b, b.pos, prev_chrom_id_b, prev_pos_b))
        throw runtime_error(
          bad_order(chrom_order, prev_chrom_b, prev_pos_b, b.chrom, b.pos));
      advance_b = site_precedes(chrom_id_b, b.pos, chrom_id_a, a.pos);
      prev_pos_b = b.pos;
    }

    if (chrom_id_a == chrom_id_b && a.pos == b.pos) {
      if (allow_uncovered || min(a.n_reads, b.n_reads) > 0) {
        const size_t meth_a = a.n_meth() + pseudocount;
        const size_t unmeth_a = a.n_unmeth() + pseudocount;
        const size_t meth_b = b.n_meth() + pseudocount;
        const size_t unmeth_b = b.n_unmeth() + pseudocount;

        const double diffscore =
          test_greater_population(meth_b, unmeth_b, meth_a, unmeth_a);

        write_methdiff_site(out, a, b, diffscore);
      }
    }
    prev_pos_a = a.pos;
  }
}

int
main_methdiff(int argc, const char **argv) {
  try {
    string outfile;
    size_t pseudocount = 1;

    // run mode flags
    bool allow_uncovered = true;
    bool VERBOSE = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "compute probability that site "
                           "has higher methylation in file A than B",
                           "<counts-a> <counts-b>");
    opt_parse.add_opt("pseudo", 'p', "pseudocount (default: 1)", false,
                      pseudocount);
    opt_parse.add_opt("nonzero-only", 'A',
                      "process only sites with coveage in both samples", false,
                      allow_uncovered);
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
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
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file_a = leftover_args[0];
    const string cpgs_file_b = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE)
      cerr << "[opening methcounts file: " << cpgs_file_a << "]" << endl;
    bgzf_file in_a(cpgs_file_a, "r");
    if (!in_a) throw runtime_error("cannot open file: " + cpgs_file_a);

    if (VERBOSE)
      cerr << "[opening methcounts file: " << cpgs_file_b << "]" << endl;
    bgzf_file in_b(cpgs_file_b, "r");
    if (!in_b) throw runtime_error("cannot open file: " + cpgs_file_b);

    if (outfile.empty() || !has_gz_ext(outfile)) {
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile);
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

      process_sites(VERBOSE, in_a, in_b, allow_uncovered, pseudocount, out);
    }
    else {
      bgzf_file out(outfile, "w");
      process_sites(VERBOSE, in_a, in_b, allow_uncovered, pseudocount, out);
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
