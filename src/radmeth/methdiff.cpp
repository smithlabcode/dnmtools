/* methdiff: Computes probability that individual CpGs have higher methylation
 *           in file A than in file B, where files A and B are specified on
 *           the command line.
 *
 * Copyright (C) 2011-2025 Andrew D Smith
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

#include <bamxx.hpp>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "MSite.hpp"

using bamxx::bgzf_file;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0)
    return q;
  if (q == 0)
    return p;
  const double larger = (p > q) ? p : q;
  const double smaller = (p > q) ? q : p;
  return larger + std::log1p(std::exp(smaller - larger));
}

static inline double
lnchoose(const unsigned int n, unsigned int m) {
  if (m == n || m == 0)
    return 0;
  if (m * 2 > n)
    m = n - m;
  return std::lgamma(n + 1.0) - std::lgamma(m + 1.0) - std::lgamma(n - m + 1.0);
}

static inline double
log_hyper_g_greater(const std::size_t meth_a, const std::size_t unmeth_a,
                    const std::size_t meth_b, const std::size_t unmeth_b,
                    std::size_t k) {
  return lnchoose(meth_b + unmeth_b - 1, k) +
         lnchoose(meth_a + unmeth_a - 1, meth_a + meth_b - 1 - k) -
         lnchoose(meth_a + unmeth_a + meth_b + unmeth_b - 2,
                  meth_a + meth_b - 1);
}

static inline double
test_greater_population(const std::size_t meth_a, const std::size_t unmeth_a,
                        const std::size_t meth_b, const std::size_t unmeth_b) {
  double p = 0;
  for (std::size_t k = (meth_b > unmeth_a) ? meth_b - unmeth_a : 0; k < meth_b;
       ++k)
    p = log_sum_log(p,
                    log_hyper_g_greater(meth_a, unmeth_a, meth_b, unmeth_b, k));
  return std::exp(p);
}

template <class T>
T &
write_methdiff_site(T &out, const MSite &a, const MSite &b,
                    const double diffscore) {
  // static constexpr auto out_fmt =
  // "%s\t%ld\t%c\t%s\t%.6g\t%ld\t%ld\t%ld\t%ld\n";
  // clang-format off
  static constexpr auto out_fmt =
    "%s"  // chrom
    "\t"
    "%ld"  // pos
    "\t"
    "%c"  // strand
    "\t"
    "%s"  // context
    "\t"
    "%.6g"  // diffscore
    "\t"
    "%ld"  // a.n_meth()
    "\t"
    "%ld"  // a.n_unmeth()
    "\t"
    "%ld"  // b.n_meth()
    "\t"
    "%ld"  // b.n_unmeth()
    "\n";
  // clang-format on
  static constexpr auto buf_size = 1024;
  static char buffer[buf_size];

  // clang-format off
  const int r = std::snprintf(buffer, buf_size, out_fmt,
                              a.chrom.data(),
                              a.pos,
                              a.strand,
                              a.context.data(),
                              diffscore,
                              a.n_meth(),
                              a.n_unmeth(),
                              b.n_meth(),
                              b.n_unmeth());
  // clang-format on
  if (r < 0)
    throw std::runtime_error("failed to write to output buffer");
  out.write(buffer, r);
  return out;
}

static inline bool
site_precedes(const std::size_t lhs_chrom_idx, const std::size_t lhs_pos,
              const std::size_t rhs_chrom_idx, const std::size_t rhs_pos) {
  return (lhs_chrom_idx < rhs_chrom_idx ||
          (lhs_chrom_idx == rhs_chrom_idx && (lhs_pos < rhs_pos)));
}

static std::size_t
get_chrom_id(std::unordered_map<std::string, std::size_t> &chrom_order,
             std::unordered_set<std::string> &chroms_seen, const MSite &s) {
  if (chroms_seen.find(s.chrom) != std::end(chroms_seen))
    throw std::runtime_error("unsorted chromosomes found");
  chroms_seen.emplace(s.chrom);
  const auto itr = chrom_order.find(s.chrom);
  if (itr == std::cend(chrom_order)) {
    const auto idx = std::size(chrom_order);
    chrom_order.emplace(s.chrom, idx);
    return idx;
  }
  return itr->second;
}

static std::string
bad_order(const std::unordered_map<std::string, std::size_t> &chrom_order,
          const std::string prev_chrom, const std::size_t prev_pos,
          const std::string chrom, const std::size_t pos) {
  std::ostringstream oss;
  const std::size_t chrom_id = chrom_order.find(chrom)->second;
  const std::size_t prev_chrom_id = chrom_order.find(prev_chrom)->second;
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

template <class T>
static void
process_sites(const bool show_progress, bgzf_file &in_a, bgzf_file &in_b,
              const bool allow_uncovered, const double pseudocount, T &out) {
  // chromosome order in the files
  std::unordered_map<std::string, std::size_t> chrom_order;
  std::unordered_set<std::string> chroms_seen_a, chroms_seen_b;

  MSite a, b;  // ADS: default "pos" might be num lim max
  a.pos = 0;
  b.pos = 0;
  std::string prev_chrom_a, prev_chrom_b;
  std::size_t chrom_id_a = 0, chrom_id_b = 0;
  std::size_t prev_chrom_id_a = 0, prev_chrom_id_b = 0;
  std::size_t prev_pos_a = 0, prev_pos_b = 0;

  bool advance_a = true;
  bool advance_b = true;

  while (true) {

    while (advance_a && read_site(in_a, a)) {
      if (prev_chrom_a.compare(a.chrom) != 0) {
        prev_chrom_id_a = chrom_id_a;
        chrom_id_a = get_chrom_id(chrom_order, chroms_seen_a, a);
        if (show_progress)
          std::cerr << "processing " << a.chrom << std::endl;
        prev_chrom_a = a.chrom;
      }
      if (site_precedes(chrom_id_a, a.pos, prev_chrom_id_a, prev_pos_a))
        throw std::runtime_error(
          bad_order(chrom_order, prev_chrom_a, prev_pos_a, a.chrom, a.pos));
      advance_a = site_precedes(chrom_id_a, a.pos, chrom_id_b, b.pos);
      prev_pos_a = a.pos;
    }

    while (advance_b && read_site(in_b, b)) {
      if (prev_chrom_b.compare(b.chrom) != 0) {
        prev_chrom_id_b = chrom_id_b;
        chrom_id_b = get_chrom_id(chrom_order, chroms_seen_b, b);
        prev_chrom_b = b.chrom;
      }
      if (site_precedes(chrom_id_b, b.pos, prev_chrom_id_b, prev_pos_b))
        throw std::runtime_error(
          bad_order(chrom_order, prev_chrom_b, prev_pos_b, b.chrom, b.pos));
      advance_b = site_precedes(chrom_id_b, b.pos, chrom_id_a, a.pos);
      prev_pos_b = b.pos;
    }

    if (!in_a || !in_b)
      break;

    if (chrom_id_a == chrom_id_b && a.pos == b.pos) {
      if (allow_uncovered || std::min(a.n_reads, b.n_reads) > 0) {
        const std::size_t meth_a = a.n_meth() + pseudocount;
        const std::size_t unmeth_a = a.n_unmeth() + pseudocount;
        const std::size_t meth_b = b.n_meth() + pseudocount;
        const std::size_t unmeth_b = b.n_unmeth() + pseudocount;

        const double diffscore =
          test_greater_population(meth_b, unmeth_b, meth_a, unmeth_a);

        write_methdiff_site(out, a, b, diffscore);
      }
      advance_a = true;
      advance_b = true;
    }
    else {
      advance_a = site_precedes(chrom_id_a, a.pos, chrom_id_b, b.pos);
      advance_b = !advance_a;
    }
  }
}

int
main_methdiff(int argc, char *argv[]) {
  try {
    std::string outfile;
    double pseudocount = 1.0;

    // run mode flags
    bool allow_uncovered = true;
    bool verbose = false;

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
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << std::endl
                << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string cpgs_file_a = leftover_args[0];
    const std::string cpgs_file_b = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS *****************/

    if (verbose)
      std::cerr << "[opening counts file: " << cpgs_file_a << "]" << std::endl;
    bgzf_file in_a(cpgs_file_a, "r");
    if (!in_a)
      throw std::runtime_error("cannot open file: " + cpgs_file_a);

    if (verbose)
      std::cerr << "[opening counts file: " << cpgs_file_b << "]" << std::endl;
    bgzf_file in_b(cpgs_file_b, "r");
    if (!in_b)
      throw std::runtime_error("cannot open file: " + cpgs_file_b);

    if (outfile.empty() || !has_gz_ext(outfile)) {
      std::ofstream out(outfile);
      process_sites(verbose, in_a, in_b, allow_uncovered, pseudocount, out);
    }
    else {
      bgzf_file out(outfile, "w");
      process_sites(verbose, in_a, in_b, allow_uncovered, pseudocount, out);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
