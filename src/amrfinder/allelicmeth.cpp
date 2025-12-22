/* Copyright (C) 2014-2022 University of Southern California,
 *                         Andrew D. Smith, and Benjamin E. Decato
 *
 * Authors: Andrew D. Smith and Benjamin E. Decato
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

#include "Epiread.hpp"
#include "MSite.hpp"

#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0)
    return q;
  if (q == 0)
    return p;
  return p > q ? p + std::log1p(std::exp(q - p))
               : q + std::log1p(std::exp(p - q));
}

static inline double
lnchoose(const unsigned int n, unsigned int m) {
  if (m == n || m == 0)
    return 0;
  if (m * 2 > n)
    m = n - m;
  using std::lgamma;
  return lgamma(n + 1.0) - lgamma(m + 1.0) - lgamma((n - m) + 1.0);
}

// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double
log_hyper_g(const std::size_t k, const std::size_t n1, const std::size_t n2,
            const std::size_t t) {
  return lnchoose(n1, k) + lnchoose(n2, t - k) - lnchoose(n1 + n2, t);
}

static double
fishers_exact(std::size_t a, std::size_t b, std::size_t c, std::size_t d) {
  const std::size_t m = a + c;  // sum of first column
  const std::size_t n = b + d;  // sum of second column
  const std::size_t k = a + b;  // sum of first row
  // ADS: want more extreme than "observed"
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (std::size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return std::exp(p);
}

static std::size_t
state_pair_to_index(const std::string &s, const std::size_t idx) {
  assert(idx < s.length() - 1);
  const char a = s[idx];
  if (a == 'C') {
    const char b = s[idx + 1];
    if (b == 'C')
      return 0;
    if (b == 'T')
      return 1;
    return 4;
  }
  if (a == 'T') {
    const char b = s[idx + 1];
    if (b == 'C')
      return 2;
    if (b == 'T')
      return 3;
    return 4;
  }
  return 4;
}

template <class T> struct PairStateCounter {
  T CC;
  T CT;
  T TC;
  T TT;

  double
  score() const {
    return (CC * TT > CT * TC) ? fishers_exact(CC, CT, TC, TT)
                               : fishers_exact(CT, CC, TT, TC);
  }
  double
  total() const {
    return CC + CT + TC + TT;
  }

  std::string
  tostring() const {
    return toa(CC) + '\t' + toa(CT) + '\t' + toa(TC) + '\t' + toa(TT);
  }

  void
  increment(const std::size_t state) {
    if (state == 0)
      ++CC;
    else if (state == 1)
      ++CT;
    else if (state == 2)
      ++TC;
    else if (state == 3)
      ++TT;
  }
};

template <typename T>
void
fit_states(const epiread &er, std::vector<PairStateCounter<T>> &counts) {
  for (std::size_t i = 0; i < er.length() - 1; ++i) {
    const std::size_t pos = er.pos + i;
    assert(pos < counts.size());
    const std::size_t curr_state = state_pair_to_index(er.seq, i);
    counts[pos].increment(curr_state);
  }
}

static void
collect_cpgs(const std::string &s,
             std::unordered_map<std::size_t, std::size_t> &cpgs) {
  std::size_t cpg_idx = 0;
  std::size_t nuc_idx = 0;
  const auto lim = end(s) - (s.size() == 0 ? 0 : 1);
  for (auto itr = begin(s); itr != lim; ++itr, ++nuc_idx)
    if (*itr == 'C' && *(itr + 1) == 'G')
      cpgs[cpg_idx++] = nuc_idx;
}

/* The "sites" in the convert_coordinates function represent pairs of
   consecutive CpG sites. So there would be one fewer of them than the
   total CpGs in a chromosome. */
static void
convert_coordinates(const std::string &chrom, std::vector<MSite> &sites) {
  std::unordered_map<std::size_t, std::size_t> cpgs;
  collect_cpgs(chrom, cpgs);
  for (std::size_t i = 0; i < sites.size(); ++i) {
    const auto pos_itr = cpgs.find(sites[i].pos);
    if (pos_itr == end(cpgs))
      throw std::runtime_error("failed converting site:\n" +
                               sites[i].tostring());
    sites[i].pos = pos_itr->second;
  }
}

template <typename T>
void
add_cytosine(const std::string &chrom_name, const std::size_t start_cpg,
             std::vector<PairStateCounter<T>> &counts,
             std::vector<MSite> &cytosines) {
  std::ostringstream s;
  s << counts[start_cpg].score() << "\t" << counts[start_cpg].total() << "\t"
    << counts[start_cpg].tostring();
  const std::string name(s.str());
  cytosines.push_back(MSite(chrom_name, start_cpg, '+', name, 0, 0));
}

template <typename T>
void
process_chrom(const std::string &chrom_name,
              const std::vector<epiread> &epireads,
              std::vector<MSite> &cytosines,
              std::vector<PairStateCounter<T>> &counts) {
  for (std::size_t i = 0; i < epireads.size(); ++i)
    fit_states(epireads[i], counts);
  for (std::size_t i = 0; i < counts.size(); ++i)
    add_cytosine(chrom_name, i, counts, cytosines);
}

static void
update_chroms_seen(const std::string &chrom_name,
                   std::unordered_set<std::string> &chroms_seen) {
  const auto chr_itr = chroms_seen.find(chrom_name);
  if (chr_itr != end(chroms_seen))
    throw std::runtime_error("chroms out of order: " + chrom_name);
  chroms_seen.insert(chrom_name);
}

static void
verify_chroms_available(
  const std::string &chrom_name,
  std::unordered_map<std::string, std::size_t> &chrom_lookup) {
  const auto chr_itr = chrom_lookup.find(chrom_name);
  if (chr_itr == end(chrom_lookup))
    throw std::runtime_error("chrom not found: " + chrom_name);
}

int
main_allelicmeth(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    static const auto description = R"(
computes probability of allele-specific methylation at each tuple of CpGs
)";
    bool VERBOSE = false;

    std::string outfile;
    std::string chroms_dir;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<epireads>");
    opt_parse.add_opt("output", 'o', "output file name", true, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory", true,
                      chroms_dir);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
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
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string epi_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::vector<std::string> chrom_names;
    std::vector<std::string> chroms;
    read_fasta_file_short_names(chroms_dir, chrom_names, chroms);
    for (auto &&i : chroms)
      transform(begin(i), end(i), begin(i),
                [](const char c) { return std::toupper(c); });

    // lookup to map chrom names to chrom sequences
    std::unordered_map<std::string, std::size_t> chrom_lookup;
    for (std::size_t i = 0; i < chrom_names.size(); ++i)
      chrom_lookup.insert(make_pair(chrom_names[i], i));

    std::unordered_map<std::string, std::size_t> chrom_sizes;
    for (std::size_t i = 0; i < chrom_names.size(); ++i) {
      std::size_t cpg_count = 0;
      for (std::size_t j = 0; j < chroms[i].size() - 1; ++j)
        cpg_count += (chroms[i][j] == 'C' && chroms[i][j + 1] == 'G');
      chrom_sizes.insert(make_pair(chrom_names[i], cpg_count));
    }

    if (VERBOSE)
      std::cerr << "number of chromosomes: " << chrom_sizes.size() << '\n';

    std::ifstream in(epi_file);
    if (!in)
      throw std::runtime_error("cannot open input file: " + epi_file);

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);

    std::unordered_set<std::string> chroms_seen;
    std::string chrom;
    epiread er;
    std::vector<epiread> epireads;
    while (in >> er) {
      if (er.chr != chrom) {
        update_chroms_seen(er.chr, chroms_seen);
        verify_chroms_available(er.chr, chrom_lookup);

        if (VERBOSE)
          std::cerr << "[processing " << er.chr << "]" << '\n';

        if (!chrom.empty()) {
          std::vector<PairStateCounter<std::uint32_t>> counts(
            chrom_sizes[chrom] - 1);
          std::vector<MSite> cytosines;
          process_chrom(chrom, epireads, cytosines, counts);
          const std::size_t chrom_idx = chrom_lookup[chrom];
          convert_coordinates(chroms[chrom_idx], cytosines);
          for (std::size_t i = 0; i < cytosines.size() - 1; ++i) {
            out << cytosines[i].chrom << "\t" << cytosines[i].pos
                << "\t+\tCpG\t" << cytosines[i].context << '\n';
          }
        }
        epireads.clear();
      }
      epireads.push_back(er);
      chrom.swap(er.chr);
    }
    if (!chrom.empty()) {
      std::vector<PairStateCounter<std::uint32_t>> counts(chrom_sizes[chrom] -
                                                          1);
      std::vector<MSite> cytosines;
      process_chrom(chrom, epireads, cytosines, counts);
      const std::size_t chrom_idx = chrom_lookup[chrom];
      convert_coordinates(chroms[chrom_idx], cytosines);
      for (std::size_t i = 0; i < cytosines.size() - 1; ++i) {
        out << cytosines[i].chrom << "\t" << cytosines[i].pos << "\t+\tCpG\t"
            << cytosines[i].context << '\n';
      }
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-avoid-magic-numbers,*-narrowing-conversions)
