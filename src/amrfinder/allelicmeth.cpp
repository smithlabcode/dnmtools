/* allelicmeth:
 *
 * Copyright (C) 2014-2022 University of Southern California,
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <utility>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <cmath>
#include <sstream>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "GenomicRegion.hpp"
#include "MSite.hpp"

#include "Epiread.hpp"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::unordered_map;
using std::unordered_set;
using std::max;
using std::min;
using std::runtime_error;

static inline double
log_sum_log(const double p, const double q) {
  if (p == 0) {return q;}
  else if (q == 0) {return p;}
  return p > q ? p + log(1.0 + exp(q - p)) : q + log(1.0 + exp(p - q));
}

static inline double
lnchoose(const unsigned int n, unsigned int m) {
  if (m == n || m == 0) return 0;
  if (m * 2 > n) m = n - m;
  using std::lgamma;
  return lgamma(n + 1.0) - lgamma(m + 1.0) - lgamma((n - m) + 1.0);
}

// p(k) =  C(n1, k) C(n2, t - k) / C(n1 + n2, t)
static double
log_hyper_g(const size_t k, const size_t n1, const size_t n2, const size_t t) {
  return lnchoose(n1, k) + lnchoose(n2, t - k) - lnchoose(n1 + n2, t);
}

static double
fishers_exact(size_t a, size_t b, size_t c, size_t d) {
  const size_t m = a + c; // sum of first column
  const size_t n = b + d; // sum of second column
  const size_t k = a + b; // sum of first row
  // ADS: want more extreme than "observed"
  const double observed = log_hyper_g(a, m, n, k);
  double p = 0.0;
  for (size_t i = (n > k ? 0ul : k - n); i <= std::min(k, m); ++i) {
    const double curr = log_hyper_g(i, m, n, k);
    if (curr <= observed)
      p = log_sum_log(p, curr);
  }
  return exp(p);
}

static size_t
state_pair_to_index(const string &s, const size_t idx) {
  assert(idx < s.length() - 1);
  const char a = s[idx];
  if (a == 'C') {
    const char b = s[idx+1];
    if (b == 'C') return 0;
    if (b == 'T') return 1;
    return 4;
  }
  if (a == 'T') {
    const char b = s[idx+1];
    if (b == 'C') return 2;
    if (b == 'T') return 3;
    return 4;
  }
  return 4;
}

template <class T>
struct PairStateCounter {
  T CC;
  T CT;
  T TC;
  T TT;

  double score() const {
    return (CC*TT > CT*TC) ?
      fishers_exact(CC, CT, TC, TT) : fishers_exact(CT, CC, TT, TC);
  }
  double total() const {return CC + CT + TC + TT;}

  string tostring() const {
    return toa(CC) + '\t' + toa(CT) + '\t' + toa(TC) + '\t' + toa(TT);
  }

  void increment(const size_t state) {
    if (state == 0) ++CC;
    else if (state == 1) ++CT;
    else if (state == 2) ++TC;
    else if (state == 3) ++TT;
  }
};


template <typename T> void
fit_states(const epiread &er, vector<PairStateCounter<T> > &counts) {
  for (size_t i = 0; i < er.length() - 1; ++i) {
    const size_t pos = er.pos + i;
    assert(pos < counts.size());
    const size_t curr_state = state_pair_to_index(er.seq, i);
    counts[pos].increment(curr_state);
  }
}

static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  size_t cpg_idx = 0;
  size_t nuc_idx = 0;
  const auto lim = end(s) - (s.size() == 0 ? 0 : 1);
  for (auto itr = begin(s); itr != lim; ++itr, ++nuc_idx)
    if (*itr == 'C' && *(itr + 1) == 'G')
      cpgs[cpg_idx++] = nuc_idx;
}

/* The "sites" in the convert_coordinates function represent pairs of
   consecutive CpG sites. So there would be one fewer of them than the
   total CpGs in a chromosome. */
static void
convert_coordinates(const string &chrom, vector<MSite> &sites) {
  unordered_map<size_t, size_t> cpgs;
  collect_cpgs(chrom, cpgs);
  for (size_t i = 0; i < sites.size(); ++i) {
    const auto pos_itr = cpgs.find(sites[i].pos);
    if (pos_itr == end(cpgs))
      throw runtime_error("failed converting site:\n" + sites[i].tostring());
    sites[i].pos = pos_itr->second;
  }
}


template <typename T>
void
add_cytosine(const string &chrom_name, const size_t start_cpg,
             vector<PairStateCounter<T>> &counts,
             vector<MSite> &cytosines) {
  std::ostringstream s;
  s << counts[start_cpg].score() << "\t"
    << counts[start_cpg].total() << "\t"
    << counts[start_cpg].tostring();
  const string name(s.str());
  cytosines.push_back(MSite(chrom_name, start_cpg, '+', name, 0, 0));
}


template<typename T>
void
process_chrom(const string &chrom_name, const vector<epiread> &epireads,
              vector<MSite> &cytosines, vector<PairStateCounter<T>> &counts) {
  for (size_t i = 0; i < epireads.size(); ++i)
    fit_states(epireads[i], counts);
  for (size_t i = 0; i < counts.size(); ++i)
    add_cytosine(chrom_name, i, counts, cytosines);
}


static void
update_chroms_seen(const string &chrom_name,
                   unordered_set<string> &chroms_seen) {
  const auto chr_itr = chroms_seen.find(chrom_name);
  if (chr_itr != end(chroms_seen))
    throw runtime_error("chroms out of order: " + chrom_name);
  chroms_seen.insert(chrom_name);
}


static void
verify_chroms_available(const string &chrom_name,
                        unordered_map<string, size_t> &chrom_lookup) {
  const auto chr_itr = chrom_lookup.find(chrom_name);
  if (chr_itr == end(chrom_lookup))
    throw runtime_error("chrom not found: " + chrom_name);
}


int
main_allelicmeth(int argc, const char **argv) {

  try {

    static const string description =
      "computes probability of allele-specific \
       methylation at each tuple of CpGs";

    static const string fasta_suffix = "fa";
    bool VERBOSE = false;

    string outfile;
    string chroms_dir;
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<epireads>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory",
                      true, chroms_dir);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string epi_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    vector<string> chrom_names;
    vector<string> chroms;
    read_fasta_file_short_names(chroms_dir, chrom_names, chroms);
    for (auto &&i: chroms)
      transform(begin(i), end(i), begin(i),
                [](const char c) { return std::toupper(c); });

    // lookup to map chrom names to chrom sequences
    unordered_map<string, size_t> chrom_lookup;
    for (size_t i = 0; i < chrom_names.size(); ++i)
      chrom_lookup.insert(make_pair(chrom_names[i], i));

    unordered_map<string, size_t> chrom_sizes;
    for (size_t i = 0; i < chrom_names.size(); ++i) {
      size_t cpg_count = 0;
      for (size_t j = 0; j < chroms[i].size() - 1; ++j)
        cpg_count += (chroms[i][j] == 'C' && chroms[i][j+1] == 'G');
      chrom_sizes.insert(make_pair(chrom_names[i], cpg_count));
    }

    if (VERBOSE)
      cerr << "number of chromosomes: " << chrom_sizes.size() << endl;

    std::ifstream in(epi_file);
    if (!in)
      throw runtime_error("cannot open input file: " + epi_file);

    std::ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());

    unordered_set<string> chroms_seen;
    string chrom;
    epiread er;
    vector<epiread> epireads;
    while (in >> er) {
      if (er.chr != chrom) {
        update_chroms_seen(er.chr, chroms_seen);
        verify_chroms_available(er.chr, chrom_lookup);

        if (VERBOSE)
          cerr << "[processing " << er.chr << "]" << endl;

        if (!chrom.empty()) {
          vector<PairStateCounter<uint32_t>> counts(chrom_sizes[chrom] - 1);
          vector<MSite> cytosines;
          process_chrom(chrom, epireads, cytosines, counts);
          const size_t chrom_idx = chrom_lookup[chrom];
          convert_coordinates(chroms[chrom_idx], cytosines);
          for (size_t i = 0; i < cytosines.size()-1; ++i) {
            out << cytosines[i].chrom << "\t"
                << cytosines[i].pos << "\t+\tCpG\t"
                << cytosines[i].context << endl;
          }
        }
        epireads.clear();
      }
      epireads.push_back(er);
      chrom.swap(er.chr);
    }
    if (!chrom.empty()) {
      vector<PairStateCounter<uint32_t>> counts(chrom_sizes[chrom] - 1);
      vector<MSite> cytosines;
      process_chrom(chrom, epireads, cytosines, counts);
      const size_t chrom_idx = chrom_lookup[chrom];
      convert_coordinates(chroms[chrom_idx], cytosines);
      for (size_t i = 0; i < cytosines.size() - 1; ++i) {
        out << cytosines[i].chrom << "\t"
            << cytosines[i].pos << "\t+\tCpG\t"
            << cytosines[i].context << endl;
      }
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
