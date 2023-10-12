/* amrfinder: program for resolving epialleles in a sliding window
 * along a chromosome.
 *
 * Copyright (C) 2014-2022 University of Southern California and
 *                         Andrew D. Smith and Benjamin E. Decato
 *
 * Authors: Fang Fang and Benjamin E. Decato and Andrew D. Smith
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
#include <atomic>

#include <bamxx.hpp>
#include <omp.h>

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "EpireadStats.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_map;
using std::runtime_error;
using std::max;
using std::min;
using std::to_string;
using std::atomic_ulong;
using std::pair;
using std::remove_if;
using std::toupper;
using std::move;
using std::size;

using bamxx::bgzf_file;

#include <bamxx.hpp>
#include <sstream>

static inline bamxx::bgzf_file &
read_epiread(bamxx::bgzf_file &f, epiread &er) {
  std::string line;
  if (getline(f, line)) er = epiread(line);
  return f;
}

static inline bool
validate_epiread_bgzf_file(const string &filename) {
  constexpr size_t max_lines_to_validate = 10000;
  bgzf_file in(filename, "r");
  if (!in) throw std::runtime_error("failed to open file: " + filename);

  string c, s, other;
  size_t p = 0;

  size_t n_lines = 0;
  string line;
  while (getline(in, line) && n_lines++ < max_lines_to_validate) {
    std::istringstream iss(line);
    if (!(iss >> c >> p >> s) || iss >> other) return false;
  }
  return true;
}

using epi_r = small_epiread;

/* merges amrs within some pre-defined distance */
static void
merge_amrs(const uint64_t gap_limit, vector<GenomicRegion> &amrs) {
  auto j = begin(amrs);
  for (auto &a : amrs)
    // check distance between two amrs is greater than gap limit
    if (j->same_chrom(a) && j->get_end() + gap_limit >= a.get_start()) {
      j->set_end(a.get_end());
      j->set_score(min(a.get_score(), j->get_score()));
    }
    else *(++j) = a;

  amrs.erase(++j, cend(amrs));
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////   CODE FOR RANDOMIZING THE READS TO GET EXPECTED NUMBER OF
/////   IDENTIFIED AMRS
/////

// static void
// set_read_states(vector<vector<char> > &state_counts, vector<epi_r> &reads) {
//   for (size_t i = 0; i < reads.size(); ++i) {
//     const size_t offset = reads[i].pos;
//     for (size_t j = 0; j < reads[i].length(); ++j) {
//       reads[i].seq[j] = state_counts[offset + j].back();
//       state_counts[offset + j].pop_back();
//     }
//   }
// }

// static void
// get_state_counts(const vector<epi_r> &reads, const size_t total_cpgs,
//               vector<vector<char> > &state_counts) {

//   state_counts = vector<vector<char> >(total_cpgs);
//   for (size_t i = 0; i < reads.size(); ++i) {
//     const size_t offset = reads[i].pos;
//     for (size_t j = 0; j < reads[i].length(); ++j)
//       state_counts[offset + j].push_back(reads[i].seq[j]);
//   }
//   for (size_t i = 0; i < state_counts.size(); ++i)
//     random_shuffle(state_counts[i].begin(), state_counts[i].end());
// }

// static void
// randomize_read_states(vector<epi_r> &reads) {
//   /**/srand(time(0) + getpid());
//   const size_t total_cpgs = get_n_cpgs(reads);
//   vector<vector<char> > state_counts;
//   get_state_counts(reads, total_cpgs, state_counts);
//   set_read_states(state_counts, reads);
// }

// Code for converting between cpg and base pair coordinates below
static inline vector<uint32_t>
collect_cpgs(const string &s) {
  vector<uint32_t> cpgs;
  for (uint32_t i = 0; i < size(s) - 1; ++i)
    if (s[i] == 'C' && s[i + 1] == 'G') cpgs.push_back(i);
  return cpgs;
}

template<typename T> static inline bool
convert_coordinates(const vector<T> &cpgs, GenomicRegion &region) {
  if (size(cpgs) <= region.get_end()) return false;
  region.set_start(cpgs[region.get_start()]);
  region.set_end(cpgs[region.get_end()]);
  return true;
}

static bool
convert_coordinates(const string &genome_file, vector<GenomicRegion> &amrs) {

  vector<string> c_name, c_seq;
  read_fasta_file_short_names(genome_file, c_name, c_seq);
  const size_t n_chroms = size(c_seq);

  for (auto &s : c_seq)
    for (auto &c : s) c = toupper(c);

  unordered_map<string, string> chrom_lookup;
  for (auto i = 0u; i < n_chroms; ++i)
    chrom_lookup.emplace(move(c_name[i]), move(c_seq[i]));

  vector<uint32_t> cpgs;
  string chrom_name;
  for (auto &a : amrs) {
    if (a.get_chrom() != chrom_name) {
      chrom_name = a.get_chrom();
      auto c_itr = chrom_lookup.find(chrom_name);
      if (c_itr == end(chrom_lookup)) return false;
      cpgs = collect_cpgs(c_itr->second);
    }
    if (!convert_coordinates(cpgs, a)) return false;
  }
  return true;
}

// make sure the current read is shortened to only include positions
// within the relevant window
static inline epi_r
clip_read(const size_t start_pos, const size_t end_pos, epi_r r) {
  if (r.pos < start_pos) {
    r.seq = r.seq.substr(start_pos - r.pos);
    r.pos = start_pos;
  }
  if (r.end() > end_pos)
    r.seq.resize(end_pos - r.pos);
  return r;
}

static vector<epi_r>
get_current_epireads(const vector<epi_r> &epireads, const size_t max_len,
                     const size_t cpg_window, const size_t start_pos,
                     size_t &read_id) {
  const auto n_epi = size(epireads);

  while (read_id < n_epi && epireads[read_id].pos + max_len <= start_pos)
    ++read_id;

  vector<epi_r> current_epireads;
  const auto end_pos = start_pos + cpg_window;
  for (auto i = read_id; i < n_epi && epireads[i].pos < end_pos; ++i)
    if (epireads[i].end() > start_pos)
      current_epireads.push_back(clip_read(start_pos, end_pos, epireads[i]));

  return current_epireads;
}

static inline uint64_t
total_states(const vector<epi_r> &epireads) {
  uint64_t tot = 0;
  for (auto &e : epireads) tot += e.length();
  return tot;
}

static inline void
add_amr(const string &chrom_name, const size_t start_cpg,
        const size_t cpg_window, const vector<epi_r> &reads, const double score,
        vector<GenomicRegion> &amrs) {
  const auto end_cpg = start_cpg + cpg_window - 1;
  const auto rds = to_string(size(reads));
  amrs.emplace_back(chrom_name, start_cpg, end_cpg, rds, score, '+');
}

static inline uint64_t
get_n_cpgs(const std::vector<epi_r> &reads) {
  uint64_t n_cpgs = 0;
  for (auto &r : reads) n_cpgs = max(n_cpgs, static_cast<uint64_t>(r.end()));
  return n_cpgs;
}

template<typename T> static inline vector<pair<T, T>>
get_block_bounds(const T start_pos, const T end_pos, const T block_size) {
  vector<pair<T, T>> blocks;
  auto block_start = start_pos;
  while (block_start < end_pos) {
    const auto block_end = min({block_start + block_size, end_pos});
    blocks.emplace_back(block_start, block_end);
    block_start = block_end;
  }
  return blocks;
}

static size_t
process_chrom(const bool verbose, const uint32_t n_threads,
              const size_t min_obs_per_cpg, const size_t window_size,
              const EpireadStats &epistat, const string &chrom_name,
              const vector<epi_r> &epireads, vector<GenomicRegion> &amrs) {
  constexpr auto blocks_per_thread = 8u;

  auto max_epiread_len = 0u;
  for (auto &e : epireads)
    max_epiread_len = max(max_epiread_len, e.length());
  const size_t min_obs_per_window = window_size*min_obs_per_cpg;

  const size_t n_cpgs = get_n_cpgs(epireads);
  if (verbose)
    cerr << "processing " << chrom_name << " "
         << "[reads: " << size(epireads) << "] "
         << "[cpgs: " << n_cpgs << "]" << endl;

  const auto n_blocks = n_threads*blocks_per_thread;

  const uint64_t lim = n_cpgs - window_size + 1;
  const auto blocks = get_block_bounds(0ul, lim, lim/n_blocks);

  atomic_ulong windows_tested = 0;

  vector<vector<GenomicRegion>> all_amrs;

#pragma omp parallel for
  for (const auto &b : blocks) {
    vector<GenomicRegion> curr_amrs;
    uint64_t start_idx = 0;
    uint64_t windows_tested_block = 0;
    for (auto i = b.first; i < b.second && start_idx < size(epireads); ++i) {
      auto curr_epireads = get_current_epireads(epireads, max_epiread_len,
                                                window_size, i, start_idx);
      if (total_states(curr_epireads) >= min_obs_per_window) {
        bool is_significant = false;
        const auto score = epistat.test_asm(curr_epireads, is_significant);
        if (is_significant)
          add_amr(chrom_name, i, window_size, curr_epireads, score, curr_amrs);
        ++windows_tested_block;
      }
    }
#pragma omp critical
    {
      all_amrs.push_back(move(curr_amrs));
    }
    windows_tested += windows_tested_block;
  }

  auto total_amrs = 0u;
  for (const auto &a : all_amrs)
    total_amrs += size(a);

  amrs.reserve(total_amrs);
  for (auto &v : all_amrs)
    for (auto &a : v) amrs.push_back(move(a));

  return windows_tested;
}

struct rename_amr {
  void operator()(GenomicRegion &r) {
    static constexpr auto label = "AMR";
    r.set_name(label + to_string(idx++) + ":" + r.get_name());
  }
  uint32_t idx{};
};

int
main_amrfinder(int argc, const char **argv) {
  try {
    const string description =
      "identify regions of allele-specific methylation";

    static const string fasta_suffix = "fa";

    bool verbose = false;
    size_t n_threads = 1;

    string outfile;
    string genome_file;

    size_t max_itr = 10;
    size_t window_size = 10;
    size_t gap_limit = 1000;

    double high_prob = 0.75, low_prob = 0.25;
    size_t min_obs_per_cpg = 4;
    double critical_value = 0.01;

    // ADS: below, for when the time comes
    // auto eng = std::default_random_engine(rng_seed);
    // std::shuffle(begin(things), end(things), eng);

    // bool RANDOMIZE_READS = false;
    bool use_bic = false;
    bool use_fdr = true;
    bool apply_correction = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<epireads>");
    opt_parse.add_opt("output", 'o', "output file", true, outfile);
    opt_parse.add_opt("chrom", 'c', "reference genome fasta file", true,
                      genome_file);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("window", 'w', "size of sliding window (in CpGs)", false,
                      window_size);
    opt_parse.add_opt("min-cov", 'm', "min coverage per cpg to test windows",
                      false, min_obs_per_cpg);
    opt_parse.add_opt("gap", 'g', "min allowed gap between amrs (in bp)", false,
                      gap_limit);
    opt_parse.add_opt("crit", 'C', "critical p-value cutoff", false,
                      critical_value);
    opt_parse.add_opt("pvals", 'h', "adjusts p-values using Hochberg step-up",
                      false, apply_correction);
    opt_parse.add_opt("nofdr", 'f', "omits FDR procedure", false, use_fdr);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.set_show_defaults();
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
    const string reads_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/
    omp_set_num_threads(n_threads);

    if (!validate_epiread_bgzf_file(reads_file))
      throw runtime_error("invalid states file: " + reads_file);

    if (verbose)
      cerr << "AMR TESTING OPTIONS: "
           << "[test=" << (use_bic ? "BIC" : "LRT") << "] "
           << "[iterations=" << max_itr << "]" << endl;

    const EpireadStats epistat(low_prob, high_prob, critical_value, max_itr,
                               use_bic);

    bamxx::bam_tpool tp(n_threads);
    bgzf_file in(reads_file, "r");
    if (!in) throw runtime_error("failed to open input file: " + reads_file);
    if (n_threads > 1) tp.set_io(in);

    vector<GenomicRegion> amrs;
    size_t windows_tested = 0;
    epiread er;
    vector<epi_r> epireads;
    string prev_chrom, curr_chrom, tmp_states;

    while (read_epiread(in, er)) {
      if (!epireads.empty() && er.chr != prev_chrom) {
        windows_tested += process_chrom(verbose, n_threads, min_obs_per_cpg,
                                        window_size, epistat, prev_chrom, epireads, amrs);
        epireads.clear();
      }
      swap(prev_chrom, er.chr);
      epireads.emplace_back(er.pos, er.seq);
    }
    if (!epireads.empty())
      windows_tested += process_chrom(verbose, n_threads, min_obs_per_cpg, window_size,
                                      epistat, prev_chrom, epireads, amrs);

    // because the threads might not add these in order
    sort(begin(amrs), end(amrs));

    // otherwise the names are not in chrom order
    for_each(begin(amrs), end(amrs), rename_amr());

    if (verbose)
      cerr << "========= POST PROCESSING =========" << endl;

    const size_t windows_accepted = amrs.size();
    if (!amrs.empty()) {
      /*
        ADS: there are several "steps" below that are done
        independently, but for which the order matters. It is not
        clear to me that this is the right order, after observing the
        behavior of the program on some newer much larger data
        sets. The steps:

        (1) Get the p-values
        (2) Get the fdr cutoff
        (3) Apply the correction if needed
        (4) collapse the AMRs that overlap in CpG sites
        (5) convert the AMRs back to base pairs
        (6) merge the AMRs based on a "gap limit"
        (7) remove the AMRs based on a score (p-value or fdr)
        (8) eliminate what remains by size using half the "gap limit"

        The score in the above is taken as the extreme from step (4),
        which interacts poorly with subsequent steps. Also, the "half
        gap limit" for eliminating AMRs at the end is not justified
        anywhere.
       */

      // get all the pvals... (but they might be BIC scores)
      vector<double> pvals(size(amrs));
      transform(cbegin(amrs), cend(amrs), begin(pvals),
                [](const GenomicRegion &r) {return r.get_score();});

      // if the pvalues are not BIC scores, get an FDR cutoff
      // corresponding to the given critical value
      const double fdr_cutoff =
        use_bic
        ? 0.0
        : smithlab::get_fdr_cutoff(windows_tested, pvals, critical_value);

      // if we are not using BIC, and if corrected p-values are
      // requested, then adjust the p-values
      if (!use_bic && apply_correction) {
        smithlab::correct_pvals(windows_tested, pvals);
        for (auto i = 0u; i < size(pvals); ++i) amrs[i].set_score(pvals[i]);
      }

      // ADS: not sure it's a good idea in this next collapse function
      // to keep the lowest among all p-values for merged regions
      merge_amrs(1, amrs);
      const size_t collapsed_amrs = amrs.size();

      if (!convert_coordinates(genome_file, amrs)) {
        cerr << "failed converting coordinates" << endl;
        return EXIT_FAILURE;
      }

      // ADS: merge AMRs if they are sufficiently close; the distance
      // should not be too large, and a distance of 1000 bp is likely
      // too large.
      merge_amrs(gap_limit, amrs);
      const size_t merged_amrs = amrs.size();

      // if BIC was not requested, then eliminate AMRs based on either
      // the p-value cutoff, or with the FDR-based cutoff, if it was
      // requested.
      if (!use_bic) {
        const auto cutoff = (apply_correction || !use_fdr) ?
          critical_value : fdr_cutoff;
        amrs.erase(remove_if(begin(amrs), end(amrs),
          [cutoff](const GenomicRegion &r) { return r.get_score() >= cutoff; }),
                   cend(amrs));
      }

      const size_t amrs_passing_fdr = amrs.size();

      // ADS: eliminating AMRs based on their size makes sense, but
      // not if that size is tied to the gap between AMRs we would
      // merge. There is no symmetry for these.
      const auto min_size = gap_limit / 2;
      amrs.erase(remove_if(begin(amrs), end(amrs),
        [min_size](const GenomicRegion &r) { return r.get_width() < min_size; }),
                 cend(amrs));

      if (verbose) {
        cerr << "WINDOWS TESTED: " << windows_tested << endl
             << "WINDOWS ACCEPTED: " << windows_accepted << endl
             << "COLLAPSED WINDOWS: " << collapsed_amrs << endl
             << "MERGED WINDOWS: " << merged_amrs << endl;
        if (use_fdr)
          cerr << "FDR CUTOFF: " << fdr_cutoff << endl
               << "WINDOWS PASSING FDR: " << amrs_passing_fdr << endl;
        cerr << "AMRS (WINDOWS PASSING MINIMUM SIZE): " << amrs.size() << endl;
      }
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile);
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
      copy(begin(amrs), end(amrs),
           std::ostream_iterator<GenomicRegion>(out, "\n"));
    }
    else if (verbose)
      cerr << "no AMRs found" << endl;
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
