/* amrfinder: program for resolving epialleles in a sliding window
 * along a chromosome.
 *
 * Copyright (C) 2014-2023 University of Southern California and
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

#include <bamxx.hpp>
#include <omp.h>

#include <atomic>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "EpireadStats.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::max;
using std::min;
using std::move;
using std::ofstream;
using std::ostream_iterator;
using std::pair;
using std::remove_if;
using std::runtime_error;
using std::size;
using std::string;
using std::to_string;
using std::vector;

using bamxx::bgzf_file;

#include <bamxx.hpp>

#include <sstream>

struct amr_summary {
  amr_summary(const vector<GenomicRegion> &amrs) {
    amr_count = size(amrs);
    amr_total_size = accumulate(cbegin(amrs), cend(amrs), 0u,
                                [](const uint64_t t, const GenomicRegion &p) {
                                  return t + p.get_width();
                                });
    amr_mean_size = static_cast<double>(amr_total_size) /
                    std::max(amr_count, static_cast<uint64_t>(1));
  }

  // amr_count is the number of identified AMRs, which are the merged
  // AMRs that are found to be significant when tested as a single
  // interval
  uint64_t amr_count{};
  // total_amr_size is the sum of the sizes of the identified AMRs
  uint64_t amr_total_size{};
  // mean_amr_size is the mean size of the identified AMRs
  double amr_mean_size{};

  auto tostring() -> string {
    static constexpr uint32_t current_precision = 2;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(current_precision);
    oss << "amr_count: " << amr_count << '\n'
        << "amr_total_size: " << amr_total_size << '\n'
        << "amr_mean_size: " << amr_mean_size;
    return oss.str();
  }
};

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
  if (!in) throw runtime_error("failed to open file: " + filename);

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
  for (const auto &a : amrs)
    // check distance between two amrs is greater than gap limit
    if (j->same_chrom(a) && j->get_end() + gap_limit >= a.get_start()) {
      j->set_end(a.get_end());
      j->set_score(min(a.get_score(), j->get_score()));
    }
    else
      *(++j) = a;

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

static vector<pair<uint32_t, uint32_t>>
get_chrom_partition(const vector<GenomicRegion> &r) {
  if (r.empty()) return {};
  vector<pair<uint32_t, uint32_t>> parts;
  string prev_chrom{r.front().get_chrom()};
  uint32_t prev_idx = 0u;
  for (auto i = 0u; i < size(r); ++i)
    if (r[i].get_chrom() != prev_chrom) {
      parts.emplace_back(prev_idx, i);
      prev_idx = i;
      prev_chrom = r[i].get_chrom();
    }
  parts.emplace_back(prev_idx, size(r));
  return parts;
}

static bool
convert_coordinates(const string &genome_file, vector<GenomicRegion> &amrs) {
  vector<string> c_name, c_seq;
  read_fasta_file_short_names(genome_file, c_name, c_seq);
  const size_t n_chroms = size(c_seq);

  for (auto &s : c_seq)
    for (auto &c : s) c = std::toupper(c);

  std::unordered_map<string, string> chrom_lookup;
  for (auto i = 0u; i < n_chroms; ++i)
    chrom_lookup.emplace(std::move(c_name[i]), std::move(c_seq[i]));

  vector<pair<uint32_t, uint32_t>> chrom_parts = get_chrom_partition(amrs);
  std::atomic_uint32_t conv_failure = 0;

#pragma omp parallel for
  for (const auto &p : chrom_parts) {
    const string chrom_name = amrs[p.first].get_chrom();
    auto c_itr = chrom_lookup.find(chrom_name);
    if (c_itr == end(chrom_lookup))
      conv_failure++;
    else {
      vector<uint32_t> cpgs = collect_cpgs(c_itr->second);
      for (uint32_t i = p.first; i < p.second; ++i)
        conv_failure += !convert_coordinates(cpgs, amrs[i]);
    }
  }
  return conv_failure == 0;
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
                     uint64_t &read_id) {

  // assert(is_sorted(cbegin(epireads), cend(epireads),
  //                  [](const epi_r &a, const epi_r &b) {
  //                    return a.pos < b.pos;
  //                  }));

  const auto n_epi = size(epireads);

  while (read_id < size(epireads) &&
         epireads[read_id].pos + max_len <= start_pos)
    ++read_id;

  const auto end_pos = start_pos + cpg_window;

  vector<epi_r> current_epireads;
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
  const auto blocks = get_block_bounds(static_cast<uint64_t>(0),
                                       lim, lim/n_blocks);

  std::atomic_ulong windows_tested = 0;

  vector<vector<GenomicRegion>> all_amrs;

#pragma omp parallel for
  for (const auto &b : blocks) {
    vector<GenomicRegion> curr_amrs;
    uint64_t start_idx = 0;
    uint64_t windows_tested_block = 0;
    for (auto i = b.first; i < b.second && start_idx < size(epireads); ++i) {
      // ADS: below, we could do binary search, but this does not seem
      // to be a bottleneck
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
      all_amrs.push_back(std::move(curr_amrs));
    }
    windows_tested += windows_tested_block;
  }

  auto total_amrs = 0u;
  for (const auto &a : all_amrs)
    total_amrs += size(a);

  amrs.reserve(total_amrs);
  for (auto &v : all_amrs)
    for (auto &a : v) amrs.push_back(std::move(a));

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

    bool verbose = false;
    uint32_t n_threads = 1;

    string outfile;
    string genome_file;
    string summary_file;

    uint32_t max_itr = 10;
    uint32_t window_size = 10;
    uint64_t gap_limit = 1000;

    double high_prob = 0.75, low_prob = 0.25;
    uint32_t min_obs_per_cpg = 4;
    double critical_value = 0.01;

    // ADS: below, for when the time comes
    // auto eng = std::default_random_engine(rng_seed);
    // std::shuffle(begin(things), end(things), eng);

    // bool RANDOMIZE_READS = false;
    bool use_bic = false;
    bool use_fdr = true;
    bool apply_correction = false;
    bool correct_for_read_count = true;

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
    opt_parse.add_opt("nordc", 'r', "turn off read count correction",
                      false, correct_for_read_count);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
    opt_parse.add_opt("summary", 'S', "write summary output here", false,
                      summary_file);
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

    const auto epistat =
      EpireadStats(low_prob, high_prob, critical_value, max_itr, use_bic);

    EpireadStats::set_read_count_correction(correct_for_read_count);

    // bamxx::bam_tpool tp(n_threads);
    bgzf_file in(reads_file, "r");
    if (!in) throw runtime_error("failed to open input file: " + reads_file);
    // if (n_threads > 1) tp.set_io(in);

    vector<GenomicRegion> amrs;

    // windows_tested is the number of sliding windows in the
    // methylome that were tested for a signal of significant
    // allele-specific methylation
    uint64_t windows_tested = 0;
    epiread er;
    vector<epi_r> epireads;
    string prev_chrom, curr_chrom, tmp_states;

    while (read_epiread(in, er)) {
      if (!epireads.empty() && er.chr != prev_chrom) {
        windows_tested +=
          process_chrom(verbose, n_threads, min_obs_per_cpg, window_size,
                        epistat, prev_chrom, epireads, amrs);
        epireads.clear();
      }
      swap(prev_chrom, er.chr);
      epireads.emplace_back(er.pos, er.seq);
    }
    if (!epireads.empty())
      windows_tested +=
        process_chrom(verbose, n_threads, min_obs_per_cpg, window_size, epistat,
                      prev_chrom, epireads, amrs);

    // because the threads might not add these in order
    sort(begin(amrs), end(amrs));

    // otherwise the names are not in chrom order
    for_each(begin(amrs), end(amrs), rename_amr());

    if (verbose)
      cerr << "========= POST PROCESSING =========" << endl;

    // windows_accepted is the number of sliding windows in the
    // methylome that were found to have a significant signal of
    // allele-specific methylation
    const auto windows_accepted = size(amrs);
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

      // collapsed_amrs is the number of intervals of consecutive CpG
      // sites that are found to reside in a window among those
      // accepted as significant
      const auto n_collapsed_amrs = size(amrs);

      if (!convert_coordinates(genome_file, amrs)) {
        cerr << "failed converting coordinates" << endl;
        return EXIT_FAILURE;
      }

      // ADS: merge AMRs if they are sufficiently close; the distance
      // should not be too large, and a distance of 1000 bp is likely
      // too large.
      merge_amrs(gap_limit, amrs);

      // merged_amrs are the number of intervals that result from
      // merging any collapsed amrs that have a distance of less than
      // gap_limit/2 from each other
      const auto n_merged_amrs = size(amrs);

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

      const auto amrs_passing_fdr = size(amrs);

      // ADS: eliminating AMRs based on their size makes sense, but
      // not if that size is tied to the gap between AMRs we would
      // merge. There is no symmetry for these.
      const auto min_size = gap_limit / 2;
      amrs.erase(remove_if(begin(amrs), end(amrs),
        [min_size](const GenomicRegion &r) { return r.get_width() < min_size; }),
                 cend(amrs));

      if (verbose) {
        cerr << "windows_tested: " << windows_tested << '\n'
             << "windows_accepted: " << windows_accepted << '\n'
             << "collapsed_amrs: " << n_collapsed_amrs << '\n'
             << "merged_amrs: " << n_merged_amrs << '\n';
        if (use_fdr)
          cerr << "fdr_cutoff: " << fdr_cutoff << '\n'
               << "amrs_passing_fdr: " << amrs_passing_fdr << '\n';
        cerr << amr_summary(amrs).tostring() << '\n';
      }
      ofstream of;
      if (!outfile.empty()) of.open(outfile);
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
      if (!out) runtime_error("failed to open output file: " + outfile);
      copy(cbegin(amrs), cend(amrs),
           ostream_iterator<GenomicRegion>(out, "\n"));
    }
    if (!summary_file.empty()) {
      ofstream summary_out(summary_file);
      if (!summary_out) throw runtime_error("failed to open: " + summary_file);
      summary_out << amr_summary(amrs).tostring() << endl;
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
