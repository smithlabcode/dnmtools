/* amrfinder: program for resolving epialleles in a sliding window
 * along a chromosome.
 *
 * Copyright (C) 2014-2025 University of Southern California and
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

#include "EpireadStats.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <atomic>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using bamxx::bgzf_file;

struct amr_summary {
  amr_summary(const std::vector<GenomicRegion> &amrs) {
    amr_count = std::size(amrs);
    amr_total_size =
      accumulate(std::cbegin(amrs), std::cend(amrs), 0ul,
                 [](const std::size_t t, const GenomicRegion &p) {
                   return t + p.get_width();
                 });
    amr_mean_size = static_cast<double>(amr_total_size) /
                    std::max(amr_count, static_cast<std::size_t>(1));
  }

  // amr_count is the number of identified AMRs, which are the merged
  // AMRs that are found to be significant when tested as a single
  // interval
  std::size_t amr_count{};
  // total_amr_size is the sum of the sizes of the identified AMRs
  std::size_t amr_total_size{};
  // mean_amr_size is the mean size of the identified AMRs
  double amr_mean_size{};

  auto
  tostring() -> std::string {
    static constexpr std::uint32_t current_precision = 2;
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
  if (getline(f, line))
    er = epiread(line);
  return f;
}

static inline bool
validate_epiread_bgzf_file(const std::string &filename) {
  constexpr std::size_t max_lines_to_validate = 10000;
  bgzf_file in(filename, "r");
  if (!in)
    throw std::runtime_error("failed to open file: " + filename);

  std::string c, s, other;
  std::size_t p = 0;

  std::size_t n_lines = 0;
  std::string line;
  while (getline(in, line) && n_lines++ < max_lines_to_validate) {
    std::istringstream iss(line);
    if (!(iss >> c >> p >> s) || iss >> other)
      return false;
  }
  return true;
}

using epi_r = small_epiread;

/* merges amrs within some pre-defined distance */
static void
merge_amrs(const std::size_t gap_limit, std::vector<GenomicRegion> &amrs) {
  auto j = std::begin(amrs);
  for (const auto &a : amrs)
    // check distance between two amrs is greater than gap limit
    if (j->same_chrom(a) && j->get_end() + gap_limit >= a.get_start()) {
      j->set_end(a.get_end());
      j->set_score(std::min(a.get_score(), j->get_score()));
    }
    else
      *(++j) = a;

  amrs.erase(++j, std::cend(amrs));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////   CODE FOR RANDOMIZING THE READS TO GET EXPECTED NUMBER OF
/////   IDENTIFIED AMRS
/////

// static void
// set_read_states(std::vector<std::vector<char> > &state_counts,
// std::vector<epi_r> &reads) {
//   for (std::size_t i = 0; i < reads.size(); ++i) {
//     const std::size_t offset = reads[i].pos;
//     for (std::size_t j = 0; j < reads[i].length(); ++j) {
//       reads[i].seq[j] = state_counts[offset + j].back();
//       state_counts[offset + j].pop_back();
//     }
//   }
// }

// static void
// get_state_counts(const std::vector<epi_r> &reads, const std::size_t
// total_cpgs,
//               std::vector<std::vector<char> > &state_counts) {

//   state_counts = std::vector<std::vector<char> >(total_cpgs);
//   for (std::size_t i = 0; i < reads.size(); ++i) {
//     const std::size_t offset = reads[i].pos;
//     for (std::size_t j = 0; j < reads[i].length(); ++j)
//       state_counts[offset + j].push_back(reads[i].seq[j]);
//   }
//   for (std::size_t i = 0; i < state_counts.size(); ++i)
//     random_shuffle(state_counts[i].std::begin(), state_counts[i].end());
// }

// static void
// randomize_read_states(std::vector<epi_r> &reads) {
//   /**/srand(time(0) + getpid());
//   const std::size_t total_cpgs = get_n_cpgs(reads);
//   std::vector<std::vector<char> > state_counts;
//   get_state_counts(reads, total_cpgs, state_counts);
//   set_read_states(state_counts, reads);
// }

// Code for converting between cpg and base pair coordinates below
static inline std::vector<std::uint32_t>
collect_cpgs(const std::string &s) {
  std::vector<std::uint32_t> cpgs;
  for (std::uint32_t i = 0; i < std::size(s) - 1; ++i)
    if (s[i] == 'C' && s[i + 1] == 'G')
      cpgs.push_back(i);
  return cpgs;
}

template <typename T>
static inline bool
convert_coordinates(const std::vector<T> &cpgs, GenomicRegion &region) {
  if (std::size(cpgs) <= region.get_end())
    return false;
  region.set_start(cpgs[region.get_start()]);
  region.set_end(cpgs[region.get_end()]);
  return true;
}

static std::vector<std::pair<std::uint32_t, std::uint32_t>>
get_chrom_partition(const std::vector<GenomicRegion> &r) {
  if (r.empty())
    return {};
  std::vector<std::pair<std::uint32_t, std::uint32_t>> parts;
  std::string prev_chrom{r.front().get_chrom()};
  std::uint32_t prev_idx = 0u;
  for (auto i = 0u; i < std::size(r); ++i)
    if (r[i].get_chrom() != prev_chrom) {
      parts.emplace_back(prev_idx, i);
      prev_idx = i;
      prev_chrom = r[i].get_chrom();
    }
  parts.emplace_back(prev_idx, std::size(r));
  return parts;
}

[[nodiscard]] static bool
convert_coordinates(const std::size_t n_threads, const std::string &genome_file,
                    std::vector<GenomicRegion> &amrs) {
  std::vector<std::string> c_name, c_seq;
  read_fasta_file_short_names(genome_file, c_name, c_seq);
  const std::size_t n_chroms = std::size(c_seq);

  for (auto &s : c_seq)
    for (auto &c : s)
      c = std::toupper(c);

  std::unordered_map<std::string, std::string> chrom_lookup;
  for (auto i = 0u; i < n_chroms; ++i)
    chrom_lookup.emplace(std::move(c_name[i]), std::move(c_seq[i]));

  const auto chrom_parts = get_chrom_partition(amrs);
  const auto n_parts = std::size(chrom_parts);
  const auto parts_beg = std::cbegin(chrom_parts);
  const std::uint32_t n_per = (n_parts + n_threads - 1) / n_threads;

  std::atomic_uint32_t conv_failure = 0;

  std::vector<std::thread> threads;
  for (auto i = 0ul; i < n_threads; ++i) {
    const auto p_beg = parts_beg + i * n_per;
    const auto p_end = parts_beg + std::min((i + 1) * n_per, n_parts);
    threads.emplace_back([&, p_beg, p_end] {
      for (auto p = p_beg; p != p_end; ++p) {
        const std::string chrom_name = amrs[p->first].get_chrom();
        auto c_itr = chrom_lookup.find(chrom_name);
        if (c_itr == std::end(chrom_lookup))
          conv_failure++;
        else {
          std::vector<std::uint32_t> cpgs = collect_cpgs(c_itr->second);
          for (auto j = p->first; j < p->second; ++j)
            conv_failure += !convert_coordinates(cpgs, amrs[j]);
        }
      }
    });
  }
  for (auto &thread : threads)
    thread.join();
  return conv_failure == 0;
}

// make sure the current read is shortened to only include positions
// within the relevant window
static inline epi_r
clip_read(const std::size_t start_pos, const std::size_t end_pos, epi_r r) {
  if (r.pos < start_pos) {
    r.seq = r.seq.substr(start_pos - r.pos);
    r.pos = start_pos;
  }
  if (r.end() > end_pos)
    r.seq.resize(end_pos - r.pos);
  return r;
}

static std::vector<epi_r>
get_current_epireads(const std::vector<epi_r> &epireads,
                     const std::size_t max_len, const std::size_t cpg_window,
                     const std::size_t start_pos, std::size_t &read_id) {
  // assert(is_sorted(std::cbegin(epireads), std::cend(epireads),
  //                  [](const epi_r &a, const epi_r &b) {
  //                    return a.pos < b.pos;
  //                  }));
  const auto n_epi = std::size(epireads);

  while (read_id < std::size(epireads) &&
         epireads[read_id].pos + max_len <= start_pos)
    ++read_id;

  const auto end_pos = start_pos + cpg_window;

  std::vector<epi_r> current_epireads;
  for (auto i = read_id; i < n_epi && epireads[i].pos < end_pos; ++i)
    if (epireads[i].end() > start_pos)
      current_epireads.push_back(clip_read(start_pos, end_pos, epireads[i]));

  return current_epireads;
}

static inline std::size_t
total_states(const std::vector<epi_r> &epireads) {
  std::size_t tot = 0;
  for (auto &e : epireads)
    tot += e.length();
  return tot;
}

static inline void
add_amr(const std::string &chrom_name, const std::size_t start_cpg,
        const std::size_t cpg_window, const std::vector<epi_r> &reads,
        const double score, std::vector<GenomicRegion> &amrs) {
  const auto end_cpg = start_cpg + cpg_window - 1;
  const auto rds = std::to_string(std::size(reads));
  amrs.emplace_back(chrom_name, start_cpg, end_cpg, rds, score, '+');
}

static inline std::size_t
get_n_cpgs(const std::vector<epi_r> &reads) {
  std::size_t n_cpgs = 0;
  for (auto &r : reads)
    n_cpgs = std::max(n_cpgs, static_cast<std::size_t>(r.end()));
  return n_cpgs;
}

template <typename T>
static inline std::vector<std::pair<T, T>>
get_block_bounds(const T start_pos, const T end_pos, const T block_size) {
  if (block_size == 0)
    return {{start_pos, end_pos}};
  std::vector<std::pair<T, T>> blocks;
  auto block_start = start_pos;
  while (block_start < end_pos) {
    const auto block_end = std::min({block_start + block_size, end_pos});
    blocks.emplace_back(block_start, block_end);
    block_start = block_end;
  }
  return blocks;
}

static std::size_t
process_chrom(const bool verbose, const std::uint32_t n_threads,
              const std::size_t min_obs_per_cpg, const std::size_t window_size,
              const EpireadStats &epistat, const std::string &chrom_name,
              const std::vector<epi_r> &epireads,
              std::vector<GenomicRegion> &amrs) {
  constexpr auto blocks_per_thread = 8u;

  auto max_epiread_len = 0u;
  for (auto &e : epireads)
    max_epiread_len = std::max(max_epiread_len, e.length());
  const std::size_t min_obs_per_window = window_size * min_obs_per_cpg;

  const std::size_t n_cpgs = get_n_cpgs(epireads);
  if (verbose)
    std::cerr << "processing " << chrom_name << " "
              << "[reads: " << std::size(epireads) << "] "
              << "[cpgs: " << n_cpgs << "]" << std::endl;

  const auto n_blocks = n_threads * blocks_per_thread;

  const std::uint32_t lim =
    n_cpgs >= window_size ? n_cpgs - window_size + 1 : 0;
  const auto blocks = get_block_bounds(0u, lim, lim / n_blocks);
  const auto blocks_beg = std::cbegin(blocks);
  const std::uint32_t n_per = (n_blocks + n_threads - 1) / n_threads;

  const auto n_epireads = std::size(epireads);
  std::atomic_ulong windows_tested = 0;

  std::vector<std::vector<GenomicRegion>> all_amrs(n_threads);

  std::vector<std::thread> threads;
  for (auto i = 0u; i < n_threads; ++i) {
    const auto b_beg = blocks_beg + i * n_per;
    const auto b_end = blocks_beg + std::min((i + 1) * n_per, n_blocks);
    threads.emplace_back([&, i, b_beg, b_end] {
      std::vector<GenomicRegion> curr_amrs;
      std::size_t windows_tested_thread = 0;
      for (auto b = b_beg; b != b_end; ++b) {
        std::size_t start_idx = 0;
        for (auto j = b->first; j < b->second && start_idx < n_epireads; ++j) {
          // ADS: below, we could do binary search, but this does not seem
          // to be a bottleneck
          auto curr_epireads = get_current_epireads(epireads, max_epiread_len,
                                                    window_size, j, start_idx);
          if (total_states(curr_epireads) >= min_obs_per_window) {
            bool is_significant = false;
            const auto score = epistat.test_asm(curr_epireads, is_significant);
            if (is_significant)
              add_amr(chrom_name, j, window_size, curr_epireads, score,
                      curr_amrs);
            ++windows_tested_thread;
          }
        }
      }
      all_amrs[i] = std::move(curr_amrs);
      windows_tested += windows_tested_thread;
    });
  }
  for (auto &thread : threads)
    thread.join();

  auto total_amrs = 0u;
  for (const auto &a : all_amrs)
    total_amrs += std::size(a);

  amrs.reserve(total_amrs);
  for (auto &v : all_amrs)
    for (auto &a : v)
      amrs.push_back(std::move(a));

  return windows_tested;
}

struct rename_amr {
  void
  operator()(GenomicRegion &r) {
    static constexpr auto label = "AMR";
    r.set_name(label + std::to_string(idx++) + ":" + r.get_name());
  }
  std::uint32_t idx{};
};

int
main_amrfinder(int argc, const char **argv) {
  try {
    const std::string description =
      "identify regions of allele-specific methylation";

    bool verbose = false;
    std::uint32_t n_threads = 1;

    std::string outfile;
    std::string genome_file;
    std::string summary_file;

    std::uint32_t max_itr = 10;
    std::uint32_t window_size = 10;
    std::size_t gap_limit = 1000;

    double high_prob = 0.75;
    double low_prob = 0.25;
    std::uint32_t min_obs_per_cpg = 4;
    double critical_value = 0.01;

    // ADS: below, for when the time comes
    // auto eng = std::default_random_engine(rng_seed);
    // std::shuffle(std::begin(things), std::end(things), eng);

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
    opt_parse.add_opt("nordc", 'r', "turn off read count correction", false,
                      correct_for_read_count);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
    opt_parse.add_opt("summary", 'S', "write summary output here", false,
                      summary_file);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.set_show_defaults();
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
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string reads_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    n_threads = std::min(std::thread::hardware_concurrency(), n_threads);

    if (!validate_epiread_bgzf_file(reads_file))
      throw std::runtime_error("invalid states file: " + reads_file);

    if (verbose)
      std::cerr << "AMR TESTING OPTIONS: "
                << "[test=" << (use_bic ? "BIC" : "LRT") << "] "
                << "[iterations=" << max_itr << "]" << std::endl;

    const EpireadStats epistat{low_prob, high_prob, critical_value,
                               max_itr,  use_bic,   correct_for_read_count};

    bamxx::bam_tpool tp(n_threads);
    bgzf_file in(reads_file, "r");
    if (!in)
      throw std::runtime_error("failed to open input file: " + reads_file);
    if (n_threads > 1 && in.is_bgzf())
      tp.set_io(in);

    std::vector<GenomicRegion> amrs;

    // windows_tested is the number of sliding windows in the
    // methylome that were tested for a signal of significant
    // allele-specific methylation
    std::size_t windows_tested = 0;
    epiread er;
    std::vector<epi_r> epireads;
    std::string prev_chrom, curr_chrom, tmp_states;

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
    std::sort(std::begin(amrs), std::end(amrs));

    // otherwise the names are not in chrom order
    std::for_each(std::begin(amrs), std::end(amrs), rename_amr());

    if (verbose)
      std::cerr << "========= POST PROCESSING =========" << std::endl;

    // windows_accepted is the number of sliding windows in the
    // methylome that were found to have a significant signal of
    // allele-specific methylation
    const auto windows_accepted = std::size(amrs);
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
      std::vector<double> pvals(std::size(amrs));
      transform(std::cbegin(amrs), std::cend(amrs), std::begin(pvals),
                [](const GenomicRegion &r) { return r.get_score(); });

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
        for (auto i = 0u; i < std::size(pvals); ++i)
          amrs[i].set_score(pvals[i]);
      }

      // ADS: not sure it's a good idea in this next collapse function
      // to keep the lowest among all p-values for merged regions
      merge_amrs(1, amrs);

      // collapsed_amrs is the number of intervals of consecutive CpG
      // sites that are found to reside in a window among those
      // accepted as significant
      const auto n_collapsed_amrs = std::size(amrs);

      if (!convert_coordinates(n_threads, genome_file, amrs)) {
        std::cerr << "failed converting coordinates" << std::endl;
        return EXIT_FAILURE;
      }

      // ADS: merge AMRs if they are sufficiently close; the distance
      // should not be too large, and a distance of 1000 bp is likely
      // too large.
      merge_amrs(gap_limit, amrs);

      // merged_amrs are the number of intervals that result from
      // merging any collapsed amrs that have a distance of less than
      // gap_limit/2 from each other
      const auto n_merged_amrs = std::size(amrs);

      // if BIC was not requested, then eliminate AMRs based on either
      // the p-value cutoff, or with the FDR-based cutoff, if it was
      // requested.
      if (!use_bic) {
        const auto cutoff =
          (apply_correction || !use_fdr) ? critical_value : fdr_cutoff;
        amrs.erase(std::remove_if(std::begin(amrs), std::end(amrs),
                                  [cutoff](const GenomicRegion &r) {
                                    return r.get_score() >= cutoff;
                                  }),
                   std::cend(amrs));
      }

      const auto amrs_passing_fdr = std::size(amrs);

      // ADS: eliminating AMRs based on their size makes sense, but
      // not if that size is tied to the gap between AMRs we would
      // merge. There is no symmetry for these.
      const auto min_size = gap_limit / 2;
      amrs.erase(
        std::remove_if(std::begin(amrs), std::end(amrs),
                       [&](const auto &r) { return r.get_width() < min_size; }),
        std::cend(amrs));

      if (verbose) {
        std::cerr << "windows_tested: " << windows_tested << '\n'
                  << "windows_accepted: " << windows_accepted << '\n'
                  << "collapsed_amrs: " << n_collapsed_amrs << '\n'
                  << "merged_amrs: " << n_merged_amrs << '\n';
        if (use_fdr)
          std::cerr << "fdr_cutoff: " << fdr_cutoff << '\n'
                    << "amrs_passing_fdr: " << amrs_passing_fdr << '\n';
        std::cerr << amr_summary(amrs).tostring() << '\n';
      }
      std::ofstream out(outfile);
      if (!out)
        std::runtime_error("failed to open output file: " + outfile);
      std::copy(std::cbegin(amrs), std::cend(amrs),
                std::ostream_iterator<GenomicRegion>(out, "\n"));
    }
    if (!summary_file.empty()) {
      std::ofstream summary_out(summary_file);
      if (!summary_out)
        throw std::runtime_error("failed to open: " + summary_file);
      summary_out << amr_summary(amrs).tostring() << std::endl;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
