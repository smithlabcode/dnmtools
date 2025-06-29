/* amrscores: produce the sliding window scores used by amrfinder.
 *
 * Copyright (C) 2025 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include "EpireadStats.hpp"
#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <array>
#include <atomic>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

[[nodiscard]]
std::tuple<std::vector<std::string>, std::vector<std::string>>
read_genome(const std::string &filename) {
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("cannot open input file " + filename);

  std::vector<std::string> names;
  std::vector<std::string> sequences;

  std::string line;
  while (getline(in, line)) {
    if (line[0] == '>') {
      const auto first_space = line.find_first_of(" \t", 1);
      if (first_space == std::string::npos)
        names.push_back(line.substr(1));
      else
        names.emplace_back(std::cbegin(line) + 1,
                           std::cbegin(line) + line.find_first_of(" \t", 1));
      sequences.emplace_back();
    }
    else
      sequences.back() += line;
  }
  return std::tuple{std::move(names), std::move(sequences)};
}

static inline bamxx::bgzf_file &
read_epiread(bamxx::bgzf_file &f, epiread &er) {
  kstring_t line = KS_INITIALIZE;
  if (getline(f, line))
    er = epiread(ks_str(&line), ks_len(&line));
  ks_free(&line);
  return f;
}

static inline bool
validate_epiread_bgzf_file(const std::string &filename) {
  constexpr std::size_t max_lines_to_validate = 10000;
  bamxx::bgzf_file in(filename, "r");
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

// Code for converting between cpg and base pair coordinates below
[[nodiscard]] static inline std::vector<std::uint32_t>
get_cpg_positions(const std::string &s) {
  const auto is_cpg = [](const auto a, const auto b) {
    return (a == 'C' || a == 'c') && (b == 'g' || b == 'G');
  };
  auto cpg_count = 0u;
  for (std::uint32_t i = 0; i + 1 < std::size(s); ++i)
    cpg_count += is_cpg(s[i], s[i + 1]);
  std::vector<std::uint32_t> cpgs(cpg_count);
  std::uint32_t j = 0;
  for (std::uint32_t i = 0; i + 1 < std::size(s); ++i)
    if (is_cpg(s[i], s[i + 1]))
      cpgs[j++] = i;
  return cpgs;
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
                     const std::size_t start_pos, std::size_t &read_idx) {
  const auto n_epi = std::size(epireads);

  if (std::size(epireads) <= read_idx)
    return {};

  while (read_idx < std::size(epireads) &&
         epireads[read_idx].pos + max_len <= start_pos)
    ++read_idx;

  const auto end_pos = start_pos + cpg_window;

  std::vector<epi_r> current_epireads;
  for (auto i = read_idx; i < n_epi && epireads[i].pos < end_pos; ++i)
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

static inline std::size_t
get_n_cpgs(const std::vector<epi_r> &reads) {
  std::size_t n_cpgs = 0;
  for (auto &r : reads)
    n_cpgs = std::max(n_cpgs, static_cast<std::size_t>(r.end()));
  return n_cpgs;
}

template <typename T>
static inline std::vector<std::pair<T, T>>
get_block_bounds(const T n_elements, const T n_blocks) {
  std::vector<std::pair<T, T>> blocks;
  const T q = n_elements / n_blocks;
  const T r = n_elements - q * n_blocks;
  T block_start{};
  for (std::size_t i = 0; i < n_blocks; ++i) {
    const auto sz = (i < r) ? q + 1 : q;
    const auto block_end = block_start + sz;
    blocks.emplace_back(block_start, block_end);
    block_start = block_end;
  }
  return blocks;
}

static std::vector<float>
process_chrom(const bool verbose, const std::uint32_t n_threads,
              const std::size_t n_cpgs, const std::size_t min_obs_per_cpg,
              const std::size_t window_size, const EpireadStats &epistat,
              const std::string &chrom_name,
              const std::vector<epi_r> &epireads) {
  constexpr auto blocks_per_thread = 8u;

  auto max_epiread_len = 0ul;
  for (auto &e : epireads)
    max_epiread_len = std::max(max_epiread_len, std::size(e));

  // minimum observations per window to be able to do the stats calcs
  const std::size_t min_obs = window_size * min_obs_per_cpg;
  const auto n_epireads = std::size(epireads);

  const std::size_t n_cpgs_from_epireads = get_n_cpgs(epireads);
  if (n_cpgs < n_cpgs_from_epireads)
    throw std::runtime_error(
      "cpgs from genome: " + std::to_string(n_cpgs) +
      ", cpgs from reads: " + std::to_string(n_cpgs_from_epireads));

  if (verbose)
    std::cerr << "processing " << chrom_name << " "
              << "[reads: " << n_epireads << "] "
              << "[cpgs: " << n_cpgs << "]" << std::endl;
  if (n_cpgs < window_size)
    return {};

  const std::size_t n_blocks = n_threads * blocks_per_thread;
  const auto blocks = get_block_bounds(n_cpgs, n_blocks);
  const auto blocks_beg = std::cbegin(blocks);
  assert(std::size(blocks) == n_blocks);

  std::vector<std::vector<float>> scores(n_threads);

  std::vector<std::thread> threads;
  for (auto thread_id = 0ul; thread_id < n_threads; ++thread_id) {
    const auto b_beg =
      blocks_beg + std::min(thread_id * blocks_per_thread, n_blocks);
    const auto b_end =
      blocks_beg + std::min((thread_id + 1) * blocks_per_thread, n_blocks);
    threads.emplace_back([&, thread_id, b_beg, b_end] {
      std::vector<float> curr_scores;
      for (auto b = b_beg; b != b_end; ++b) {
        // ADS: below, we could do binary search for the first read_idx of the
        // set of blocks, but this does not seem to be a bottleneck
        std::size_t read_idx = 0;
        for (auto start_pos = b->first; start_pos < b->second; ++start_pos) {
          if (n_epireads == read_idx) {
            curr_scores.push_back(0.0);
            continue;
          }
          auto curr_epireads = get_current_epireads(
            epireads, max_epiread_len, window_size, start_pos, read_idx);
          const auto n_states = total_states(curr_epireads);
          bool is_sig{};
          const auto score =
            n_states < min_obs ? 0.0 : epistat.test_asm(curr_epireads, is_sig);
          curr_scores.push_back(score);
        }
      }
      scores[thread_id] = std::move(curr_scores);
    });
  }
  for (auto &thread : threads)
    thread.join();

  const auto n_scores = std::transform_reduce(
    std::cbegin(scores), std::cend(scores), 0ul, std::plus<>(),
    [&](const auto &x) { return std::size(x); });

  std::vector<float> all_scores(n_scores);
  auto scr_itr = std::begin(all_scores);
  for (const auto &s : scores)
    scr_itr = std::copy(std::cbegin(s), std::cend(s), scr_itr);

  assert(scr_itr == std::cend(all_scores));
  assert(std::size(all_scores) == n_cpgs);

  // make scores associated with the CpG site in the middle of the window
  if (std::size(all_scores) > window_size) {
    std::copy_backward(std::cbegin(all_scores),
                       std::cend(all_scores) - window_size / 2,
                       std::end(all_scores));
    std::fill_n(std::begin(all_scores), window_size / 2, 0.0);
  }
  return all_scores;
}

static inline void
write_output(std::ofstream &out, const std::string &chrom_name,
             const std::vector<std::uint32_t> &cpgs,
             const std::vector<float> &scores) {
  static constexpr auto buf_size = 128ul;
  static const auto excpt = std::runtime_error("failed to write output line");
  static std::array<char, buf_size> buf{};
  const int m = std::sprintf(buf.data(), "%s\t", chrom_name.data());
  if (m <= 0)
    throw excpt;
  auto j = std::cbegin(cpgs);
  auto cursor = buf.data() + m;
  for (const auto scr : scores) {
    const int n = std::sprintf(cursor, "%d\t%.6g\n", *j++, scr);
    if (n <= 0)
      throw excpt;
    out.write(buf.data(), m + n);
  }
}

int
main_amrscores(int argc, const char **argv) {
  try {
    const std::string description =
      "produce the sliding window scores used by amrfinder";

    bool verbose = false;
    std::uint32_t n_threads = 1;

    std::string outfile;
    std::string genome_file;

    std::uint32_t max_itr = 10;
    std::uint32_t window_size = 10;

    double high_prob = 0.75;
    double low_prob = 0.25;
    std::uint32_t min_obs_per_cpg = 4;

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
    opt_parse.add_opt("pvals", 'h', "adjusts p-values using Hochberg step-up",
                      false, apply_correction);
    opt_parse.add_opt("nofdr", 'f', "omits FDR procedure", false, use_fdr);
    opt_parse.add_opt("nordc", 'r', "turn off read count correction", false,
                      correct_for_read_count);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
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

    std::vector<std::vector<std::uint32_t>> cpgs;
    auto [names, seqs] = read_genome(genome_file);
    for (auto i = 0u; i < std::size(names); ++i) {
      cpgs.emplace_back(get_cpg_positions(seqs[i]));
      std::string().swap(seqs[i]);
    }

    std::unordered_map<std::string, std::uint32_t> name_to_idx;
    for (auto i = 0u; i < std::size(names); ++i)
      name_to_idx.emplace(names[i], i);

    // clang-format off
    const EpireadStats epistat{
      low_prob,
      high_prob,
      0.001, // critical_value,
      max_itr,
      use_bic,
      correct_for_read_count,
    };
    // clang-format on

    bamxx::bam_tpool tp(n_threads);
    bamxx::bgzf_file in(reads_file, "r");
    if (!in)
      throw std::runtime_error("failed to open input file: " + reads_file);
    if (n_threads > 1 && in.is_bgzf())
      tp.set_io(in);

    std::ofstream out(outfile);
    if (!out)
      std::runtime_error("failed to open output file: " + outfile);

    epiread er;
    std::vector<epi_r> epireads;
    std::string prev_chrom;
    std::vector<std::string> chrom_names;

    while (read_epiread(in, er)) {
      if (er.chr != prev_chrom) {
        if (!epireads.empty()) {
          const auto chrom_itr = name_to_idx.find(prev_chrom);
          if (chrom_itr == std::cend(name_to_idx))
            throw std::runtime_error("failed to find chrom: " + prev_chrom);
          const auto scores = process_chrom(
            verbose, n_threads, std::size(cpgs[chrom_itr->second]),
            min_obs_per_cpg, window_size, epistat, prev_chrom, epireads);
          assert(std::size(scores) == std::size(cpgs[chrom_itr->second]));
          std::vector<epi_r>().swap(epireads);
          write_output(out, prev_chrom, cpgs[chrom_itr->second], scores);
        }
        chrom_names.emplace_back(er.chr);
        prev_chrom = er.chr;
      }
      epireads.emplace_back(er.pos, er.seq);
    }
    if (!epireads.empty()) {
      const auto chrom_itr = name_to_idx.find(prev_chrom);
      if (chrom_itr == std::cend(name_to_idx))
        throw std::runtime_error("failed to find chrom: " + prev_chrom);
      const auto scores = process_chrom(
        verbose, n_threads, std::size(cpgs[chrom_itr->second]), min_obs_per_cpg,
        window_size, epistat, prev_chrom, epireads);
      assert(std::size(scores) == std::size(cpgs[chrom_itr->second]));
      std::vector<epi_r>().swap(epireads);
      write_output(out, prev_chrom, cpgs[chrom_itr->second], scores);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
