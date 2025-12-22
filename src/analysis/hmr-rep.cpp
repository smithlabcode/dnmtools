/* Copyright (C) 2019-2025 Andrew D. Smith, Song Qiang, Jenny Qu
 *
 * This is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "Interval6.hpp"
#include "MSite.hpp"
#include "TwoStateHMM.hpp"

#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-avoid-magic-numbers,*-narrowing-conversions)

static Interval6
as_gen_rgn(const MSite &s) {
  return Interval6(s.chrom, s.pos, s.pos + 1, std::string{}, 0.0, '+');
}

static double
get_stepup_cutoff(std::vector<double> scores, const double cutoff) {
  if (cutoff <= 0)
    return std::numeric_limits<double>::max();
  else if (cutoff > 1)
    return std::numeric_limits<double>::min();

  const size_t n = scores.size();
  std::sort(std::begin(scores), std::end(scores));
  size_t i = 1;
  while (i < n && scores[i - 1] < (cutoff * i) / n)
    ++i;
  return scores[i - 1];
}

template <class T>
[[nodiscard]] T
pair_sum(const std::pair<T, T> &t) {
  return t.first + t.second;
}

static void
get_domain_scores_rep(
  const std::vector<bool> &state_ids,
  const std::vector<std::vector<std::pair<double, double>>> &meth,
  const std::vector<size_t> &reset_points, std::vector<double> &scores) {
  const size_t n_reps = meth.size();
  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0.0;
  for (size_t i = 0; i < state_ids.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        scores.push_back(score);
        score = 0;
      }
      ++reset_idx;
    }
    if (state_ids[i]) {
      in_domain = true;
      for (size_t r = 0; r < n_reps; ++r)
        if (pair_sum(meth[r][i]) >= 1)
          score += 1.0 - meth[r][i].first / pair_sum(meth[r][i]);
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }
}

static void
build_domains(const std::vector<MSite> &cpgs,
              const std::vector<size_t> &reset_points,
              const std::vector<bool> &state_ids,
              std::vector<Interval6> &domains) {
  size_t n_cpgs = 0, n_domains = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  for (size_t i = 0; i < state_ids.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        domains.back().stop = prev_end;
        domains.back().score = n_cpgs;
        n_cpgs = 0;
      }
      ++reset_idx;
    }
    if (state_ids[i]) {
      if (!in_domain) {
        in_domain = true;
        domains.push_back(as_gen_rgn(cpgs[i]));
        domains.back().name = "HYPO" + std::to_string(n_domains++);
      }
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      domains.back().stop = prev_end;
      domains.back().score = n_cpgs;
      n_cpgs = 0;
    }
    prev_end = cpgs[i].pos + 1;
  }
}

template <class T, class U>
static void
separate_regions(const bool VERBOSE, const size_t desert_size,
                 std::vector<MSite> &cpgs, std::vector<std::vector<T>> &meth,
                 std::vector<std::vector<U>> &reads,
                 std::vector<size_t> &reset_points) {
  if (VERBOSE)
    std::cerr << "[separating by cpg desert]\n";

  // eliminate the zero-read cpg sites if no coverage in any replicates
  const size_t n_reps = meth.size();

  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    bool has_data = true;
    for (size_t rep_idx = 0; rep_idx < n_reps; ++rep_idx)
      has_data = (has_data && reads[rep_idx][i] > 0);

    if (has_data) {
      cpgs[j] = cpgs[i];
      for (size_t r = 0; r < n_reps; ++r) {
        meth[r][j] = meth[r][i];
        reads[r][j] = reads[r][i];
      }
      ++j;
    }
  }

  cpgs.erase(std::begin(cpgs) + j, std::end(cpgs));
  for (size_t r = 0; r < n_reps; ++r) {
    meth[r].erase(std::begin(meth[r]) + j, std::end(meth[r]));
    reads[r].erase(std::begin(reads[r]) + j, std::end(reads[r]));
  }

  // segregate cpgs
  size_t prev_cpg = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom)
                          ? cpgs[i].pos - prev_cpg
                          : std::numeric_limits<size_t>::max();
    if (dist > desert_size)
      reset_points.push_back(i);
    prev_cpg = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());

  if (VERBOSE)
    std::cerr << "[cpgs retained: " << cpgs.size() << "]\n"
              << "[deserts removed: " << reset_points.size() - 2 << "]\n";
}

static void
shuffle_cpgs_rep(const size_t rng_seed, const TwoStateHMM &hmm,
                 std::vector<std::vector<std::pair<double, double>>>
                   meth,  // cppcheck-suppress passedByValue
                 const std::vector<size_t> &reset_points,
                 const double f_to_b_trans, const double b_to_f_trans,
                 const std::vector<double> &fg_alpha,
                 const std::vector<double> &fg_beta,
                 const std::vector<double> &bg_alpha,
                 const std::vector<double> &bg_beta,
                 std::vector<double> &domain_scores) {
  auto eng = std::default_random_engine(rng_seed);
  for (size_t r = 0; r < meth.size(); ++r)
    std::shuffle(std::begin(meth[r]), std::end(meth[r]), eng);

  std::vector<bool> state_ids;
  std::vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, f_to_b_trans, b_to_f_trans,
                        fg_alpha, fg_beta, bg_alpha, bg_beta, state_ids,
                        scores);

  get_domain_scores_rep(state_ids, meth, reset_points, domain_scores);
  std::sort(std::begin(domain_scores), std::end(domain_scores));
}

static void
assign_p_values(const std::vector<double> &random_scores,
                const std::vector<double> &observed_scores,
                std::vector<double> &p_values) {
  const double n_randoms = std::max(random_scores.size(), 1ul);
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back(
      (std::end(random_scores) - upper_bound(std::begin(random_scores),
                                             std::end(random_scores),
                                             observed_scores[i])) /
      n_randoms);
}

static void
read_params_file(const bool VERBOSE, const std::string &params_file,
                 double &fg_alpha, double &fg_beta, double &bg_alpha,
                 double &bg_beta, double &f_to_b_trans, double &b_to_f_trans,
                 double &fdr_cutoff) {
  std::string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw std::runtime_error("failed to parse params file: " + params_file);

  in >> jnk >> fg_alpha >> jnk >> fg_beta >> jnk >> bg_alpha >> jnk >>
    bg_beta >> jnk >> f_to_b_trans >> jnk >> b_to_f_trans >> jnk >> fdr_cutoff;

  if (VERBOSE)
    std::cerr << "read in params from " << params_file << '\n'
              << "FG_ALPHA\t" << fg_alpha << '\n'
              << "FG_BETA\t" << fg_beta << '\n'
              << "BG_ALPHA\t" << bg_alpha << '\n'
              << "BG_BETA\t" << bg_beta << '\n'
              << "F_B\t" << f_to_b_trans << '\n'
              << "B_F\t" << b_to_f_trans << '\n'
              << "FDR_CUTOFF\t" << fdr_cutoff << '\n';
}

static void
write_params_file(const std::string &outfile,
                  const std::vector<double> &fg_alpha,
                  const std::vector<double> &fg_beta,
                  const std::vector<double> &bg_alpha,
                  const std::vector<double> &bg_beta, const double f_to_b_trans,
                  const double b_to_f_trans, const double fdr_cutoff) {
  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile);
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  for (size_t r = 0; r < fg_alpha.size(); ++r)
    out << "FG_ALPHA_" << r + 1 << '\t' << fg_alpha[r] << '\t' << "FG_BETA_"
        << r + 1 << '\t' << fg_beta[r] << '\t' << "BG_ALPHA_" << r + 1 << '\t'
        << bg_alpha[r] << '\t' << "BG_BETA_" << r + 1 << '\t' << bg_beta[r]
        << '\n';

  out << "F_B\t" << f_to_b_trans << '\n'
      << "B_F\t" << b_to_f_trans << '\n'
      << "FDR_CUTOFF\t" << fdr_cutoff << '\n';
}

static void
load_cpgs(const std::string &cpgs_file, std::vector<MSite> &cpgs,
          std::vector<std::pair<double, double>> &meth,
          std::vector<uint32_t> &reads) {
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  while (read_site(in, the_site)) {
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    meth.push_back(std::make_pair(the_site.n_meth(), the_site.n_unmeth()));
  }
}

static void
check_consistent_sites(const std::string &expected_filename,
                       const std::vector<MSite> &expected,
                       const std::string &observed_filename,
                       const std::vector<MSite> &observed) {
  if (expected.size() != observed.size()) {
    std::ostringstream err_msg;
    err_msg << "inconsistent number of sites\n"
            << "file=" << expected_filename << ","
            << "sites=" << expected.size() << '\n'
            << "file=" << observed_filename << ","
            << "sites=" << observed.size() << '\n';
    throw std::runtime_error(err_msg.str());
  }
}

template <class InputIterator>
[[nodiscard]] double
get_mean(InputIterator first, InputIterator last) {
  return std::accumulate(first, last, 0.0) / std::distance(first, last);
}

static std::vector<std::string>
split_comma(const std::string &orig) {
  std::string tmp(orig);
  std::replace(std::begin(tmp), std::end(tmp), ',', ' ');
  std::istringstream iss(tmp);
  std::vector<std::string> parts;
  while (iss >> tmp)
    parts.push_back(tmp);
  return parts;
}

int
main_hmr_rep(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    std::string outfile;
    std::string hypo_post_outfile;
    std::string meth_post_outfile;

    size_t desert_size = 1000;
    size_t max_iterations = 10;
    size_t rng_seed = 408;

    // run mode flags
    bool VERBOSE = false;

    const double tolerance = 1e-10;  // corrections for small values

    std::string params_in_files;
    std::string params_out_file;

    const std::string description =
      R"(Identify HMRs in a set of replicate methylomes. Methylation must be
      provided in the methcounts format (chrom, position, strand, context,
      methylation, reads). See the methcounts documentation for details
      for details. This program assumes only data at CpG sites and that
      strands are collapsed so only the positive site appears in the file.)";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description,
                           "<methcount-file-1> <methcount-file-2> ...");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn cpgs with reads in HMR",
                      false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("post-hypo", '\0',
                      "output file for single-CpG posteiror "
                      "hypomethylation probability (default: NULL)",
                      false, hypo_post_outfile);
    opt_parse.add_opt("post-meth", '\0',
                      "output file for single-CpG posteiror "
                      "methylation probability (default: NULL)",
                      false, meth_post_outfile);
    opt_parse.add_opt("params-in", 'P',
                      "HMM parameter files for "
                      "individual methylomes (separated with comma)",
                      false, params_in_files);
    opt_parse.add_opt("params-out", 'p', "write HMM parameters to this file",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed", false, rng_seed);
    opt_parse.set_show_defaults();
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
    if (leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::vector<std::string> cpgs_files(leftover_args);
    /****************** END COMMAND LINE OPTIONS *****************/

    for (const auto &filename : cpgs_files)
      // cppcheck-suppress useStlAlgorithm
      if (!is_msite_file(filename))
        throw std::runtime_error("malformed counts file: " + filename);

    std::vector<std::string> params_in_file;
    if (!params_in_files.empty()) {
      params_in_file = split_comma(params_in_files);
      assert(cpgs_files.size() == params_in_file.size());
    }

    const size_t n_reps = cpgs_files.size();

    std::vector<MSite> cpgs;
    std::vector<std::vector<std::pair<double, double>>> meth(n_reps);
    std::vector<std::vector<uint32_t>> reads(n_reps);

    if (VERBOSE)
      std::cerr << "[reading methylation levels]\n";
    for (size_t i = 0; i < n_reps; ++i) {
      if (VERBOSE)
        std::cerr << "[filename=" << cpgs_files[i] << "]\n";
      std::vector<MSite> curr_rep;
      load_cpgs(cpgs_files[i], curr_rep, meth[i], reads[i]);
      if (VERBOSE)
        std::cerr << "[total_cpgs=" << curr_rep.size() << "]\n"
                  << "[mean_coverage="
                  << get_mean(std::begin(reads[i]), std::end(reads[i]))
                  << "]\n";
      if (i > 0)
        check_consistent_sites(cpgs_files[0], cpgs, cpgs_files[i], curr_rep);
      else
        swap(cpgs, curr_rep);
    }

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    std::vector<size_t> reset_points;
    separate_regions(VERBOSE, desert_size, cpgs, meth, reads, reset_points);

    /****************** initalize params *****************/
    const TwoStateHMM hmm(tolerance, max_iterations, VERBOSE);
    std::vector<double> fg_alpha(n_reps), fg_beta(n_reps);
    std::vector<double> bg_alpha(n_reps), bg_beta(n_reps);
    double fdr_cutoff = std::numeric_limits<double>::max();

    double f_to_b_trans = 0.25;
    double b_to_f_trans = 0.25;

    if (!params_in_file.empty()) {  // read parameters files
      double fdr_cutoff_rep{};      // ignore this cutoff
      for (size_t i = 0; i < n_reps; ++i)
        read_params_file(VERBOSE, params_in_file[i], fg_alpha[i], fg_beta[i],
                         bg_alpha[i], bg_beta[i], f_to_b_trans, b_to_f_trans,
                         fdr_cutoff_rep);
      max_iterations = 0;
    }
    else {
      for (size_t i = 0; i < n_reps; ++i) {
        // JQU: there are many 0s in reads[r], but the parameter start
        // points don't need to be perfect
        const double mean_reads =
          get_mean(std::begin(reads[i]), std::end(reads[i]));
        fg_alpha[i] = 0.33 * mean_reads;
        fg_beta[i] = 0.67 * mean_reads;
        bg_alpha[i] = 0.67 * mean_reads;
        bg_beta[i] = 0.33 * mean_reads;
      }
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, f_to_b_trans, b_to_f_trans,
                            fg_alpha, fg_beta, bg_alpha, bg_beta);

    std::vector<bool> state_ids;
    std::vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, f_to_b_trans, b_to_f_trans,
                          fg_alpha, fg_beta, bg_alpha, bg_beta, state_ids,
                          posteriors);

    std::vector<double> domain_scores;
    get_domain_scores_rep(state_ids, meth, reset_points, domain_scores);

    std::vector<double> random_scores;
    shuffle_cpgs_rep(rng_seed, hmm, meth, reset_points, f_to_b_trans,
                     b_to_f_trans, fg_alpha, fg_beta, bg_alpha, bg_beta,
                     random_scores);

    std::vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (fdr_cutoff == std::numeric_limits<double>::max())
      fdr_cutoff = get_stepup_cutoff(p_values, 0.01);

    std::vector<Interval6> domains;
    build_domains(cpgs, reset_points, state_ids, domains);

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (size_t i = 0; i < domains.size(); ++i)
      if (p_values[i] < fdr_cutoff) {
        domains[i].name = "HYPO" + std::to_string(good_hmr_count++);
        out << domains[i] << '\n';
      }

    // write all the hmm parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
                        f_to_b_trans, b_to_f_trans, fdr_cutoff);

    if (!hypo_post_outfile.empty()) {
      if (VERBOSE)
        std::cerr << "[writing=" << hypo_post_outfile << "]\n";
      std::ofstream out_post(hypo_post_outfile);
      for (size_t i = 0; i < cpgs.size(); ++i) {
        size_t m_reads = 0, u_reads = 0;
        for (size_t j = 0; j < n_reps; ++j) {
          m_reads += meth[j][i].first;
          u_reads += meth[j][i].second;
        }
        Interval6 cpg(as_gen_rgn(cpgs[i]));
        cpg.name =
          "CpG:" + std::to_string(m_reads) + ":" + std::to_string(u_reads);
        cpg.score = posteriors[i];
        out_post << cpg << '\n';
      }
    }

    if (!meth_post_outfile.empty()) {
      std::ofstream out_post(meth_post_outfile);
      if (VERBOSE)
        std::cerr << "[writing=" << meth_post_outfile << "]\n";
      for (size_t i = 0; i < cpgs.size(); ++i) {
        size_t m_reads = 0, u_reads = 0;
        for (size_t j = 0; j < n_reps; ++j) {
          m_reads += meth[j][i].first;
          u_reads += meth[j][i].second;
        }
        Interval6 cpg(as_gen_rgn(cpgs[i]));
        cpg.name =
          "CpG:" + std::to_string(m_reads) + ":" + std::to_string(u_reads);
        cpg.score = 1.0 - posteriors[i];
        out_post << cpg << '\n';
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
