/* Copyright (C) 2009-2025 Song Qiang and Andrew D Smith
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
 *
 * You should have received a copy of the GNU General Public License along
 * with this software; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 */

[[maybe_unused]] static constexpr auto about = R"(
Identify regions of elevated methylation. Designed for methylomes of plants,
especially A. thaliana, and other organisms with a low background methylation
level, but with increased methylation associated with biological activity at a
specific regions of the genome.
)";

#include "BetaBin.hpp"
#include "Interval6.hpp"
#include "MSite.hpp"
#include "ThreeStateHMM.hpp"

#include "OptionParser.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-narrowing-conversions)

static Interval6
as_gen_rgn(const MSite &s) {
  return Interval6(s.chrom, s.pos, s.pos + 1, std::string{}, 0.0, '+');
}

static void
load_cpgs(const std::string &cpgs_file, std::vector<MSite> &cpgs,
          std::vector<std::pair<double, double>> &meth) {
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  std::string line;
  while (getline(in, line))
    cpgs.push_back(MSite(line));

  meth.resize(std::size(cpgs));
  for (std::size_t i = 0; i < std::size(cpgs); ++i)
    meth[i] = std::make_pair(cpgs[i].n_meth(), cpgs[i].n_unmeth());
}

template <class T>
static void
separate_regions(const std::size_t desert_size, std::vector<MSite> &cpgs,
                 std::vector<T> &meth, std::vector<std::size_t> &reset_points,
                 std::size_t &total_bases, std::size_t &bases_in_deserts) {
  // eliminate the zero-read cpgs
  std::size_t j = 0;
  for (std::size_t i = 0; i < std::size(cpgs); ++i)
    if (cpgs[i].n_reads > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      ++j;
    }
  cpgs.erase(std::begin(cpgs) + j, std::end(cpgs));
  meth.erase(std::begin(meth) + j, std::end(meth));

  total_bases = 0;
  bases_in_deserts = 0;
  std::size_t prev_pos = 0;
  for (std::size_t i = 0; i < std::size(cpgs); ++i) {
    const std::size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom)
                               ? cpgs[i].pos - prev_pos
                               : std::numeric_limits<std::size_t>::max();
    if (dist > desert_size) {
      reset_points.push_back(i);
      if (dist < std::numeric_limits<std::size_t>::max())
        bases_in_deserts += dist;
    }
    if (dist < std::numeric_limits<std::size_t>::max())
      total_bases += dist;

    prev_pos = cpgs[i].pos;
  }
  reset_points.push_back(std::size(cpgs));
}

// ADS (!!!) this function seems to not be working at all right now
static void
read_params_file(const std::string &params_file,
                 // betabin &hypo_emission,
                 // betabin &HYPER_emission,
                 // betabin &HYPO_emission,
                 std::vector<std::vector<double>> &trans) {
  std::ifstream in(params_file);
  if (!in)
    throw std::runtime_error("failed to read param file: " + params_file);

  std::string hypo_emission_str;
  getline(in, hypo_emission_str);

  std::string HYPER_emission_str;
  getline(in, HYPER_emission_str);

  std::string HYPO_emission_str;
  getline(in, HYPO_emission_str);

  trans.resize(3, std::vector<double>(3, 0.0));
  for (std::size_t i = 0; i < std::size(trans); ++i)
    for (std::size_t j = 0; j < std::size(trans[i]); ++j)
      in >> trans[i][j];
}

static void
write_params_file(const std::string &params_file, const betabin &hypo_emission,
                  const betabin &HYPER_emission, const betabin &HYPO_emission,
                  const std::vector<std::vector<double>> &trans) {
  std::ofstream out(params_file);
  out << hypo_emission.tostring() << '\n'
      << HYPER_emission.tostring() << '\n'
      << HYPO_emission.tostring() << '\n';
  for (const auto &i : trans) {
    std::copy(std::cbegin(i), std::cend(i),
              std::ostream_iterator<double>(out, "\t"));
    out << '\n';
  }
}

static void
build_domains(const std::vector<MSite> &cpgs,
              const std::vector<std::pair<double, double>> &meth,
              const std::vector<std::size_t> &reset_points,
              const std::vector<STATE_LABELS> &classes,
              std::vector<Interval6> &domains) {
  domains.clear();

  for (std::size_t i = 0; i + 1 < std::size(reset_points); ++i) {
    const std::size_t start = reset_points[i];
    const std::size_t end = reset_points[i + 1];

    Interval6 domain(as_gen_rgn(cpgs[start]));
    STATE_LABELS prev_state = classes[start];
    std::size_t n = 1;
    double meth_sum =
      meth[start].first / (meth[start].first + meth[start].second);

    for (std::size_t j = start + 1; j < end; ++j) {
      if ((prev_state == hypo && classes[j] == hypo) ||
          (prev_state != hypo && classes[j] != hypo)) {
        ++n;
        meth_sum += meth[j].first / (meth[j].first + meth[j].second);
      }
      else {
        domain.stop = cpgs[j - 1].pos + 1;
        const std::string label = (prev_state == hypo ? "hypo" : "hyper");
        domain.name = label + ":" + std::to_string(n);
        domain.score = meth_sum;
        domain.strand = '+';
        if (prev_state == HYPER || prev_state == HYPO)
          domains.push_back(domain);

        domain = Interval6(as_gen_rgn(cpgs[j]));
        n = 1;
        prev_state = classes[j];
        meth_sum = meth[j].first / (meth[j].first + meth[j].second);
      }
    }
    domain.stop = cpgs[end - 1].pos + 1;
    const std::string label = (prev_state == hypo ? "hypo" : "hyper");
    domain.name = label + ":" + std::to_string(n);
    domain.score = meth_sum;
    domain.strand = '+';
    if (prev_state == HYPER || prev_state == HYPO)
      domains.push_back(domain);
  }
}

static void
filter_domains(const double min_cumulative_meth,
               std::vector<Interval6> &domains) {
  std::size_t j = 0;
  for (std::size_t i = 0; i < std::size(domains); ++i)
    if (domains[i].score >= min_cumulative_meth)
      domains[j++] = domains[i];
  domains.erase(std::begin(domains) + j, std::end(domains));
}

static void
initialize_transitions(std::vector<std::vector<double>> &trans) {
  // "hypo" -> "hypo" is high, and goes only to "HYPER";
  // "HYPER" -> "HYPER" or "HYPO" equally;
  // "HYPO" -> "HYPER" usually; never to "hypo"

  // NOLINTBEGIN(*-avoid-magic-numbers)
  trans = {
    {0.990, 0.010, 0.000},
    {0.050, 0.475, 0.475},
    {0.000, 0.666, 0.334},
  };
  // NOLINTEND(*-avoid-magic-numbers)
}

int
main_hypermr(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    // NOLINTBEGIN(*-avoid-magic-numbers)
    std::size_t desert_size = 1000;
    std::size_t max_iterations = 10;
    // corrections for small values (not parameters):
    double tolerance = 1e-10;
    double min_cumulative_meth = 4.0;

    // NOLINTEND(*-avoid-magic-numbers)

    // run mode flags
    bool VERBOSE = false;
    bool USE_VITERBI_DECODING = false;

    std::string outfile;
    std::string scores_file;
    std::string params_in_file;
    std::string params_out_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           about, "<counts-file>");
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
    opt_parse.add_opt("scores", 's', "output file for posterior scores", false,
                      scores_file);
    opt_parse.add_opt("tolerance", 't', "numerical tolerance", false,
                      tolerance);
    opt_parse.add_opt("desert", 'd', "desert size", false, desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("viterbi", 'V', "Use Viterbi decoding", false,
                      USE_VITERBI_DECODING);
    opt_parse.add_opt("min-meth", 'M',
                      "min cumulative methylation level in a hypemr", false,
                      min_cumulative_meth);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("params-in", 'P', "parameters input file", false,
                      params_in_file);
    opt_parse.add_opt("params-out", 'p', "parameters ouptut file", false,
                      params_out_file);
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
    if (std::size(leftover_args) != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (!is_msite_file(cpgs_file))
      throw std::runtime_error("malformed counts file: " + cpgs_file);

    if (VERBOSE)
      std::cerr << "[loading_data]" << '\n';
    std::vector<MSite> cpgs;
    std::vector<std::pair<double, double>> meth;
    load_cpgs(cpgs_file, cpgs, meth);
    const std::size_t n_sites = std::size(cpgs);
    const double mean_cov =
      std::accumulate(
        std::cbegin(cpgs), std::cend(cpgs), 0.0,
        [](const auto a, const auto &c) { return a + c.n_reads; }) /
      static_cast<double>(n_sites);
    std::cerr << "[n_sites=" << n_sites << "]" << '\n'
              << "[mean_coverage=" << mean_cov << "]" << '\n';

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    std::vector<std::size_t> reset_points;
    std::size_t total_bases = 0, bases_in_deserts = 0;
    separate_regions(desert_size, cpgs, meth, reset_points, total_bases,
                     bases_in_deserts);

    if (VERBOSE) {
      const auto des_frac = static_cast<double>(bases_in_deserts) / total_bases;
      std::cerr << "[n_sites_retained=" << std::size(cpgs) << "]" << '\n'
                << "[deserts_removed=" << std::size(reset_points) - 2 << "]\n"
                << "[remaining_genome_fraction=" << 1.0 - des_frac << "]"
                << '\n';
    }

    ThreeStateHMM hmm(meth, reset_points, tolerance, max_iterations, VERBOSE);

    std::vector<std::vector<double>> trans;
    initialize_transitions(trans);

    betabin hypo_emission, HYPER_emission, HYPO_emission;

    if (!params_in_file.empty())
      read_params_file(params_in_file,
                       // hypo_emission,
                       // HYPER_emission,
                       // HYPO_emission,
                       trans);
    else {
      // NOLINTBEGIN(*-avoid-magic-numbers)
      const double n_reads = mean_cov;
      const double fg_alpha = 0.33 * n_reads;
      const double fg_beta = 0.67 * n_reads;
      const double bg_alpha = 0.67 * n_reads;
      const double bg_beta = 0.33 * n_reads;
      // NOLINTEND(*-avoid-magic-numbers)
      hypo_emission = betabin(fg_alpha, fg_beta);
      HYPER_emission = betabin(bg_alpha, bg_beta);
      HYPO_emission = hypo_emission;
    }

    hmm.set_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);
    if (max_iterations > 0)
      hmm.BaumWelchTraining();
    hmm.get_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);

    if (!params_out_file.empty())
      write_params_file(params_out_file, hypo_emission, HYPER_emission,
                        HYPO_emission, trans);

    // decode the states
    std::vector<STATE_LABELS> classes;
    if (USE_VITERBI_DECODING)
      hmm.ViterbiDecoding();
    else
      hmm.PosteriorDecoding();
    hmm.get_classes(classes);

    // identify the domains of hypermethylation
    std::vector<Interval6> domains;
    build_domains(cpgs, hmm.observations, reset_points, classes, domains);
    filter_domains(min_cumulative_meth, domains);

    // write the results
    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open outfile: " + outfile);
    std::copy(std::cbegin(domains), std::cend(domains),
              std::ostream_iterator<Interval6>(out, "\n"));

    // if requested, write the posterior scores
    if (!scores_file.empty()) {
      if (USE_VITERBI_DECODING)
        hmm.PosteriorDecoding();
      std::vector<Triplet> scores;
      hmm.get_state_posteriors(scores);
      std::ofstream score_out(scores_file);
      for (std::size_t i = 0; i < std::size(cpgs); ++i) {
        score_out << cpgs[i] << "\t";
        if (classes[i] == hypo)
          score_out << "hypo\n" << scores[i].hypo << '\n';
        else if (classes[i] == HYPER)
          score_out << "HYPER\n" << scores[i].HYPER << '\n';
        else  // if (classes[i] == HYPO)
          score_out << "HYPO\n" << scores[i].HYPO << '\n';
      }
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-narrowing-conversions)
