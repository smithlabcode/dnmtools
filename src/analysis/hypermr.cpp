/* Copyright (C) 2009-2023 University of Southern California
 *                         Song Qiang and Andrew D Smith
 * Author: Song Qiang, Andrew D. Smith
 *
 * This file is part of dnmtools.
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
 * 02110-1301 USA
 */

#include <bamxx.hpp>

#include <cmath>
#include <fstream>
#include <numeric>
#include <stdexcept>

#include "GenomicRegion.hpp"
#include "MSite.hpp"
#include "OptionParser.hpp"
#include "ThreeStateHMM.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

using std::begin;
using std::cerr;
using std::cout;
using std::end;
using std::endl;
using std::istream_iterator;
using std::make_pair;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::runtime_error;
using std::string;
using std::to_string;
using std::vector;

using std::ofstream;
using std::ostream_iterator;

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static void
load_cpgs(const string &cpgs_file, vector<MSite> &cpgs,
          vector<pair<double, double>> &meth) {
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in) throw runtime_error("failed opening file: " + cpgs_file);

  MSite the_site;
  string line;
  while (getline(in, line)) cpgs.push_back(MSite(line));

  meth.resize(cpgs.size());
  for (size_t i = 0; i < cpgs.size(); ++i)
    meth[i] = make_pair(cpgs[i].n_meth(), cpgs[i].n_unmeth());
}

template<class T> static void
separate_regions(const size_t desert_size, vector<MSite> &cpgs, vector<T> &meth,
                 vector<size_t> &reset_points, size_t &total_bases,
                 size_t &bases_in_deserts) {
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (cpgs[i].n_reads > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      ++j;
    }
  cpgs.erase(begin(cpgs) + j, end(cpgs));
  meth.erase(begin(meth) + j, end(meth));

  total_bases = 0;
  bases_in_deserts = 0;
  size_t prev_pos = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom)
                          ? cpgs[i].pos - prev_pos
                          : numeric_limits<size_t>::max();
    if (dist > desert_size) {
      reset_points.push_back(i);
      if (dist < numeric_limits<size_t>::max()) bases_in_deserts += dist;
    }
    if (dist < numeric_limits<size_t>::max()) total_bases += dist;

    prev_pos = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());
}

static void
read_params_file(const string &params_file, betabin &hypo_emission,
                 betabin &HYPER_emission, betabin &HYPO_emission,
                 vector<vector<double>> &trans) {
  std::ifstream in(params_file);
  if (!in) throw runtime_error("failed to read param file: " + params_file);

  string hypo_emission_str;
  getline(in, hypo_emission_str);

  string HYPER_emission_str;
  getline(in, HYPER_emission_str);

  string HYPO_emission_str;
  getline(in, HYPO_emission_str);

  trans.resize(3, vector<double>(3, 0.0));
  for (size_t i = 0; i < trans.size(); ++i)
    for (size_t j = 0; j < trans[i].size(); ++j) in >> trans[i][j];
}

static void
write_params_file(const string &params_file, const betabin &hypo_emission,
                  const betabin &HYPER_emission, const betabin &HYPO_emission,
                  const vector<vector<double>> &trans) {
  ofstream out(params_file);
  out << hypo_emission.tostring() << endl
      << HYPER_emission.tostring() << endl
      << HYPO_emission.tostring() << endl;
  for (auto &i : trans) {
    copy(begin(i), end(i), ostream_iterator<double>(out, "\t"));
    out << endl;
  }
}

static void
build_domains(const bool VERBOSE, const vector<MSite> &cpgs,
              const vector<pair<double, double>> &meth,
              const vector<size_t> &reset_points,
              const vector<STATE_LABELS> &classes,
              vector<GenomicRegion> &domains) {
  domains.clear();

  for (size_t i = 0; i < reset_points.size() - 1; ++i) {
    const size_t start = reset_points[i];
    const size_t end = reset_points[i + 1];

    GenomicRegion domain(as_gen_rgn(cpgs[start]));
    STATE_LABELS prev_state = classes[start];
    size_t n = 1;
    double meth_sum =
      meth[start].first / (meth[start].first + meth[start].second);

    for (size_t j = start + 1; j < end; ++j) {
      if ((prev_state == hypo && classes[j] == hypo) ||
          (prev_state != hypo && classes[j] != hypo)) {
        ++n;
        meth_sum += meth[j].first / (meth[j].first + meth[j].second);
      }
      else {
        domain.set_end(cpgs[j - 1].pos + 1);
        const string label = (prev_state == hypo ? "hypo" : "hyper");
        domain.set_name(label + ":" + std::to_string(n));
        domain.set_score(meth_sum);
        domain.set_strand('+');
        if (prev_state == HYPER || prev_state == HYPO)
          domains.push_back(domain);

        domain = GenomicRegion(as_gen_rgn(cpgs[j]));
        n = 1;
        prev_state = classes[j];
        meth_sum = meth[j].first / (meth[j].first + meth[j].second);
      }
    }
    domain.set_end(cpgs[end - 1].pos + 1);
    const string label = (prev_state == hypo ? "hypo" : "hyper");
    domain.set_name(label + ":" + std::to_string(n));
    domain.set_score(meth_sum);
    domain.set_strand('+');
    if (prev_state == HYPER || prev_state == HYPO) domains.push_back(domain);
  }
}

static void
filter_domains(const bool VERBOSE, const double min_cumulative_meth,
               vector<GenomicRegion> &domains) {
  size_t j = 0;
  for (size_t i = 0; i < domains.size(); ++i)
    if (domains[i].get_score() >= min_cumulative_meth)
      domains[j++] = domains[i];
  domains.erase(begin(domains) + j, end(domains));
}

static void
initialize_transitions(vector<vector<double>> &trans) {
  // "hypo" -> "hypo" is high, and goes only to "HYPER";
  // "HYPER" -> "HYPER" or "HYPO" equally;
  // "HYPO" -> "HYPER" usually; never to "hypo"
  trans = {{0.990, 0.010, 0.000}, {0.050, 0.475, 0.475}, {0.000, 0.666, 0.334}};
}

int
main_hypermr(int argc, const char **argv) {
  try {
    static const string description =
      "Identify regions of elevated methylation. Designed for "
      "methylomes of plants, especially A. thaliana, and other "
      "organisms with a low background methylation level, but "
      "with increased methylation associated with biological "
      "activity at a specific regions of the genome.";

    string outfile;
    string scores_file;
    string trans_file;

    size_t desert_size = 1000;
    size_t max_iterations = 10;

    // run mode flags
    bool VERBOSE = false;

    // corrections for small values (not parameters):
    double tolerance = 1e-10;

    double min_cumulative_meth = 4.0;
    bool USE_VITERBI_DECODING = false;

    string params_in_file;
    string params_out_file;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description, "<counts-file>");
    opt_parse.add_opt("out", 'o', "output file (BED format)", false, outfile);
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
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (VERBOSE) cerr << "[loading_data]" << endl;
    vector<MSite> cpgs;
    vector<pair<double, double>> meth;
    load_cpgs(cpgs_file, cpgs, meth);
    const size_t n_sites = cpgs.size();

    double mean_cov = 0.0;
    for (auto &&c : cpgs) mean_cov += c.n_reads;
    mean_cov /= n_sites;

    if (VERBOSE)
      cerr << "[n_sites=" << n_sites << "]" << endl
           << "[mean_coverage=" << mean_cov << "]" << endl;

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    size_t total_bases = 0, bases_in_deserts = 0;
    separate_regions(desert_size, cpgs, meth, reset_points, total_bases,
                     bases_in_deserts);

    if (VERBOSE) {
      const auto des_frac = static_cast<double>(bases_in_deserts) / total_bases;
      cerr << "[n_sites_retained=" << cpgs.size() << "]" << endl
           << "[deserts_removed=" << reset_points.size() - 2 << "]" << endl
           << "[remaining_genome_fraction=" << 1.0 - des_frac << "]" << endl;
    }

    ThreeStateHMM hmm(meth, reset_points, tolerance, max_iterations, VERBOSE);

    vector<vector<double>> trans;
    initialize_transitions(trans);

    betabin hypo_emission, HYPER_emission, HYPO_emission;

    if (!params_in_file.empty())
      read_params_file(params_in_file, hypo_emission, HYPER_emission,
                       HYPO_emission, trans);
    else {
      const double n_reads = mean_cov;
      const double fg_alpha = 0.33 * n_reads;
      const double fg_beta = 0.67 * n_reads;
      const double bg_alpha = 0.67 * n_reads;
      const double bg_beta = 0.33 * n_reads;

      hypo_emission = betabin(fg_alpha, fg_beta);
      HYPER_emission = betabin(bg_alpha, bg_beta);
      HYPO_emission = hypo_emission;
    }

    hmm.set_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);
    if (max_iterations > 0) hmm.BaumWelchTraining();
    hmm.get_parameters(hypo_emission, HYPER_emission, HYPO_emission, trans);

    if (!params_out_file.empty())
      write_params_file(params_out_file, hypo_emission, HYPER_emission,
                        HYPO_emission, trans);

    // decode the states
    vector<STATE_LABELS> classes;
    if (USE_VITERBI_DECODING)
      hmm.ViterbiDecoding();
    else
      hmm.PosteriorDecoding();
    hmm.get_classes(classes);

    // identify the domains of hypermethylation
    vector<GenomicRegion> domains;
    build_domains(VERBOSE, cpgs, hmm.observations, reset_points, classes,
                  domains);
    filter_domains(VERBOSE, min_cumulative_meth, domains);

    // write the results
    ofstream of;
    if (!outfile.empty()) of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
    copy(begin(domains), end(domains),
         ostream_iterator<GenomicRegion>(out, "\n"));

    // if requested, write the posterior scores
    if (!scores_file.empty()) {
      if (USE_VITERBI_DECODING) hmm.PosteriorDecoding();
      vector<Triplet> scores;
      hmm.get_state_posteriors(scores);
      ofstream score_out(scores_file);
      for (size_t i = 0; i < cpgs.size(); ++i) {
        score_out << cpgs[i] << "\t";
        if (classes[i] == hypo)
          score_out << "hypo\n" << scores[i].hypo << endl;
        else if (classes[i] == HYPER)
          score_out << "HYPER\n" << scores[i].HYPER << endl;
        else  // if (classes[i] == HYPO)
          score_out << "HYPO\n" << scores[i].HYPO << endl;
      }
    }
  }
  catch (const std::exception &e) {
    cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
