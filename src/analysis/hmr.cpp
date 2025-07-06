/* Copyright (C) 2009-2023 University of Southern California
 *                         Andrew D Smith
 *
 * Author: Andrew D. Smith, Song Qiang
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
 */

#include <cmath>
#include <cstdint>  // for [u]int[0-9]+_t
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include <bamxx.hpp>

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include "MSite.hpp"
#include "TwoStateHMM.hpp"
#include "counts_header.hpp"

using std::accumulate;
using std::cerr;
using std::cout;
using std::endl;
using std::make_pair;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::reduce;
using std::runtime_error;
using std::string;
using std::to_string;
using std::unordered_set;
using std::vector;

using bamxx::bgzf_file;

struct hmr_summary {
  hmr_summary(const vector<GenomicRegion> &hmrs) {
    hmr_count = hmrs.size();
    hmr_total_size = accumulate(cbegin(hmrs), cend(hmrs), 0ul,
                                [](const uint64_t t, const GenomicRegion &p) {
                                  return t + p.get_width();
                                });
    hmr_mean_size = static_cast<double>(hmr_total_size) /
                    std::max(hmr_count, static_cast<uint64_t>(1));
  }
  // hmr_count is the number of identified HMRs.
  uint64_t hmr_count{};
  // total_hmr_size is the sum of the sizes of the identified HMRs
  uint64_t hmr_total_size{};
  // mean_hmr_size is the mean size of the identified HMRs
  double hmr_mean_size{};

  string
  tostring() {
    std::ostringstream oss;
    oss << "hmr_count: " << hmr_count << endl
        << "hmr_total_size: " << hmr_total_size << endl
        << "hmr_mean_size: " << std::fixed << std::setprecision(2)
        << hmr_mean_size;

    return oss.str();
  }
};

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static string
format_cpg_meth_tag(const pair<double, double> &m) {
  return "CpG:" + to_string(static_cast<size_t>(m.first)) + ":" +
         to_string(static_cast<size_t>(m.second));
}

static double
get_stepup_cutoff(vector<double> scores, const double cutoff) {
  if (cutoff <= 0)
    return numeric_limits<double>::max();
  else if (cutoff > 1)
    return numeric_limits<double>::min();

  const size_t n = scores.size();
  std::sort(begin(scores), end(scores));
  size_t i = 1;
  while (i < n && scores[i - 1] < (cutoff * i) / n)
    ++i;
  return scores[i - 1];
}

static void
get_domain_scores(const vector<bool> &state_ids,
                  const vector<pair<double, double>> &meth,
                  const vector<size_t> &reset_points, vector<double> &scores) {

  size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
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
      score += 1.0 - (meth[i].first / (meth[i].first + meth[i].second));
    }
    else if (in_domain) {
      in_domain = false;
      scores.push_back(score);
      score = 0;
    }
  }

  if (in_domain)
    scores.push_back(score);
}

static void
build_domains(const vector<MSite> &cpgs, const vector<size_t> &reset_points,
              const vector<bool> &state_ids, vector<GenomicRegion> &domains) {

  size_t n_cpgs = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  for (size_t i = 0; i < state_ids.size(); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        domains.back().set_end(prev_end);
        domains.back().set_score(n_cpgs);
        n_cpgs = 0;
      }
      ++reset_idx;
    }
    if (state_ids[i]) {  // currently in an hmr
      if (!in_domain) {
        in_domain = true;
        domains.push_back(as_gen_rgn(cpgs[i]));
        domains.back().set_name("HYPO" + to_string(domains.size()));
      }
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      domains.back().set_end(prev_end);
      domains.back().set_score(n_cpgs);
      n_cpgs = 0;
    }
    prev_end = cpgs[i].pos + 1;
  }
  if (in_domain) {
    domains.back().set_end(prev_end);
    domains.back().set_score(n_cpgs);
  }
}

template <class T, class U>
static void
separate_regions(const bool verbose, const size_t desert_size,
                 vector<MSite> &cpgs, vector<T> &meth, vector<U> &reads,
                 vector<size_t> &reset_points) {
  if (verbose)
    cerr << "[separating by cpg desert]" << endl;
  // eliminate the zero-read cpgs
  size_t j = 0;
  for (size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(begin(cpgs) + j, end(cpgs));
  meth.erase(begin(meth) + j, end(meth));
  reads.erase(begin(reads) + j, end(reads));

  double total_bases = 0;
  double bases_in_deserts = 0;
  // segregate cpgs
  size_t prev_pos = 0;
  for (size_t i = 0; i < cpgs.size(); ++i) {
    const size_t dist = (i > 0 && cpgs[i].chrom == cpgs[i - 1].chrom)
                          ? cpgs[i].pos - prev_pos
                          : numeric_limits<size_t>::max();
    if (dist > desert_size) {
      reset_points.push_back(i);
      if (dist < numeric_limits<size_t>::max())
        bases_in_deserts += dist;
    }
    if (dist < numeric_limits<size_t>::max())
      total_bases += dist;

    prev_pos = cpgs[i].pos;
  }
  reset_points.push_back(cpgs.size());

  if (verbose)
    cerr << "[cpgs retained: " << cpgs.size() << "]" << endl
         << "[deserts removed: " << reset_points.size() - 2 << "]" << endl
         << "[genome fraction covered: "
         << 1.0 - (bases_in_deserts / total_bases) << "]" << endl;
}

/* function to "fold" the methylation profile so that the middle
 * methylation becomes lower methylation, and both the low and high
 * methylation become high. this method actually seems to work.
 */
static void
make_partial_meth(const vector<uint32_t> &reads,
                  vector<pair<double, double>> &meth) {
  for (size_t i = 0; i < reads.size(); ++i) {
    double m = meth[i].first / reads[i];
    m = (m <= 0.5) ? (1.0 - 2 * m) : (1.0 - 2 * (1.0 - m));
    meth[i].first = reads[i] * m;
    meth[i].second = (reads[i] - meth[i].first);
  }
}

static void
shuffle_cpgs(const size_t rng_seed, const TwoStateHMM &hmm,
             vector<pair<double, double>> meth, vector<size_t> reset_points,
             const double p_fb, const double p_bf, const double fg_alpha,
             const double fg_beta, const double bg_alpha, const double bg_beta,
             vector<double> &domain_scores) {

  auto eng = std::default_random_engine(rng_seed);
  std::shuffle(begin(meth), end(meth), eng);

  vector<bool> state_ids;
  vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                        bg_alpha, bg_beta, state_ids, scores);
  get_domain_scores(state_ids, meth, reset_points, domain_scores);
  sort(begin(domain_scores), end(domain_scores));
}

static void
assign_p_values(const vector<double> &random_scores,
                const vector<double> &observed_scores,
                vector<double> &p_values) {
  const double n_randoms = random_scores.empty() ? 1 : random_scores.size();
  for (size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back((end(random_scores) - upper_bound(begin(random_scores),
                                                         end(random_scores),
                                                         observed_scores[i])) /
                       n_randoms);
}

static void
read_params_file(const bool verbose, const string &params_file,
                 double &fg_alpha, double &fg_beta, double &bg_alpha,
                 double &bg_beta, double &p_fb, double &p_bf,
                 double &domain_score_cutoff) {
  string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw runtime_error("failed to parse params file: " + params_file);
  in >> jnk >> fg_alpha >> jnk >> fg_beta >> jnk >> bg_alpha >> jnk >>
    bg_beta >> jnk >> p_fb >> jnk >> p_bf >> jnk >> domain_score_cutoff;
  if (verbose)
    cerr << "FG_ALPHA\t" << fg_alpha << endl
         << "FG_BETA\t" << fg_beta << endl
         << "BG_ALPHA\t" << bg_alpha << endl
         << "BG_BETA\t" << bg_beta << endl
         << "F_B\t" << p_fb << endl
         << "B_F\t" << p_bf << endl
         << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << endl;
}

static void
write_params_file(const string &outfile, const double fg_alpha,
                  const double fg_beta, const double bg_alpha,
                  const double bg_beta, const double p_fb, const double p_bf,
                  const double domain_score_cutoff) {

  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "FG_ALPHA\t" << fg_alpha << endl
      << "FG_BETA\t" << fg_beta << endl
      << "BG_ALPHA\t" << bg_alpha << endl
      << "BG_BETA\t" << bg_beta << endl
      << "F_B\t" << p_fb << endl
      << "B_F\t" << p_bf << endl
      << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << endl;
}

[[nodiscard]]
static inline double
get_n_meth(const MSite &s) {
  return s.meth * s.n_reads;
}

[[nodiscard]]
static inline double
get_n_unmeth(const MSite &s) {
  return (1.0 - s.meth) * s.n_reads;
}

static void
load_cpgs(const string &cpgs_file, vector<MSite> &cpgs,
          vector<pair<double, double>> &meth, vector<uint32_t> &reads) {

  bgzf_file in(cpgs_file, "r");
  if (!in)
    throw runtime_error("failed opening file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  MSite prev_site, the_site;
  while (read_site(in, the_site)) {
    if (!the_site.is_cpg() || distance(prev_site, the_site) < 2)
      throw runtime_error("error: input is not symmetric-CpGs: " + cpgs_file);
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    meth.push_back(make_pair(get_n_meth(the_site), get_n_unmeth(the_site)));
    prev_site = the_site;
  }
}

template <typename InputIterator, typename T = double>
static auto
get_mean(InputIterator first, InputIterator last) -> T {
  return accumulate(first, last, static_cast<T>(0)) /
         std::distance(first, last);
}

template <class T>
static void
check_sorted_within_chroms(T first, const T last) {

  // empty, or a single element
  if (first == last || first + 1 == last)
    return;

  unordered_set<string> chroms_seen;
  string cur_chrom = "";

  for (auto fast = first + 1; fast != last; ++fast, ++first) {
    if (fast->chrom != cur_chrom) {
      if (chroms_seen.find(fast->chrom) != end(chroms_seen)) {
        throw runtime_error("input not grouped by chromosomes. "
                            "Error in the following line:\n" +
                            fast->tostring());
      }
      cur_chrom = fast->chrom;
      chroms_seen.insert(cur_chrom);
      // don't check ordering if chroms have changed
    }
    else {  // first and fast are in the same chrom
      if (first->pos >= fast->pos) {
        throw runtime_error("input file not sorted properly. "
                            "Error in the following lines:\n" +
                            first->tostring() + "\n" + fast->tostring());
      }
    }
  }
}

int
main_hmr(int argc, char *argv[]) {

  try {

    constexpr double min_coverage = 1.0;

    string outfile;
    string hypo_post_outfile;
    string meth_post_outfile;

    size_t desert_size = 1000;
    size_t max_iterations = 10;
    size_t rng_seed = 408;

    // run mode flags
    bool verbose = false;
    bool PARTIAL_METH = false;
    bool allow_extra_fields = false;

    // corrections for small values
    const double tolerance = 1e-10;

    string summary_file;

    string params_in_file;
    string params_out_file;

    const string description =
      "Identify HMRs in methylomes. Methylation must be provided in the \
      methcounts format (chrom, position, strand, context,              \
      methylation, reads). See the methcounts documentation for         \
      details. This program assumes only data at CpG sites and that     \
      strands are collapsed so only the positive site appears in the    \
      file.";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<methylation-file>");
    opt_parse.add_opt("out", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn covered cpgs in HMR", false,
                      desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("partial", '\0', "identify PMRs instead of HMRs", false,
                      PARTIAL_METH);
    opt_parse.add_opt("post-hypo", '\0',
                      "output file for single-CpG posterior "
                      "hypomethylation probability (default: none)",
                      false, hypo_post_outfile);
    opt_parse.add_opt("post-meth", '\0',
                      "output file for single-CpG posteiror "
                      "methylation probability (default: none)",
                      false, meth_post_outfile);
    opt_parse.add_opt("params-in", 'P',
                      "HMM parameter file "
                      "(override training)",
                      false, params_in_file);
    opt_parse.add_opt("params-out", 'p',
                      "write HMM parameters to this "
                      "file (default: none)",
                      false, params_out_file);
    opt_parse.add_opt("seed", 's', "specify random seed", false, rng_seed);
    opt_parse.add_opt("summary", 'S', "write summary output here", false,
                      summary_file);
    opt_parse.add_opt("relaxed", '\0',
                      "input has extra fields (used for nanopore)", false,
                      allow_extra_fields);
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
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    MSite::no_extra_fields = (allow_extra_fields == false);

    if (!is_msite_file(cpgs_file))
      throw runtime_error("malformed counts file: " + cpgs_file);

    // separate the regions by chrom and by desert
    vector<MSite> cpgs;
    vector<pair<double, double>> meth;
    vector<uint32_t> reads;
    if (verbose)
      cerr << "[reading methylation levels]" << endl;
    load_cpgs(cpgs_file, cpgs, meth, reads);

    if (verbose)
      cerr << "[checking if input is properly formatted]" << endl;
    check_sorted_within_chroms(begin(cpgs), end(cpgs));
    if (PARTIAL_METH)
      make_partial_meth(reads, meth);

    const auto mean_coverage = get_mean(cbegin(reads), cend(reads));

    if (verbose)
      cerr << "[total_cpgs=" << size(cpgs) << "]" << '\n'
           << "[mean_coverage=" << mean_coverage << "]" << endl;

    // check for sufficient data
    if (mean_coverage < static_cast<double>(min_coverage)) {
      if (verbose)
        cerr << "error: insufficient data"
             << " mean_coverage=" << mean_coverage
             << " min_coverage=" << min_coverage << endl;
      if (!summary_file.empty()) {
        std::ofstream summary_out(summary_file);
        if (!summary_out)
          throw runtime_error("failed to open: " + summary_file);
        summary_out << hmr_summary({}).tostring() << endl;
      }
      return EXIT_SUCCESS;
    }

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    vector<size_t> reset_points;
    separate_regions(verbose, desert_size, cpgs, meth, reads, reset_points);

    const TwoStateHMM hmm(tolerance, max_iterations, verbose);

    double p_fb = 0.25;
    double p_bf = 0.25;

    double fg_alpha = 0, fg_beta = 0;
    double bg_alpha = 0, bg_beta = 0;
    double domain_score_cutoff = std::numeric_limits<double>::max();

    if (!params_in_file.empty()) {  // read parameters file
      read_params_file(verbose, params_in_file, fg_alpha, fg_beta, bg_alpha,
                       bg_beta, p_fb, p_bf, domain_score_cutoff);
      max_iterations = 0;
    }
    else {
      const double n_reads = get_mean(begin(reads), end(reads));
      fg_alpha = 0.33 * n_reads;
      fg_beta = 0.67 * n_reads;
      bg_alpha = 0.67 * n_reads;
      bg_beta = 0.33 * n_reads;
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                            bg_alpha, bg_beta);
    // DECODE THE DOMAINS
    vector<bool> state_ids;
    vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                          bg_alpha, bg_beta, state_ids, posteriors);

    vector<double> domain_scores;
    get_domain_scores(state_ids, meth, reset_points, domain_scores);

    vector<double> random_scores;
    shuffle_cpgs(rng_seed, hmm, meth, reset_points, p_fb, p_bf, fg_alpha,
                 fg_beta, bg_alpha, bg_beta, random_scores);

    vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (domain_score_cutoff == numeric_limits<double>::max() &&
        !domain_scores.empty())
      domain_score_cutoff = get_stepup_cutoff(p_values, 0.01);

    // write parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
                        p_fb, p_bf, domain_score_cutoff);

    vector<GenomicRegion> domains;
    if (!domain_scores.empty())
      build_domains(cpgs, reset_points, state_ids, domains);

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    size_t good_hmr_count = 0;
    for (auto i = 0u; i < size(domains); ++i)
      if (p_values[i] < domain_score_cutoff) {
        domains[good_hmr_count] = domains[i];
        domains[good_hmr_count].set_name("HYPO" + to_string(good_hmr_count));
        ++good_hmr_count;
      }
    domains.resize(good_hmr_count);

    for (const auto &d : domains)
      out << d << '\n';

    if (!hypo_post_outfile.empty()) {
      if (verbose)
        cerr << "[writing=" << hypo_post_outfile << "]" << endl;
      std::ofstream out(hypo_post_outfile);
      for (size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(posteriors[i]);
        out << cpg << '\n';
      }
    }

    if (!meth_post_outfile.empty()) {
      std::ofstream out(meth_post_outfile);
      if (verbose)
        cerr << "[writing=" << meth_post_outfile << "]" << endl;
      for (size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(1.0 - posteriors[i]);
        out << cpg << '\n';
      }
    }
    if (!summary_file.empty()) {
      std::ofstream summary_out(summary_file);
      if (!summary_out)
        throw runtime_error("failed to open: " + summary_file);
      summary_out << hmr_summary(domains).tostring() << endl;
    }
  }
  catch (runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
