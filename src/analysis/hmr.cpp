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

#include "MSite.hpp"
#include "TwoStateHMM.hpp"
#include "counts_header.hpp"

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <bamxx.hpp>

#include <cmath>
#include <cstdint>  // for [u]int[0-9]+_t
#include <fstream>
#include <iomanip>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_set>

struct hmr_summary {
  hmr_summary(const std::vector<GenomicRegion> &hmrs) {
    hmr_count = hmrs.size();
    hmr_total_size =
      std::accumulate(std::cbegin(hmrs), std::cend(hmrs), 0ul,
                      [](const std::uint64_t t, const GenomicRegion &p) {
                        return t + p.get_width();
                      });
    hmr_mean_size = static_cast<double>(hmr_total_size) /
                    std::max(hmr_count, static_cast<std::uint64_t>(1));
  }
  // hmr_count is the number of identified HMRs.
  std::uint64_t hmr_count{};
  // total_hmr_size is the sum of the sizes of the identified HMRs
  std::uint64_t hmr_total_size{};
  // mean_hmr_size is the mean size of the identified HMRs
  double hmr_mean_size{};

  std::string
  tostring() {
    std::ostringstream oss;
    oss << "hmr_count: " << hmr_count << '\n'
        << "hmr_total_size: " << hmr_total_size << '\n'
        << "hmr_mean_size: " << std::fixed << std::setprecision(2)
        << hmr_mean_size;

    return oss.str();
  }
};

static GenomicRegion
as_gen_rgn(const MSite &s) {
  return GenomicRegion(s.chrom, s.pos, s.pos + 1);
}

static std::string
format_cpg_meth_tag(const std::pair<double, double> &m) {
  return "CpG:" + std::to_string(static_cast<std::size_t>(m.first)) + ":" +
         std::to_string(static_cast<std::size_t>(m.second));
}

static double
get_stepup_cutoff(std::vector<double> scores, const double cutoff) {
  if (cutoff <= 0)
    return std::numeric_limits<double>::max();
  else if (cutoff > 1)
    return std::numeric_limits<double>::min();

  const std::size_t n = scores.size();
  std::sort(begin(scores), std::end(scores));
  std::size_t i = 1;
  while (i < n && scores[i - 1] < (cutoff * i) / n)
    ++i;
  return scores[i - 1];
}

static void
get_domain_scores(const std::vector<bool> &state_ids,
                  const std::vector<std::pair<double, double>> &meth,
                  const std::vector<std::size_t> &reset_points,
                  std::vector<double> &scores) {

  std::size_t reset_idx = 1;
  bool in_domain = false;
  double score = 0;
  for (std::size_t i = 0; i < state_ids.size(); ++i) {
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
build_domains(const std::vector<MSite> &cpgs,
              const std::vector<std::size_t> &reset_points,
              const std::vector<bool> &state_ids,
              std::vector<GenomicRegion> &domains) {

  std::size_t n_cpgs = 0, reset_idx = 1, prev_end = 0;
  bool in_domain = false;
  for (std::size_t i = 0; i < state_ids.size(); ++i) {
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
        domains.back().set_name("HYPO" + std::to_string(domains.size()));
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
separate_regions(const bool verbose, const std::size_t desert_size,
                 std::vector<MSite> &cpgs, std::vector<T> &meth,
                 std::vector<U> &reads,
                 std::vector<std::size_t> &reset_points) {
  if (verbose)
    std::cerr << "[separating by cpg desert]" << '\n';
  // eliminate the zero-read cpgs
  std::size_t j = 0;
  for (std::size_t i = 0; i < cpgs.size(); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(std::begin(cpgs) + j, std::end(cpgs));
  meth.erase(std::begin(meth) + j, std::end(meth));
  reads.erase(std::begin(reads) + j, std::end(reads));

  double total_bases = 0;
  double bases_in_deserts = 0;
  // segregate cpgs
  std::size_t prev_pos = 0;
  for (std::size_t i = 0; i < cpgs.size(); ++i) {
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
  reset_points.push_back(cpgs.size());

  if (verbose)
    std::cerr << "[cpgs retained: " << cpgs.size() << "]" << '\n'
              << "[deserts removed: " << reset_points.size() - 2 << "]" << '\n'
              << "[genome fraction covered: "
              << 1.0 - (bases_in_deserts / total_bases) << "]" << '\n';
}

/* function to "fold" the methylation profile so that the middle
 * methylation becomes lower methylation, and both the low and high
 * methylation become high. this method actually seems to work.
 */
static void
make_partial_meth(const std::vector<std::uint32_t> &reads,
                  std::vector<std::pair<double, double>> &meth) {
  for (std::size_t i = 0; i < reads.size(); ++i) {
    double m = meth[i].first / reads[i];
    m = (m <= 0.5) ? (1.0 - 2 * m) : (1.0 - 2 * (1.0 - m));
    meth[i].first = reads[i] * m;
    meth[i].second = (reads[i] - meth[i].first);
  }
}

static void
shuffle_cpgs(const std::size_t rng_seed, const TwoStateHMM &hmm,
             std::vector<std::pair<double, double>> meth,
             std::vector<std::size_t> reset_points, const double p_fb,
             const double p_bf, const double fg_alpha, const double fg_beta,
             const double bg_alpha, const double bg_beta,
             std::vector<double> &domain_scores) {

  auto eng = std::default_random_engine(rng_seed);
  std::shuffle(std::begin(meth), std::end(meth), eng);

  std::vector<bool> state_ids;
  std::vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                        bg_alpha, bg_beta, state_ids, scores);
  get_domain_scores(state_ids, meth, reset_points, domain_scores);
  sort(std::begin(domain_scores), std::end(domain_scores));
}

static void
assign_p_values(const std::vector<double> &random_scores,
                const std::vector<double> &observed_scores,
                std::vector<double> &p_values) {
  const double n_randoms = random_scores.empty() ? 1 : random_scores.size();
  for (std::size_t i = 0; i < observed_scores.size(); ++i)
    p_values.push_back(
      (std::end(random_scores) - upper_bound(std::begin(random_scores),
                                             std::end(random_scores),
                                             observed_scores[i])) /
      n_randoms);
}

static void
read_params_file(const bool verbose, const std::string &params_file,
                 double &fg_alpha, double &fg_beta, double &bg_alpha,
                 double &bg_beta, double &p_fb, double &p_bf,
                 double &domain_score_cutoff) {
  std::string jnk;
  std::ifstream in(params_file);
  if (!in)
    throw std::runtime_error("failed to parse params file: " + params_file);
  in >> jnk >> fg_alpha >> jnk >> fg_beta >> jnk >> bg_alpha >> jnk >>
    bg_beta >> jnk >> p_fb >> jnk >> p_bf >> jnk >> domain_score_cutoff;
  if (verbose)
    std::cerr << "FG_ALPHA\t" << fg_alpha << '\n'
              << "FG_BETA\t" << fg_beta << '\n'
              << "BG_ALPHA\t" << bg_alpha << '\n'
              << "BG_BETA\t" << bg_beta << '\n'
              << "F_B\t" << p_fb << '\n'
              << "B_F\t" << p_bf << '\n'
              << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << '\n';
}

static void
write_params_file(const std::string &outfile, const double fg_alpha,
                  const double fg_beta, const double bg_alpha,
                  const double bg_beta, const double p_fb, const double p_bf,
                  const double domain_score_cutoff) {

  std::ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str());
  std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

  out.precision(30);
  out << "FG_ALPHA\t" << fg_alpha << '\n'
      << "FG_BETA\t" << fg_beta << '\n'
      << "BG_ALPHA\t" << bg_alpha << '\n'
      << "BG_BETA\t" << bg_beta << '\n'
      << "F_B\t" << p_fb << '\n'
      << "B_F\t" << p_bf << '\n'
      << "DOMAIN_SCORE_CUTOFF\t" << domain_score_cutoff << '\n';
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
load_cpgs(const std::string &cpgs_file, std::vector<MSite> &cpgs,
          std::vector<std::pair<double, double>> &meth,
          std::vector<std::uint32_t> &reads) {

  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("failed opening file: " + cpgs_file);

  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);

  MSite prev_site, the_site;
  while (read_site(in, the_site)) {
    if (!the_site.is_cpg() || distance(prev_site, the_site) < 2)
      throw std::runtime_error("error: input is not symmetric-CpGs: " +
                               cpgs_file);
    cpgs.push_back(the_site);
    reads.push_back(the_site.n_reads);
    meth.push_back(
      std::make_pair(get_n_meth(the_site), get_n_unmeth(the_site)));
    prev_site = the_site;
  }
}

template <typename InputIterator, typename T = double>
static auto
get_mean(InputIterator first, InputIterator last) -> T {
  return std::accumulate(first, last, static_cast<T>(0)) /
         std::distance(first, last);
}

template <class T>
static void
check_sorted_within_chroms(T first, const T last) {

  // empty, or a single element
  if (first == last || first + 1 == last)
    return;

  std::unordered_set<std::string> chroms_seen;
  std::string cur_chrom = "";

  for (auto fast = first + 1; fast != last; ++fast, ++first) {
    if (fast->chrom != cur_chrom) {
      if (chroms_seen.find(fast->chrom) != std::end(chroms_seen)) {
        throw std::runtime_error("input not grouped by chromosomes. "
                                 "Error in the following line:\n" +
                                 fast->tostring());
      }
      cur_chrom = fast->chrom;
      chroms_seen.insert(cur_chrom);
      // don't check ordering if chroms have changed
    }
    else {  // first and fast are in the same chrom
      if (first->pos >= fast->pos) {
        throw std::runtime_error("input file not sorted properly. "
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

    std::string outfile;
    std::string hypo_post_outfile;
    std::string meth_post_outfile;

    std::size_t desert_size = 1000;
    std::size_t max_iterations = 10;
    std::size_t rng_seed = 408;

    // run mode flags
    bool verbose = false;
    bool PARTIAL_METH = false;
    bool allow_extra_fields = false;

    // corrections for small values
    const double tolerance = 1e-10;

    std::string summary_file;

    std::string params_in_file;
    std::string params_out_file;

    const std::string description =
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
    const std::string cpgs_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    MSite::no_extra_fields = (allow_extra_fields == false);

    if (!is_msite_file(cpgs_file))
      throw std::runtime_error("malformed counts file: " + cpgs_file);

    // separate the regions by chrom and by desert
    std::vector<MSite> cpgs;
    std::vector<std::pair<double, double>> meth;
    std::vector<std::uint32_t> reads;
    if (verbose)
      std::cerr << "[reading methylation levels]\n";
    load_cpgs(cpgs_file, cpgs, meth, reads);

    if (verbose)
      std::cerr << "[checking if input is properly formatted]\n";
    check_sorted_within_chroms(std::begin(cpgs), std::end(cpgs));
    if (PARTIAL_METH)
      make_partial_meth(reads, meth);

    const auto mean_coverage = get_mean(std::cbegin(reads), std::cend(reads));

    if (verbose)
      std::cerr << "[total_cpgs=" << size(cpgs) << "]" << '\n'
                << "[mean_coverage=" << mean_coverage << "]" << '\n';

    // check for sufficient data
    if (mean_coverage < static_cast<double>(min_coverage)) {
      if (verbose)
        std::cerr << "error: insufficient data"
                  << " mean_coverage=" << mean_coverage
                  << " min_coverage=" << min_coverage << '\n';
      if (!summary_file.empty()) {
        std::ofstream summary_out(summary_file);
        if (!summary_out)
          throw std::runtime_error("failed to open: " + summary_file);
        summary_out << hmr_summary({}).tostring() << '\n';
      }
      return EXIT_SUCCESS;
    }

    // separate the regions by chrom and by desert, and eliminate
    // those isolated CpGs
    std::vector<std::size_t> reset_points;
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
      const double n_reads = get_mean(std::cbegin(reads), std::cend(reads));
      fg_alpha = 0.33 * n_reads;
      fg_beta = 0.67 * n_reads;
      bg_alpha = 0.67 * n_reads;
      bg_beta = 0.33 * n_reads;
    }

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                            bg_alpha, bg_beta);
    // DECODE THE DOMAINS
    std::vector<bool> state_ids;
    std::vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                          bg_alpha, bg_beta, state_ids, posteriors);

    std::vector<double> domain_scores;
    get_domain_scores(state_ids, meth, reset_points, domain_scores);

    std::vector<double> random_scores;
    shuffle_cpgs(rng_seed, hmm, meth, reset_points, p_fb, p_bf, fg_alpha,
                 fg_beta, bg_alpha, bg_beta, random_scores);

    std::vector<double> p_values;
    assign_p_values(random_scores, domain_scores, p_values);

    if (domain_score_cutoff == std::numeric_limits<double>::max() &&
        !domain_scores.empty())
      domain_score_cutoff = get_stepup_cutoff(p_values, 0.01);

    // write parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, fg_alpha, fg_beta, bg_alpha, bg_beta,
                        p_fb, p_bf, domain_score_cutoff);

    std::vector<GenomicRegion> domains;
    if (!domain_scores.empty())
      build_domains(cpgs, reset_points, state_ids, domains);

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    std::size_t good_hmr_count = 0;
    for (auto i = 0u; i < size(domains); ++i)
      if (p_values[i] < domain_score_cutoff) {
        domains[good_hmr_count] = domains[i];
        domains[good_hmr_count].set_name("HYPO" +
                                         std::to_string(good_hmr_count));
        ++good_hmr_count;
      }
    domains.resize(good_hmr_count);

    for (const auto &d : domains)
      out << d << '\n';

    if (!hypo_post_outfile.empty()) {
      if (verbose)
        std::cerr << "[writing=" << hypo_post_outfile << "]\n";
      std::ofstream out(hypo_post_outfile);
      for (std::size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(posteriors[i]);
        out << cpg << '\n';
      }
    }

    if (!meth_post_outfile.empty()) {
      std::ofstream out(meth_post_outfile);
      if (verbose)
        std::cerr << "[writing=" << meth_post_outfile << "]\n";
      for (std::size_t i = 0; i < cpgs.size(); ++i) {
        GenomicRegion cpg(as_gen_rgn(cpgs[i]));
        cpg.set_name(format_cpg_meth_tag(meth[i]));
        cpg.set_score(1.0 - posteriors[i]);
        out << cpg << '\n';
      }
    }
    if (!summary_file.empty()) {
      std::ofstream summary_out(summary_file);
      if (!summary_out)
        throw std::runtime_error("failed to open: " + summary_file);
      summary_out << hmr_summary(domains).tostring() << '\n';
    }
  }
  catch (std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
