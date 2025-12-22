/* Copyright (C) 2009-2025 Andrew D Smith and Song Qiang
 *
 * Author: Andrew D. Smith and Song Qiang
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
#include "OptionParser.hpp"
#include "TwoStateHMM.hpp"
#include "counts_header.hpp"

#include <bamxx.hpp>

#include "nlohmann/json.hpp"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

template <typename InputIterator, typename T = double>
[[nodiscard]] static auto
get_mean(InputIterator first, InputIterator last) -> T {
  return std::accumulate(first, last, static_cast<T>(0)) /
         std::distance(first, last);
}

struct hmr_params {
  static constexpr auto frac_init = 0.33;
  static constexpr auto trans_init = 0.25;
  double fg_alpha{};
  double fg_beta{};
  double bg_alpha{};
  double bg_beta{};
  double p_fb{};
  double p_bf{};
  double domain_score_cutoff{};

  [[nodiscard]] auto
  fg_meth() const {
    return fg_alpha / (fg_alpha + fg_beta);
  }

  [[nodiscard]] auto
  bg_meth() const {
    return bg_alpha / (bg_alpha + bg_beta);
  }

  auto
  init(const double n_reads) {
    fg_alpha = frac_init * n_reads;
    fg_beta = (1.0 - frac_init) * n_reads;
    bg_alpha = (1.0 - frac_init) * n_reads;
    bg_beta = frac_init * n_reads;
    p_fb = trans_init;
    p_bf = trans_init;
    domain_score_cutoff = std::numeric_limits<double>::max();
  }

  /// Estimate methylation level using param values and the posterior
  /// probability on HMR status
  [[nodiscard]] auto
  est_meth(const double posterior) -> double {
    return fg_meth() * posterior + bg_meth() * (1.0 - posterior);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(hmr_params, fg_alpha, fg_beta, bg_alpha,
                                 bg_beta, p_fb, p_bf, domain_score_cutoff);
};

struct hmr_summary {
  explicit hmr_summary(const std::vector<Interval6> &hmrs) :
    hmr_count{std::size(hmrs)},
    hmr_total_size{std::accumulate(
      std::cbegin(hmrs), std::cend(hmrs), 0ul,
      [](const std::uint64_t t, const Interval6 &p) { return t + size(p); })} {
    hmr_mean_size = static_cast<double>(hmr_total_size) /
                    std::max(static_cast<std::uint32_t>(hmr_count), 1u);
  }
  // hmr_count is the number of identified HMRs.
  std::uint64_t hmr_count{};
  // total_hmr_size is the sum of the sizes of the identified HMRs
  std::uint64_t hmr_total_size{};
  // mean_hmr_size is the mean size of the identified HMRs
  double hmr_mean_size{};

  auto
  tostring() -> std::string {
    std::ostringstream oss;
    oss << "hmr_count: " << hmr_count << '\n'
        << "hmr_total_size: " << hmr_total_size << '\n'
        << "hmr_mean_size: " << std::fixed << std::setprecision(2)
        << hmr_mean_size;
    return oss.str();
  }
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(hmr_summary, hmr_count, hmr_total_size,
                                 hmr_mean_size);
};

static auto
to_interval(const MSite &s) -> Interval6 {
  const auto p = static_cast<std::uint32_t>(s.pos);
  return {s.chrom, p, p + 1, {}, {}, '+'};
}

static auto
get_stepup_cutoff(std::vector<double> scores, const double cutoff) -> double {
  if (cutoff <= 0)
    return std::numeric_limits<double>::max();
  if (cutoff > 1)
    return std::numeric_limits<double>::min();

  const std::uint32_t n = std::size(scores);
  std::sort(begin(scores), std::end(scores));
  std::uint32_t i = 1;
  while (i < n && scores[i - 1] < (cutoff * i) / n)
    ++i;
  return scores[i - 1];
}

static auto
get_domain_scores(const std::vector<bool> &state_ids,
                  const std::vector<std::pair<double, double>> &meth,
                  const std::vector<std::size_t> &reset_points)
  -> std::vector<double> {
  std::size_t reset_idx = 1;
  bool in_domain{};
  double score{};
  std::vector<double> domain_scores;
  for (std::size_t i = 0; i < std::size(state_ids); ++i) {
    if (reset_points[reset_idx] == i) {
      if (in_domain) {
        in_domain = false;
        domain_scores.push_back(score);
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
      domain_scores.push_back(score);
      score = 0;
    }
  }
  if (in_domain)
    domain_scores.push_back(score);
  return domain_scores;
}

static void
build_domains(const std::vector<MSite> &cpgs,
              const std::vector<std::size_t> &reset_points,
              const std::vector<bool> &state_ids,
              std::vector<Interval6> &domains) {
  std::uint32_t n_cpgs{};
  std::uint32_t prev_stop{};
  bool in_domain = false;
  // assume reset points include 0 and size(cpgs)
  auto r_itr = std::cbegin(reset_points) + 1;
  for (std::uint32_t i = 0; i < std::size(state_ids); ++i) {
    if (i == *r_itr) {
      if (in_domain) {
        in_domain = false;
        domains.back().stop = prev_stop;
        domains.back().score = n_cpgs;
        n_cpgs = 0;
      }
      ++r_itr;
    }
    if (state_ids[i]) {  // currently in an hmr
      if (!in_domain) {
        in_domain = true;
        domains.push_back(to_interval(cpgs[i]));
        domains.back().name = "HYPO" + std::to_string(std::size(domains));
      }
      ++n_cpgs;
    }
    else if (in_domain) {
      in_domain = false;
      domains.back().stop = prev_stop;
      domains.back().score = n_cpgs;
      n_cpgs = 0;
    }
    prev_stop = cpgs[i].pos + 1;
  }
  if (in_domain) {
    domains.back().stop = prev_stop;
    domains.back().score = n_cpgs;
  }
}

template <class T, class U>
static auto
separate_regions(const bool verbose, const std::size_t desert_size,
                 std::vector<MSite> &cpgs, std::vector<T> &meth,
                 std::vector<U> &reads) -> std::vector<std::size_t> {
  if (verbose)
    std::cerr << "[separating by cpg desert]\n";
  // eliminate the zero-read cpgs
  std::int64_t j = 0;
  for (std::size_t i = 0; i < std::size(cpgs); ++i)
    if (reads[i] > 0) {
      cpgs[j] = cpgs[i];
      meth[j] = meth[i];
      reads[j] = reads[i];
      ++j;
    }
  cpgs.erase(std::begin(cpgs) + j, std::end(cpgs));
  meth.erase(std::begin(meth) + j, std::end(meth));
  reads.erase(std::begin(reads) + j, std::end(reads));

  std::vector<std::size_t> reset_points;

  std::size_t total_bases{};
  std::size_t bases_in_deserts{};
  // segregate cpgs
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

  if (verbose)
    std::cerr << "[cpgs retained: " << std::size(cpgs) << "]\n"
              << "[deserts removed: " << std::size(reset_points) - 2 << "]\n"
              << "[genome fraction covered: "
              << 1.0 - (static_cast<double>(bases_in_deserts) /
                        static_cast<double>(total_bases))
              << "]\n";
  return reset_points;
}

/* function to "fold" the methylation profile so that the middle methylation
 * becomes lower methylation, and both the low and high methylation become
 * high. this method actually seems to work.
 */
static void
make_partial_meth(const std::vector<std::uint32_t> &reads,
                  std::vector<std::pair<double, double>> &meth) {
  static constexpr auto one_half = 0.5;
  for (std::size_t i = 0; i < std::size(reads); ++i) {
    const auto r = reads[i];
    double m = meth[i].first / r;
    m = (m <= one_half) ? (1.0 - 2 * m) : (1.0 - 2 * (1.0 - m));
    meth[i] = {r * m, r * (1.0 - m)};
  }
}

[[nodiscard]] static auto
shuffle_cpgs(const std::size_t rng_seed, const TwoStateHMM &hmm,
             std::vector<std::pair<double, double>> meth,
             const std::vector<std::size_t> &reset_points, const double p_fb,
             const double p_bf, const double fg_alpha, const double fg_beta,
             const double bg_alpha,
             const double bg_beta) -> std::vector<double> {
  auto eng = std::default_random_engine(rng_seed);
  std::shuffle(std::begin(meth), std::end(meth), eng);
  std::vector<bool> state_ids;
  std::vector<double> scores;
  hmm.PosteriorDecoding(meth, reset_points, p_fb, p_bf, fg_alpha, fg_beta,
                        bg_alpha, bg_beta, state_ids, scores);
  auto domain_scores = get_domain_scores(state_ids, meth, reset_points);
  std::sort(std::begin(domain_scores), std::end(domain_scores));
  return domain_scores;
}

[[nodiscard]] static auto
assign_p_values(const std::vector<double> &random,
                const std::vector<double> &observed) -> std::vector<double> {
  const auto n_rand =
    random.empty() ? 1.0 : static_cast<double>(std::size(random));
  const auto n_scores = std::size(observed);
  const auto r_end = std::cend(random);
  std::vector<double> p_values(n_scores);
  for (auto i = 0u; i < n_scores; ++i) {
    const auto ub = std::upper_bound(std::cbegin(random), r_end, observed[i]);
    p_values[i] = static_cast<double>(std::distance(ub, r_end)) / n_rand;
  }
  return p_values;
}

[[nodiscard]] static auto
read_params_file(const std::string &params_file) -> hmr_params {
  std::ifstream in(params_file);
  if (!in)
    throw std::runtime_error("failed to open params file: " + params_file);
  const nlohmann::json data = nlohmann::json::parse(in);
  try {
    hmr_params p(data);
    return p;
  }
  catch (const std::exception &e) {
    const auto msg = "failed to parse parameters file: " + params_file;
    throw std::runtime_error(std::string(e.what()) + "\n" + msg);
  }
}

static void
write_params_file(const std::string &outfile, const hmr_params &params) {
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open params output file: " + outfile);
  out << nlohmann::json(params);
}

[[nodiscard]]
static inline auto
get_n_meth(const MSite &s) -> double {
  return s.meth * static_cast<double>(s.n_reads);
}

[[nodiscard]]
static inline auto
get_n_unmeth(const MSite &s) -> double {
  return (1.0 - s.meth) * static_cast<double>(s.n_reads);
}

static auto
load_cpgs(const std::string &cpgs_file) -> std::vector<MSite> {
  bamxx::bgzf_file in(cpgs_file, "r");
  if (!in)
    throw std::runtime_error("failed opening file: " + cpgs_file);
  if (get_has_counts_header(cpgs_file))
    skip_counts_header(in);
  std::vector<MSite> cpgs;
  MSite prev_site, the_site;
  while (read_site(in, the_site)) {
    if (!the_site.is_cpg() || distance(prev_site, the_site) < 2)
      throw std::runtime_error("input is not symmetric-CpGs: " + cpgs_file);
    cpgs.push_back(the_site);
    prev_site = the_site;
  }
  return cpgs;
}

template <class T>
static void
check_sorted_within_chroms(T first, const T last) {
  // empty, or a single element
  if (last <= first + 1)
    return;
  std::unordered_set<std::string> chroms_seen;
  std::string cur_chrom;
  for (auto next = first + 1; next != last; ++next, ++first) {
    if (next->chrom != cur_chrom) {
      if (chroms_seen.find(next->chrom) != std::end(chroms_seen))
        throw std::runtime_error("input not grouped by chromosomes. "
                                 "Error in the following line:\n" +
                                 next->tostring());
      cur_chrom = next->chrom;
      chroms_seen.insert(cur_chrom);
      // don't check ordering if chroms have changed
    }
    else {  // first and next are in the same chrom
      if (first->pos >= next->pos) {
        throw std::runtime_error("input file not sorted properly. "
                                 "Error in the following lines:\n" +
                                 first->tostring() + "\n" + next->tostring());
      }
    }
  }
}

auto
main_hmr(int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    static constexpr auto min_coverage = 1.0;
    static constexpr auto tolerance = 1e-10;  // corrections for small values
    static constexpr auto p_value_cutoff = 0.01;

    std::string outfile;
    std::string post_outfile;
    std::string meth_outfile;

    std::size_t desert_size{1000};   // NOLINT(*-avoid-magic-numbers)
    std::size_t max_iterations{10};  // NOLINT(*-avoid-magic-numbers)
    std::size_t rng_seed{408};       // NOLINT(*-avoid-magic-numbers)

    // run mode flags
    bool verbose{false};
    bool PARTIAL_METH{false};
    bool allow_extra_fields{false};

    std::string summary_file;
    std::string params_in_file;
    std::string params_out_file;

    constexpr auto description =
      R"(Identify HMRs in methylomes. Methylation must be provided
      in the methcounts format (chrom, position, strand, context,
      methylation, reads). See the methcounts documentation for
      details. This program assumes only data at CpG sites and that
      strands are collapsed so only the positive site appears in the
      file.)";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<methylation-file>");
    opt_parse.add_opt("out", 'o', "output file", true, outfile);
    opt_parse.add_opt("desert", 'd', "max dist btwn covered cpgs in HMR", false,
                      desert_size);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_iterations);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    opt_parse.add_opt("partial", '\0', "identify PMRs instead of HMRs", false,
                      PARTIAL_METH);
    opt_parse.add_opt("post-hypo", '\0',
                      "output file for single-CpG posterior "
                      "hypomethylation probability (default: none)",
                      false, post_outfile);
    opt_parse.add_opt("post-meth", '\0',
                      "output file for single-CpG posteiror "
                      "methylation probability (default: none)",
                      false, meth_outfile);
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

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open output file: " + outfile);

    // separate the regions by chrom and by desert
    if (verbose)
      std::cerr << "[reading methylation levels]\n";
    const auto orig_cpgs = load_cpgs(cpgs_file);
    if (verbose)
      std::cerr << "[checking if input is properly formatted]\n";
    check_sorted_within_chroms(std::cbegin(orig_cpgs), std::cend(orig_cpgs));

    auto cpgs = orig_cpgs;

    std::vector<std::pair<double, double>> meth;
    std::vector<std::uint32_t> reads;
    for (const auto &cpg : cpgs) {
      reads.push_back(cpg.n_reads);
      meth.emplace_back(get_n_meth(cpg), get_n_unmeth(cpg));
    }

    if (PARTIAL_METH)
      make_partial_meth(reads, meth);

    const auto mean_coverage = get_mean(std::cbegin(reads), std::cend(reads));

    if (verbose)
      std::cerr << "[total_cpgs=" << std::size(cpgs) << "]\n"
                << "[mean_coverage=" << mean_coverage << "]\n";

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

    // separate the regions by chrom and by desert, and eliminate those
    // isolated CpGs
    const auto reset_points =
      separate_regions(verbose, desert_size, cpgs, meth, reads);

    const TwoStateHMM hmm(tolerance, max_iterations, verbose);

    hmr_params params;

    if (!params_in_file.empty()) {  // read parameters file
      params = read_params_file(params_in_file);
      max_iterations = 0;
    }
    else
      params.init(mean_coverage);

    if (max_iterations > 0)
      hmm.BaumWelchTraining(meth, reset_points, params.p_fb, params.p_bf,
                            params.fg_alpha, params.fg_beta, params.bg_alpha,
                            params.bg_beta);
    // decode the domains
    std::vector<bool> state_ids;
    std::vector<double> posteriors;
    hmm.PosteriorDecoding(meth, reset_points, params.p_fb, params.p_bf,
                          params.fg_alpha, params.fg_beta, params.bg_alpha,
                          params.bg_beta, state_ids, posteriors);

    const auto domain_scores = get_domain_scores(state_ids, meth, reset_points);
    const auto random_scores = shuffle_cpgs(
      rng_seed, hmm, meth, reset_points, params.p_fb, params.p_bf,
      params.fg_alpha, params.fg_beta, params.bg_alpha, params.bg_beta);
    const auto p_values = assign_p_values(random_scores, domain_scores);

    if (params.domain_score_cutoff == std::numeric_limits<double>::max() &&
        !domain_scores.empty())
      params.domain_score_cutoff = get_stepup_cutoff(p_values, p_value_cutoff);

    // write parameters if requested
    if (!params_out_file.empty())
      write_params_file(params_out_file, params);

    std::vector<Interval6> domains;
    if (!domain_scores.empty())
      build_domains(cpgs, reset_points, state_ids, domains);

    {
      std::size_t good_hmr_count = 0;
      for (auto i = 0u; i < std::size(domains); ++i)
        if (p_values[i] < params.domain_score_cutoff) {
          domains[good_hmr_count] = domains[i];
          domains[good_hmr_count].name =
            "HYPO" + std::to_string(good_hmr_count);
          ++good_hmr_count;
        }
      domains.resize(good_hmr_count);
      for (const auto &d : domains)
        out << to_string(d) << '\n';
    }

    if (!post_outfile.empty()) {
      if (verbose)
        std::cerr << "[writing=" << post_outfile << "]\n";
      std::ofstream out_post(post_outfile);
      if (!out_post)
        throw std::runtime_error("failed to open output file: " + post_outfile);
      auto j = 0u;
      for (auto cpg : orig_cpgs) {
        cpg.meth = (cpg.n_reads > 0) ? 1.0 - posteriors[j++] : 0.0;
        out_post << cpg << '\n';
      }
    }

    if (!meth_outfile.empty()) {
      std::ofstream out_meth(meth_outfile);
      if (!out_meth)
        throw std::runtime_error("failed to open output file: " + meth_outfile);
      if (verbose)
        std::cerr << "[writing=" << meth_outfile << "]\n";
      auto j = 0u;
      for (auto cpg : orig_cpgs) {
        cpg.meth = cpg.n_reads > 0 ? params.est_meth(posteriors[j++]) : 0.0;
        out_meth << cpg << '\n';
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
