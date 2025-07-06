/* Copyright (C) 2013-2025 University of Southern California
 *                         Andrew D Smith and Jenny Qu
 *
 * Author: Jenny Qu and Andrew D. Smith
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

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

[[nodiscard]] static auto
read_fasta(const std::string &fasta_file,
           std::unordered_map<std::string, std::uint32_t> &chrom_names,
           std::vector<std::string> &chroms) {
  std::ifstream in(fasta_file);
  if (!in)
    throw std::runtime_error("Failed to read fasta file: " + fasta_file);

  std::string line;
  std::string seq;
  std::uint32_t index{};
  while (std::getline(in, line)) {
    if (line[0] == '>') {
      const auto name = line.substr(1);
      const auto itr = chrom_names.find(name);
      if (itr != std::cend(chrom_names))
        throw std::runtime_error("Duplicate chrom found: " + name);
      if (!seq.empty()) {
        chroms.push_back(std::move(seq));
        seq.clear();
      }
      chrom_names.emplace(name, index++);
    }
    else
      seq += line;
  }
  if (!seq.empty())
    chroms.push_back(std::move(seq));
}

inline static bool
is_cpg(const std::string &s, const std::size_t idx) {
  return std::toupper(s[idx]) == 'C' && std::toupper(s[idx + 1]) == 'G';
}

static void
build_coordinate_converter(
  const std::string &chrom,
  std::unordered_map<std::size_t, std::size_t> &cpg_lookup) {
  cpg_lookup.clear();
  const auto lim = std::size(chrom) - 1u;
  std::size_t cpg_count = 0;
  for (std::size_t i = 0; i < lim; ++i)
    if (is_cpg(chrom, i))
      cpg_lookup[cpg_count++] = i;
}

static std::size_t
convert_coordinates(const std::unordered_map<std::size_t, std::size_t> &cpgs,
                    const std::size_t position) {
  const std::unordered_map<std::size_t, std::size_t>::const_iterator i(
    cpgs.find(position));
  if (i == cpgs.end())
    throw std::runtime_error("could not convert:\n" + toa(position));
  return i->second;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR WORKING WITH EPIREADS BELOW HERE
///////

struct epiread {
  bool
  operator<(const epiread &other) const {
    return (chr < other.chr || (chr == other.chr && pos < other.pos));
  }
  std::size_t
  get_end() const {
    return pos + std::size(seq);
  }
  std::size_t
  size() const {
    return std::size(seq);
  }
  bool
  flip_states();

  std::string chr;
  std::size_t pos{};
  std::string seq;
};

bool
epiread::flip_states() {
  const std::size_t meth_states_count = std::count(seq.begin(), seq.end(), 'C');
  const auto sz = std::size(seq);
  if (meth_states_count < 0.5 * sz) {
    for (std::size_t i = 0; i < sz; ++i)
      seq[i] = (seq[i] == 'T') ? 'C' : ((seq[i] == 'C') ? 'T' : seq[i]);
    return true;
  }
  return false;
}

static std::istream &
operator>>(std::istream &in, epiread &er) {
  std::string buffer;
  if (getline(in, buffer)) {
    std::istringstream is(buffer);
    if (!(is >> er.chr >> er.pos >> er.seq))
      throw std::runtime_error("malformed epiread line:\n" + buffer);
  }
  return in;
}

static bool
check_sorted(const std::vector<epiread> &epireads) {
  for (std::size_t i = 1; i < std::size(epireads); ++i)
    if (epireads[i] < epireads[i - 1])
      return false;
  return true;
}

// code for entropy calculations below here

static double
compute_prob_read_has_state(const std::vector<double> &site_probs,
                            const std::size_t start_cpg,
                            const std::size_t end_cpg, const std::size_t state,
                            const epiread &er) {
  double prob = 1.0;
  std::size_t cpg_pos = start_cpg;
  const std::size_t lim = end_cpg - start_cpg;
  for (std::size_t i = 0; i < lim; ++i) {
    const char curr_state = ((state & (1ul << i)) > 0) ? 'C' : 'T';
    const std::size_t er_idx = (cpg_pos >= er.pos)
                                 ? cpg_pos - er.pos
                                 : std::numeric_limits<std::size_t>::max();
    if (er_idx < std::size(er) && er.seq[er_idx] != 'N')
      prob *= static_cast<double>(curr_state == er.seq[er_idx]);
    else
      prob *=
        (curr_state == 'C') ? site_probs[cpg_pos] : 1.0 - site_probs[cpg_pos];
    ++cpg_pos;
  }
  return prob;
}

static double
compute_entropy_for_window(const std::vector<double> &site_probs,
                           const std::vector<epiread> &epireads,
                           const std::size_t start_idx,
                           const std::size_t end_idx,
                           const std::size_t start_cpg,
                           const std::size_t end_cpg,
                           std::size_t &reads_in_window) {

  const std::size_t n_states = 1ul << (end_cpg - start_cpg);

  double entropy = 0.0;
  for (std::size_t i = 0; i < n_states; ++i) {

    double state_prob = 0.0;
    reads_in_window = 0;
    for (std::size_t j = start_idx; j < end_idx; ++j)
      if (epireads[j].get_end() > start_cpg && epireads[j].pos < end_cpg) {
        state_prob += compute_prob_read_has_state(site_probs, start_cpg,
                                                  end_cpg, i, epireads[j]);
        ++reads_in_window;
      }
    if (reads_in_window > 0) {
      state_prob /= reads_in_window;
      entropy += (state_prob > 0.0) ? state_prob * log2(state_prob) : 0.0;
    }
  }
  return entropy;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
///////
/////// CODE FOR SLIDING THE WINDOW ALONG THE CHROMOSOME BELOW HERE
///////

/* This function just basically computes the same thing as methcounts
   output, so that unobserved states can be imputed  */
static std::vector<double>
compute_site_probs(const std::size_t n_cpgs,
                   const std::vector<epiread> &epireads) {

  std::vector<double> site_probs(n_cpgs);
  std::vector<std::uint32_t> totals(n_cpgs);

  for (std::size_t i = 0; i < epireads.size(); ++i) {
    const std::size_t len = std::size(epireads[i]);
    std::size_t idx = epireads[i].pos;
    for (std::size_t j = 0; j < len; ++j, ++idx) {
      site_probs[idx] += (epireads[i].seq[j] == 'C');
      totals[idx] += (epireads[i].seq[j] != 'N');
    }
  }

  for (auto i = 0u; i < std::size(site_probs); ++i)
    if (totals[i] > 0.0)
      site_probs[i] /= totals[i];

  return site_probs;
}

static void
move_start_index(const std::size_t max_epiread_len,
                 const std::vector<epiread> &epireads,
                 const std::size_t start_cpg, std::size_t &idx) {
  while (idx < epireads.size() &&
         epireads[idx].get_end() + max_epiread_len <= start_cpg)
    ++idx;
}

static void
move_end_index(const std::vector<epiread> &epireads,
               const std::size_t start_cpg, const std::size_t cpg_window,
               std::size_t &idx) {
  while (idx < epireads.size() && epireads[idx].pos < start_cpg + cpg_window)
    ++idx;
}

static void
process_chrom(const bool VERBOSE, const std::size_t cpg_window,
              const std::vector<epiread> &epireads,
              const std::unordered_map<std::size_t, std::size_t> &cpg_lookup,
              std::ostream &out) {

  const std::string chrom(epireads.front().chr);
  if (!check_sorted(epireads))
    throw std::runtime_error("epireads not sorted in chrom: " + chrom);

  const std::size_t n_cpgs = cpg_lookup.size();
  if (VERBOSE)
    std::cerr << "processing " << chrom << " (cpgs = " << n_cpgs << ")\n";

  const auto site_probs = compute_site_probs(n_cpgs, epireads);

  std::size_t max_epiread_len = 0;
  for (const auto &er : epireads)
    max_epiread_len = std::max(max_epiread_len, std::size(er));

  std::size_t start_cpg = 0;
  std::size_t start_idx = 0, end_idx = 0;
  while (start_cpg + cpg_window < n_cpgs) {

    move_start_index(max_epiread_len, epireads, start_cpg, start_idx);
    move_end_index(epireads, start_cpg, cpg_window, end_idx);

    std::size_t reads_used = 0;
    const double entropy =
      compute_entropy_for_window(site_probs, epireads, start_idx, end_idx,
                                 start_cpg, start_cpg + cpg_window, reads_used);

    out << chrom << '\t'
        << convert_coordinates(cpg_lookup, start_cpg + cpg_window / 2) << '\t'
        << "+\tCpG\t" << entropy << '\t' << reads_used << '\n';

    ++start_cpg;
  }
}

int
main_methentropy(int argc, char *argv[]) {

  try {

    bool VERBOSE = false;
    bool FLIP_MAJORITY_STATE = false;

    std::size_t cpg_window = 4;
    std::string outfile;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]),
                           "compute methylation entropy in sliding window",
                           "<fasta-genome-file> <epireads-file>");
    opt_parse.add_opt("window", 'w', "number of CpGs in sliding window", false,
                      cpg_window);
    opt_parse.add_opt("flip", 'F', "flip read majority state to meth", false,
                      FLIP_MAJORITY_STATE);
    opt_parse.add_opt("output", 'o', "output file (default: stdout)", false,
                      outfile);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
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
    if (leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string genome_file = leftover_args.front();
    const std::string epi_file = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    std::unordered_map<std::string, std::uint32_t> chrom_names;
    std::vector<std::string> chroms;
    read_fasta(genome_file, chrom_names, chroms);

    std::ifstream in(epi_file);
    if (!in)
      throw std::runtime_error("cannot open input file: " + epi_file);

    std::ofstream of;
    if (!outfile.empty())
      of.open(outfile);
    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

    std::unordered_map<std::size_t, std::size_t> cpg_lookup;

    std::vector<epiread> epireads;
    epiread tmp_er;
    const std::string chrom;
    while (in >> tmp_er) {
      const auto &name = tmp_er.chr;
      if (!epireads.empty() && name != epireads.back().chr) {
        const auto chrom = epireads.back().chr;
        const auto itr = chrom_names.find(chrom);
        if (itr == std::cend(chrom_names))
          throw std::runtime_error("chrom not found: " + chrom);
        build_coordinate_converter(chroms[itr->second], cpg_lookup);
        process_chrom(VERBOSE, cpg_window, epireads, cpg_lookup, out);
        epireads.clear();
      }
      if (FLIP_MAJORITY_STATE)
        tmp_er.flip_states();
      epireads.push_back(tmp_er);
    }
    if (!epireads.empty()) {
      const auto chrom = epireads.back().chr;
      const auto itr = chrom_names.find(chrom);
      if (itr == std::cend(chrom_names))
        throw std::runtime_error("chrom not found: " + chrom);
      build_coordinate_converter(chroms[itr->second], cpg_lookup);
      process_chrom(VERBOSE, cpg_window, epireads, cpg_lookup, out);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
