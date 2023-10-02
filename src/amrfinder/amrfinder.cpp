/* amrfinder: program for resolving epialleles in a sliding window
 * along a chromosome.
 *
 * Copyright (C) 2014-2022 University of Southern California and
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

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include <bamxx.hpp>

#include "GenomicRegion.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "EpireadStats.hpp"
#include "GenomicRegion.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::unordered_map;
using std::runtime_error;

using bamxx::bgzf_file;

#include <bamxx.hpp>
#include <sstream>

static inline bamxx::bgzf_file &
read_epiread(bamxx::bgzf_file &f, epiread &er) {
  std::string line;
  if (getline(f, line))
    er = epiread(line);
  return f;
}

static inline bool
validate_epiread_bgzf_file(const string &filename) {
  constexpr size_t max_lines_to_validate = 10000;
  bgzf_file in(filename, "r");
  if (!in) throw std::runtime_error("failed to open file: " + filename);

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

////////////////////////////////////////////////////////////////////////
/// CODE BELOW HERE IS FOR FILTERING THE AMR WINDOWS AND MERGING THEM
/// TO OBTAIN FINAL AMRS

static void
eliminate_amrs_by_cutoff(const double cutoff, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 0; i < amrs.size(); ++i)
    if (amrs[i].get_score() < cutoff) {
      amrs[j] = amrs[i];
      ++j;
    }
  amrs.erase(begin(amrs) + j, end(amrs));
}


static void
eliminate_amrs_by_size(const size_t min_size, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 0; i < amrs.size(); ++i)
    if (amrs[i].get_width() >= min_size) {
      amrs[j] = amrs[i];
      ++j;
    }
  amrs.erase(begin(amrs) + j, end(amrs));
}

/* merges amrs that overlap, and uses the minimum score among them */
static void
collapse_amrs(vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    // +1 below is because intervals in terms of CpGs are inclusive
    if (amrs[j].same_chrom(amrs[i]) &&
        amrs[j].get_end() + 1 >= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
      amrs[j].set_score(std::min(amrs[i].get_score(), amrs[j].get_score()));
    }
    else {
      ++j;
      amrs[j] = amrs[i];
    }
  ++j;
  amrs.erase(begin(amrs) + j, end(amrs));
}


/* merges amrs within some pre-defined distance */
static void
merge_amrs(const size_t gap_limit, vector<GenomicRegion> &amrs) {
  size_t j = 0;
  for (size_t i = 1; i < amrs.size(); ++i)
    // check distance between two amrs is greater than gap limit
    if (amrs[j].same_chrom(amrs[i]) &&
        amrs[j].get_end() + gap_limit >= amrs[i].get_start()) {
      amrs[j].set_end(amrs[i].get_end());
      amrs[j].set_score(std::min(amrs[i].get_score(), amrs[j].get_score()));
    }
    else {
      ++j;
      amrs[j] = amrs[i];
    }
  ++j;
  amrs.erase(begin(amrs) + j, end(amrs));
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////
/////   CODE FOR RANDOMIZING THE READS TO GET EXPECTED NUMBER OF
/////   IDENTIFIED AMRS
/////

// static void
// set_read_states(vector<vector<char> > &state_counts, vector<epiread> &reads) {
//   for (size_t i = 0; i < reads.size(); ++i) {
//     const size_t offset = reads[i].pos;
//     for (size_t j = 0; j < reads[i].length(); ++j) {
//       reads[i].seq[j] = state_counts[offset + j].back();
//       state_counts[offset + j].pop_back();
//     }
//   }
// }

// static void
// get_state_counts(const vector<epiread> &reads, const size_t total_cpgs,
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
// randomize_read_states(vector<epiread> &reads) {
//   /**/srand(time(0) + getpid());
//   const size_t total_cpgs = get_n_cpgs(reads);
//   vector<vector<char> > state_counts;
//   get_state_counts(reads, total_cpgs, state_counts);
//   set_read_states(state_counts, reads);
// }


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////////////
//////////////  CODE FOR CONVERTING BETWEEN CPG AND BASE PAIR
//////////////  COORDINATES BELOW HERE
//////////////

// ADS: need to check that the use if this `is_cpg` function will not
// overrun the string; why is this used and not `is_cpg` function used
// in other programs?
inline static bool
is_cpg(const string &s, const size_t idx) {
  return toupper(s[idx]) == 'C' && toupper(s[idx + 1]) == 'G';
}


static void
collect_cpgs(const string &s, unordered_map<size_t, size_t> &cpgs) {
  const size_t lim = s.length() - 1;
  size_t cpg_count = 0;
  for (size_t i = 0; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs[cpg_count++] = i;
}


static void
convert_coordinates(const unordered_map<size_t, size_t> &cpgs,
                    GenomicRegion &region)  {
  const auto start_itr = cpgs.find(region.get_start());
  const auto end_itr = cpgs.find(region.get_end());
  if (start_itr == end(cpgs) || end_itr == end(cpgs))
    throw runtime_error("could not convert:\n" + region.tostring());
  region.set_start(start_itr->second);
  region.set_end(end_itr->second);
}


static void
get_chrom(const string &chrom_name,
          const unordered_map<string, string> &chrom_files,
          GenomicRegion &chrom_region, string &chrom) {
  // get the right filename
  auto fn = chrom_files.find(chrom_name);
  if (fn == end(chrom_files))
    throw runtime_error("could not find chrom: " + chrom_name);
  chrom.clear();
  read_fasta_file(fn->second, chrom_name, chrom);
  if (chrom.empty())
    throw runtime_error("could not find chrom: " + chrom_name);

  chrom_region.set_chrom(chrom_name);
}

static void
convert_coordinates(const bool VERBOSE, const string chroms_dir,
                    const string fasta_suffix, vector<GenomicRegion> &amrs) {

  unordered_map<string, string> chrom_files;
  identify_and_read_chromosomes(chroms_dir, fasta_suffix, chrom_files);
  if (VERBOSE)
    cerr << "CHROMS:\t" << chrom_files.size() << endl;

  unordered_map<size_t, size_t> cpgs;
  string chrom;
  // ADS: this should be done with a chrom_name and not using the
  // dummy "region"
  GenomicRegion chrom_region("chr0", 0, 0);
  for (size_t i = 0; i < amrs.size(); ++i) {
    if (!amrs[i].same_chrom(chrom_region)) {
      get_chrom(amrs[i].get_chrom(), chrom_files, chrom_region, chrom);
      collect_cpgs(chrom, cpgs);
      if (VERBOSE)
        cerr << "CONVERTING: " << chrom_region.get_chrom() << endl;
    }
    convert_coordinates(cpgs, amrs[i]);
  }
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
/////////
/////////  CODE FOR DOING THE SLIDING WINDOW STUFF BELOW HERE
/////////


/* make sure the current read is shortened to only include positions
   within the relevant window */
static void
clip_read(const size_t start_pos, const size_t end_pos, epiread &r) {
  if (r.pos < start_pos) {
    r.seq = r.seq.substr(start_pos - r.pos);
    r.pos = start_pos;
  }
  if (r.end() > end_pos)
    r.seq = r.seq.substr(0, end_pos - r.pos);
}


static void
get_current_epireads(const vector<epiread> &epireads,
                     const size_t max_epiread_len,
                     const size_t cpg_window, const size_t start_pos,
                     size_t &read_id, vector<epiread> &current_epireads) {
  while (read_id < epireads.size() &&
         epireads[read_id].pos + max_epiread_len <=start_pos)
    ++read_id;

  const size_t end_pos = start_pos + cpg_window;
  for (size_t i = read_id; (i < epireads.size() &&
                            epireads[i].pos < end_pos); ++i) {
    if (epireads[i].end() > start_pos) {
      current_epireads.push_back(epireads[i]);
      clip_read(start_pos, end_pos, current_epireads.back());
    }
  }
}


static size_t
total_states(const vector<epiread> &epireads) {
  size_t total = 0;
  for (size_t i = 0; i < epireads.size(); ++i)
    total += epireads[i].length();
  return total;
}


static void
add_amr(const string &chrom_name, const size_t start_cpg,
        const size_t cpg_window, const vector<epiread> &reads,
        const double score, vector<GenomicRegion> &amrs) {
  static const string name_label("AMR");
  const size_t end_cpg = start_cpg + cpg_window - 1;
  const string amr_name(name_label + toa(amrs.size()) + ":" + toa(reads.size()));
  amrs.push_back(GenomicRegion(chrom_name, start_cpg, end_cpg,
                               amr_name, score, '+'));
}


static size_t
process_chrom(const bool VERBOSE, const bool PROGRESS,
              const size_t min_obs_per_cpg, const size_t window_size,
              const EpireadStats &epistat, const string &chrom_name,
              const vector<epiread> &epireads, vector<GenomicRegion> &amrs) {
  size_t max_epiread_len = 0;
  for (size_t i = 0; i < epireads.size(); ++i)
    max_epiread_len = std::max(max_epiread_len, epireads[i].length());
  const size_t min_obs_per_window = window_size*min_obs_per_cpg;

  const size_t chrom_cpgs = get_n_cpgs(epireads);
  if (VERBOSE)
    cerr << "processing " << chrom_name << " "
         << "[reads: " << epireads.size() << "] "
         << "[cpgs: " << chrom_cpgs << "]" << endl;

  ProgressBar progress(epireads.size());

  size_t windows_tested = 0;
  size_t start_idx = 0;
  const size_t lim = chrom_cpgs - window_size + 1;
  vector<epiread> current_epireads;
  for (size_t i = 0; i < lim && start_idx < epireads.size(); ++i) {

    if (PROGRESS && progress.time_to_report(start_idx))
      progress.report(cerr, start_idx);

    current_epireads.clear();
    get_current_epireads(epireads, max_epiread_len,
                         window_size, i, start_idx, current_epireads);

    if (total_states(current_epireads) >= min_obs_per_window) {
      bool is_significant = false;
      const double score = epistat.test_asm(current_epireads, is_significant);
      if (is_significant)
        add_amr(chrom_name, i, window_size, current_epireads, score, amrs);
      ++windows_tested;
    }
  }
  if (PROGRESS)
    progress.report(cerr, epireads.size());
  return windows_tested;
}


int
main_amrfinder(int argc, const char **argv) {
  try {

    const string description =
      "identify regions of allele-specific methylation";

    static const string fasta_suffix = "fa";

    bool VERBOSE = false;
    bool PROGRESS = false;

    string outfile;
    string chroms_dir;

    size_t max_itr = 10;
    size_t window_size = 10;
    size_t gap_limit = 1000;

    double high_prob = 0.75, low_prob = 0.25;
    size_t min_obs_per_cpg = 4;
    double critical_value = 0.01;

    // ADS: below, for when the time comes
    // auto eng = std::default_random_engine(rng_seed);
    // std::shuffle(begin(things), end(things), eng);

    // bool RANDOMIZE_READS = false;
    bool use_bic = false;
    bool use_fdr = true;
    bool apply_correction = false;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description, "<epireads>");
    opt_parse.add_opt("output", 'o', "output file", true, outfile);
    opt_parse.add_opt("chrom", 'c', "genome sequence file/directory",
                      true, chroms_dir);
    opt_parse.add_opt("itr", 'i', "max iterations", false, max_itr);
    opt_parse.add_opt("window", 'w', "size of sliding window (in CpGs)",
                      false, window_size);
    opt_parse.add_opt("min-cov", 'm', "min coverage per cpg to test windows",
                      false, min_obs_per_cpg);
    opt_parse.add_opt("gap", 'g', "min allowed gap between amrs (in bp)",
                      false, gap_limit);
    opt_parse.add_opt("crit", 'C', "critical p-value cutoff",
                      false, critical_value);
    opt_parse.add_opt("pvals", 'h', "adjusts p-values using Hochberg step-up",
                      false, apply_correction);
    opt_parse.add_opt("nofdr", 'f', "omits FDR procedure", false, use_fdr);
    opt_parse.add_opt("bic", 'b', "use BIC to compare models", false, use_bic);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    opt_parse.add_opt("progress", 'P', "print progress info", false, PROGRESS);
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

    if (!validate_epiread_bgzf_file(reads_file))
      throw runtime_error("invalid states file: " + reads_file);

    if (VERBOSE)
      cerr << "AMR TESTING OPTIONS: "
           << "[test=" << (use_bic ? "BIC" : "LRT") << "] "
           << "[iterations=" << max_itr << "]" << endl;

    const EpireadStats epistat(low_prob, high_prob,
                               critical_value, max_itr, use_bic);

    bgzf_file in(reads_file, "r");
    if (!in) throw runtime_error("failed to open input file: " + reads_file);

    vector<GenomicRegion> amrs;
    size_t windows_tested = 0;
    epiread er;
    vector<epiread> epireads;
    string prev_chrom, curr_chrom, tmp_states;

    while (read_epiread(in, er)) {
      curr_chrom = er.chr;
      if (!epireads.empty() && curr_chrom != prev_chrom) {
        windows_tested +=
          process_chrom(VERBOSE, PROGRESS, min_obs_per_cpg, window_size,
                        epistat, prev_chrom, epireads, amrs);
        epireads.clear();
      }
      epireads.push_back(er);
      prev_chrom = curr_chrom;
    }

    if (!epireads.empty())
      windows_tested +=
        process_chrom(VERBOSE, PROGRESS, min_obs_per_cpg, window_size,
                      epistat, prev_chrom, epireads, amrs);

    if (VERBOSE)
      cerr << "========= POST PROCESSING =========" << endl;

    const size_t windows_accepted = amrs.size();
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
      vector<double> pvals;
      for (size_t i = 0; i < amrs.size(); ++i)
        pvals.push_back(amrs[i].get_score());

      // if the pvalues are not BIC scores, get an FDR cutoff
      // corresponding to the given critical value
      const double fdr_cutoff = use_bic ? 0.0 :
        smithlab::get_fdr_cutoff(windows_tested, pvals, critical_value);

      // if we are not using BIC, and if corrected p-values are
      // requested, then adjust the p-values
      if (!use_bic && apply_correction) {
        smithlab::correct_pvals(windows_tested, pvals);
        for (size_t i = 0; i < pvals.size(); ++i)
          amrs[i].set_score(pvals[i]);
      }

      // ADS: not sure it's a good idea in this next collapse function
      // to keep the lowest among all p-values for merged regions
      collapse_amrs(amrs);
      const size_t collapsed_amrs = amrs.size();

      convert_coordinates(VERBOSE, chroms_dir, fasta_suffix, amrs);

      // ADS: merge AMRs if they are sufficiently close; the distance
      // should not be too large, and a distance of 1000 bp is likely
      // too large.
      merge_amrs(gap_limit, amrs);
      const size_t merged_amrs = amrs.size();

      // if BIC was not requested, then eliminate AMRs based on either
      // the p-value cutoff, or with the FDR-based cutoff, if it was
      // requested.
      if (!use_bic) {
        if (apply_correction || !use_fdr)
          eliminate_amrs_by_cutoff(critical_value, amrs);
        else eliminate_amrs_by_cutoff(fdr_cutoff, amrs);
      }

      const size_t amrs_passing_fdr = amrs.size();

      // ADS: eliminating AMRs based on their size makes sense, but
      // not if that size is tied to the gap between AMRs we would
      // merge. There is no symmetry for these.
      eliminate_amrs_by_size(gap_limit/2, amrs);

      if (VERBOSE) {
        cerr << "WINDOWS TESTED: " << windows_tested << endl
             << "WINDOWS ACCEPTED: " << windows_accepted << endl
             << "COLLAPSED WINDOWS: " << collapsed_amrs << endl
             << "MERGED WINDOWS: " << merged_amrs << endl;
        if (use_fdr)
          cerr  << "FDR CUTOFF: " << fdr_cutoff << endl
                << "WINDOWS PASSING FDR: " << amrs_passing_fdr << endl;
        cerr << "AMRS (WINDOWS PASSING MINIMUM SIZE): " << amrs.size() << endl;
      }
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? cout.rdbuf() : of.rdbuf());
      copy(begin(amrs), end(amrs),
           std::ostream_iterator<GenomicRegion>(out, "\n"));
    }
    else if (VERBOSE)
      cerr << "no AMRs found" << endl;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
