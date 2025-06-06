/* Copyright (C) 2022-2024 Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Guilherme Sena
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

#include <config.h>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::begin;
using std::cout;
using std::end;
using std::endl;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

static const string PROGRAM_NAME = "dnmtools";

struct dnmtools_command {
  string tag;
  string description;
  std::function<int(const int, const char **)> fun;

  auto
  operator()(const int argc, const char **argv) const -> int {
    return fun(argc - 1, argv + 1);
  }
};

auto
operator<<(std::ostream &out, const dnmtools_command &cmd) -> std::ostream & {
  static const size_t pad_size = 4;
  static const size_t offset = 15;
  static const string pad(pad_size, ' ');
  return out << pad << std::left << std::setw(offset) << (cmd.tag + ":")
             << cmd.description;
}

// ADS: not sure of best way to acquire these below beyond simply
// declaring them here
int
abismal(int argc, const char **argv);
int
abismalidx(int argc, const char **argv);
int
simreads(int argc, const char **argv);
int
main_counts(int argc, const char **argv);
int
main_nanopore(int argc, const char **argv);
int
main_allelicmeth(int argc, const char **argv);
int
main_amrfinder(int argc, const char **argv);
int
main_amrtester(int argc, const char **argv);
int
main_bsrate(int argc, const char **argv);
int
main_hmr(int argc, const char **argv);
int
main_hmr_rep(int argc, const char **argv);
int
main_hypermr(int argc, const char **argv);
int
main_levels(int argc, const char **argv);
int
main_methentropy(int argc, const char **argv);
int
main_methstates(int argc, const char **argv);
int
main_multimethstat(int argc, const char **argv);
int
main_pmd(int argc, const char **argv);
int
main_roimethstat(int argc, const char **argv);
int
main_cpgbins(int argc, const char **argv);
int
main_mlml(int argc, const char **argv);
int
main_dmr(int argc, const char **argv);
int
main_methdiff(int argc, const char **argv);
int
main_radmeth_adjust(int argc, const char **argv);
int
main_radmeth(int argc, const char **argv);
int
main_radmeth_merge(int argc, const char **argv);
int
main_clean_hairpins(int argc, const char **argv);
int
main_uniq(int argc, const char **argv);
int
main_fast_liftover(int argc, const char **argv);
int
main_format(int argc, const char **argv);
int
main_guessprotocol(int argc, const char **argv);
int
main_lift_filter(int argc, const char **argv);
int
main_merge_bsrate(int argc, const char **argv);
int
main_merge_methcounts(int argc, const char **argv);
int
main_selectsites(int argc, const char **argv);
int
main_symmetric_cpgs(int argc, const char **argv);
int
metagene(int argc, const char **argv);
int
main_covered(int argc, const char **argv);
int
main_xcounts(int argc, const char **argv);
int
main_unxcounts(int argc, const char **argv);
int
main_recovered(int argc, const char **argv);
int
kmersites(int argc, const char **argv);

void
print_help(
  const vector<pair<string, vector<dnmtools_command>>> &command_groups) {
  cout << "Program: " << PROGRAM_NAME << "\n"
       << "Version: " << VERSION << "\n"
       << "Usage: " << PROGRAM_NAME << " <command> [options]\n"
       << "Commands:" << endl;
  for (auto &&g : command_groups) {
    cout << "  " << g.first << ":" << endl;
    for (auto &&c : g.second)
      cout << c << endl;
    cout << endl;
  }
}

int
main(int argc, const char **argv) {
  try {
    vector<pair<string, vector<dnmtools_command>>> command_groups = {
      // clang-format off
{{"mapping",
 {{{"abismal",    "map FASTQ reads to a FASTA reference genome or an index", abismal},
   {"abismalidx", "convert a FASTA reference genome to an abismal index",    abismalidx},
   {"simreads",   "simulate reads in a FASTA reference genome",              simreads}}}},

{"methylome construction",
 {{{"format",    "convert SAM/BAM mapped bs-seq reads to standard dnmtools format", main_format},
   {"uniq",      "remove duplicate reads from sorted mapped reads",           main_uniq},
   {"bsrate",    "compute the BS conversion rate from BS-seq reads mapped to a genome",  main_bsrate},
   {"counts",    "get methylation levels from mapped WGBS reads",             main_counts},
   {"nano",      "get methylation levels from mapped nanopore reads",         main_nanopore},
   {"sym",       "get CpG sites and make methylation levels symmetric",       main_symmetric_cpgs},
   {"levels",    "compute methylation summary statistics from a counts file", main_levels}}}},

{"methylome analysis",
 {{{"hmr",       "identify hypomethylated regions", main_hmr},
   {"hmr-rep",   "identify hypomethylated regions in a set of replicate methylomes", main_hmr_rep},
   {"hypermr",   "identify hypermethylated regions in plant methylomes", main_hypermr},
   {"entropy",   "compute methylation entropy in sliding window", main_methentropy},
   {"pmd",       "identify partially methylated domains", main_pmd},
   {"roi",       "get average CpG methylation in each of a set of genomic interval", main_roimethstat},
   {"cpgbins",    "get average CpG methylation in genomic bins", main_cpgbins},
   {"multistat", "same as roi except for multiple samples/input files", main_multimethstat},
   {"mlml",      "program to estimate hydroxymethylation levels", main_mlml}}}},

{"allele-specific methylation (ASM)",
 {{{"states",    "convert reads in SAM format into methylation states at CpGs", main_methstates},
   {"allelic",   "get probability of ASM for each pair of neighboring CpGs", main_allelicmeth},
   {"amrfinder", "identify regions of ASM in the genome", main_amrfinder},
   {"amrtester", "test a set of genomic intervals for ASM", main_amrtester}}}},

{"differential methylation (DM)",
 {{{"dmr",       "identify DMRs from genomic intervals and single-CpG DM probabilities",  main_dmr},
   {"diff",      "compute single-CpG DM probability between two methylomes", main_methdiff},
   {"radmeth",   "compute DM probabilities for each CpG using multiple methylomes", main_radmeth},
   {"radadjust", "adjust p-values from radmeth output", main_radmeth_adjust},
   {"radmerge",  "merge significant CpGs in radmeth output", main_radmeth_merge}}}},

{"methylation visualization",
 {{{"fastlift",   "liftover methylation levels between species", main_fast_liftover},
   {"metagene",   "summarize methylation around genomic features", metagene},
   {"liftfilter", "filter CpGs that are not CpGs in the target genome", main_lift_filter}}}},

{"utilities",
 {{{"cleanhp",       "fix and stat invdup/hairping reads", main_clean_hairpins},
   {"guessprotocol", "guess whether protocol is ordinary, pbat or random", main_guessprotocol},
   {"merge-bsrate",  "merge bisulfite conversion rates files from bsrate", main_merge_bsrate},
   {"merge",         "merge multiple counts files into a counts file or a table", main_merge_methcounts},
   {"covered",       "filter a counts file for only covered sites", main_covered},
   {"recovered",     "replace missing sites in a counts file", main_recovered},
   {"xcounts",       "compress counts files by removing information", main_xcounts},
   {"unxcounts",     "reverse the xcounts process yielding a counts file", main_unxcounts},
   {"selectsites",   "sites inside a set of genomic intervals", main_selectsites},
   {"kmersites",     "make track file for sites matching kmer", kmersites}}}}}};
    // clang-format on

    if (argc < 2) {
      print_help(command_groups);
      return EXIT_SUCCESS;
    }

    const auto has_tag = [&](const dnmtools_command &a) {
      return a.tag == argv[1];
    };

    for (auto &g : command_groups) {
      const auto the_cmd = find_if(begin(g.second), end(g.second), has_tag);
      if (the_cmd != end(g.second))
        return (*the_cmd)(argc, argv);
    }

    std::cerr << "ERROR: invalid command " << argv[1] << std::endl;
  }
  catch (const std::exception &e) {
    std::cerr << "ERROR:\t" << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
