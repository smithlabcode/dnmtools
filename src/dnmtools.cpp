/* Copyright (C) 2022 University of Southern California and
 *                    Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Song Qiang and Guilherme Sena
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

#include <iostream>
#include <string>
#include <cstring>

using std::string;
using std::to_string;
using std::cerr;
using std::cout;
using std::endl;

#define PROGRAM_NAME "dnmtools"
#define PROGRAM_VERSION "1.2.2"

int abismal(int argc, const char **argv);
int abismalidx(int argc, const char **argv);
int simreads(int argc, const char **argv);
int main_methcounts(int argc, const char **argv);
int main_allelicmeth(int argc, const char **argv);
int main_amrfinder(int argc, const char **argv);
int main_amrtester(int argc, const char **argv);
int main_bsrate(int argc, const char **argv);
int main_hmr(int argc, const char **argv);
int main_hmr_rep(int argc, const char **argv);
int main_hypermr(int argc, const char **argv);
int main_levels(int argc, const char **argv);
int main_methentropy(int argc, const char **argv);
int main_methstates(int argc, const char **argv);
int main_multimethstat(int argc, const char **argv);
int main_pmd(int argc, const char **argv);
int main_roimethstat(int argc, const char **argv);
int main_mlml(int argc, const char **argv);
int main_dmr(int argc, const char **argv);
int main_methdiff(int argc, const char **argv);
int main_radmeth_adjust(int argc, const char **argv);
int main_radmeth(int argc, const char **argv);
int main_radmeth_merge(int argc, const char **argv);
int main_clean_hairpins(int argc, const char **argv);
int main_duplicate_remover(int argc, const char **argv);
int main_fast_liftover(int argc, const char **argv);
int main_format_reads(int argc, const char **argv);
int main_guessprotocol(int argc, const char **argv);
int main_lc_approx(int argc, const char **argv);
int main_lift_filter(int argc, const char **argv);
int main_merge_bsrate(int argc, const char **argv);
int main_merge_methcounts(int argc, const char **argv);
int main_selectsites(int argc, const char **argv);
int main_symmetric_cpgs(int argc, const char **argv);

void
print_help() {
  static const string sep = "  ";
  cout <<  "Program: " << PROGRAM_NAME << "\n";
  cout <<  "Version: " << PROGRAM_VERSION << "\n";
  cout <<  "Usage: " << PROGRAM_NAME << " <command> [options]\n";
  cout <<  "Commands:";

  cout <<  "\n" << sep << "read mapping:\n";
  cout <<  sep+sep << "abismal:       map FASTQ reads to a FASTA reference genome or an index\n";
  cout <<  sep+sep << "abismalidx:    convert a FASTA reference genome to an abismal index\n";
  cout <<  sep+sep << "simreads:      simulate reads in a FASTA reference genome\n";

  cout <<  "\n" << sep << "methylome construction:\n";
  cout <<  sep+sep << "format:        convert SAM/BAM mapped bs-seq reads to standard dnmtools format\n";
  cout <<  sep+sep << "uniq:          remove duplicate reads from sorted mapped reads\n";
  cout <<  sep+sep << "bsrate:        compute the BS conversion rate from BS-seq reads mapped to a genome\n";
  cout <<  sep+sep << "counts:        get methylation levels from mapped WGBS reads\n";
  cout <<  sep+sep << "sym:           get CpG sites and make methylation levels symmetric.\n";
  cout <<  sep+sep << "levels:        compute methylation summary statistics from a counts file\n";

  cout <<  "\n" << sep << "methylome analysis:\n";
  cout <<  sep+sep << "hmr:           identify hypomethylated regions.\n";
  cout <<  sep+sep << "hmr-rep:       identify hypomethylated regions in a set of replicate methylomes\n";
  cout <<  sep+sep << "entropy:       compute methylation entropy in sliding window\n";
  cout <<  sep+sep << "multistat:     summarize methylation from to genomic intervals in a BED file.\n";
  cout <<  sep+sep << "pmd:           identify partially methylated domains\n";
  cout <<  sep+sep << "roi:           compute average CpG methylation in each of a set of genomic interval\n";
  cout <<  sep+sep << "mlml:          program to estimate hydroxymethylation levels\n";

  cout <<  "\n" << sep << "allele-specific methylation:\n";
  cout <<  sep+sep << "states:        convert read sequences in SAM format to methylation states at CpGs covered by those reads\n";
  cout <<  sep+sep << "allelic:       computes probability of allele-specific methylation at each tuple of CpGs\n";
  cout <<  sep+sep << "amrfinder:     identify regions of allele-specific methylation\n";
  cout <<  sep+sep << "amrtester:     resolve epi-alleles\n";

  cout <<  "\n" << sep << "differential methylation:\n";
  cout <<  sep+sep << "dmr:           computes DMRs based on HMRs and probability of differences at individual CpGs\n";
  cout <<  sep+sep << "diff:          compute probability that site has higher methylation in file A than B\n";
  cout <<  sep+sep << "radmeth:       computes differentially methylated CpGs\n";
  cout <<  sep+sep << "radadjust:     adjust p-values of radmeth output\n";
  cout <<  sep+sep << "radmerge:      merge significant CpGs in radmeth output\n";

  cout <<  "\n" << sep << "methylation visualization:\n";
  cout <<  sep+sep << "fastlift:      liftover methylation levels between species\n";
  cout <<  sep+sep << "liftfilter:    filter CpG regions that do not exist in resulting genome\n";

  cout <<  "\n" << sep << "general-purpose tools:\n";
  cout <<  sep+sep << "cleanhp:       fix and stat invdup/hairping reads\n";
  cout <<  sep+sep << "guessprotocol: guess whether protocol is ordinary, pbat or random\n";
  cout <<  sep+sep << "lc:            approximate line counts in a file\n";
  cout <<  sep+sep << "merge-bsrate:  merge the BS conversion rate from two sets of BS-seq reads mapped to a genome\n";
  cout <<  sep+sep << "merge:         merge multiple methcounts files\n";
  cout <<  sep+sep << "selectsites:   sites inside a set of genomic intervals\n";

  cout <<  "\n";
}

int
main(int argc, const char **argv) {
  int ret = 0;
  if (argc < 2) { print_help(); return ret; }

  if (strcmp(argv[1], "abismal") == 0) ret = abismal(argc - 1, argv + 1);
  else if (strcmp(argv[1], "abismalidx") == 0) ret = abismalidx(argc - 1, argv + 1);
  else if (strcmp(argv[1], "simreads") == 0) ret = simreads(argc - 1, argv + 1);
  else if (strcmp(argv[1], "counts") == 0) ret = main_methcounts(argc - 1, argv + 1);
  else if (strcmp(argv[1], "allelic") == 0) ret = main_allelicmeth(argc - 1, argv + 1);
  else if (strcmp(argv[1], "amrfinder") == 0) ret = main_amrfinder(argc - 1, argv + 1);
  else if (strcmp(argv[1], "amrtester") == 0) ret = main_amrtester(argc - 1, argv + 1);
  else if (strcmp(argv[1], "bsrate") == 0) ret = main_bsrate(argc - 1, argv + 1);
  else if (strcmp(argv[1], "hmr") == 0) ret = main_hmr(argc - 1, argv + 1);
  else if (strcmp(argv[1], "hmr-rep") == 0) ret = main_hmr_rep(argc - 1, argv + 1);
  else if (strcmp(argv[1], "hypermr") == 0) ret = main_hypermr(argc - 1, argv + 1);
  else if (strcmp(argv[1], "levels") == 0) ret = main_levels(argc - 1, argv + 1);
  else if (strcmp(argv[1], "entropy") == 0) ret = main_methentropy(argc - 1, argv + 1);
  else if (strcmp(argv[1], "states") == 0) ret = main_methstates(argc - 1, argv + 1);
  else if (strcmp(argv[1], "multistat") == 0) ret = main_multimethstat(argc - 1, argv + 1);
  else if (strcmp(argv[1], "pmd") == 0) ret = main_pmd(argc - 1, argv + 1);
  else if (strcmp(argv[1], "roi") == 0) ret = main_roimethstat(argc - 1, argv + 1);
  else if (strcmp(argv[1], "mlml") == 0) ret = main_mlml(argc - 1, argv + 1);
  else if (strcmp(argv[1], "dmr") == 0) ret = main_dmr(argc - 1, argv + 1);
  else if (strcmp(argv[1], "diff") == 0) ret = main_methdiff(argc - 1, argv + 1);
  else if (strcmp(argv[1], "radmeth") == 0) ret = main_radmeth(argc - 1, argv + 1);
  else if (strcmp(argv[1], "radadjust") == 0) ret = main_radmeth_adjust(argc - 1, argv + 1);
  else if (strcmp(argv[1], "radmerge") == 0) ret = main_radmeth_merge(argc - 1, argv + 1);
  else if (strcmp(argv[1], "cleanhp") == 0) ret = main_clean_hairpins(argc - 1, argv + 1);
  else if (strcmp(argv[1], "uniq") == 0) ret = main_duplicate_remover(argc - 1, argv + 1);
  else if (strcmp(argv[1], "fastlift") == 0) ret = main_fast_liftover(argc - 1, argv + 1);
  else if (strcmp(argv[1], "format") == 0) ret = main_format_reads(argc - 1, argv + 1);
  else if (strcmp(argv[1], "guessprotocol") == 0) ret = main_guessprotocol(argc - 1, argv + 1);
  else if (strcmp(argv[1], "lc") == 0) ret = main_lc_approx(argc - 1, argv + 1);
  else if (strcmp(argv[1], "liftfilter") == 0) ret = main_lift_filter(argc - 1, argv + 1);
  else if (strcmp(argv[1], "merge-bsrate") == 0) ret = main_merge_bsrate(argc - 1, argv + 1);
  else if (strcmp(argv[1], "merge") == 0) ret = main_merge_methcounts(argc - 1, argv + 1);
  else if (strcmp(argv[1], "selectsites") == 0) ret = main_selectsites(argc - 1, argv + 1);
  else if (strcmp(argv[1], "sym") == 0) ret = main_symmetric_cpgs(argc - 1, argv + 1);

  else {
    cerr << "command not found: " << argv[1] << endl;
    ret = 1;
  }

  return ret;
}
