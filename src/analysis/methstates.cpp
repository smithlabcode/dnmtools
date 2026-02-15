/* methstates: a program for converting read sequences in SAM format
 * files into methylation states at CpGs covered by those reads
 *
 * Copyright (C) 2011-2022 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith and Masaru Nakajima
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

#include "OptionParser.hpp"
#include "bam_record_utils.hpp"
#include "smithlab_os.hpp"

#include <bamxx.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static constexpr auto b2c = "TNGNNNCNNNNNNNNNNNNA";  // NOLINT(*-avoid-c-arrays)

template <class BidirIt, class OutputIt>
// constexpr // since C++20
OutputIt
revcomp_copy(BidirIt first, BidirIt last, OutputIt d_first) {
  for (; first != last; ++d_first)
    *d_first = b2c[*(--last) - 'A'];  // NOLINT(*-constant-array-index)
  return d_first;
}

inline static bool
is_cpg(const std::string &s, const std::uint64_t idx) {
  return s[idx] == 'C' && s[idx + 1] == 'G';
}

static void
collect_cpgs(const std::string &s, std::vector<std::uint64_t> &cpgs) {
  cpgs.clear();
  const std::uint64_t lim = std::size(s) - 1;
  for (auto i = 0u; i < lim; ++i)
    if (is_cpg(s, i))
      cpgs.push_back(i);
}

static bool
convert_meth_states_pos(const std::vector<std::uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bamxx::bam_rec &aln,
                        std::uint64_t &first_cpg_index, std::string &states) {
  states.clear();

  const std::uint64_t seq_start = get_pos(aln);
  const std::uint64_t width = rlen_from_cigar(aln);
  const std::uint64_t seq_end = seq_start + width;

  std::string seq_str;
  get_seq_str(aln, seq_str);
  apply_cigar(aln, seq_str, 'N');

  if (std::size(seq_str) != width)
    throw std::runtime_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start
  auto cpg_itr =
    std::lower_bound(std::cbegin(cpgs), std::cend(cpgs), seq_start);
  auto first_cpg_itr = std::cend(cpgs);

  if (cpg_itr == std::cend(cpgs))
    return false;

  for (; cpg_itr != std::cend(cpgs) && *cpg_itr < seq_end; ++cpg_itr) {
    const char x = seq_str[*cpg_itr - seq_start];
    states += x == 'T' ? 'T' : (x == 'C' ? 'C' : 'N');
    if (first_cpg_itr == std::cend(cpgs))
      first_cpg_itr = cpg_itr;
  }

  if (first_cpg_itr != std::cend(cpgs))
    first_cpg_index = std::distance(std::cbegin(cpgs), first_cpg_itr);

  return states.find_first_of("CT") != std::string::npos;
}

static bool
convert_meth_states_neg(const std::vector<std::uint64_t> &cpgs,
                        const bamxx::bam_header &hdr, const bamxx::bam_rec &aln,
                        std::uint64_t &first_cpg_index, std::string &states) {
  /* ADS: the "revcomp" on the read sequence is needed for the cigar
     to be applied, since the cigar is relative to the genome
     coordinates and not the read's sequence. But the read sequence
     may is assumed to have been T-rich to begin with, so it becomes
     A-rich. And the position of the C in the CpG becomes the G
     position.
   */

  states.clear();

  const std::uint64_t seq_start = get_pos(aln);
  const std::uint64_t width = rlen_from_cigar(aln);
  const std::uint64_t seq_end = seq_start + width;

  std::string orig_seq;
  get_seq_str(aln, orig_seq);

  std::string seq_str;
  seq_str.resize(std::size(orig_seq));
  revcomp_copy(std::cbegin(orig_seq), std::cend(orig_seq), std::begin(seq_str));
  apply_cigar(aln, seq_str, 'N');

  if (std::size(seq_str) != width)
    throw std::runtime_error("bad sam record format: " + to_string(hdr, aln));

  // get the first cpg site equal to or large than seq_start - 1
  // the -1 is because we look for G in the read corresponding to a
  // CpG in chromosome, which are indexed in cpgs based on the position of C
  auto cpg_itr = std::lower_bound(std::cbegin(cpgs), std::cend(cpgs),
                                  seq_start > 0 ? seq_start - 1 : 0);
  auto first_cpg_itr = std::cend(cpgs);

  if (cpg_itr == std::cend(cpgs))
    return false;

  for (; cpg_itr != std::cend(cpgs) && *cpg_itr < seq_end - 1; cpg_itr++) {
    const char x = seq_str[*cpg_itr - seq_start + 1];
    states += (x == 'G') ? 'C' : ((x == 'A') ? 'T' : 'N');
    if (first_cpg_itr == std::cend(cpgs))
      first_cpg_itr = cpg_itr;
  }

  if (first_cpg_itr != std::cend(cpgs))
    first_cpg_index = std::distance(std::cbegin(cpgs), first_cpg_itr);

  return states.find_first_of("CT") != std::string::npos;
}

static void
get_chrom(const std::string &chrom_name,
          const std::vector<std::string> &all_chroms,
          const std::unordered_map<std::string, std::uint64_t> &chrom_lookup,
          std::string &chrom) {
  auto the_chrom = chrom_lookup.find(chrom_name);
  if (the_chrom == std::cend(chrom_lookup))
    throw std::runtime_error("could not find chrom: " + chrom_name);

  chrom = all_chroms[the_chrom->second];
  if (chrom.empty())
    throw std::runtime_error("problem with chrom: " + chrom_name);
}

int
main_methstates(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  static constexpr std::int64_t output_buffer_size{4096};
  std::array<char, output_buffer_size> buf{};
  static constexpr auto output_format = "%s\t%lu\t%s\n";

  try {
    // clang-format off
    const auto description =
R"(Convert mapped reads in SAM format into a format that indicates binary
sequences of methylation states in each read, indexed by the identity of the
CpG they cover, along with the chromosome. Only reads that cover a CpG site
are included in the output. All output is relative to the positive reference
strand. This format is used as input to other tools, and is not intended to be
human-interpretable. All chromosome sequences are loaded at once.)";
    // clang-format on

    bool verbose{};
    bool compress_output{};

    std::string chrom_file;
    std::string outfile("-");
    int n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description, "<sam-file>");
    opt_parse.add_opt("output", 'o', "output file name", false, outfile);
    opt_parse.add_opt("chrom", 'c', "fasta format reference genome file", true,
                      chrom_file);
    opt_parse.add_opt("threads", 't', "threads to use for reading input", false,
                      n_threads);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);

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
    const std::string mapped_reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (n_threads < 1)
      throw std::runtime_error("thread count must be at least 1");

    /* first load in all the chromosome sequences and names, and make a map
       from chromosome name to the location of the chromosome itself */
    std::vector<std::string> all_chroms, chrom_names;
    read_fasta_file_short_names(chrom_file, chrom_names, all_chroms);
    for (auto &i : all_chroms)
      std::transform(std::cbegin(i), std::cend(i), std::begin(i),
                     [](const char c) { return std::toupper(c); });

    std::unordered_map<std::string, std::uint64_t> chrom_lookup;
    for (std::uint64_t i = 0; i < std::size(chrom_names); ++i)
      chrom_lookup[chrom_names[i]] = i;

    if (verbose)
      std::cerr << "n_chroms: " << std::size(all_chroms) << '\n';

    bamxx::bam_tpool tp(n_threads);  // declared first; destroyed last

    bamxx::bam_in in(mapped_reads_file);
    if (!in)
      throw std::runtime_error("cannot open input file " + mapped_reads_file);
    bamxx::bam_header hdr(in);
    if (!hdr)
      throw std::runtime_error("cannot read heade" + mapped_reads_file);

    // open the output file
    const std::string output_mode = compress_output ? "w" : "wu";
    bamxx::bgzf_file out(outfile, output_mode);
    if (!out)
      throw std::runtime_error("error opening output file: " + outfile);

    /* set the threads for the input file decompression */
    if (n_threads > 1) {
      tp.set_io(in);
      tp.set_io(out);
    }

    std::vector<std::uint64_t> cpgs;
    std::unordered_set<std::int32_t> chroms_seen;
    std::int32_t chrom_idx{-1};

    // iterate over records/reads in the SAM file, sequentially processing
    // each before considering the next
    bamxx::bam_rec aln;
    while (in.read(hdr, aln)) {
      if (get_tid(aln) != chrom_idx) {  // get correct chrom if it has changed
        chrom_idx = get_tid(aln);
        // make sure all reads from same chrom are contiguous in the file
        if (chroms_seen.find(chrom_idx) != std::cend(chroms_seen))
          throw std::runtime_error("wrong chrom order (check SAM/BAM sorting)");
        const auto chrom_name = sam_hdr_tid2name(hdr, chrom_idx);
        if (verbose)
          std::cerr << "processing " << chrom_name << '\n';

        std::string chrom;
        get_chrom(chrom_name, all_chroms, chrom_lookup, chrom);
        collect_cpgs(chrom, cpgs);
      }

      std::uint64_t first_cpg_index = std::numeric_limits<std::uint64_t>::max();
      std::string seq;

      const bool has_cpgs =
        bam_is_rev(aln)
          ? convert_meth_states_neg(cpgs, hdr, aln, first_cpg_index, seq)
          : convert_meth_states_pos(cpgs, hdr, aln, first_cpg_index, seq);

      if (has_cpgs) {
        const auto n =
          std::snprintf(std::data(buf), output_buffer_size, output_format,
                        sam_hdr_tid2name_ptr(hdr, chrom_idx), first_cpg_index,
                        std::data(seq));
        if (n < 0 || n >= output_buffer_size || !out.write(std::data(buf), n))
          throw std::runtime_error("failure writing output");
      }
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
