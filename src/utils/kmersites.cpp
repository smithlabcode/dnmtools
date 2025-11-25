/* kmersites: a program to generate a wiggle format file (using the UCSC
 * Genome Browser wiggle format) to indicate the location of sites matching a
 * specific k-mer
 *
 * Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
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

#include <bamxx.hpp>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

static inline auto
process_chrom_wig(const std::string &kmer, const int offset,
                  const std::string &name, const std::string &chrom,
                  bamxx::bgzf_file &out) -> void {
  static const auto variable_step_chrom_header = "variableStep chrom=";

  out.write(variable_step_chrom_header + name + "\n");

  const auto kmer_size = size(kmer);
  const auto chrom_size = size(chrom);
  if (kmer_size > chrom_size)
    throw std::runtime_error("kmer size " + std::to_string(kmer_size) +
                             " larger than chrom size " +
                             std::to_string(chrom_size));

  const auto beg_kmer = std::cbegin(kmer);
  const auto end_kmer = std::cend(kmer);

  const auto end_chrom = std::cend(chrom);
  auto chrom_itr = std::cbegin(chrom);
  auto chrom_itr_k = chrom_itr + kmer_size;  // NOLINT(*-narrowing-conversions)

  auto pos = 0;
  while (chrom_itr_k != end_chrom) {
    if (std::equal(beg_kmer, end_kmer, chrom_itr++, chrom_itr_k++))
      out.write(std::to_string(pos + offset) + "\t1\n");
    ++pos;
  }
}

[[nodiscard]] static auto
read_fasta_file(const std::string &filename)
  -> std::tuple<std::vector<std::string>, std::vector<std::string>> {
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("cannot open input file " + filename);

  std::vector<std::string> names;
  std::vector<std::string> sequences;

  std::string line;
  while (std::getline(in, line)) {
    if (line[0] == '>') {
      const auto first_space = line.find_first_of(" \t", 1);
      if (first_space == std::string::npos)
        names.push_back(line.substr(1));
      else
        names.push_back(line.substr(1, first_space - 1));
      sequences.emplace_back();
    }
    else
      sequences.back() += line;
  }
  return {names, sequences};
}

static inline auto
process_chrom_with_named_lines(const std::string &kmer, const int offset,
                               const std::string &name,
                               const std::string &chrom,
                               bamxx::bgzf_file &out) {
  const auto kmer_size = size(kmer);
  const auto chrom_size = size(chrom);
  if (kmer_size > chrom_size)
    throw std::runtime_error("kmer size " + std::to_string(kmer_size) +
                             " larger than chrom size " +
                             std::to_string(chrom_size));

  const auto beg_kmer = std::cbegin(kmer);
  const auto end_kmer = std::cend(kmer);

  const auto end_chrom = std::cend(chrom);
  auto chrom_itr = std::cbegin(chrom);
  auto chrom_itr_k = chrom_itr + kmer_size;  // NOLINT(*-narrowing-conversions)

  auto pos = 0;
  while (chrom_itr_k != end_chrom) {
    if (std::equal(beg_kmer, end_kmer, chrom_itr++, chrom_itr_k++))
      out.write(name + "\t" + std::to_string(pos + offset) + "\t1\n");
    ++pos;
  }
}

[[nodiscard]] static inline auto
bad_dna_kmer(const std::string &kmer) -> bool {
  const auto x =
    std::find_if(std::cbegin(kmer), std::cend(kmer), [](const auto c) {
      return c != 'A' && c != 'C' && c != 'G' && c != 'T';
    });
  return x != std::cend(kmer);
}

auto
kmersites(const int argc, char *argv[]) -> int {  // NOLINT(*-avoid-c-arrays)
  try {
    bool verbose{false};
    bool show_progress{false};
    bool compress_output{false};
    bool name_each_line{false};

    std::string kmer = "CG";
    std::string outfile;
    int offset = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("dnmtools kmersites", "get sites matching kmer",
                           "<fasta-file>");
    opt_parse.add_opt("output", 'o', "output file name (default: stdout)",
                      false, outfile);
    opt_parse.add_opt("offset", 'O', "offset within kmer to report", false,
                      offset);
    opt_parse.add_opt("kmer", 'k', "kmer to report", false, kmer);
    opt_parse.add_opt("zip", 'z', "output gzip format", false, compress_output);
    opt_parse.add_opt("name-each-line", '\0', "name each line with chrom",
                      false, name_each_line);
    opt_parse.add_opt("progress", '\0', "show progress", false, show_progress);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.about_requested() || opt_parse.help_requested() ||
        leftover_args.empty()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string chroms_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    if (bad_dna_kmer(kmer)) {
      std::cerr << "invalid DNA kmer: " << kmer << "\n";
      return EXIT_FAILURE;
    }

    if (offset < 0)
      throw std::runtime_error("offset must be non-negative (specified=" +
                               std::to_string(offset) + ")");

    std::ostringstream cmd;
    std::copy(argv, argv + argc, std::ostream_iterator<const char *>(cmd, " "));

    // file types from HTSlib use "-" for the filename to go to stdout
    if (outfile.empty())
      outfile = "-";

    if (verbose)
      std::cerr << "[input fastq file: " << chroms_file << "]\n"
                << "[output file: " << outfile << "]\n"
                << "[output format: " << (compress_output ? "bgzf" : "text")
                << "]\n"
                << "[k-mer sequence to report: " << kmer << "]\n"
                << "[command line: " << cmd.str() << "]\n";

    auto [names, chroms] = read_fasta_file(chroms_file);
    for (auto &chrom : chroms)  // cppcheck-suppress constVariableReference
      std::transform(std::cbegin(chrom), std::cend(chrom), std::begin(chrom),
                     [](const char c) { return std::toupper(c); });

    // open the output file
    bamxx::bgzf_file out(outfile, compress_output ? "w" : "wu");
    if (!out)
      throw std::runtime_error("error opening output file: " + outfile);

    auto chrom_itr = std::cbegin(chroms);
    for (const auto &name : names) {
      if (show_progress)
        std::cerr << "processing: " << name << '\n';
      if (name_each_line)
        process_chrom_with_named_lines(kmer, offset, name, *chrom_itr++, out);
      else
        process_chrom_wig(kmer, offset, name, *chrom_itr++, out);
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
