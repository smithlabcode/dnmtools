/* xcounts: reformat counts so they only give the m and u counts in a dynamic
 * step wig format
 *
 * Copyright (C) 2023 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include "MSite.hpp"
#include "OptionParser.hpp"
#include "counts_header.hpp"
#include "dnmt_error.hpp"

#include <bamxx.hpp>

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <string>
#include <system_error>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// NOLINTBEGIN(*-avoid-c-arrays,*-avoid-magic-numbers,*-avoid-non-const-global-variables,*-narrowing-conversions,*-constant-array-index,*-pointer-arithmetic)

enum class xcounts_err : std::uint8_t {
  // clang-format off
  ok                            = 0,
  open_failure                  = 1,
  header_failure                = 2,
  chromosome_not_found          = 3,
  inconsistent_chromosome_order = 4,
  incorrect_chromosome_size     = 5,
  failed_to_write_file          = 6,
  failed_to_parse_site          = 7,
  // clang-format on
};

// clang-format off
struct xcounts_err_cat : std::error_category {
  auto name() const noexcept -> const char * override {
    return "xcounts error";
  }
  auto message(const int condition) const -> std::string override {
    using std::string_literals::operator""s;
    switch (condition) {
    case 0: return "ok"s;
    case 1: return "failed to open methylome file"s;
    case 2: return "failed to parse xcounts header"s;
    case 3: return "failed to find chromosome in xcounts header"s;
    case 4: return "inconsistent chromosome order"s;
    case 5: return "incorrect chromosome size"s;
    case 6: return "failed to write to file"s;
    case 7: return "failed to parse site"s;
    }
    std::abort();  // unreacheable
  }
};
// clang-format on

template <>
struct std::is_error_code_enum<xcounts_err> : public std::true_type {};

std::error_code
make_error_code(xcounts_err e) {
  static auto category = xcounts_err_cat{};
  return std::error_code(static_cast<std::underlying_type_t<xcounts_err>>(e),
                         category);
}

template <typename T>
static inline std::uint32_t
fill_output_buffer(const std::uint32_t offset, const MSite &s, T &buf) {
  auto buf_end = buf.data() + buf.size();
  auto res = std::to_chars(buf.data(), buf_end, s.pos - offset);
  *res.ptr++ = '\t';
  res = std::to_chars(res.ptr, buf_end, s.n_meth());
  *res.ptr++ = '\t';
  res = std::to_chars(res.ptr, buf_end, s.n_unmeth());
  *res.ptr++ = '\n';
  return std::distance(buf.data(), res.ptr);
}

int
main_xcounts(int argc, char *argv[]) {  // NOLINT(*-avoid-c-arrays)
  try {
    // ADS: It might happen that a "chromosome" has no CpG sites (like
    // Scaffold113377 in strPur2). Therefore, we can't assume each chrom will
    // appear in the input when attempting to check consistency.
    bool verbose = false;
    bool gzip_output = false;
    bool keep_all_sites = false;
    bool require_coverage = false;
    std::size_t n_threads = 1;
    std::string genome_file;
    std::string header_file;

    std::string outfile{"-"};
    const std::string description =
      "compress counts files by removing context information";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           description,
                           "<counts-file> (\"-\" for standard input)", 1);
    opt_parse.add_opt("output", 'o', "output file (default is standard out)",
                      false, outfile);
    opt_parse.add_opt("chroms", 'c', "make header from this reference", false,
                      genome_file);
    opt_parse.add_opt("reads", 'r',
                      "output only sites with reads (on by default)", false,
                      require_coverage);
    opt_parse.add_opt("keep", 'k', "keep sites that are not covered", false,
                      keep_all_sites);
    opt_parse.add_opt("header", 'h', "use this file to generate header", false,
                      header_file);
    opt_parse.add_opt("threads", 't', "threads for compression (use few)",
                      false, n_threads);
    opt_parse.add_opt("zip", 'z',
                      "gzip compress output (automatic if input is gzip)",
                      false, gzip_output);
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
    if (leftover_args.size() != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string filename(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    std::vector<std::string> chrom_names;
    std::vector<std::uint64_t> chrom_sizes;
    if (!genome_file.empty()) {
      const int ret = get_chrom_sizes_for_counts_header(
        n_threads, genome_file, chrom_names, chrom_sizes);
      if (ret)
        throw dnmt_error{"failed to get chrom sizes from: " + genome_file};
    }

    bamxx::bam_tpool tpool(n_threads);
    bamxx::bgzf_file in(filename, "r");
    if (!in)
      throw dnmt_error{"could not open file: " + filename};

    const auto outfile_mode = (gzip_output || in.is_compressed()) ? "w" : "wu";

    bamxx::bgzf_file out(outfile, outfile_mode);
    if (!out)
      throw dnmt_error{"error opening output file: " + outfile};

    if (n_threads > 1) {
      if (in.is_bgzf())
        tpool.set_io(in);
      tpool.set_io(out);
    }

    std::unordered_map<std::string, std::uint32_t> chrom_order;
    std::unordered_set<std::uint32_t> chroms_seen;
    if (!header_file.empty())
      chrom_order = write_counts_header_from_file(header_file, out);
    else if (!genome_file.empty())
      chrom_order =
        write_counts_header_from_chrom_sizes(chrom_names, chrom_sizes, out);

    // use the kstring_t type to more directly use the BGZF file
    kstring_t line{0, 0, nullptr};
    const int ret = ks_resize(&line, 1024);
    if (ret)
      throw dnmt_error("failed to acquire buffer");

    std::vector<char> buf(128);

    std::uint32_t offset = 0;
    std::string prev_chrom;
    bool found_header = (!genome_file.empty() || !header_file.empty());

    std::error_code ec{};

    std::uint32_t chrom_counter = 0;

    MSite site;
    while (ec == std::errc{} && bamxx::getline(in, line)) {
      if (is_counts_header_line(line.s)) {
        if (!genome_file.empty() || !header_file.empty())
          continue;
        found_header = true;
        const std::string header_line{line.s};
        write_counts_header_line(header_line, out);
        continue;
      }

      if (!site.initialize(line.s, line.s + line.l)) {
        ec = xcounts_err::failed_to_parse_site;
        break;
      }
      if (!found_header) {
        ec = xcounts_err::header_failure;
        break;
      }

      if (site.chrom != prev_chrom) {
        if (verbose)
          std::cerr << "processing: " << site.chrom << '\n';

        if (!chrom_order.empty()) {
          const auto expected_chrom_order = chrom_order.find(site.chrom);
          if (expected_chrom_order == std::cend(chrom_order)) {
            ec = xcounts_err::chromosome_not_found;
            break;
          }
          if (expected_chrom_order->second < chrom_counter) {
            ec = xcounts_err::inconsistent_chromosome_order;
            break;
          }
          if (chroms_seen.count(expected_chrom_order->second) > 0) {
            ec = xcounts_err::inconsistent_chromosome_order;
            break;
          }
          chroms_seen.emplace(expected_chrom_order->second);
        }

        ++chrom_counter;

        prev_chrom = site.chrom;
        offset = 0;

        site.chrom += '\n';
        const std::int64_t sz = size(site.chrom);
        if (bgzf_write(out.f, site.chrom.data(), sz) != sz)
          ec = xcounts_err::failed_to_write_file;
      }
      if (keep_all_sites || site.n_reads > 0) {
        const std::int64_t sz = fill_output_buffer(offset, site, buf);
        if (bgzf_write(out.f, buf.data(), sz) != sz)
          ec = xcounts_err::failed_to_write_file;
        offset = site.pos;
      }
    }
    ks_free(&line);

    if (ec) {
      std::cerr << "failed converting " << filename << " to " << outfile << '\n'
                << ec.message() << '\n';
      return EXIT_FAILURE;
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

// NOLINTEND(*-avoid-c-arrays,*-avoid-magic-numbers,*-avoid-non-const-global-variables,*-narrowing-conversions,*-constant-array-index,*-pointer-arithmetic)
