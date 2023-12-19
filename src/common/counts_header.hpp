/* xcounts_utils: code for doing things with xcounts format and some
 * for counts format that is common to several tools.
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

#ifndef COUNTS_HEADER_HPP
#define COUNTS_HEADER_HPP

#include <string>
#include <vector>
#include <cstdint>

#include "bamxx.hpp"

void
write_counts_header_from_chrom_sizes(const std::vector<std::string> &chrom_names,
                                     const std::vector<uint64_t> &chrom_sizes,
                                     bamxx::bgzf_file &out);

void
write_counts_header_from_file(const std::string &header_file, 
                              bamxx::bgzf_file &out);

// returns -1 on failure, 0 on success
int
get_chrom_sizes_for_counts_header(const uint32_t n_threads,
                                  const std::string &filename,
                                  std::vector<std::string> &chrom_names,
                                  std::vector<uint64_t> &chrom_sizes);

void
write_counts_header_from_bam_header(const bamxx::bam_header &hdr,
                                    bamxx::bgzf_file &out);

bool
write_counts_header_line(std::string line, bamxx::bgzf_file &out);

bamxx::bgzf_file &
skip_counts_header(bamxx::bgzf_file &in);

bool
has_counts_header(const std::string &filename);

inline bool
is_counts_header_version_line(const std::string &line) {
  const auto version_line = "#DNMTOOLS";
  return line.compare(0, 9, version_line) == 0;
}

template<typename T>
inline bool
is_counts_header_line(T &line) {
  return line[0] == '#';
}

#endif
