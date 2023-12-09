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

#include "counts_header.hpp"

#include <string>
#include <vector>
#include <cassert>
#include <cstdint>
#include <sstream>

#include "bam_record_utils.hpp"

// generated by autotools
#include <config.h>

#include "bamxx.hpp"

using std::vector;
using std::string;
using std::to_string;

using bamxx::bgzf_file;

void
write_counts_header_from_chom_sizes(const vector<string> &chrom_names,
                                    const vector<uint64_t> &chrom_sizes,
                                    bgzf_file &out) {
  const auto version = "#DNMTOOLS " + string(VERSION) + "\n";
  out.write(version.c_str());
  for (auto i = 0u; i < size(chrom_sizes); ++i) {
    const string tmp =
      "#" + chrom_names[i] + " " + to_string(chrom_sizes[i]) + "\n";
    out.write(tmp.c_str());
  }
  out.write("#\n");
}


inline bgzf_file &
getline(bgzf_file &file, kstring_t &line) {
  if (file.f == nullptr) return file;
  const int x = bgzf_getline(file.f, '\n', &line);
  if (x == -1) {
    file.destroy();
    free(line.s);
    line = {0, 0, nullptr};
  }
  if (x < -1) {
    // ADS: this is an error condition and should be handled
    // differently from the EOF above.
    file.destroy();
    free(line.s);
    line = {0, 0, nullptr};
  }
  return file;
}


bamxx::bgzf_file &
skip_counts_header(bamxx::bgzf_file &in) {

  // use the kstring_t type to more directly use the BGZF file
  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, 1024);
  if (ret) return in;

  while (getline(in, line) && line.s[0] == '#') {
    if (line.s[0] == '#' && line.l == 1)
      return in;
  }
  // otherwise we have missed the final line of the header
  assert(line.s[0] != '#');
  return in;
}


int
get_chrom_sizes_for_counts_header(const uint32_t n_threads,
                                  const string &filename,
                                  vector<string> &chrom_names,
                                  vector<uint64_t> &chrom_sizes) {

  bamxx::bam_tpool tpool(n_threads);

  bgzf_file in(filename, "r");
  if (!in) return -1;
  if (n_threads > 1 && in.is_bgzf())
    tpool.set_io(in);

  chrom_names.clear();
  chrom_sizes.clear();

  // use the kstring_t type to more directly use the BGZF file
  kstring_t line{0, 0, nullptr};
  const int ret = ks_resize(&line, 1024);
  if (ret) return -1;

  uint64_t chrom_size = 0;
  while (getline(in, line)) {
    if (line.s[0] == '>') {
      if (!chrom_names.empty()) chrom_sizes.push_back(chrom_size);
      chrom_names.emplace_back(line.s + 1);
      chrom_size = 0;
    }
    else chrom_size += line.l;
  }
  if (!chrom_names.empty()) chrom_sizes.push_back(chrom_size);

  ks_free(&line);

  assert(size(chrom_names) == size(chrom_sizes));

  return 0;
}


void
write_counts_header_from_bam_header(const bamxx::bam_header &hdr,
                                    bgzf_file &out) {
  const auto version = "#DNMTOOLS " + string(VERSION) + "\n";
  out.write(version.c_str());
  for (auto i = 0; i < hdr.h->n_targets; ++i) {
    const size_t tid_size = sam_hdr_tid2len(hdr, i);
    const string tid_name = sam_hdr_tid2name(hdr, i);
    std::ostringstream oss;
    oss << "#" << tid_name << ' ' << tid_size << '\n';
    out.write(oss.str());
  }
  out.write("#\n");
}


bool
write_counts_header_line(string line, bgzf_file &out) {
  line += '\n';
  const int64_t ret = bgzf_write(out.f, line.data(), size(line));
  return ret == static_cast<int64_t>(size(line));
}

bool
has_counts_header(const string &filename) {
  bgzf_file in(filename, "r");
  if (!in) return false;
  string line;
  if (!getline(in, line)) return false;
  return line[0] == '#';
}
