/*
  Copyright (C) 2015-2025 Andrew D Smith

  Authors: Andrew D. Smith

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.
*/

#ifndef MSITE_HPP
#define MSITE_HPP

#include <bamxx.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iterator>  // IWYU pragma: keep
#include <sstream>
#include <string>
#include <unordered_map>

struct MSite {
  // This defaults to true, and when true MSite lines cannot be parsed if they
  // have additional columns in the file.
  static bool no_extra_fields;

  std::string chrom{};
  size_t pos{};
  char strand{};
  std::string context{};
  double meth{};
  size_t n_reads{};

  MSite() = default;
  MSite(const std::string &chrom, const size_t pos, const char strand,
        const std::string &context, const double meth, const size_t n_reads) :
    chrom{chrom}, pos{pos}, strand{strand}, context{context}, meth{meth},
    n_reads{n_reads} {}
  explicit MSite(const std::string &line);
  explicit MSite(const char *line, const int n);

  [[nodiscard]] bool
  initialize(const char *c, const char *c_end);

  [[nodiscard]] bool
  operator<(const MSite &other) const {
    const int r = chrom.compare(other.chrom);
    return (r < 0 || (r == 0 && (pos < other.pos ||
                                 (pos == other.pos && strand < other.strand))));
  }

  [[nodiscard]] size_t
  n_meth() const {
    return std::round(meth * n_reads);
  }

  [[nodiscard]] size_t
  n_unmeth() const {
    return n_reads - n_meth();
  }

  [[nodiscard]] double
  n_meth_f() const {
    return meth * n_reads;
  }

  [[nodiscard]] double
  n_unmeth_f() const {
    return n_reads - n_meth_f();
  }

  // functions below are for manipulating symmetric CpG sites
  void
  add(const MSite &other) {
    // ADS: possible that this function has specific behavior that
    // should be placed in the 'sym' source, since we might not want
    // this line below to operate generally.
    if (!is_mutated() && other.is_mutated())
      context += 'x';
    // ADS: order matters below as n_reads update invalidates n_meth()
    // function until meth has been updated
    const double total_c_reads = (meth * n_reads + other.meth * other.n_reads);
    n_reads += other.n_reads;
    meth = total_c_reads / std::max(1ul, n_reads);
  }

  // ADS: function below has redundant check for is_cpg, which is
  // expensive and might be ok to remove
  [[nodiscard]] bool
  is_mate_of(const MSite &first) {
    return (first.pos + 1 == pos && first.is_cpg() && is_cpg() &&
            first.strand == '+' && strand == '-');
  }

  // Functions below test the type of site. These are CpG, CHH, and CHG
  // divided into two kinds: CCG and CXG, the former including a CpG
  // within. Also included is a function that tests if a site has a mutation.

  [[nodiscard]] bool
  is_cpg() const {
    return std::size(context) >= 3 &&
           (context[0] == 'C' && context[1] == 'p' && context[2] == 'G');
  }

  [[nodiscard]] bool
  is_chh() const {
    return std::size(context) >= 3 &&
           (context[0] == 'C' && context[1] == 'H' && context[2] == 'H');
  }

  [[nodiscard]] bool
  is_ccg() const {
    return std::size(context) >= 3 &&
           (context[0] == 'C' && context[1] == 'C' && context[2] == 'G');
  }

  [[nodiscard]] bool
  is_cxg() const {
    return std::size(context) >= 3 &&
           (context[0] == 'C' && context[1] == 'X' && context[2] == 'G');
  }

  [[nodiscard]] bool
  is_mutated() const {
    return std::size(context) == 4 && context[3] == 'x';
  }

  void
  set_mutated() {
    if (!is_mutated())
      context += 'x';
  }

  void
  set_unmutated() {
    if (is_mutated())
      context.resize(std::size(context) - 1);
  }

  [[nodiscard]] std::string
  tostring() const;
};

template <class T>
T &
operator>>(T &in, MSite &s) {
  std::string line;
  if (getline(in, line))
    s = MSite(line);
  return in;
}

template <class T>
T &
operator<<(T &out, const MSite &s) {
  out << s.tostring();  // seems to be an issue returning this directly
  return out;
}

[[nodiscard]] size_t
distance(const MSite &a, const MSite &b);

// find the byte offset within the given file of the first site in the
// file (for the specified stream) that does not precede the specified
// position given by the (chrom, start_pos) pair.
void
find_offset_for_msite(const std::string &chrom, const size_t start_pos,
                      std::ifstream &site_in);

void
find_offset_for_msite(
  const std::unordered_map<std::string, uint32_t> &chrom_order,
  const std::string &chr, const size_t idx, std::ifstream &site_in);

// ADS: using the functions below for now for static polymorphism
// because I don't want to have insertion and extraction operators
// hanging around for MSite with bamxx::bgzf_file. The reason is that
// if we are going to use them, we need to be aware of which functions
// we are calling, and not just do it unconsciously.

inline bamxx::bgzf_file &
write_site(bamxx::bgzf_file &f, const MSite &s) {
  // ADS: too slow??
  f.write(s.tostring() + "\n");
  return f;
}

inline std::ostream &
write_site(std::ostream &out, const MSite &s) {
  return out << s << '\n';
}

inline bamxx::bgzf_file &
read_site(bamxx::bgzf_file &f, MSite &s) {
  std::string line;
  if (getline(f, line))
    s = MSite(line);
  return f;
}

[[nodiscard]] bool
is_msite_file(const std::string &file);

#endif
