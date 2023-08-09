/*
  Copyright (C) 2015-2022 University of Southern California
  Authors: Andrew D. Smith

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include "MSite.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <regex>

#include "smithlab_utils.hpp"

using std::string;
using std::runtime_error;
using std::regex_match;

MSite::MSite(const string &line) {
  /* GS: this is faster but seems to be genenerating issues when
   * compiled with clang
  std::istringstream iss;
  iss.rdbuf()->pubsetbuf(const_cast<char*>(line.c_str()), line.length());
  */
  std::istringstream iss(line);
  string strand_tmp;
  if (!(iss >> chrom >> pos >> strand_tmp >> context >> meth >> n_reads))
    throw std::runtime_error("bad line: \"" + line + "\"");
  strand = strand_tmp[0];
  if (strand != '-' && strand != '+')
    throw std::runtime_error("bad line: \"" + line + "\"");
}


string
MSite::tostring() const {
  std::ostringstream oss;
  oss << chrom << '\t'
      << pos << '\t'
      << strand << '\t'
      << context << '\t'
      << meth << '\t'
      << n_reads;
  return oss.str();
}


size_t
distance(const MSite &a, const MSite &b) {
  return a.chrom == b.chrom ?
    std::max(a.pos, b.pos) - std::min(a.pos, b.pos) :
    std::numeric_limits<size_t>::max();
}


using std::ifstream;
using std::ios_base;

static void
move_to_start_of_line(ifstream &in) {
  char next;
  // move backwards by: one step forward, two steps back
  while (in.good() && in.get(next) && next != '\n') {
    in.unget();
    in.unget();
  }
  // ADS: need to reconsider this below
  if (in.bad())
    // hope this only happens when hitting the start of the file
    in.clear();
}


void
find_offset_for_msite(const std::string &chr,
                      const size_t idx,
                      std::ifstream &site_in) {

  site_in.seekg(0, ios_base::beg);
  const size_t begin_pos = site_in.tellg();
  site_in.seekg(0, ios_base::end);
  const size_t end_pos = site_in.tellg();

  if (end_pos - begin_pos < 2)
    throw runtime_error("empty counts file");

  size_t step_size = (end_pos - begin_pos)/2;

  site_in.seekg(0, ios_base::beg);
  string low_chr;
  size_t low_idx = 0;
  if (!(site_in >> low_chr >> low_idx))
    throw runtime_error("failed search in counts file");

  // MAGIC: need the -2 here to get past the EOF and possibly a '\n'
  site_in.seekg(-2, ios_base::end);
  move_to_start_of_line(site_in);
  string high_chr;
  size_t high_idx;
  if (!(site_in >> high_chr >> high_idx))
    throw runtime_error("failed search in counts file");

  size_t pos = step_size;
  site_in.seekg(pos, ios_base::beg);
  move_to_start_of_line(site_in);
  size_t prev_pos = 0; // keep track of previous position in file

  // binary search inside sorted file on disk
  // iterate until step size is 0 or positions are identical
  while (step_size > 0 && prev_pos != pos) {

    // track (mid) position in file to make sure it keeps moving
    prev_pos = pos;

    string mid_chr; // chromosome name at mid point
    size_t mid_idx = 0; // position at mid point
    if (!(site_in >> mid_chr >> mid_idx))
      throw runtime_error("failed navigating inside file");

    // this check will never give a false indication of unsorted, but
    // might catch an unsorted file
    if (mid_chr < low_chr || mid_chr > high_chr)
      throw runtime_error("chromosomes unsorted inside file: "
                          "low=" + low_chr + ",mid=" + mid_chr +
                          ",high=" + high_chr);

    step_size /= 2; // cut the range in half

    if (chr < mid_chr || (chr == mid_chr && idx <= mid_idx)) {
      // move to the left
      std::swap(mid_chr, high_chr);
      std::swap(mid_idx, high_idx);
      pos -= step_size;
    }
    else {
      // move to the left
      std::swap(mid_chr, low_chr);
      std::swap(mid_idx, low_idx);
      pos += step_size;
    }
    site_in.seekg(pos, ios_base::beg);
    move_to_start_of_line(site_in);
  }
}



bool
is_msite_file(const string &file) {
  ifstream in(file);
  if (!in)
    throw runtime_error("cannot open file: " + file);

  string line;
  if(!getline(in, line)) return false;

  std::istringstream iss(line);

  string chrom;
  if (!(iss >> chrom)) return false;
  
  long int pos = 0;
  if (!(iss >> pos)) return false;

  string strand;
  if (!(iss >> strand) || 
      (strand.size() != 1) || 
      ((strand != "+") && (strand != "-")) ) 
    return false;

  string context;
  std::regex pattern("^C[pHWX][GH]$");
  if (!(iss >> context) || !regex_match(context, pattern)) return false;

  double level = 0.0;
  if (!(iss >> level) || level < 0 || level > 1) return false;

  long int n_reads = 0;
  if (!(iss >> n_reads)) return false;

  string temp;
  if (iss >> temp) return false;
  else return true;

}

