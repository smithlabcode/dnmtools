/* Copyright (C) 2020-2023 Masaru Nakajima and Andrew D. Smith
 *
 * Authors: Masaru Nakajima and Andrew D. Smith
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

/* ADS: need to control all the macros from HTSlib pollution. For
   functions maybe:

   $ gcc -dM -E sam.h | grep "define [a-z]" | awk '{print $2}' |\
       grep "[(]" | awk -v FS="(" '{print "#undef",$1}'

   This gives about 65 symbols that need to be deleted. For the others
   I don't know what to do because some of them have "#define _" which
   means they should be system symbols.
*/

#ifndef BAM_RECORD_HPP
#define BAM_RECORD_HPP

#include <stdexcept>
#include <string>

#include <htslib/sam.h>  // from HTSlib
#include <htslib/thread_pool.h>  // from HTSlib

struct bam_rec {
  bam_rec() : b{bam_init1()} {}
  bam_rec(const bam_rec &other) : b{bam_copy1(bam_init1(), other.b)} {}
  bam_rec &operator=(bam_rec rhs) { std::swap(b, rhs.b); return *this; }
  ~bam_rec() { if (b != nullptr) bam_destroy1(b); }
  bam1_t *b;
};

struct bam_infile {
  bam_infile(const std::string &fn) : f{hts_open(fn.c_str(), "r")} {}
  ~bam_infile() { if (f != nullptr) hts_close(f); }
  operator bool() const { return f != nullptr; }
  template <typename T> bool read(T &h, bam_rec &b) {
    const int x = sam_read1(f, h.h, b.b);  // -1 on EOF; args non-const
    if (x < -1) throw std::runtime_error("failed reading bam record");
    return x >= 0;
  }
  bool is_mapped_reads_file() const {
    const htsFormat *fmt = hts_get_format(f);
    return fmt->category == sequence_data &&
      (fmt->format == bam || fmt->format == sam);
  }
  samFile *f{};
};

struct bam_header {
  bam_header() = default;
  bam_header(const bam_header &rhs) : h{bam_hdr_dup(rhs.h)} {}
  bam_header(bam_infile &in) : h{sam_hdr_read(in.f)} {}
  ~bam_header() { if (h != nullptr) bam_hdr_destroy(h); }
  operator bool() const { return h != nullptr; }
  bool add_pg_line(const std::string cmd, const std::string id,
                   const std::string vn) {
    return sam_hdr_add_line(h, "PG", "ID", id.c_str(),  "VN",
                            vn.c_str(), "CL", cmd.c_str(), nullptr) == 0;
  }
  std::string tostring() const { return sam_hdr_str(h); }
  sam_hdr_t *h{};
};

struct bam_outfile {
  bam_outfile(const std::string &fn, const bool fmt = false) :
    f{hts_open(fn.c_str(), fmt ? "bw" : "w")} {}
  ~bam_outfile() { if (f != nullptr) hts_close(f); }
  operator bool() const { return f != nullptr; }
  bool write(const bam_header &h, const bam_rec &b) {
    return sam_write1(f, h.h, b.b) >= 0;
  }
  bool write(const bam_header &h) { return sam_hdr_write(f, h.h) == 0; }
  htsFile *f{};
};

struct bam_tpool{
  bam_tpool(const size_t n_threads) : tpool{hts_tpool_init(n_threads), 0} {}
  ~bam_tpool() { hts_tpool_destroy(tpool.pool); }
  template<class T> void set_io(const T &bam_file) {
    const int ret = hts_set_thread_pool(bam_file.f, &tpool);
    if (ret < 0) throw std::runtime_error("failed to set thread pool");
  }
  htsThreadPool tpool{};
};

#endif
