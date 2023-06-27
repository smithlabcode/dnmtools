/* format: a command within dnmtools to ensure SAM and BAM format
 * reads are conforming to expectations of dnmtools software
 *
 * Copyright (C) 2020-2023 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew Smith and Guilherme Sena
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

/* The output of this program should include mapped reads that are
 * T-rich, which might require reverse-complementing sequences, and
 * switching their strand, along with any tag that indicates T-rich
 * vs. A-rich.
 *
 * Focusing only on single end reads, for each supported read mapping
 * tool, we require a means of determining whether or not the read is
 * A-rich and then changing that format to indicate T-rich.
 */

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <stdexcept>
#include <sstream>
#include <cstring>
#include <cstdint> // for [u]int[0-9]+_t

#include <config.h>

// from HTSlib
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

// from smithlab
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;


struct fr_expt: public std::exception {
  int64_t err;        // error possibly from HTSlib
  int the_errno;      // ERRNO at time of construction
  string msg;         // the message
  string the_what;    // to report
  fr_expt(const int64_t _err, const string &_msg) :
    err{_err}, the_errno{errno}, msg{_msg} {
    std::ostringstream oss;
    oss << "[error: " << err << "][" << "ERRNO: " << the_errno << "]"
        << "[" << strerror(the_errno) << "][" << msg << "]";
    the_what = oss.str();
  }
  fr_expt(const string &_msg) : fr_expt(0, _msg) {}
  const char*
  what() const noexcept override {return the_what.c_str();}
};


static inline bool
eats_ref(const uint32_t c) {return bam_cigar_type(bam_cigar_op(c)) & 2;}

static inline bool
eats_query(const uint32_t c) {return bam_cigar_type(bam_cigar_op(c)) & 1;}

static inline size_t
bam_get_n_cigar(const bam1_t *b) {return b->core.n_cigar;}


static inline uint32_t
to_insertion(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CINS;
}


static void
fix_internal_softclip(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 3) return;
  // find first non-softclip
  auto c_beg = cigar;
  auto c_end = cigar + n_cigar;

  while (!eats_ref(*c_beg) && ++c_beg != c_end);
  if (c_beg == c_end) throw fr_expt("cigar eats no ref");

  while (!eats_ref(*(c_end-1)) && --c_end != c_beg);
  if (c_beg == c_end) throw fr_expt("cigar eats no ref");

  for (auto c_itr = c_beg; c_itr != c_end; ++c_itr)
    if (bam_cigar_op(*c_itr) == BAM_CSOFT_CLIP)
      *c_itr = to_insertion(*c_itr);
}


static inline uint32_t
to_softclip(const uint32_t x) {
  return (x & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
}


static void
fix_external_insertion(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return;

  auto c_itr = cigar;
  const auto c_end = c_itr + n_cigar;

  for (; !eats_ref(*c_itr) && c_itr != c_end; ++c_itr)
    *c_itr = to_softclip(*c_itr);

  if (c_itr == c_end) throw fr_expt("cigar eats no ref");

  c_itr = cigar + n_cigar - 1;
  for (; !eats_ref(*c_itr) && c_itr != cigar; --c_itr)
    *c_itr = to_softclip(*c_itr);
}


static size_t
merge_cigar_ops(const size_t n_cigar, uint32_t *cigar) {
  if (n_cigar < 2) return n_cigar;
  auto c_itr1 = cigar;
  auto c_end = c_itr1 + n_cigar;
  auto c_itr2 = c_itr1 + 1;
  auto op1 = bam_cigar_op(*c_itr1);
  while (c_itr2 != c_end) {
    auto op2 = bam_cigar_op(*c_itr2);
    if (op1 == op2) {
      *c_itr1 = bam_cigar_gen(bam_cigar_oplen(*c_itr1) +
                              bam_cigar_oplen(*c_itr2), op1);
    }
    else {
      *(++c_itr1) = *c_itr2;
      op1 = op2;
    }
    ++c_itr2;
  }
  // another increment to move past final "active" element for c_itr1
  ++c_itr1;
  return std::distance(cigar, c_itr1);
}


static size_t
correct_cigar(bam1_t *b) {
  /* This function will change external insertions into soft clip
     operations. Not sure why those would be present. It will also
     change internal soft-clip operations into insertions. This could
     be needed if soft-clipped ends of reads were moved to the middle
     of a merged fragment. Finally, it will collapse adjacent
     identical operations. None of this impacts the seq/qual/aux which
     get moved as a block */

  uint32_t *cigar = bam_get_cigar(b);
  size_t n_cigar = b->core.n_cigar;
  fix_external_insertion(n_cigar, cigar);
  fix_internal_softclip(n_cigar, cigar);

  // merge identical adjacent cigar ops and get new number of ops
  n_cigar = merge_cigar_ops(n_cigar, cigar);
  // difference in bytes to shift the internal data
  const size_t delta = (b->core.n_cigar - n_cigar)*sizeof(uint32_t);
  if (delta > 0) { // if there is a difference; do the shift
    auto data_end = bam_get_aux(b) + bam_get_l_aux(b);
    std::copy(bam_get_seq(b), data_end, bam_get_seq(b) - delta);
    b->core.n_cigar = n_cigar; // and update number of cigar ops
  }
  return delta;
}


static inline size_t
get_rlen(const bam1_t *b) { // less tedious
  return bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
}


static inline size_t
get_qlen(const bam1_t *b) { // less tedious
  return b->core.l_qseq;
}


static inline void
complement_seq(char *first, char *last) {
  for (; first != last; ++first) {
    assert(valid_base(*first));
    *first = complement(*first);
  }
}


static inline void
reverse(char *a, char *b) {
  char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
    assert(valid_base(*p1) && valid_base(*p2));
  }
}

static inline void
reverse(unsigned char *a, unsigned char *b) {
  unsigned char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
}



// return value is the number of cigar ops that are fully consumed in
// order to read n_ref, while "partial_oplen" is the number of bases
// that could be taken from the next operation, which might be merged
// with the other read.
static uint32_t
get_full_and_partial_ops(const uint32_t *cig_in, const uint32_t in_ops,
                         const uint32_t n_ref_full, uint32_t *partial_oplen) {
  // assume: n_ops <= size(cig_in) <= size(cig_out)
  size_t rlen = 0;
  uint32_t i = 0;
  for (i = 0; i < in_ops; ++i) {
    if (eats_ref(cig_in[i])) {
      if (rlen + bam_cigar_oplen(cig_in[i]) > n_ref_full)
        break;
      rlen += bam_cigar_oplen(cig_in[i]);
    }
  }
  *partial_oplen = n_ref_full - rlen;
  return i;
}


template <class T> void
revcomp_inplace(T first, T last) {
  std::transform(first, last, first, complement);
  std::reverse(first, last);
}


static void
revcomp_seq(bam1_t *aln) {
  // generate the sequence in ascii
  const auto seq = bam_get_seq(aln);
  const size_t l_qseq = get_qlen(aln);
  unsigned char *buf = (unsigned char *)malloc(l_qseq*sizeof(unsigned char));
  for (size_t i = 0; i < l_qseq; ++i)
    buf[i] = seq_nt16_str[bam_seqi(seq, i)];

  revcomp_inplace(buf, buf + l_qseq); // point of this function

  // copy it back...
  for (size_t i = 0; i < l_qseq; ++i)
    bam_set_seqi(seq, i, seq_nt16_table[buf[i]]);
}


const uint8_t byte_revcom_table[] ={
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  8, 136, 72, 0, 40, 0, 0, 0, 24, 0, 0, 0, 0, 0, 0, 248,
  4, 132, 68, 0, 36, 0, 0, 0, 20, 0, 0, 0, 0, 0, 0, 244,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  2, 130, 66, 0, 34, 0, 0, 0, 18, 0, 0, 0, 0, 0, 0, 242,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  1, 129, 65, 0, 33, 0, 0, 0, 17, 0, 0, 0, 0, 0, 0, 241,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  15, 143, 79, 0, 47, 0, 0, 0, 31, 0, 0, 0, 0, 0, 0, 255
};


static void
revcomp_seq_by_byte(bam1_t *aln) {
  const size_t l_qseq = get_qlen(aln);
  auto seq = bam_get_seq(aln);
  size_t num_bytes = ceil(l_qseq / 2.0);
  auto seq_end = seq + num_bytes;
  for (size_t i = 0; i < num_bytes; i++) {
    seq[i] = byte_revcom_table[seq[i]];
  }
  reverse(seq, seq_end);
  if (l_qseq % 2 == 1){
    for (size_t i = 0; i < num_bytes - 1; i++ ) {
      seq[i] = (seq[i] << 4) | (seq[i+1] >> 4 );
    }
    seq[num_bytes-1] <<= 4;
  }  
}


static inline bool
is_a_rich(const bam1_t *b) {return bam_aux2A(bam_aux_get(b, "CV")) == 'A';}


static inline bool
format_is_bam_or_sam(htsFile *hts) {
  const htsFormat *fmt = hts_get_format(hts);
  return fmt->category == sequence_data &&
    (fmt->format == bam || fmt->format == sam);
}


static void
flip_conversion(bam1_t *aln) {
  if (aln->core.flag & BAM_FREVERSE)
    aln->core.flag = aln->core.flag & (~BAM_FREVERSE);
  else
    aln->core.flag = aln->core.flag | BAM_FREVERSE;

  // generate the sequence in ascii
  auto a_seq = bam_get_seq(aln);
  const size_t a_len = aln->core.l_qseq;
  char *buf = (char *)malloc(a_len*sizeof(char));
  for (size_t i = 0; i < a_len; ++i)
    buf[i] = seq_nt16_str[bam_seqi(a_seq, i)];

  revcomp_inplace(buf, buf + a_len);

  for (size_t i = 0; i < a_len; ++i)
    bam_set_seqi(a_seq, i, seq_nt16_table[(unsigned char)buf[i]]);

  // ADS: don't like *(cv + 1) below, but no HTSlib function for it?
  uint8_t *cv = bam_aux_get(aln, "CV");
  if (!cv) throw fr_expt("bam_aux_get failed for CV");
  *(cv + 1) = 'T';
}


static bool
are_mates(const bam1_t *one, const bam1_t *two) {
  return one->core.mtid == two->core.tid &&
    one->core.mpos == two->core.pos &&
    bam_is_rev(one) != bam_is_rev(two);
  // below is a consistency check and should not be necessary
  /* &&
     two->core.mtid == one->core.tid &&
     two->core.mpos == one->core.pos; */
}


static int
truncate_overlap(const bam1_t *a, const uint32_t overlap, bam1_t *c) {

  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  uint32_t part_op = 0;
  const uint32_t c_cur =
    get_full_and_partial_ops(a_cig, a_ops, overlap, &part_op);

  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a->core.n_cigar && part_op > 0);

  const uint32_t c_ops = c_cur + use_partial;
  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));

  // ADS: replace this with a std::copy
  memcpy(c_cig, a_cig, c_cur*sizeof(uint32_t));
  // ADS: warning, if !use_partial, the amount of part_op used below
  // would make no sense.
  if (use_partial)
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
  /* after this point the cigar is set and should decide everything */

  const uint32_t c_seq_len = bam_cigar2qlen(c_ops, c_cig);
  char *c_seq = (char *)calloc(c_seq_len + 1, sizeof(char));
  if (!c_seq) throw fr_expt("allocating sequence");

  // copy the prefix of a into c; must be easier
  for (size_t i = 0; i < c_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  const uint16_t flag = a->core.flag & (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);

  int ret = bam_set1(c,
                     a->core.l_qname - (a->core.l_extranul + 1),
                     bam_get_qname(a),
                     flag,       // flags (SR and revcomp info)
                     a->core.tid,
                     a->core.pos,
                     a->core.qual,
                     c_ops,      // merged cigar ops
                     c_cig,      // merged cigar
                     -1,         // (no mate)
                     -1,         // (no mate)
                     isize,      // rlen from new cigar
                     c_seq_len,  // truncated seq length
                     c_seq,      // truncated sequence
                     NULL,       // no qual info
                     8);         // enough for the 2 tags?
  if (ret < 0) throw fr_expt(ret, "bam_set1");

  /* add the tags */
  const int64_t nm = bam_aux2i(bam_aux_get(a, "NM")); // ADS: do better here!
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw fr_expt(ret, "bam_aux_update_int");

  const uint8_t conversion = bam_aux2A(bam_aux_get(a, "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &conversion);
  if (ret < 0) throw fr_expt(ret, "bam_aux_append");

  return ret;
}


static int
merge_overlap(const bam1_t *a, const bam1_t *b,
              const uint32_t head, bam1_t *c) {
  assert(head > 0);

  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;

  uint32_t part_op = 0;
  uint32_t c_cur = get_full_and_partial_ops(a_cig, a_ops, head, &part_op);
  // ADS: hack here because the get_full_and_partial_ops doesn't do
  // exactly what is needed for this.
  const bool use_partial = (c_cur < a->core.n_cigar && part_op > 0);

  // check if the middle op would be the same
  const bool merge_mid =
    (use_partial > 0 ?
     bam_cigar_op(a_cig[c_cur]) == bam_cigar_op(b_cig[0]) :
     bam_cigar_op(a_cig[c_cur-1]) == bam_cigar_op(b_cig[0]));

  // c_ops: include the prefix of a_cig we need; then add for the
  // partial op; subtract for the identical op in the middle; finally
  // add the rest of b_cig.
  const uint32_t c_ops = c_cur + use_partial - merge_mid + b_ops;

  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));
  // std::fill(c_cig, c_cig + c_ops, std::numeric_limits<uint32_t>::max());
  memcpy(c_cig, a_cig, c_cur*sizeof(uint32_t));

  if (use_partial) {
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
    c_cur++; // index of dest for copying b_cig; faciltates corner case
  }
  // Here we get the length of a's sequence part contribution to c's
  // sequence before the possibility of merging the last entry with
  // the first entry in b's cigar. This is done with the cigar, so
  // everything depends on the "use_partial"
  const size_t a_seq_len = bam_cigar2qlen(c_cur, c_cig);
  /* ADS: above the return type of bam_cigar2qlen is uint64_t, but
     according to the source as of 05/2023 it cannot become
     negative; no possible error code returned */

  if (merge_mid) // update the middle op if it's the same
    c_cig[c_cur-1] = bam_cigar_gen(bam_cigar_oplen(c_cig[c_cur-1]) +
                                   bam_cigar_oplen(b_cig[0]),
                                   bam_cigar_op(b_cig[0]));
  // copy the cigar from b into c
  memcpy(c_cig + c_cur, b_cig + merge_mid, (b_ops-merge_mid)*sizeof(uint32_t));
  /* done with cigar string here */

  /* now deal with sequence */
  // allocate the right size for c's seq using the query-consuming ops
  // corresponding to the prefix of a copied into c.
  const uint32_t c_seq_len = a_seq_len + b->core.l_qseq;
  char *c_seq = (char *)calloc(c_seq_len + 1, sizeof(char));

  // copy the prefix of a into c
  for (size_t i = 0; i < a_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];
  // copy all of b into c
  for (size_t i = 0; i < get_qlen(b); ++i)
    c_seq[a_seq_len + i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];

  // reverse and complement the part of c corresponding to b
  revcomp_inplace(c_seq + a_seq_len, c_seq + c_seq_len);
  // reverse(c_seq + a_seq_len, c_seq + c_seq_len);
  // complement_seq(c_seq + a_seq_len, c_seq + c_seq_len);

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  uint16_t flag = (a->core.flag & (BAM_FREAD1 |
                                   BAM_FREAD2 |
                                   BAM_FREVERSE));

  int ret = bam_set1(c,
                     a->core.l_qname - (a->core.l_extranul + 1),
                     bam_get_qname(a),
                     flag,             // (no PE; revcomp info)
                     a->core.tid,
                     a->core.pos,
                     a->core.qual,     // mapq from "a" (consider update)
                     c_ops,            // merged cigar ops
                     c_cig,            // merged cigar
                     -1,               // (no mate)
                     -1,               // (no mate)
                     isize,            // updated
                     c_seq_len,        // merged sequence length
                     c_seq,            // merged sequence
                     NULL,             // no qual info
                     8);               // enough for 2 tags?
  if (ret < 0) throw fr_expt(ret, "bam_set1 in merge_overlap");

  // add the tag for mismatches
  const int64_t nm = (bam_aux2i(bam_aux_get(a, "NM")) +
                      bam_aux2i(bam_aux_get(b, "NM")));
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw fr_expt(ret, "bam_aux_update_int in merge_overlap");

  // add the tag for conversion
  const uint8_t cv = bam_aux2A(bam_aux_get(a, "CV"));
  ret = bam_aux_append(c, "CV", 'A', 1, &cv);
  if (ret < 0) throw fr_expt(ret, "bam_aux_append in merge_overlap");

  return ret;
}


static int
merge_non_overlap(const bam1_t *a, const bam1_t *b,
                  const uint32_t spacer, bam1_t *c) {

  /* make the cigar string */
  // collect info about the cigar strings
  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;
  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;
  // allocate the new cigar string
  const uint32_t c_ops = a_ops + b_ops + 1;
  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));
  // concatenate the new cigar strings with a "skip" in the middle
  memcpy(c_cig, a_cig, a_ops*sizeof(uint32_t));
  c_cig[a_ops] = bam_cigar_gen(spacer, BAM_CREF_SKIP);
  memcpy(c_cig + a_ops + 1, b_cig, b_ops*sizeof(uint32_t));
  /* done with cigars */

  /* now make the sequence */
  // get info about the lengths
  const size_t a_seq_len = get_qlen(a);
  const size_t b_seq_len = get_qlen(b);
  const size_t c_seq_len = a_seq_len + b_seq_len;
  // allocate and fill the new one as a char array
  char *c_seq = (char *)calloc(c_seq_len + 1, sizeof(char));
  for (size_t i = 0; i < a_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];
  for (size_t i = 0; i < b_seq_len; ++i)
    c_seq[a_seq_len + i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
  // reverse and complement the part corresponding to "b"
  revcomp_inplace(c_seq + a_seq_len, c_seq + c_seq_len);
  // reverse(c_seq + a_seq_len, c_seq + c_seq_len);
  // complement_seq(c_seq + a_seq_len, c_seq + c_seq_len);

  // get the template length from the cigar
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag: only need to keep strand and single-end info
  const uint16_t flag = a->core.flag & (BAM_FREAD1 |
                                        BAM_FREAD2 |
                                        BAM_FREVERSE);

  int ret =
    bam_set1(c,
             a->core.l_qname - (a->core.l_extranul + 1),
             bam_get_qname(a),
             flag,             // flags (no PE; revcomp info)
             a->core.tid,
             a->core.pos,
             a->core.qual,     // mapq from a (consider update)
             c_ops,            // merged cigar ops
             c_cig,            // merged cigar
             -1,               // (no mate)
             -1,               // (no mate)
             isize,            // TLEN (relative to reference; SAM docs)
             c_seq_len,        // merged sequence length
             c_seq,            // merged sequence
             NULL,             // no qual info
             8);               // enough for 2 tags of 1 byte value?
  if (ret < 0) throw fr_expt(ret, "bam_set1 in merge_non_overlap");

  /* add the tags */
  const int64_t nm = (bam_aux2i(bam_aux_get(a, "NM")) +
                      bam_aux2i(bam_aux_get(b, "NM")));
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) throw fr_expt(ret, "merge_non_overlap:bam_aux_update_int");

  const uint8_t cv = bam_aux2A(bam_aux_get(a, "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &cv);
  if (ret < 0) throw fr_expt(ret, "merge_non_overlap:bam_aux_append");

  return ret;
}


static int
keep_better_end(const bam1_t *a, const bam1_t *b, bam1_t *c) {
  c = bam_copy1(c, get_rlen(a) >= get_rlen(b) ? a : b);
  c->core.mtid = -1;
  c->core.mpos = -1;
  c->core.isize = get_rlen(c);
  c->core.flag &= (BAM_FREAD1 | BAM_FREAD2 | BAM_FREVERSE);
  return 0;
}


static size_t
merge_mates(const size_t range,
            bam1_t *one, bam1_t *two, bam1_t *merged) {

  if (!are_mates(one, two))
    return -std::numeric_limits<int>::max();

  // arithmetic easier using base 0 so subtracting 1 from pos
  const int one_s = one->core.pos;
  const int one_e = bam_endpos(one);
  const int two_s = two->core.pos;
  const int two_e = bam_endpos(two);
  assert(one_s >= 0 && two_s >= 0);

  const int spacer = two_s - one_e;
  if (spacer >= 0) {
    /* fragments longer enough that there is space between them: this
     * size of the spacer ("_") is determined based on the reference
     * positions of the two ends, and here we assume "one" maps to
     * positive genome strand.
     *                               spacer
     *                              <======>
     * left                                                         right
     * one_s                    one_e      two_s                    two_e
     * [------------end1------------]______[------------end2------------]
     */
    merge_non_overlap(one, two, spacer, merged);
  }
  else {
    const int head = two_s - one_s;
    if (head >= 0) {
      /* (Even if "head == 0" we will deal with it here.)
       *
       * CASE 1: head > 0
       *
       * fragment longer than or equal to the length of the left-most
       * read, but shorter than twice the read length (hence spacer
       * was < 0): this is determined by obtaining the size of the
       * "head" in the diagram below: the portion of end1 that is not
       * within [=]. If the read maps to the positive strand, this
       * depends on the reference start of end2 minus the reference
       * start of end1. For negative strand, this is reference start
       * of end1 minus reference start of end2.
       *
       * <======= head =========>
       *
       * left                                             right
       * one_s              two_s      one_e              two_e
       * [------------end1------[======]------end2------------]
       */
      if (head > 0) {
        merge_overlap(one, two, head, merged);
      }
      /* CASE 2: head == 0
       *
       * CASE 2A: one_e < two_e
       * left                                             right
       * one_s/two_s               one_e                  two_e
       * [=========== end1/end2========]------ end2 ----------]
       * keep "two"
       *
       * CASE 2B: one_e >= two_e
       * left                                             right
       * one_s/two_s               two_e                  one_e
       * [=========== end1/end2========]------ end1 ----------]
       * keep "one"
       *
       */
      // *** ELSE ***
      if (head == 0) { // keep the end with more ref bases
        keep_better_end(one, two, merged);
      }
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       *                 <=== overlap ==>
       * left                                       right
       * two_s           one_s      two_e           one_e
       * [--end2---------[==============]---------end1--]
       */
      const int overlap = two_e - one_s;
      if (overlap > 0) {
        truncate_overlap(one, overlap, merged);
      }
    }
  }

  // if merging two ends caused strange things in the cigar, fix them.
  correct_cigar(merged);

  return two_e - one_s;
}
/********Above are functions for merging pair-end reads********/


// ADS: will move to using this function once it is written
static void
standardize_format(const string &input_format, bam1_t *aln) {
  int err_code = 0;

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    // A/T rich; get the ZS tag value
    auto zs_tag = bam_aux_get(aln, "ZS");
    if (!zs_tag) throw fr_expt("bam_aux_get for ZS (invalid bsmap)");
    // ADS: test for errors on the line below
    const uint8_t cv = string(bam_aux2Z(zs_tag))[1] == '-' ? 'A' : 'T';
    // get the "mismatches" tag
    auto nm_tag = bam_aux_get(aln, "NM");
    if (!nm_tag) throw fr_expt("bam_aux_get for NM (invalid bsmap)");
    const int64_t nm = bam_aux2i(nm_tag);

    aln->l_data = bam_get_aux(aln) - aln->data; // del aux (no data resize)

    /* add the tags we want */
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here.
    err_code = bam_aux_update_int(aln, "NM", nm);
    if (err_code < 0) throw fr_expt(err_code, "bam_aux_update_int");
    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln, "CV", 'A', 1, &cv);
    if (err_code < 0) throw fr_expt(err_code, "bam_aux_append");

    if (bam_is_rev(aln)) revcomp_seq(aln); // reverse complement if needed
  }
  if (input_format == "bismark") {
    // ADS: Previously we modified the read names at the first
    // underscore. Even if the names are still that way, it should no
    // longer be needed since we compare names up to a learned suffix.

    // A/T rich; get the XR tag value
    auto xr_tag = bam_aux_get(aln, "XR");
    if (!xr_tag) throw fr_expt("bam_aux_get for XR (invalid bismark)");
    const uint8_t cv = string(bam_aux2Z(xr_tag)) == "GA" ? 'A' : 'T';
    // get the "mismatches" tag
    auto nm_tag = bam_aux_get(aln, "NM");
    if (!nm_tag) throw fr_expt("bam_aux_get for NM (invalid bismark)");
    const int64_t nm = bam_aux2i(nm_tag);

    aln->l_data = bam_get_aux(aln) - aln->data; // del aux (no data resize)

    /* add the tags we want */
    // "udpate" for "int" because it determines the right size; even
    // though we just deleted all tags, it will add it back here.
    err_code = bam_aux_update_int(aln, "NM", nm);
    if (err_code < 0) throw fr_expt(err_code, "bam_aux_update_int");
    // "append" for "char" because there is no corresponding update
    err_code = bam_aux_append(aln, "CV", 'A', 1, &cv);
    if (err_code < 0) throw fr_expt(err_code, "bam_aux_append");

    if (bam_is_rev(aln)) revcomp_seq(aln); // reverse complement if needed
  }

  // Be sure this doesn't depend on mapper! Removes the "qual" part of
  // the data in a bam1_t struct but does not change its uncompressed
  // size.
  const auto qs = bam_get_qual(aln);
  std::fill(qs, qs + aln->core.l_qseq, '\xff'); // deletes qseq
}


static vector<string>
load_read_names(const string &inputfile, const size_t n_reads) {
  samFile *hts = hts_open(inputfile.c_str(), "r");
  if (!hts) throw fr_expt("failed to open file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr) throw fr_expt("failed to read header: " + inputfile);

  bam1_t *aln = bam_init1();
  vector<string> names;
  size_t count = 0;
  int err_code = 0;

  while ((err_code = sam_read1(hts, hdr, aln)) >= 0 && count++ < n_reads)
    names.push_back(string(bam_get_qname(aln)));
  // err_core == -1 means EOF
  if (err_code < -1) fr_expt(err_code, "load_read_names:sam_read1");

  bam_destroy1(aln);

  return names;
}


static size_t
get_max_repeat_count(const vector<string> &names, const size_t suff_len) {
  // assume "suff_len" is shorter than the shortest entry in "names"
  size_t repeat_count = 0;
  size_t tmp_repeat_count = 0;
  // allow the repeat_count to go to 2, which might not be the "max"
  // but still would indicate that this suffix length is too long and
  // would result in more that two reads identified mutually as mates.
  for (size_t i = 1; i < names.size() && repeat_count < 2; ++i) {
    if (names[i-1].size() == names[i].size() &&
        equal(begin(names[i-1]), end(names[i-1]) - suff_len, begin(names[i])))
      ++tmp_repeat_count;
    else tmp_repeat_count = 0;
    repeat_count = std::max(repeat_count, tmp_repeat_count);
  }
  return repeat_count;
}


static bool
check_suff_len(const string &inputfile, const size_t suff_len,
               const size_t n_names_to_check) {
  /* thus function will indicate if the given suff_len would result in
     more than two reads being mutually considered mates */
  auto names(load_read_names(inputfile, n_names_to_check));
  // get the minimum read name length
  size_t min_name_len = std::numeric_limits<size_t>::max();
  for (auto &&i: names)
    min_name_len = std::min(min_name_len, i.size());
  if (min_name_len <= suff_len)
    throw fr_expt("given suffix length exceeds min read name length");
  sort(begin(names), end(names));
  return get_max_repeat_count(names, suff_len) < 2;
}


static size_t
guess_suff_len(const string &inputfile, const size_t n_names_to_check,
               size_t &repeat_count) {

  // ADS: assuming copy elision but should test it
  auto names(load_read_names(inputfile, n_names_to_check));

  // get the minimum read name length
  size_t min_name_len = std::numeric_limits<size_t>::max();
  for (auto &&i: names)
    min_name_len = std::min(min_name_len, i.size());
  assert(min_name_len > 0);

  sort(begin(names), end(names));

  // these will be returned
  size_t suff_len = 0;
  repeat_count = 0;

  // check the possible read name suffix lengths; if any causes a
  // repeat count of more than 1 (here this means == 2), all greater
  // suffix lengths will also
  const size_t max_suff_len = min_name_len - 1;
  while (suff_len < max_suff_len && repeat_count == 0) {
    // check current suffix length guess
    repeat_count = get_max_repeat_count(names, suff_len);
    // we want to lag by one iteration
    if (repeat_count == 0)
      ++suff_len;
  }
  // repeat_count should be equal to the greatest value of
  // tmp_repeat_count over all inner iterations above. If this value
  // is not 1, it will indicate whether we have exactly hit a max
  // repeat count of 2, indicating mates, or exceeded it, indicating
  // there seems not to be a good suffix length to remove for
  // identifying mates

  return suff_len;
}


static string
remove_suff(const string &s, const size_t suff_len) {
  return s.size() > suff_len ? s.substr(0, s.size() - suff_len) : s;
}


static bool
check_sorted(const string &inputfile, const size_t suff_len, size_t n_reads) {
  // In order to check if mates are consecutive we need to check if a
  // given end has a mate and that mate is not adjacent. This requires
  // storing previous reads, not simply checking for adjacent pairs.
  auto names(load_read_names(inputfile, n_reads));
  for (auto &&i : names)
    i = remove_suff(i, suff_len);

  std::unordered_map<string, size_t> mate_lookup;
  for (size_t i = 0; i < names.size(); ++i) {
    auto the_mate = mate_lookup.find(names[i]);
    if (the_mate == end(mate_lookup)) // 1st time seeing this one
      mate_lookup[names[i]] = i;
    else if (the_mate->second != i - 1)
      return false;
  }
  // made it here: all reads with mates are consecutive
  return true;
}


static bool
check_input_file(const string &infile) {
  samFile* hts = hts_open(infile.c_str(), "r");
  if (!hts || errno) throw fr_expt("error opening: " + infile);
  const htsFormat *fmt = hts_get_format(hts);
  if (fmt->category != sequence_data)
    throw fr_expt("not sequence data: " + infile);
  if (fmt->format != bam && fmt->format != sam)
    throw fr_expt("not SAM/BAM format: " + infile);
  return true;
}


static bool
check_format_in_header(const string &input_format, const string &inputfile) {
  samFile* hts = hts_open(inputfile.c_str(), "r");
  if (!hts) throw fr_expt("error opening file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr) throw fr_expt("failed to read header: " + inputfile);

  auto begin_hdr = sam_hdr_str(hdr);
  auto end_hdr = begin_hdr + std::strlen(begin_hdr);
  auto it = std::search(begin_hdr, end_hdr,
                        begin(input_format), end(input_format),
                        [](const unsigned char a, const unsigned char b) {
                          return std::toupper(a) == std::toupper(b);
                        });
  return it != end_hdr;
}


static bool
same_name(const bam1_t *a, const bam1_t *b, const size_t suff_len) {
  // "+ 1" below: extranul counts *extras*; we don't want *any* nulls
  const uint16_t a_l = a->core.l_qname - (a->core.l_extranul + 1);
  const uint16_t b_l = b->core.l_qname - (b->core.l_extranul + 1);
  if (a_l != b_l) return false;
  assert(a_l > suff_len);
  return !std::strncmp(bam_get_qname(a), bam_get_qname(b), a_l - suff_len);
}


static void
add_pg_line(const string &cmd, sam_hdr_t *hdr) {
  int err_code =
    sam_hdr_add_line(hdr, "PG", "ID", "DNMTOOLS", "VN",
                     VERSION, "CL", cmd.c_str(), NULL);
  if (err_code) throw fr_expt(err_code, "failed to add pg header line");
}


static void
format(const string &cmd, const size_t n_threads,
       const string &inputfile, const string &outfile,
       const bool bam_format, const string &input_format,
       const size_t suff_len, const size_t max_frag_len) {

  int err_code = 0;

  // open the hts files; assume already checked
  samFile *hts = hts_open(inputfile.c_str(), "r");
  samFile *out = hts_open(outfile.c_str(), bam_format ? "wb" : "w");

  // set the threads
  htsThreadPool the_thread_pool{hts_tpool_init(n_threads), 0};
  err_code = hts_set_thread_pool(hts, &the_thread_pool);
  if (err_code < 0) throw fr_expt("error setting threads");
  err_code = hts_set_thread_pool(out, &the_thread_pool);
  if (err_code < 0) throw fr_expt("error setting threads");

  // headers: load the input file's header, and then update to the
  // output file's header, then write it and destroy; we will only use
  // the input file header.
  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr) throw fr_expt("failed to read header");
  sam_hdr_t *hdr_out = bam_hdr_dup(hdr);
  if (!hdr_out) throw fr_expt("failed create header");
  add_pg_line(cmd, hdr_out);
  err_code = sam_hdr_write(out, hdr_out);
  if (err_code) throw fr_expt(err_code, "failed to output header");

  // now process the reads
  bam1_t *aln = bam_init1();
  bam1_t *prev_aln = bam_init1();
  bam1_t *merged = bam_init1();
  bool previous_was_merged = false;

  err_code = sam_read1(hts, hdr, aln); // for EOF, err_code == -1
  if (err_code < -1) throw fr_expt(err_code, "format:sam_read1");

  std::swap(aln, prev_aln); // start with prev_aln being first read

  while ((err_code = sam_read1(hts, hdr, aln)) >= 0) {
    standardize_format(input_format, aln);
    if (same_name(prev_aln, aln, suff_len)) {
      // below: essentially check for dovetail
      if (!bam_is_rev(aln)) std::swap(prev_aln, aln);
      const size_t frag_len = merge_mates(max_frag_len, prev_aln, aln, merged);
      if (frag_len > 0 && frag_len < max_frag_len) {
        if (is_a_rich(merged)) flip_conversion(merged);
        err_code = sam_write1(out, hdr, merged);
        if (err_code < 0) throw fr_expt(err_code, "format:sam_write1");
      }
      else {
        if (is_a_rich(prev_aln)) flip_conversion(prev_aln);
        err_code = sam_write1(out, hdr, prev_aln);
        if (err_code < 0) throw fr_expt(err_code, "format:sam_write1");
        if (is_a_rich(aln)) flip_conversion(aln);
        err_code = sam_write1(out, hdr, aln);
        if (err_code < 0) throw fr_expt(err_code, "format:sam_write1");
      }
      previous_was_merged = true;
    }
    else {
      if (!previous_was_merged) {
        if (is_a_rich(prev_aln)) flip_conversion(prev_aln);
        err_code = sam_write1(out, hdr, prev_aln);
        if (err_code < 0) throw fr_expt(err_code, "format:sam_write1");
      }
      previous_was_merged = false;
    }
    std::swap(prev_aln, aln);
  }
  if (err_code < -1) throw fr_expt(err_code, "format:sam_read1");

  if (!previous_was_merged) {
    if (is_a_rich(prev_aln)) flip_conversion(prev_aln);
    err_code = sam_write1(out, hdr, prev_aln);
    if (err_code < 0) throw fr_expt(err_code, "format:sam_write1");
  }

  // turn off the lights
  bam_destroy1(prev_aln);
  bam_destroy1(aln);
  bam_destroy1(merged);
  bam_hdr_destroy(hdr);
  bam_hdr_destroy(hdr_out);
  err_code = hts_close(hts);
  if (err_code < 0) throw fr_expt(err_code, "format:hts_close");
  err_code = hts_close(out);
  if (err_code < 0) throw fr_expt(err_code, "format:hts_close");
  // do this after the files have been closed
  hts_tpool_destroy(the_thread_pool.pool);
}


int
main_format(int argc, const char **argv) {

  try {

    size_t n_reads_to_check = 1000000;

    bool bam_format = false;
    bool use_stdout = false;

    string input_format;
    string outfile;
    int max_frag_len = std::numeric_limits<int>::max();
    size_t suff_len = 0;
    bool single_end = false;
    bool VERBOSE = false;
    size_t n_threads = 1;

    const string description = "convert SAM/BAM mapped bs-seq reads "
      "to standard dnmtools format";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<sam/bam-file> [out-file]", 2);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("bam", 'B', "output in BAM format", false, bam_format);
    opt_parse.add_opt("stdout", '\0',
                      "write to standard output", false, use_stdout);
    opt_parse.add_opt("format", 'f', "input format {abismal, bsmap, bismark}",
                      false, input_format);
    opt_parse.add_opt("suff", 's', "read name suffix length",
                      false, suff_len);
    opt_parse.add_opt("single-end", '\0',
                      "assume single-end [do not use with -suff]",
                      false, single_end);
    opt_parse.add_opt("max-frag", 'L', "maximum allowed insert size",
                      false, max_frag_len);
    opt_parse.add_opt("check", 'c',
                      "check this many reads to validate read name suffix",
                      false, n_reads_to_check);
    opt_parse.add_opt("verbose", 'v', "print more information",
                      false, VERBOSE);
    opt_parse.set_show_defaults();
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (opt_parse.about_requested() || opt_parse.help_requested() ||
        leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_FAILURE;
    }
    if (suff_len != 0 && single_end) {
      cerr << "incompatible arguments specified" << endl
           << opt_parse.help_message() << endl;
      return EXIT_FAILURE;
    }
    if ((leftover_args.size() == 1 && !use_stdout) ||
        (leftover_args.size() == 2 && use_stdout)) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_FAILURE;
    }
    const string infile(leftover_args.front());
    if (leftover_args.size() == 2 && !use_stdout)
      outfile = leftover_args.back();
    else
      outfile = string("-"); // so htslib can write to stdout
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char*>(cmd, " "));

    check_input_file(infile);

    if (VERBOSE)
      if (!check_format_in_header(input_format, infile))
        cerr << "[warning: input format not found in header "
             << "(" << input_format << ", " << infile << ")]" << endl;

    if (!single_end) {
      if (suff_len == 0) {
        size_t repeat_count = 0;
        suff_len = guess_suff_len(infile, n_reads_to_check, repeat_count);
        if (repeat_count > 1)
          throw fr_expt("failed to identify read name suffix length\n"
                        "verify reads are not single-end\n"
                        "specify read name suffix length directly");
        if (VERBOSE)
          cerr << "[read name suffix length guess: " << suff_len << "]" << endl;
      }
      else if (!check_suff_len(infile, suff_len, n_reads_to_check))
        throw fr_expt("wrong read name suffix length [" +
                      std::to_string(suff_len) + "] in: " + infile);
      if (!check_sorted(infile, suff_len, n_reads_to_check))
        throw fr_expt("mates not consecutive in: " + infile);
    }

    if (VERBOSE)
      cerr << "[input file: " << infile << "]" << endl
           << "[mapper: " << input_format << "]" << endl
           << "[configuration: " << (single_end ? "SE" : "PE") << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: " << (bam_format ? "B" : "S") << "AM]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl
           << "[readname suffix length: " << suff_len << "]" << endl;

    format(cmd.str(), n_threads, infile, outfile,
           bam_format, input_format, suff_len, max_frag_len);
  }
  catch (const std::exception &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
