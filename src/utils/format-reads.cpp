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
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::max;
using std::min;
using std::runtime_error;
using std::unordered_map;
using std::swap;
using std::to_string;
using std::ostringstream;


void
complement_seq(char *first, char *last) {
  while (first != last) {
    *first = complement(*first);
    assert(valid_base(*first));
  }
}


void
reverse(char *a, char *b) {
  char *p1, *p2;
  for (p1 = a, p2 = b - 1; p2 > p1; ++p1, --p2) {
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
    assert(*p1 == 'A' ||
           *p1 == 'C' ||
           *p1 == 'G' ||
           *p1 == 'T');
    assert(*p2 == 'A' ||
           *p2 == 'C' ||
           *p2 == 'G' ||
           *p2 == 'T');
  }
}


void
error_and_exit(const int64_t e, const char *msg) {
  fprintf(stderr, "ERROR: %ld %d %s %s\n",
          e, errno, strerror(errno), msg);
  exit(e);
}


static inline bool
consumes_reference(uint32_t op) {
  return bam_cigar_type(bam_cigar_op(op)) & 2;
}


static uint32_t
get_full_and_partial_ops(const uint32_t *cig_in, const uint32_t in_ops,
                         const uint32_t n_ref, uint32_t *partial_oplen) {
  // assume: n_ops <= size(cig_in) <= size(cig_out)
  size_t bases = 0;
  int i = 0;
  for (i = 0; i < in_ops; ++i) {
    if (consumes_reference(cig_in[i])) {
      if (bases + bam_cigar_oplen(cig_in[i]) > n_ref)
        break;
      bases += bam_cigar_oplen(cig_in[i]);
    }
  }
  *partial_oplen = n_ref - bases;
  return i;
}


void
swap_bams(bam1_t **a, bam1_t **b) {
  bam1_t *c = *a;
  *a = *b;
  *b = c;
}


template <class T> void
revcomp_inplace(T first, T last) {
  std::transform(first, last, first, complement);
  std::reverse(first, last);
}


inline bool
is_a_rich(const bam1_t *aln) {
  return bam_aux2A(bam_aux_get(aln, "CV")) == 'A';
}


inline bool
format_is_bam_or_sam(htsFile *hts) {
  const htsFormat *fmt = hts_get_format(hts);
  return fmt->category == sequence_data &&
    (fmt->format == bam || fmt->format == sam);
}


static bam1_t *
get_read(samFile *hts, sam_hdr_t *hdr) {
  bam1_t *b = bam_init1();
  const int result = sam_read1(hts, hdr, b);
  if (result >= 0) return b;

  if (result < -1)
    throw runtime_error("error reading file: " + string(hts->fn));
  else // -1 should mean EOF, so we free this read
    bam_destroy1(b);
  return 0;
}


inline bool
is_rc(const bam1_t *aln) {
  return bam_is_rev(aln);
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
  for (int i = 0; i < a_len; ++i)
    buf[i] = seq_nt16_str[bam_seqi(a_seq, i)];

  revcomp_inplace(buf, buf + a_len);

  for (int i = 0; i < a_len; ++i)
    bam_set_seqi(a_seq, i, seq_nt16_table[(unsigned char)buf[i]]);

  uint8_t *cv = bam_aux_get(aln, "CV");
  *(cv + 1) = 'T';

  // uint8_t updated_cv = bam_aux2A(original_cv) == 'A' ? 'T' : 'A';
  // int e = bam_aux_del(aln, original_cv);
  // e = bam_aux_append(aln, "CV", 'A', 1, &updated_cv);
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
truncate_overlapping(const bam1_t *a, const uint32_t overlap, bam1_t *c) {

  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  uint32_t part_op = 0;
  const uint32_t c_cur =
    get_full_and_partial_ops(a_cig, a_ops, overlap, &part_op);

  const uint32_t c_ops = c_cur + (part_op > 0);
  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));

  memcpy(c_cig, a_cig, c_cur*sizeof(uint32_t));
  if (part_op > 0)
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));

  const uint32_t c_seq_len = bam_cigar2qlen(c_ops, c_cig);
  char *c_seq = (char *)calloc(c_seq_len + 1, sizeof(char));
  if (!c_seq) error_and_exit((int64_t)c_seq, "allocating c_seq");

  // copy the prefix of a into c
  for (int i = 0; i < c_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  uint16_t flag = (a->core.flag & (BAM_FREAD1 |
                                   BAM_FREAD2 |
                                   BAM_FREVERSE));

  int ret = bam_set1(c,
                     strlen(bam_get_qname(a)), // name length from "a"
                     bam_get_qname(a), // get name from "a"
                     flag,             // flags (no PE; revcomp info)
                     a->core.tid,      // tid
                     a->core.pos,      // pos
                     a->core.qual,     // mapq from a (consider update)
                     c_ops,            // merged cigar ops
                     c_cig,            // merged cigar
                     -1,               // (no mate)
                     -1,               // (no mate)
                     isize,            // TLEN (relative to reference; SAM docs)
                     c_seq_len,        // merged sequence length
                     c_seq,            // merged sequence
                     NULL,             // no qual info
                     8);               // enough for 2 tags?
  if (ret < 0) error_and_exit(ret, "bam_set1 in truncate_overlapping");

  /* add the tags */
  const int64_t nm = bam_aux2i(bam_aux_get(a, "NM")); // ADS: do better here!
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) error_and_exit(ret, "bam_aux_update_int in truncate_overlapping");

  const uint8_t conversion = bam_aux2A(bam_aux_get(a, "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &conversion);
  if (ret < 0) error_and_exit(ret, "bam_aux_append in truncate_overlapping");

  return ret;
}



static int
merge_overlapping(const bam1_t *a, const bam1_t *b,
                  const uint32_t head, bam1_t *c) {

  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;

  uint32_t part_op = 0;
  uint32_t c_cur = get_full_and_partial_ops(a_cig, a_ops, head, &part_op);

  // check if the middle op would be the same
  const bool merge_mid = head != 0 &&
    (part_op > 0 ?
     bam_cigar_op(a_cig[c_cur]) == bam_cigar_op(b_cig[0]) :
     bam_cigar_op(a_cig[c_cur-1]) == bam_cigar_op(b_cig[0]));

  // c_ops: include the prefix of a_cig we need; then add for the
  // partial op; subtract for the identical op in the middle; finally
  // add the rest of b_cig.
  const uint32_t c_ops = c_cur + (head != 0 || part_op > 0) - merge_mid + b_ops;
  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));

  memcpy(c_cig, a_cig, c_cur*sizeof(uint32_t));
  if (head != 0 && part_op > 0) {
    c_cig[c_cur] = bam_cigar_gen(part_op, bam_cigar_op(a_cig[c_cur]));
    c_cur++; // index of dest for copying b_cig; faciltates corner case
  }
  // here we get the length of a's sequence part contribution to c's
  // sequence before the possibility of merging the last entry with
  // the first entry in b's cigar
  const uint32_t a_seq_len = bam_cigar2qlen(c_cur, c_cig);

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
  for (int i = 0; i < a_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];
  // copy all of b into c
  for (int i = 0; i < b->core.l_qseq; ++i)
    c_seq[a_seq_len + i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];

  // reverse and complement the part of c corresponding to b
  reverse(c_seq + a_seq_len, c_seq + c_seq_len);
  complement_seq(c_seq + a_seq_len, c_seq + c_seq_len);

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // check for internal soft-clip?

  // flag only needs to worry about strand and single-end stuff
  uint16_t flag = (a->core.flag & (BAM_FREAD1 |
                                   BAM_FREAD2 |
                                   BAM_FREVERSE));

  int ret = bam_set1(c,
                     strlen(bam_get_qname(a)), // name length from "a"
                     bam_get_qname(a), // get name from "a"
                     flag,             // flags (no PE; revcomp info)
                     a->core.tid,      // tid
                     a->core.pos,      // pos
                     a->core.qual,     // mapq from a (consider update)
                     c_ops,            // merged cigar ops
                     c_cig,            // merged cigar
                     -1,               // (no mate)
                     -1,               // (no mate)
                     isize,            // TLEN (relative to reference; SAM docs)
                     c_seq_len,        // merged sequence length
                     c_seq,            // merged sequence
                     NULL,             // no qual info
                     8);               // enough for 2 tags?
  if (ret < 0) error_and_exit(ret, "bam_set1 in merge_overlapping");

  // add the tag for mismatches
  const int64_t n_mismatches =
    bam_aux2i(bam_aux_get(a, "NM")) + bam_aux2i(bam_aux_get(b, "NM"));
  ret = bam_aux_update_int(c, "NM", n_mismatches);
  if (ret < 0) error_and_exit(ret, "set mismatches in merge_overlapping");

  // add the tag for conversion
  const uint8_t conversion = bam_aux2A(bam_aux_get(a, "CV"));
  ret = bam_aux_append(c, "CV", 'A', 1, &conversion);
  if (ret < 0) error_and_exit(ret, "set cv in merge_overlapping");

  return ret;
}


static int
merge_non_overlapping(const bam1_t *a, const bam1_t *b,
                      const uint32_t spacer, bam1_t *c) {
  // ADS: convert internal softclip to insertion (BAM_CINS consumes query)

  /* make the cigar string */
  const uint32_t *a_cig = bam_get_cigar(a);
  const uint32_t a_ops = a->core.n_cigar;

  const uint32_t *b_cig = bam_get_cigar(b);
  const uint32_t b_ops = b->core.n_cigar;

  const uint32_t c_ops = a_ops + b_ops + 1;
  uint32_t *c_cig = (uint32_t*)calloc(c_ops, sizeof(uint32_t));

  // copy the cigars into c
  memcpy(c_cig, a_cig, a_ops*sizeof(uint32_t));
  c_cig[a_ops] = bam_cigar_gen(spacer, BAM_CREF_SKIP);
  memcpy(c_cig + a_ops + 1, b_cig, b_ops*sizeof(uint32_t));
  /* done with cigars */

  /* now make the sequence */
  const uint32_t a_seq_len = a->core.l_qseq;
  const uint32_t b_seq_len = b->core.l_qseq;
  const uint32_t c_seq_len = a_seq_len + b_seq_len;
  char *c_seq = (char *)calloc(c_seq_len + 1, sizeof(char));
  for (int i = 0; i < a_seq_len; ++i)
    c_seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(a), i)];
  for (int i = 0; i < b_seq_len; ++i)
    c_seq[a_seq_len + i] = seq_nt16_str[bam_seqi(bam_get_seq(b), i)];
  // reverse and complement the part of c corresponding to b
  reverse(c_seq + a_seq_len, c_seq + c_seq_len);
  complement_seq(c_seq + a_seq_len, c_seq + c_seq_len);

  // get the template length
  const hts_pos_t isize = bam_cigar2rlen(c_ops, c_cig);

  // flag only needs to worry about strand and single-end stuff
  uint16_t flag = (a->core.flag & (BAM_FREAD1 |
                                   BAM_FREAD2 |
                                   BAM_FREVERSE));

  int ret = bam_set1(c,
                     strlen(bam_get_qname(a)), // name length from "a"
                     bam_get_qname(a), // get name from "a"
                     flag,             // flags (no PE; revcomp info)
                     a->core.tid,      // tid
                     a->core.pos,      // pos
                     a->core.qual,     // mapq from a (consider update)
                     c_ops,            // merged cigar ops
                     c_cig,            // merged cigar
                     -1,               // (no mate)
                     -1,               // (no mate)
                     isize,            // TLEN (relative to reference; SAM docs)
                     c_seq_len,        // merged sequence length
                     c_seq,            // merged sequence
                     NULL,             // no qual info
                     8);               // enough for 2 tags?
  if (ret < 0) error_and_exit(ret, "bam_set1 in merge_non_overlapping");

  /* add the tags */
  const int64_t nm = (bam_aux2i(bam_aux_get(a, "NM")) +
                      bam_aux2i(bam_aux_get(b, "NM")));
  // "udpate" for "int" because it determines the right size
  ret = bam_aux_update_int(c, "NM", nm);
  if (ret < 0) error_and_exit(ret, "set n_mismatches in merge_non_overlapping");

  const uint8_t conversion = bam_aux2A(bam_aux_get(a, "CV"));
  // "append" for "char" because there is no corresponding update
  ret = bam_aux_append(c, "CV", 'A', 1, &conversion);
  if (ret < 0) error_and_exit(ret, "set cv in merge_non_overlapping");

  return ret;
}


static size_t
merge_mates(const size_t range,
            const bam1_t *one, const bam1_t *two, bam1_t *merged) {

  if (!are_mates(one, two))
    return -std::numeric_limits<int>::max();

  // ADS: tests show that revcomp is now different from before. Not
  // sure if this has to do with a revcomp deleted recently, or if
  // it's because of something else.

  // ADS: not sure this can be consistent across mappers
  // GS: true after standardization

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
     *
     * left                                                         right
     * one_s                    one_e      two_s                    two_e
     * [------------end1------------]______[------------end2------------]
     */
    int res = merge_non_overlapping(one, two, spacer, merged);
    // ADS: need to take care of soft clipping in between;
  }
  else {
    const int head = two_s - one_s;
    if (head >= 0) {
      /* fragment longer than or equal to the read length, but shorter
       * than twice the read length: this is determined by obtaining
       * the size of the "head" in the diagram below: the portion of
       * end1 that is not within [=]. If the read maps to the positive
       * strand, this depends on the reference start of end2 minus the
       * reference start of end1. For negative strand, this is
       * reference start of end1 minus reference start of end2.
       *
       * <======= head =========>
       *
       * left                                             right
       * one_s              two_s      one_e              two_e
       * [------------end1------[======]------end2------------]
       */
      int res = merge_overlapping(one, two, head, merged);
      // ADS: need to take care of soft clipping in between
    }
    else {
      /* dovetail fragments shorter than read length: this is
       * identified if the above conditions are not satisfied, but
       * there is still some overlap. The overlap will be at the 5'
       * ends of reads, which in theory shouldn't happen unless the
       * two ends are covering identical genomic intervals.
       *
       * left                                       right
       * two_s            one_s    two_e            one_e
       * [--end2----------[============]----------end1--]
       */
      const int overlap = two_e - one_s;
      if (overlap >= 0) {
        int res = truncate_overlapping(one, overlap, merged);
      }
    }
  }

  return two_e - one_s;
}
/********Above are functions for merging pair-end reads********/

// ADS: there is a bug somewhere when a value of 0 is given for
// suffix_len
static string
remove_suff(const string &x, const size_t suffix_len) {
  return x.size() > suffix_len ? x.substr(0, x.size() - suffix_len) : x;
}

bool
bsmap_get_rc(const string &strand_tag) {
  return strand_tag.size() > 5 && strand_tag[5] == '-';
}

bool
bsmap_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 6 && richness_tag[6] == '-';
}

bool
bismark_get_a_rich(const string &richness_tag) {
  return richness_tag.size() > 5 && richness_tag[5] == 'G';
}


// ADS: will move to using this function once it is written
static void
standardize_format(const string &input_format, bam1_t *aln) {

  if (input_format == "abismal" || input_format == "walt") return;

  if (input_format == "bsmap") {
    int res = 0; // result of htslib calls
    kstring_t s = KS_INITIALIZE;
    res = bam_aux_get_str(aln, "ZS", &s);
    if (res != 1)
      throw runtime_error("record appears to be invalid for bsmap");
    const bool a_rich = s.l > 6 && s.s[6] == '-';

    uint8_t *zs_tag = bam_aux_get(aln, "ZS");
    if (!zs_tag)
      throw runtime_error("error interpreting bam record");
    res = bam_aux_del(aln, zs_tag);
    if (res < 0)
      throw runtime_error("error updating tag");

    const uint8_t cv_val = a_rich ? 'A' : 'T';
    res = bam_aux_append(aln, "CV", 'A', sizeof(char), &cv_val);
    if (res == -1)
      throw runtime_error("error updating tag");

    if (bam_is_rev(aln)) {
      uint8_t *begin_seq = bam_get_seq(aln);
      uint8_t *end_seq = begin_seq + aln->core.l_qseq;
      revcomp_inplace(begin_seq, end_seq);
    }
  }
  // doesn't depend on mapper
  // aln.qual = "*";
}


static vector<string>
load_read_names(const string &inputfile, const size_t n_reads) {
  samFile *hts = hts_open(inputfile.c_str(), "r");
  if (!hts)
    throw runtime_error("error opening file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr)
    throw runtime_error("failed to read header: " + inputfile);

  bam1_t *aln = 0;
  vector<string> names;
  size_t count = 0;
  while ((aln = get_read(hts, hdr)) && count++ < n_reads)
    names.push_back(string(bam_get_qname(aln)));
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
    repeat_count = max(repeat_count, tmp_repeat_count);
  }
  return repeat_count;
}


static bool
check_suff_len(const string &inputfile, const size_t suff_len,
               const size_t n_names_to_check) {
  auto names(load_read_names(inputfile, n_names_to_check));
  // get the minimum read name length
  size_t min_name_len = std::numeric_limits<size_t>::max();
  for (auto &&i: names)
    min_name_len = min(min_name_len, i.size());
  if (min_name_len <= suff_len)
    throw runtime_error("given suffix length exceeds min read name length" );
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
    min_name_len = min(min_name_len, i.size());
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


static bool
check_sorted(const string &inputfile, const size_t suff_len, size_t n_reads) {
  // In order to check if mates are consecutive we need to check if a
  // given end has a mate and that mate is not adjacent. This requires
  // storing previous reads, not simply checking for adjacent pairs.

  auto names(load_read_names(inputfile, n_reads));
  for (auto &&i : names)
    i = remove_suff(i, suff_len);

  unordered_map<string, size_t> mate_lookup;
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
check_input_file(const string &input_filename) {
  samFile* hts = hts_open(input_filename.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("error opening file: " + input_filename);
  const htsFormat *fmt = hts_get_format(hts);
  if (fmt->category != sequence_data)
    throw runtime_error("not sequence data: " + input_filename);
  if (fmt->format != bam && fmt->format != sam)
    throw runtime_error("not SAM/BAM file: " + input_filename);
  return true;
}


static bool
check_format_in_header(const string &input_format, const string &inputfile) {
  samFile* hts = hts_open(inputfile.c_str(), "r");
  if (!hts || errno)
    throw runtime_error("error opening file: " + inputfile);

  sam_hdr_t *hdr = sam_hdr_read(hts);
  if (!hdr)
    throw runtime_error("failed to read header: " + inputfile);

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
  const uint16_t a_l = a->core.l_qname - (a->core.l_extranul + 1);
  const uint16_t b_l = b->core.l_qname - (b->core.l_extranul + 1);

  if (a_l != b_l) return false;
  assert(a_l > suff_len);

  return !std::strncmp(bam_get_qname(a), bam_get_qname(b), a_l - suff_len);
}


static void
format(const string &cmd, const size_t n_threads,
       const string &inputfile, const string &outfile,
       const bool bam_format, const string &input_format,
       const size_t suff_len, const size_t max_frag_len) {

  // open the hts files
  samFile* hts = hts_open(inputfile.c_str(), "r");
  samFile *out = hts_open(outfile.c_str(), bam_format ? "wb" : "w");

  // set the threads
  htsThreadPool the_thread_pool{hts_tpool_init(n_threads), 0};
  if (hts_set_thread_pool(hts, &the_thread_pool) < 0)
    throw runtime_error("error setting threads");
  if (hts_set_thread_pool(out, &the_thread_pool) < 0)
    throw runtime_error("error setting threads");

  // headers: load the input file's header, and then update to the
  // output file's header, then write it and destroy; we will only use
  // the input file header.
  sam_hdr_t *hdr = sam_hdr_read(hts);
  sam_hdr_t *hdr_out = bam_hdr_dup(hdr);
  if (sam_hdr_add_line(hdr_out, "PG", "ID", "DNMTOOLS",
                       "VN", VERSION, "CL", cmd.c_str(), NULL))
    throw runtime_error("failed to format header");
  if (sam_hdr_write(out, hdr_out))
    throw runtime_error("failed to output header");

  // now process the reads
  // sam_rec aln, prev_aln;
  bam1_t *aln = bam_init1();
  bam1_t *prev_aln = bam_init1();
  bam1_t *merged = bam_init1();
  // string read_name, prev_name;
  bool previous_was_merged = false;

  int res = sam_read1(hts, hdr, aln);
  if (res < 0) error_and_exit((int64_t)aln, "K");

  swap_bams(&aln, &prev_aln);

  // bam1_t *bam_rec;
  while ((res = sam_read1(hts, hdr, aln)) >= 0) {

    // get_sam_record(hdr, bam_rec, aln);
    standardize_format(input_format, aln);
    // read_name = remove_suff(aln.qname, suff_len);

    if (same_name(prev_aln, aln, suff_len)) {
      if (!is_rc(aln)) // essentially check for dovetail
        swap_bams(&prev_aln, &aln); // swap(prev_aln, aln);
      const int frag_len = merge_mates(max_frag_len, prev_aln, aln, merged);
      if (frag_len > 0 && frag_len < max_frag_len) {
        if (is_a_rich(merged))
          flip_conversion(merged);
        if (sam_write1(out, hdr, merged) < 0)
          throw runtime_error("failed writing bam record");
      }
      else {
        if (is_a_rich(aln))
          flip_conversion(aln);
        if (is_a_rich(prev_aln))
          flip_conversion(prev_aln);
        if (sam_write1(out, hdr, prev_aln) < 0 ||
            sam_write1(out, hdr, aln) < 0)
          throw runtime_error("failed writing bam record");
      }
      previous_was_merged = true;
    }
    else {
      if (!previous_was_merged) {
        if (is_a_rich(prev_aln))
          flip_conversion(prev_aln);
        if (sam_write1(out, hdr, prev_aln) < 0)
          throw runtime_error("failed writing bam record");
      }
      previous_was_merged = false;
    }
    swap_bams(&prev_aln, &aln);
  }

  if (!previous_was_merged) {
    if (is_a_rich(prev_aln))
      flip_conversion(prev_aln);
    if (sam_write1(out, hdr, prev_aln) < 0)
      throw runtime_error("failed writing bam record");
  }

  // turn off the lights
  bam_destroy1(prev_aln);
  bam_destroy1(aln);
  bam_hdr_destroy(hdr);
  hts_close(hts);
  bam_hdr_destroy(hdr_out);
  hts_close(out);
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
    const string inputfile(leftover_args.front());
    if (leftover_args.size() == 2 && !use_stdout)
      outfile = leftover_args.back();
    else
      outfile = string("-"); // so htslib can write to stdout
    /****************** END COMMAND LINE OPTIONS *****************/

    std::ostringstream cmd;
    copy(argv, argv + argc, std::ostream_iterator<const char*>(cmd, " "));

    check_input_file(inputfile);

    if (VERBOSE)
      if (!check_format_in_header(input_format, inputfile))
        cerr << "[warning: input format not found in header "
             << "(" << input_format << ", " << inputfile << ")]" << endl;

    if (!single_end) {
      if (suff_len == 0) {
        size_t repeat_count = 0;
        suff_len = guess_suff_len(inputfile, n_reads_to_check, repeat_count);
        if (repeat_count > 1)
          throw runtime_error("failed to identify read name suffix length\n"
                              "verify reads are not single-end\n"
                              "specify read name suffix length directly");
        if (VERBOSE)
          cerr << "[read name suffix length guess: " << suff_len << "]" << endl;
      }
      else if (!check_suff_len(inputfile, suff_len, n_reads_to_check))
        throw runtime_error("wrong read name suffix length [" +
                            to_string(suff_len) + "] in: " + inputfile);
      if (!check_sorted(inputfile, suff_len, n_reads_to_check))
        throw runtime_error("mates not consecutive in: " + inputfile);
    }

    if (VERBOSE)
      cerr << "[input file: " << inputfile << "]" << endl
           << "[mapper: " << input_format << "]" << endl
           << "[configuration: " << (single_end ? "SE" : "PE") << "]" << endl
           << "[output file: " << outfile << "]" << endl
           << "[output format: " << (bam_format ? "B" : "S") << "AM]" << endl
           << "[threads requested: " << n_threads << "]" << endl
           << "[command line: \"" << cmd.str() << "\"]" << endl
           << "[readname suffix length: " << suff_len << "]" << endl;

    format(cmd.str(), n_threads, inputfile, outfile,
           bam_format, input_format, suff_len, max_frag_len);

  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
