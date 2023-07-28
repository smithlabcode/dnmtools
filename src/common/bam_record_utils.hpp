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

#ifndef BAM_RECORD_UTILS_HPP
#define BAM_RECORD_UTILS_HPP

#include "bam_record.hpp"

#include <string>

int
truncate_overlap(const bam_rec &a, const uint32_t overlap, bam_rec &c);

int
merge_overlap(const bam_rec &a, const bam_rec &b,
              const uint32_t head, bam_rec &c);

int
merge_non_overlap(const bam_rec &a, const bam_rec &b,
                  const uint32_t spacer, bam_rec &c);

int
keep_better_end(const bam_rec &a, const bam_rec &b, bam_rec &c);

size_t correct_cigar(bam_rec &b);

void flip_conversion(bam_rec &aln);

inline bool
is_a_rich(const bam_rec &b) { return bam_aux2A(bam_aux_get(b.b, "CV")) == 'A'; }

inline bool
is_rev(const bam_rec &a) { return bam_is_rev(a.b); }

void
standardize_format(const std::string &input_format, bam_rec &aln);

inline bool
are_mates(const bam_rec &one, const bam_rec &two) {
  return one.b->core.mtid == two.b->core.tid &&
    one.b->core.mpos == two.b->core.pos &&
    bam_is_rev(one.b) != bam_is_rev(two.b);
  // below is a consistency check and should not be necessary
  /* &&
     two->core.mtid == one->core.tid &&
     two->core.mpos == one->core.pos; */
}


inline size_t
get_l_qseq(const bam_rec &b) { return b.b->core.l_qseq; }

inline size_t
get_n_targets(const bam_header &bh) { return bh.h->n_targets; }

inline std::string
get_qname(const bam_rec &b) { return bam_get_qname(b.b); }

inline int32_t
get_tid(const bam_rec &b) { return b.b->core.tid; }

inline hts_pos_t
get_pos(const bam_rec &b) { return b.b->core.pos; }

inline hts_pos_t
get_endpos(const bam_rec &b) { return bam_endpos(b.b); }

inline bool
precedes_by_start(const bam_rec &a, const bam_rec &b) {
  // assumes a.get_tid() <= b.get_tid()
  return get_tid(a) == get_tid(b) && get_pos(a) < get_pos(b);
}

inline bool
precedes_by_end_and_strand(const bam_rec &a, const bam_rec &b) {
  const auto end_a = bam_endpos(a.b);
  const auto end_b = bam_endpos(b.b);
  return end_a < end_b || (end_a == end_b && bam_is_rev(a.b) < bam_is_rev(b.b));
}

inline bool
equivalent_chrom_and_start(const bam_rec &a, const bam_rec &b) {
  return a.b->core.pos == b.b->core.pos && a.b->core.tid == b.b->core.tid;
}

inline bool
equivalent_end_and_strand(const bam_rec &a, const bam_rec &b) {
  return bam_endpos(a.b) == bam_endpos(b.b) && bam_is_rev(a.b) == bam_is_rev(b.b);
}

#endif
