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
#include "bam_header.hpp"
#include "bam_infile.hpp"
#include "bam_outfile.hpp"
#include "bam_tpool.hpp"

#include <string>

bool are_mates(const bam_rec &one, const bam_rec &two);

int truncate_overlap(const bam_rec &a, const uint32_t overlap, bam_rec &c);

int merge_overlap(const bam_rec &a, const bam_rec &b,
                  const uint32_t head, bam_rec &c);

int merge_non_overlap(const bam_rec &a, const bam_rec &b,
                      const uint32_t spacer, bam_rec &c);

int keep_better_end(const bam_rec &a, const bam_rec &b, bam_rec &c);

size_t correct_cigar(bam_rec &b);

void revcomp_seq_by_byte(bam_rec &aln);

bool is_a_rich(const bam_rec &b);

void flip_conversion(bam_rec &aln);

void standardize_format(const std::string &input_format, bam_rec &aln);

bool is_rev(const bam_rec &aln);

#endif
