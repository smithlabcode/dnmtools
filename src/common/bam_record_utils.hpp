#ifndef BAM_RECORD_UTILS_HPP
#define BAM_RECORD_UTILS_HPP

#include "bam_record.hpp"

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
