/*    Copyright (C) 2013 University of Southern California and
 *                       Egor Dolzhenko
 *                       Andrew D Smith
 *
 *    Authors: Andrew D. Smith and Egor Dolzhenko
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */
 
#ifndef COMBINE_PVALS_HPP_
#define COMBINE_PVALS_HPP_

#include <vector>

#include "locus.hpp"
#include "bin_for_distance.hpp"

void combine_pvals(std::vector<LocusIterator> &loci_iterators, 
                    const BinForDistance &bin_for_distance);

#endif //COMBINE_PVALS_HPP_
