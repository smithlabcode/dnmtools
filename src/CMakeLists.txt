# Copyright (C) 2025 Andrew D Smith
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

if(NOT TARGET smithlab_cpp)
  add_subdirectory(bamxx)
endif()
if(NOT TARGET smithlab_cpp)
  add_subdirectory(smithlab_cpp)
endif()
add_subdirectory(common)
add_subdirectory(radmeth)
add_subdirectory(utils)
add_subdirectory(analysis)
add_subdirectory(amrfinder)
add_subdirectory(abismal)
add_subdirectory(mlml)

add_executable(dnmtools dnmtools.cpp)
# ADS: below, for config.h
target_include_directories(dnmtools PUBLIC ${CMAKE_BINARY_DIR})
target_link_libraries(dnmtools PUBLIC
  dnmtools_objs
  smithlab_cpp
  bamxx
  abismal_objs
  radmeth
  format-reads
  uniq
  bsrate
  methcounts
  nanopore
  symmetric-cpgs
  levels
  hmr
  hmr-rep
  hypermr
  methentropy
  pmd
  roimethstat
  autocorr
  cpgbins
  multimethstat
  mlml
  methstates
  allelicmeth
  amrfinder
  amrtester
  fast-liftover
  metagene
  lift-filter
  clean-hairpins
  guessprotocol
  merge-bsrate
  merge-methcounts
  covered
  recovered
  xcounts
  unxcounts
  selectsites
  kmersites
)
