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

# Packages
if(USE_LIBDEFLATE)
  find_package(LIBDEFLATE REQUIRED)
endif()
if(BUILD_NANOPORE)
  find_package(HTSLIB 1.20 REQUIRED)
  add_compile_definitions(BUILD_NANOPORE)
else()
  find_package(HTSLIB REQUIRED)
endif()
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# Subdirs
if(NOT TARGET smithlab_cpp)
  add_subdirectory(smithlab_cpp)
endif()
if(NOT TARGET bamxx)
  add_subdirectory(bamxx)
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
  abismal_objs
  dnmtools_analysis_objs
  dnmtools_utils_objs
  dnmtools_radmeth_objs
  dnmtools_mlml_objs
  dnmtools_amrfinder_objs
)
