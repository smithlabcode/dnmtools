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

find_package(HTSLIB REQUIRED)

add_library(radmeth OBJECT
  dmr.cpp
  methdiff.cpp
  radmeth-adjust.cpp
  radmeth-merge.cpp
  radmeth.cpp
  radmeth_design.cpp
  radmeth_nano.cpp
  radmeth_optimize_gamma.cpp
  radmeth_optimize_series.cpp
  radmeth_utils.cpp
)

target_include_directories(radmeth PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_BINARY_DIR}
)
target_link_libraries(radmeth PUBLIC
  dnmtools_objs
  bamxx
  smithlab_cpp
  HTSLIB::HTSLIB
)
