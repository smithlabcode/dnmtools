# This file is part of dnmtools
#
# Copyright (C) 2025 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.28)
project(
  dnmtools
  VERSION 1.4.4
  DESCRIPTION
  "Tools for analyzing DNA methylation data"
  HOMEPAGE_URL https://github.com/smithlabcode/dnmtools
  LANGUAGES CXX)

# Set language version used
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_EXTENSIONS off)  # prevents std=gnu++17
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

include(CheckIncludeFileCXX)
include(CheckFunctionExists)
include(CheckCXXCompilerFlag)

include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

configure_file(data/config.h.in config.h)

# Collect any global linker options as needed and then apply them
# individually to targets
set(GLOBAL_COMPILE_OPTIONS "")
set(GLOBAL_LINKER_OPTIONS "")

if(USE_STATIC_LIBS)
  # This needs to come before finding any libraries so that the static
  # versions are identified
  message(STATUS "Enabling static linkage for all non-system libraries")
  message(STATUS "Configuring to clone ZLib")
  include(ExternalProject)
  set(ZLIB_CMAKE_ARGS
    -DZLIB_BUILD_EXAMPLES=off
    -DSKIP_INSTALL_FILES=on
    -DCMAKE_POSITION_INDEPENDENT_CODE=on
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/src/zlib
    -DCMAKE_BUILD_TYPE=Release
  )
  if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
    list(APPEND ZLIB_CMAKE_ARGS
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
    )
  endif()
  ExternalProject_Add(
    ZLIB
    GIT_REPOSITORY https://github.com/madler/zlib.git
    GIT_TAG master
    CMAKE_ARGS ${ZLIB_CMAKE_ARGS}
  )
  # Include the built zlib headers and link against the built zlib library
  set(ZLIB_INCLUDE_DIR "${PROJECT_BINARY_DIR}/src/zlib/include")
  set(ZLIB_LIBRARY "${PROJECT_BINARY_DIR}/src/zlib/lib/libz.a")

  # Create the CMake target for the built zlib
  add_library(ZLIB_IMPORTED INTERFACE)
  set_target_properties(ZLIB_IMPORTED PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${ZLIB_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES ${ZLIB_LIBRARY}
  )

  # This alias means we don't need to worry where zlib came from
  add_library(ZLIB::ZLIB ALIAS ZLIB_IMPORTED)

  ExternalProject_Add(
    HTSLIB
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG master
    CONFIGURE_COMMAND ""
    # "autoreconf -i <SOURCE_DIR>; <SOURCE_DIR>/configure"
    BUILD_COMMAND make -C <SOURCE_DIR> lib-static -j
    INSTALL_COMMAND ""
    # make install
    TEST_COMMAND ""
    BUILD_BYPRODUCTS libhts.a
  )
  # ExternalProject_Add_Step(HTSLIB
  #   bootstrap
  #   COMMAND autoreconf -i <SOURCE_DIR>
  #   DEPENDEES download
  #   DEPENDERS configure
  # )
  ExternalProject_Get_Property(HTSLIB SOURCE_DIR)
  message(STATUS "SOURCE_DIR: ${SOURCE_DIR}")

  # Include the built zlib headers and link against the built zlib library
  set(HTSLIB_INCLUDE_DIR "${SOURCE_DIR}/src/htslib/htslib")
  set(HTSLIB_LIBRARY "${SOURCE_DIR}/htslib.a")

  # Create the CMake target for the built zlib
  add_library(HTSLIB_IMPORTED INTERFACE)
  set_target_properties(HTSLIB_IMPORTED PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${HTSLIB_INCLUDE_DIR}
    INTERFACE_LINK_LIBRARIES ${HTSLIB_LIBRARY}
  )
  set_target_properties(HTSLIB_IMPORTED PROPERTIES IMPORTED_LOCATION ${HTSLIB_LIBRARY})

  # This alias means we don't need to worry where htslib came from
  add_library(HTSLIB::HTSLIB ALIAS HTSLIB_IMPORTED)

  # Set static for the linker so the compiler's libraries will be static
  ## ADS: using this instead of forcing -static for everything avoids the
  ## static linkage that Aiso warns against, but also means it's not 100%
  ## static linked ADS: can't do this if the compiler is AppleClang because
  ## they don't have the libc++.a and the libgcc wouldn't make sense anyway.
  if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")  # macOS
    list(APPEND GLOBAL_LINKER_OPTIONS -static-libgcc -static-libstdc++)
  endif()

endif()

if(ENABLE_LTO)
  # Turn on LTO if we are building for distribution
  include(CheckIPOSupported)
  check_ipo_supported(RESULT result OUTPUT output)
  if(result)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(FATAL_ERROR "IPO is not supported: ${output}")
  endif()
endif()

if(STATIC_ANALYSIS)
  include(cmake/static_analysis.cmake)
endif()

# ADS: set the most stringent warnings we can
list(APPEND GLOBAL_COMPILE_OPTIONS
  -Wall -Wextra -Wpedantic -Werror -Wfatal-errors
)

if(STRIP_PATHS_FROM_BINARIES)
  # This is set if we have configured to distribute
  set(PREFIX_MAP_ARG -ffile-prefix-map=)
  list(TRANSFORM STRIP_SUB_LIST PREPEND ${PREFIX_MAP_ARG})
  list(APPEND GLOBAL_COMPILE_OPTIONS ${STRIP_SUB_LIST})
endif()

message(STATUS "Finished global compile and linker configuration")

add_subdirectory(src)
