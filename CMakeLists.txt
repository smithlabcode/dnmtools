# This file is part of dnmtools
#
# Copyright (C) 2025 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This software is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
# more details.

# to find the version of cmake do
# $ cmake --version
cmake_minimum_required(VERSION 3.28)
project(
  dnmtools
  VERSION 1.5.1
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
add_compile_options(
  -Wall
  -Wextra
  -Wpedantic
  -Werror
  -Wfatal-errors
)

add_subdirectory(src)
