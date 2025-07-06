[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/dnmtools/total?style=social)](https://github.com/smithlabcode/dnmtools/releases)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/version.svg)](https://anaconda.org/bioconda/dnmtools)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/platforms.svg)](https://anaconda.org/bioconda/dnmtools)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/downloads.svg)](https://anaconda.org/bioconda/dnmtools)
[![Documentation Status](https://readthedocs.org/projects/dnmtools/badge/?version=latest)](https://dnmtools.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

DNMTools is a set of tools for analyzing DNA methylation data from
high-throughput sequencing experiments, especially whole genome bisulfite
sequencing (WGBS), but also reduced representation bisulfite sequencing
(RRBS). These tools focus on overcoming the computing challenges imposed by
the scale of genome-wide DNA methylation data, which is usually the early
parts of data analysis.

**Nanopore** As of v1.5.0, DNMTools has funcionality to start analysis with a
BAM file from Nanopore sequencing with 5mC and 5hmC calls at CpG sites.

## Usage

The documentation for DNMTools can be found
[here](https://dnmtools.readthedocs.io).

## Installation

- **Linux**
  [binary](https://github.com/smithlabcode/dnmtools/releases/download/v1.5.0/dnmtools-1.5.0-Linux.tar.gz).
  Should work on any Linux distribution since roughly 2017.

- **Mac**
  [binary](https://github.com/smithlabcode/dnmtools/releases/download/v1.5.0/dnmtools-1.5.0-macOS.tar.gz).
  Should work on any Mac hardware and macOS-13 (Ventura) or newer.

- **Conda**
  ```console
  conda install -c bioconda dnmtools
  ```

- **Source**
  [dnmtools-1.5.0.tar.gz](https://github.com/smithlabcode/dnmtools/releases/download/v1.5.0/dnmtools-1.5.0.tar.gz). Dependencies:
  [GSL](http://www.gnu.org/software/gsl),
  [HTSlib](https://github.com/samtools/htslib),
  [libdeflate](https://github.com/ebiggers/libdeflate) and
  [ZLib](https://github.com/madler/zlib). Installing HTSlib as a package
  should also give you ZLib and libdeflate.  System-specific details below.

  Build DNMTools like this:
  ```console
  tar -xf dnmtools-1.5.0.tar.gz
  cd dnmtools-1.5.0
  ./configure --prefix=$HOME
  make
  make install
  ```

  To get dependencies and a compiler on (these might with OS/package updates):

  Ubuntu/Debian
  ```console
  apt-get install build-essential htslib-dev libgsl-dev
  ```

  RedHat/Fedora
  ```console
  dnf install @c-development @development-tools htslib-devel gsl-devel awk
  ```

  Homebrew (see notes below)
  ```console
  brew install gcc htslib gsl
  ```

  Conda (see notes below)
  ```console
  conda create -n build-env -c conda-forge -c bioconda \
      gcc gxx make autoconf automake htslib gsl zlib binutils && \
  conda activate build-env
  ```

  Notes: If you use only Homebrew or only Conda to setup your environment, you
  could need additional dependencies, and some of what I listed you might
  already have. You might need to set additional environment variables or run
  configure differently. For example with Homebrew:
  ```console
  ./configure CPPFLAGS="-I$(brew --prefix)/include" LDFLAGS="-L$(brew --prefix)/lib"
  ```

## Contact

Andrew D. Smith
andrewds@usc.edu

## Copyright and License Information

Copyright (C) 2022-2025

Andrew D. Smith and Guilherme de Sena Brandine

Authors of DNMTools: Andrew D. Smith and Guilherme de Sena Brandine

Essential contributors: Ben Decato, Meng Zhou, Liz Ji, Terence Li, Jenny Qu,
Qiang Song, Fang Fang and Masaru Nakajima

This is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This software is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details.
