[![GitHub Downloads](https://img.shields.io/github/downloads/smithlabcode/dnmtools/total?style=social)](https://github.com/smithlabcode/dnmtools/releases)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/version.svg)](https://anaconda.org/bioconda/dnmtools)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/platforms.svg)](https://anaconda.org/bioconda/dnmtools)
[![Install with Conda](https://anaconda.org/bioconda/dnmtools/badges/downloads.svg)](https://anaconda.org/bioconda/dnmtools)
[![Documentation Status](https://readthedocs.org/projects/dnmtools/badge/?version=latest)](https://dnmtools.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

DNMTools is a set of tools for analyzing DNA methylation data from
high-throughput sequencing experiments, especially whole genome
bisulfite sequencing (WGBS), but also reduced representation bisulfite
sequencing (RRBS). These tools focus on overcoming the computing
challenges imposed by the scale of genome-wide DNA methylation data,
which is usually the early parts of data analysis.

## Installing release 1.3.0

The documentation for DNMTools can be found
[here](https://dnmtools.readthedocs.io). But if you want to install
from source and you are reading this on GitHub or in a source tree you
unpacked, then keep reading. And if you are in a terminal, sorry for
all the formatting.

### Required libraries

* A recent compiler. Most users will be building and installing this
  software with GCC. We require a compiler that fully supports C++11,
  so we recommend using at least GCC 5.8. There are still many systems
  that install a very old version of GCC by default, so if you have
  problems with building this software, that might be the first thing
  to check.
* The GNU Scientific Library. It can be installed using apt on Linux
  (Ubuntu, Debian), using brew on macOS, or from source available
  [here](http://www.gnu.org/software/gsl).
* The Zlib compression library. Most likely you already have this
  installed on your system. If not, it can be installed using apt on
  Linux (Ubuntu, Debian) through the package `zlib1g-dev`. On macOS,
  Zlib can be installed with brew (I have not needed to install this
  separately on any macOS).
* The HTSlib library. This can be installed through brew on macOS,
  through apt on Linux (Ubuntu, Debian), or from source downloadable
  [here](https://github.com/samtools/htslib).

All the above can also be installed using conda. If you use conda for
these dependencies, even if you are building dnmtools from the source
repo, it is easiest if all dependencies are available through conda.

### Configuration

* Download [dnmtools-1.3.0.tar.gz](https://github.com/smithlabcode/dnmtools/releases/download/v1.3.0/dnmtools-1.3.0.tar.gz).
* Unpack the archive:
```console
tar -zxvf dnmtools-1.3.0.tar.gz
```
* Move into the dnmtools directory and create a build directory:
```console
cd dnmtools-1.3.0 && mkdir build && cd build
```
* Run the configuration script:
```console
../configure
```
If you do not want to install DNMTools system-wide, or if you do
not have admin privileges, specify a prefix directory:
```console
../configure --prefix=/some/reasonable/place
```
If you installed HTSlib yourself in some non-standard directory,
you must specify the location like this:
```console
../configure CPPFLAGS='-I /path/to/htslib/headers' \
             LDFLAGS='-L/path/to/htslib/lib'
```
Depending on how you obtained HTSlib, the headers may not be
in a directory at the same depth as the library file.

### Building and installing the tools

If you are still in the `build` directory, run `make` to compile the
tools, and then `make install` to install them:
```console
make && make install
```
If your HTSlib (or some other library) is not installed system-wide,
then you might need to udpate your library path:
```console
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/htslib/lib
```

### Testing the program

To test if everything was successful, simply run `dnmtools` without
any arguments and you should see the list of available commands:
```console
dnmtools
```

### Using a clone of the repo

Not recommended, but if you want to do it this way, we assume you know
what you are doing. We strongly recommend using DNMTools through the
latest stable release under the releases section on GitHub. Developers
who wish to work on the latest commits, which are unstable, can
compile the source using a `Makefile` left in the root of the source
tree. If HTSLib and other libraries are available system-wide,
compile by running:
```console
make
```
This functionality will probably be removed soon, and if you want to
build the code this way, you should know what you are doing any be
able to make it work yourself.

## Usage

Read the [documentation](https://dnmtools.readthedocs.io) for usage of
individual tools within DNMTools.

## Contacts and bug reports

Andrew D. Smith
andrewds@usc.edu

Guilherme de Sena Brandine
desenabr@usc.edu

## Copyright and License Information

Copyright (C) 2022-2023
Andrew D. Smith and Guilherme de Sena Brandine

Authors of DNMTools: Andrew D. Smith and Guilherme de Sena Brandine

Essential contributors: Ben Decato, Meng Zhou, Liz Ji, Terence Li,
Jenny Qu, Qiang Song and Fang Fang

This is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
