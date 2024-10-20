Installation
============

## Installation via conda

*Note (9/26/2023)* Although conda currently works, it seems that using
`mamba` is more reliable. If you have `conda` installed, you can
install mamba very easily by following instructions here:

https://anaconda.org/conda-forge/mamba

In the instructions below, replacing `conda` with `mamba` should work
the same (and, in some cases, more reliably), but with mamba.

Right now conda is probably the easiest way to install
dnmtools. Dnmtools is among the bioconda recipes. If you know how to
use conda and are setup to use bioconda, then you might simply be able
to do:
```console
$ conda install dnmtools
```

If you have problems with this, then you can try:
```console
$ conda install -c conda-forge -c bioconda --strict-channel-priority dnmtools
```

The dependencies for dnmtools might themselves have dependencies that
are in the `conda-forge` channel. Using `--strict-channel-priority` is
recommended by bioconda.

* Conda: you can find installation instructions
  [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

* Bioconda: additional helpful setup instructions can be found
  [here](https://bioconda.github.io).

If you encounter further problems, try creating a new environment for
dnmtools within conda:
```console
$ conda create --name dnmtools_env
$ conda activate dnmtools_env
$ conda install -c conda-forge -c bioconda --strict-channel-priority dnmtools
```

Above, you would have a new environment called `dnmtools_env` that
would need to be activated when you want to use dnmtools.

## Building dnmtools from source release

### Required libraries

* A recent compiler. Most users will be building and installing this
  software with GCC. We require a compiler that fully supports C++17,
  so we recommend using at least GCC 9 (released in 2019). There are
  still many systems that install a very old version of GCC by
  default, so if you have problems with building this software, that
  might be the first thing to check. The clang LLVM compiler can
  also be used with a recent enough version.
* The GNU Scientific Library. It can be installed using apt on Linux
  (Ubuntu, Debian), using brew on macOS, or from source available
  [here](http://www.gnu.org/software/gsl).
* The HTSlib library. This can be installed through brew on macOS,
  through apt on Linux (Ubuntu, Debian), or from source downloadable
  [here](https://github.com/samtools/htslib).
* The Zlib compression library: Most likely you already have this
  installed on your system. If not, it can be installed using apt on
  Linux (Ubuntu, Debian) through the package `zlib1g-dev`. On macOS,
  Zlib can be installed with brew.

### Configuration

* Download [dnmtools-1.4.4.tar.gz](https://github.com/smithlabcode/dnmtools/releases/download/v1.4.4/dnmtools-1.4.4.tar.gz).
* Unpack the archive:
```console
$ tar -zxvf dnmtools-1.4.4.tar.gz
```
* Move into the dnmtools directory and create a build directory:
```console
$ cd dnmtools-1.4.4
$ mkdir build && cd build
```
* Run the configuration script:
```console
$ ../configure
```
If you do not want to install dnmtools system-wide, or if you do
not have admin privileges, specify a prefix directory:
```console
$ ../configure --prefix=/some/reasonable/place
```
If you installed HTSlib yourself in some non-standard directory,
you must specify the location like this:
```console
$ ../configure CPPFLAGS='-I /path/to/htslib/headers' \
               LDFLAGS='-L/path/to/htslib/lib'
```
Depending on how you obtained HTSlib, the headers may not be
in a directory at the same depth as the library file.

### Building and installing the tools

If you are still in the `build` directory, run `make` to compile the
tools, and then `make install` to install them:
```console
$ make
$ make install
```
If your HTSlib (or some other library) is not installed system-wide,
then you might need to udpate your library path:
```console
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/htslib/lib
```

### Testing the program

To test if everything was successful, simply run `dnmtools` without
any arguments and you should see the list of available commands:
```console
dnmtools
```
There is a test suite for `dnmtools` and these test can be performed
as follows:
```console
make check
```
This must be done from the build directory. Note that the tests
performed with `make check` are mostly regression tests that cover
prior issues rather than coverage tests to test all the functionality
of `dnmtools`.

## Using a clone of the repo

We strongly recommend using dnmtools through the latest stable release
under the releases section on GitHub. Developers who wish to work on
the latest commits, which are unstable, can start the process using
the `autogen.sh` script and then using the usual autotools approach.
