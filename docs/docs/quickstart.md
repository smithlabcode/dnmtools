Installation
============

## Building dnmtools from source release

### Required libraries

* A recent compiler: most users will be building and installing this
  software with GCC. We require a compiler that fully supports C++11,
  so we recommend using at least GCC 5.8. There are still many systems
  that install a very old version of GCC by default, so if you have
  problems with building this software, that might be the first thing
  to check.
* The GNU Scientific Library: this has always been required. It can be
  installed using `apt` on Linux (Ubuntu, Debian), using `brew` on
  macOS, or from source available
  [here](http://www.gnu.org/software/gsl).
* The Zlib compression library. Most likely you already have this
  installed on your system. If not, it can be installed using `apt` on
  Linux (Ubuntu, Debian) through the package `zlib1g-dev`. On macOS,
  Zlib can be installed with `brew`.
* The HTSlib library, which can be installed through `brew` on macOS,
  through `apt` on Linux (Ubuntu, Debian), or from source downloadable
  [here](https://github.com/samtools/htslib).

### Configuration

1. Download [dnmtools-1.0.0.tar.gz](https://github.com/smithlabcode/dnmtools/releases/download/v1.0.0/dnmtools-1.0.0.tar.gz).
2. Unpack the archive:

```
$ tar -zxvf dnmtools-1.0.0.tar.gz
```
3. Move into the dnmtools directory and create a build directory:
```
$ cd dnmtools-1.0.0
$ mkdir build && cd build
```
4. Run the configuration script:
```
$ ../configure
```

If you do not want to install `dnmtools` system-wide, or if you do
not have admin privileges, specify a prefix directory:
```
$ ../configure --prefix=/some/reasonable/place
```
If you installed HTSlib yourself in some non-standard directory,
you must specify the location like this:
```
$ ../configure CPPFLAGS='-I /path/to/htslib/headers' \
               LDFLAGS='-L/path/to/htslib/lib'
```

### Building and installing the tools

If you are still in the `build` directory, run `make` to compile the
tools, and then `make install` to install them:

```
$ make
$ make install
```

If your HTSlib (or some other library) is not installed system-wide,
then you might need to udpate your library path:

```
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/htslib/lib
```

### Testing the program

To test if everything was successful, simply run `dnmtools` without
any arguments and you should see the list of available commands
```
$ dnmtools
```


## Using a clone of the repo

Not recommended, but if you want to do it this way, we assume you know
what you are doing. We strongly recommend using dnmtools through the
latest stable release under the releases section on GitHub. Developers
who wish to work on the latest commits, which are unstable, can
compile the source using the `Makefile` available in the root of the
source tree. If HTSLib and other libraries are available system-wide,
compile by running:

```
$ make
```
