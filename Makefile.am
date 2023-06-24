# This file is part of dnmtools
#
# Copyright (C) 2023 Andrew D. Smith
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

EXTRA_DIST = \
             README.md \
	     LICENSE \
             docs/README.md \
             docs/mkdocs.yml \
             docs/content/abismal.md \
             docs/content/allelic.md \
             docs/content/amrfinder.md \
             docs/content/amrtester.md \
             docs/content/bsrate.md \
             docs/content/cleanhp.md \
             docs/content/counts.md \
             docs/content/diff.md \
             docs/content/dmr.md \
             docs/content/entropy.md \
             docs/content/fastlift.md \
             docs/content/format.md \
             docs/content/guessprotocol.md \
             docs/content/hmr.md \
             docs/content/hmr-rep.md \
             docs/content/hypermr.md \
             docs/content/index.md \
             docs/content/lc.md \
             docs/content/levels.md \
             docs/content/liftfilter.md \
             docs/content/merge-bsrate.md \
             docs/content/merge.md \
             docs/content/mlml.md \
             docs/content/multistat.md \
             docs/content/pmd.md \
             docs/content/quickstart.md \
             docs/content/radadjust.md \
             docs/content/radmerge.md \
             docs/content/radmeth.md \
             docs/content/roi.md \
             docs/content/selectsites.md \
             docs/content/states.md \
             docs/content/sym.md \
             docs/content/uniq.md \
             pipeline/config.yaml \
             pipeline/runconfig.yaml \
             pipeline/Snakefile \
             test/test.fa \
             test/test_abismalidx.test \
             test/test_simreads.test \
             test/test_abismal.test \
             test/test_format.test \
             test/test_bsrate.test

ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp src/abismal
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp src/abismal, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/common -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 $(OPENMP_CXXFLAGS) # default has optimization on

TESTS = test/test_abismalidx.test \
	test/test_simreads.test \
	test/test_abismal.test \
	test/test_format.test \
	test/test_bsrate.test

TEST_EXTENSIONS = .test

test_simreads.log: test_abismalidx.log
test_abismal.log: test_simreads.log
test_format.log: test_abismal.log
test_bsrate.log: test_format.log

noinst_LIBRARIES = libdnmtools.a
libdnmtools_a_SOURCES = \
        src/common/BetaBin.cpp \
        src/common/Distro.cpp \
        src/common/EmissionDistribution.cpp \
        src/common/Epiread.cpp \
        src/common/EpireadStats.cpp \
        src/common/LevelsCounter.cpp \
        src/common/MSite.cpp \
        src/common/Smoothing.cpp \
        src/common/ThreeStateHMM.cpp \
        src/common/TwoStateHMM.cpp \
        src/common/TwoStateHMM_PMD.cpp \
        src/common/bsutils.cpp \
        src/common/numerical_utils.cpp

libdnmtools_a_SOURCES += \
        src/common/BetaBin.hpp \
        src/common/Distro.hpp \
        src/common/EmissionDistribution.hpp \
        src/common/Epiread.hpp \
        src/common/EpireadStats.hpp \
        src/common/LevelsCounter.hpp \
        src/common/MSite.hpp \
        src/common/Smoothing.hpp \
        src/common/ThreeStateHMM.hpp \
        src/common/TwoStateHMM.hpp \
        src/common/TwoStateHMM_PMD.hpp \
        src/common/bsutils.hpp \
        src/common/numerical_utils.hpp

LDADD = libdnmtools.a  src/abismal/libabismal.a src/smithlab_cpp/libsmithlab_cpp.a

bin_PROGRAMS = dnmtools

dnmtools_SOURCES = src/dnmtools.cpp

dnmtools_SOURCES += src/analysis/pmd.cpp
dnmtools_SOURCES += src/analysis/methstates.cpp
dnmtools_SOURCES += src/analysis/bsrate.cpp
dnmtools_SOURCES += src/analysis/methentropy.cpp
dnmtools_SOURCES += src/analysis/methcounts.cpp
dnmtools_SOURCES += src/analysis/roimethstat.cpp
dnmtools_SOURCES += src/analysis/multimethstat.cpp
dnmtools_SOURCES += src/analysis/hmr.cpp
dnmtools_SOURCES += src/analysis/hmr-rep.cpp
dnmtools_SOURCES += src/analysis/levels.cpp
dnmtools_SOURCES += src/analysis/hypermr.cpp

dnmtools_SOURCES += src/utils/clean-hairpins.cpp
dnmtools_SOURCES += src/utils/guessprotocol.cpp
dnmtools_SOURCES += src/utils/uniq.cpp
dnmtools_SOURCES += src/utils/merge-bsrate.cpp
dnmtools_SOURCES += src/utils/format-reads.cpp
dnmtools_SOURCES += src/utils/lc-approx.cpp
dnmtools_SOURCES += src/utils/selectsites.cpp
dnmtools_SOURCES += src/utils/symmetric-cpgs.cpp
dnmtools_SOURCES += src/utils/merge-methcounts.cpp
dnmtools_SOURCES += src/utils/lift-filter.cpp
dnmtools_SOURCES += src/utils/fast-liftover.cpp

dnmtools_SOURCES += src/amrfinder/allelicmeth.cpp
dnmtools_SOURCES += src/amrfinder/amrfinder.cpp
dnmtools_SOURCES += src/amrfinder/amrtester.cpp

dnmtools_SOURCES += src/radmeth/dmr.cpp
dnmtools_SOURCES += src/radmeth/methdiff.cpp
dnmtools_SOURCES += src/radmeth/radmeth.cpp
dnmtools_SOURCES += src/radmeth/radmeth-adjust.cpp
dnmtools_SOURCES += src/radmeth/radmeth-merge.cpp

dnmtools_SOURCES += src/mlml/mlml.cpp
