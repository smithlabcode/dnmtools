# This file is part of dnmtools
#
# Copyright (C) 2022 Andrew D. Smith
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

EXTRA_DIST = README.md
ACLOCAL_AMFLAGS = -I m4

SUBDIRS := src/smithlab_cpp
install installdirs: SUBDIRS := $(filter-out src/smithlab_cpp, $(SUBDIRS))
AM_CPPFLAGS = -I $(top_srcdir)/src/common -I $(top_srcdir)/src/smithlab_cpp

CXXFLAGS = -O3 # default has optimization on

noinst_LIBRARIES = libdnmtools.a
libdnmtools_a_SOURCES = \
        src/common/BetaBin.cpp \
        src/common/Distro.cpp \
        src/common/EmissionDistribution.cpp \
        src/common/Epiread.cpp \
        src/common/EpireadStats.cpp \
        src/common/LevelsCounter.cpp \
        src/common/MethpipeSite.cpp \
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
        src/common/MethpipeFiles.hpp \
        src/common/MethpipeSite.hpp \
        src/common/Smoothing.hpp \
        src/common/ThreeStateHMM.hpp \
        src/common/TwoStateHMM.hpp \
        src/common/TwoStateHMM_PMD.hpp \
        src/common/bsutils.hpp \
        src/common/numerical_utils.hpp

LDADD = libdnmtools.a src/smithlab_cpp/libsmithlab_cpp.a

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
dnmtools_SOURCES += src/analysis/hmr_rep.cpp
dnmtools_SOURCES += src/analysis/levels.cpp
dnmtools_SOURCES += src/analysis/hypermr.cpp

dnmtools_SOURCES += src/utils/clean-hairpins.cpp
dnmtools_SOURCES += src/utils/guessprotocol.cpp
dnmtools_SOURCES += src/utils/duplicate-remover.cpp
dnmtools_SOURCES += src/utils/merge-bsrate.cpp
dnmtools_SOURCES += src/utils/format_reads.cpp
dnmtools_SOURCES += src/utils/lc_approx.cpp
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
