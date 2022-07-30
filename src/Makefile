# Makefile from DNMTools software
#
# Copyright (C) 2022 University of Southern California and
#                    Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# check if recursive clone was done
SMITHLAB_CPP=$(abspath $(dir $(MAKEFILE_LIST)))/smithlab_cpp
ABISMAL=$(abspath $(dir $(MAKEFILE_LIST)))/abismal

ifeq (,$(wildcard $(SMITHLAB_CPP)/Makefile))
$(error src/smithlab_cpp does not have a Makefile. \
Did you use --recursive when running git clone?)
endif

ifeq (,$(wildcard $(ABISMAL)/Makefile))
$(error src/abismal does not have a Makefile. \
Did you use --recursive when running git clone?)
endif

ifndef SRC_ROOT
SRC_ROOT=..
endif

CXX = g++
CXXFLAGS = -Wall -std=c++11 -fopenmp
OPTFLAGS = -O3
DEBUGFLAGS = -g

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

LDLIBS = -lgsl -lgslcblas -lz -lhts
SMITHLAB_CPP_LIB = $(addprefix $(SMITHLAB_CPP)/, libsmithlab_cpp.a)
ABISMAL_LIB = $(addprefix $(ABISMAL)/src/, libabismal.a)
LDFLAGS = $(SMITHLAB_CPP_LIB) $(ABISMAL_LIB)

COMMON_DIR = common
INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_DIR)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))
override CPPFLAGS += $(INCLUDEARGS)

PROGS = dnmtools
OBJS = amrfinder/allelicmeth.o amrfinder/amrfinder.o amrfinder/amrtester.o \
analysis/bsrate.o analysis/hmr.o analysis/hmr-rep.o analysis/hypermr.o \
analysis/levels.o analysis/methcounts.o analysis/methentropy.o \
analysis/methstates.o analysis/multimethstat.o analysis/pmd.o \
analysis/roimethstat.o mlml/mlml.o radmeth/dmr.o radmeth/methdiff.o \
radmeth/radmeth-adjust.o radmeth/radmeth-merge.o radmeth/radmeth.o \
utils/clean-hairpins.o utils/duplicate-remover.o utils/fast-liftover.o \
utils/format-reads.o utils/guessprotocol.o utils/lc-approx.o \
utils/lift-filter.o utils/merge-bsrate.o utils/merge-methcounts.o \
utils/symmetric-cpgs.o utils/selectsites.o abismal/src/abismal.o \
abismal/src/abismalidx.o abismal/src/simreads.o abismal/src/AbismalIndex.o

COMMON_OBJS = $(addprefix $(COMMON_DIR)/, \
BetaBin.o bsutils.o Distro.o EmissionDistribution.o Epiread.o \
EpireadStats.o LevelsCounter.o MethpipeSite.o numerical_utils.o \
Smoothing.o ThreeStateHMM.o TwoStateHMM.o TwoStateHMM_PMD.o)

all: $(PROGS)

%.o : %.cpp %.hpp

$(SMITHLAB_CPP_LIB):
	@$(MAKE) -C $(SMITHLAB_CPP) HAVE_HTSLIB=1

$(ABISMAL_LIB):
	@$(MAKE) -C $(ABISMAL) HAVE_HTSLIB=1

$(PROGS): dnmtools.cpp $(LDFLAGS) $(OBJS) $(COMMON_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJS) $(COMMON_OBJS) \
	$(LDFLAGS) $(CPPFLAGS) $(LDLIBS)

install:
	@mkdir -p $(SRC_ROOT)/bin
	@install -m 755 $(PROGS) $(SRC_ROOT)/bin

clean:
	$(MAKE) -C $(SMITHLAB_CPP) clean
	$(MAKE) -C $(ABISMAL) clean
	@-rm -f $(PROGS) */*.o *.o *.so */*.a *.a *~
