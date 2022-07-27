#    This file is part of dnmtools
#
#    Copyright (C) 2010-2022 University of Southern California and
#                            Andrew D. Smith
#
#    Authors: Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

DNMTOOLS_ROOT = $(shell pwd)

all:
	@$(MAKE) -C src DNMTOOLS_ROOT=$(DNMTOOLS_ROOT)

install:
	@$(MAKE) -C src DNMTOOLS_ROOT=$(DNMTOOLS_ROOT) install

clean:
	@$(MAKE) -C src DNMTOOLS_ROOT=$(DNMTOOLS_ROOT) clean
.PHONY: clean

distclean: clean
	@rm -rf $(DNMTOOLS_ROOT)/bin
.PHONY: distclean
