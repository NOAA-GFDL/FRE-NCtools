# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
# ------------------------------------------------------------------------------

# Bring in settings for the site
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

CXX      := icpc
CXXFLAGS := -O3 -g -traceback $(CFLAGS_SITE)

EZBAR=$(SRCDIR)/extern/ezProgressBar-2.1.0
EZNC=$(SRCDIR)/extern/eznc-0.2.1/src
EZODO=$(SRCDIR)/extern/ezOdometer-0.1.1
EZOPT=$(SRCDIR)/extern/ezOptionParser-0.1.2
EZSLICE=$(SRCDIR)/extern/ezSlice-0.0.0
EZSTR=$(SRCDIR)/extern/ezStringUtil-0.1.0
EZTEST=$(SRCDIR)/extern/ezTest-0.0.0

INCLUDES := -I${NETCDF_HOME}/include -I$(SRCDIR) -I$(EZBAR) -I$(EZNC) -I$(EZODO) -I$(EZOPT) -I$(EZSLICE) -I$(EZSTR)
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz  -lirc $(CLIBS_SITE) $(STATIC)

TARGETS  := ncx

SOURCES  := ncx.cpp

OBJECTS  := $(SOURCES:cpp=o)

all: $(TARGETS)

ncx: $(OBJECTS)
	$(CXX) -o $@ $^ $(CLIBS)

ncx.o: ncx.cpp ncx.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

install: all
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)

docs:

install-docs: docs

clean:
	-rm -f *.o $(TARGETS)
