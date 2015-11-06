#
# $Id: fre-nctools.mk,v 1.1.2.1.2.2.4.2 2012/11/05 15:25:06 keo Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build Regridding Executables
# ------------------------------------------------------------------------------
# afy    Ver   1.00  Initial version (Makefile, ver 17.0.4.2)       June 10
# afy    Ver   1.01  Add rules to build MPI-based executable        June 10
# afy    Ver   1.02  Simplified according to fre-nctools standards  June 10
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
#
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
#
include ./env.$(SITE)

CXX      := icpc
CXXFLAGS := -O3 -g -traceback

EZBAR=./extern/ezProgressBar-2.1.0
EZNC=./extern/eznc-0.2.1/src
EZODO=./extern/ezOdometer-0.1.1
EZOPT=./extern/ezOptionParser-0.1.2
EZSLICE=./extern/ezSlice-0.0.0
EZSTR=./extern/ezStringUtil-0.1.0
EZTEST=./extern/ezTest-0.0.0

INCLUDES := -I${NETCDF_HOME}/include -I. -I$(EZBAR) -I$(EZNC) -I$(EZODO) -I$(EZOPT) -I$(EZSLICE) -I$(EZSTR)
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz  -lirc $(CLIBS2) $(STATIC)

TARGETS  := ncx

SOURCES  := ncx.cpp

OBJECTS  := $(SOURCES:cpp=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

ncx: $(OBJECTS)
	$(CXX) -o $@ $^ $(CLIBS)

ncx.o: ncx.cpp ncx.hpp $(HEADERS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
