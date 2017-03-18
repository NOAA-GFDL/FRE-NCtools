# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to build make_solo_mosaic
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011, 2017
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst, Aleksey Yakovlev and
# Seth Underwood
# ------------------------------------------------------------------------------

# MPICC and CC default settings
MPICC    := mpicc
CC       := icc

# Get the site specific settings, CC and MPICC may be overridden
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

CFLAGS   := -O3 -traceback $(CFLAGS_SITE)
CFLAGS_O2:= -O2 -traceback $(CFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include -I$(SRCDIR)/ -I$(SRCDIR)/../shared -I$(SRCDIR)/../../shared/mosaic
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := make_solo_mosaic

SOURCES  := make_solo_mosaic.c get_contact.c
SOURCES  += mpp.c mpp_domain.c mpp_io.c tool_util.c
SOURCES  += create_xgrid.c interp.c mosaic_util.c read_mosaic.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = ../shared/mpp.h  ../shared/mpp_domain.h  ../shared/mpp_io.h ../shared/tool_util.h   \
	  ../../shared/mosaic/constant.h ../../shared/mosaic/create_xgrid.h  \
	  ../../shared/mosaic/interp.h  ../../shared/mosaic/mosaic_util.h  \
	  ./get_contact.h ../../shared/mosaic/read_mosaic.h

all: $(TARGETS)

make_solo_mosaic: $(OBJECTS)
	$(CC) -o $@ $^ $(CLIBS)

make_solo_mosaic.o: make_solo_mosaic.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

get_contact.o: get_contact.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

mosaic_util.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

read_mosaic.o: ../../shared/mosaic/read_mosaic.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

interp.o: ../../shared/mosaic/interp.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

mpp_io.o: ../shared/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp_domain.o: ../shared/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

mpp.o: ../shared/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

tool_util.o: ../shared/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $< 

create_xgrid.o: ../../shared/mosaic/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS_O2) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

install: all
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)

docs:

install-docs: docs

clean:
	-rm -f *.o $(TARGETS)
