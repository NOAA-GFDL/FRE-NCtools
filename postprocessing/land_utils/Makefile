# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build land combine/decompress tools
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
# ------------------------------------------------------------------------------

# Default compilers
CC       := icc
FC       := ifort

# Bring in setting for the build
# Compilers may be overridden
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

CFLAGS   := -O3 -g -traceback $(CFLAGS_SITE)
CLIBS    := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 64 -traceback $(FFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz $(LIBS_SITE) $(STATIC)

TARGETS  := scatter-ncc combine-ncc decompress-ncc is-compressed

SOURCES  := nfu.F90 nfu_compress.F90

OBJECTS  := $(SOURCES:F90=o)

all: $(TARGETS)

scatter-ncc: scatter-ncc.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

scatter-ncc.o: scatter-ncc.F90 $(OBJECTS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

combine-ncc: combine-ncc.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

combine-ncc.o: combine-ncc.F90 $(OBJECTS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

decompress-ncc: decompress-ncc.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

decompress-ncc.o: decompress-ncc.F90 $(OBJECTS)
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

nfu_compress.o: nfu_compress.F90 nfu.o
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

nfu.o: nfu.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

is-compressed: is-compressed.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(CLIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

install: $(TARGETS)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)

docs:

install-docs: docs

clean:
	-rm -f *.o *.mod $(TARGETS)
