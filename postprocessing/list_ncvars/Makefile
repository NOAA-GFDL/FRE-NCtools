# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build list_ncvars
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
# ------------------------------------------------------------------------------

# Default compilers
FC       := ifort

# Bring in the site build settings
# Compiler may be overridden
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -r8 -traceback $(FFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include
LIBS     :=  -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz $(LIBS_SITE) $(STATIC)

TARGETS  := list_ncvars.exe
TARGETS_SHARE := list_ncvars.csh

SOURCES  := list_ncvars.f90

OBJECTS  := $(SOURCES:f90=o)

all: $(TARGETS)

list_ncvars.exe: list_ncvars.o
	$(FC) -o $@ $(OBJECTS) $(LIBS)

list_ncvars.o: list_ncvars.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

install: all $(TARGETS_SHARE)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 -d $(PREFIX_SHARE)
	install -m 755 -d $(PREFIX_SHARE_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)
	install -m 755 $(TARGETS_SHARE:%=$(SRCDIR)/%) $(PREFIX_SHARE_EXEC)

docs:

install-docs: docs

clean:
	-rm -f *.o $(TARGETS)
