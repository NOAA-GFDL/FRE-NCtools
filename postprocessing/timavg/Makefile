# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build the timeaverage tool
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
# ------------------------------------------------------------------------------

# Default compilers
FC       := ifort

# Bring in make settings for the site
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

FFLAGS   := -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 64 -traceback $(FFLAGS_SITE)
FFLAGS_r4:= -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 32 -traceback $(FFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz $(LIBS_SITE) $(STATIC)

TARGETS  := TAVG.exe TAVG.r4.exe
TARGETS_SHARE := timavg.csh

OBJECTS  := time_average.o
OBJECTS_r4= time_average.r4.o

all: $(TARGETS)

TAVG.exe: time_average.o
	$(FC) -o $@ $(OBJECTS) $(LIBS)

time_average.o: time_average.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $< 

TAVG.r4.exe: time_average.r4.o
	$(FC) -o $@ $(OBJECTS_r4) $(LIBS)

time_average.r4.o: time_average.f90
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $< 

%.o: %.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

%.r4..o: %.f90
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

docs:

install: all $(TARGETS_SHARE)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 -d $(PREFIX_SHARE)
	install -m 755 -d $(PREFIX_SHARE_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)
	install -m 755 $(TARGETS_SHARE:%=$(SRCDIR)/%) $(PREFIX_SHARE_EXEC)

install-docs: docs

clean:
	-rm -f *.o *.mod $(TARGETS)
