# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build the plevel tool
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2011
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by V. Balaji, Amy Langenhorst and Aleksey Yakovlev
# ------------------------------------------------------------------------------

# Bring in make settings for the platform
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin


VPATH = $(SRCDIR):.

FC       := ifort
FFLAGS   := -module ./modules.r8 -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 64 -traceback $(FFLAGS_SITE)
FFLAGS_r4:= -module ./modules.r4 -fltconsistency -fno-alias -stack_temps -safe_cray_ptr -ftz -assume byterecl -g -O2 -i4 -real_size 32 -traceback $(FFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include
LIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lz $(LIBS_SITE) $(STATIC)

TARGETS  := PLEV.exe PLEV.r4.exe
TARGETS_SHARE := plevel.sh

OBJECTS  := plev_constants.o moisture_convert.o pressure_interp.o pinterp_utilities.o
OBJECTS_r4= plev_constants.r4.o moisture_convert.r4.o pressure_interp.r4.o pinterp_utilities.r4.o

MODULES := modules.r8
MODULES_r4 := modules.r4

all: $(TARGETS)

PLEV.exe: run_pressure_interp.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

run_pressure_interp.o: run_pressure_interp.F90 $(OBJECTS) $(MODULES)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $< 

PLEV.r4.exe: run_pressure_interp.r4.o
	$(FC) -o $@ $^ $(OBJECTS) $(LIBS)

run_pressure_interp.r4.o: run_pressure_interp.F90 $(OBJECTS_r4) $(MODULES_r4)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $< 

pinterp_utilities.o: pinterp_utilities.F90 $(MODULES)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $< 

pinterp_utilities.r4.o: pinterp_utilities.F90 $(MODULES_r4)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $< 

pressure_interp.o: pressure_interp.F90 moisture_convert.o plev_constants.o $(MODULES)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

pressure_interp.r4.o: pressure_interp.F90 moisture_convert.r4.o plev_constants.r4.o $(MODULES_r4)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

moisture_convert.o: moisture_convert.F90  plev_constants.o $(MODULES)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

moisture_convert.r4.o: moisture_convert.F90  plev_constants.r4.o $(MODULES_r4)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

plev_constants.o: plev_constants.F90 $(MODULES)
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

plev_constants.r4.o: plev_constants.F90 $(MODULES_r4) $(MODULES_r4)
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

modules.r8:
	mkdir -p ./modules.r8

modules.r4:
	mkdir -p ./modules.r4

%.o: %.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $<

%.r4..o: %.F90
	$(FC) $(FFLAGS_r4) $(INCLUDES) -c -o $@ $<

install: all $(TARGETS_SAVE)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_EXEC)
	install -m 755 -d $(PREFIX_SHARE)
	install -m 755 -d $(PREFIX_SHARE_EXEC)
	install -m 755 $(TARGETS) $(PREFIX_EXEC)
	install -m 755 $(TARGETS_SHARE:%=$(SRCDIR)/%) $(PREFIX_SHARE_EXEC)

docs:

install-docs:

clean:
	-rm -rf *.o *.mod $(TARGETS) $(MODULES) $(MODULES_r4)
