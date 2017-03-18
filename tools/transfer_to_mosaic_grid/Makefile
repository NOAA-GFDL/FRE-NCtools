# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to build transfer_to_mosaic_grid
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2009-2010, 2017
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

CFLAGS   := -O3 -g -traceback $(CFLAGS_SITE)
CFLAGS_O2:= -O2 -g -traceback $(CFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include -I$(SRCDIR)/ -I$(SRCDIR)/../shared -I$(SRCDIR)/../../shared/mosaic
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := transfer_to_mosaic_grid

SOURCES  := transfer_to_mosaic.c
SOURCES  += create_xgrid.c gradient_c2l.c interp.c read_mosaic.c
SOURCES  += mpp_domain.c mpp_io.c tool_util.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = ../shared/mpp.h  ../shared/mpp_domain.h  ../shared/mpp_io.h ../shared/tool_util.h   \
	  ../../shared/mosaic/constant.h ../../shared/mosaic/create_xgrid.h  \
	  ../../shared/mosaic/gradient_c2l.h ../../shared/mosaic/interp.h  \
	  ../../shared/mosaic/mosaic_util.h  ../../shared/mosaic/read_mosaic.h 

all: $(TARGETS)

transfer_to_mosaic_grid: $(OBJECTS) mosaic_util.o mpp.o
	$(CC) -o $@ $^ $(CLIBS)

mosaic_util.o: ../../shared/mosaic/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

read_mosaic.o: ../../shared/mosaic/read_mosaic.c $(HEADERS)
	$(CC) -Duse_netCDF $(CFLAGS) $(INCLUDES) -c $< 

gradient_c2l.o: ../../shared/mosaic/gradient_c2l.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

interp.o: ../../shared/mosaic/interp.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

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
