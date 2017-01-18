# $Id: fre-nctools.mk,v 1.1.2.3 2012/09/12 19:40:57 Kyle.Olivo Exp $
# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build the 'combine_blobs' Program
# ------------------------------------------------------------------------------
# keo    Ver   1.00  Initial version                                September 12
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2012
# Written by Kyle Olivo
include ./env.$(SITE)

CC       := icc
CFLAGS   := -O3 -g -traceback $(CFLAGS_SITE)
CFLAGS_O2:= -O2 -g -traceback $(CFLAGS_SITE)
INCLUDES := -I${NETCDF_HOME}/include
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := combine_blobs

SOURCES  := combine_blobs.c

OBJECTS  := $(SOURCES:c=o)

HEADERS = fre-nctools.mk

all: $(TARGETS)

combine_blobs: $(OBJECTS)
	$(CC) -o $@ $^ $(CLIBS)

combine_blobs.o: combine_blobs.c $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

%.o: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(TARGETS)
