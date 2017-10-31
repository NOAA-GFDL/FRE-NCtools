# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile template to Build the 'combine_blobs' Program
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2012, 2017
# Written by Kyle Olivo and Seth Underwood
# ------------------------------------------------------------------------------

# Default compilers
CC       := icc

# Include the build settings for the site
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
INCLUDES := -I${NETCDF_HOME}/include
CLIBS     := -L${NETCDF_HOME}/lib -L${HDF5_HOME}/lib -lnetcdf -lhdf5_hl -lhdf5 -lz $(CLIBS_SITE) $(STATIC)

TARGETS  := combine_blobs

SOURCES  := combine_blobs.c

OBJECTS  := $(SOURCES:c=o)

all: $(TARGETS)

combine_blobs: $(OBJECTS)
	$(CC) -o $@ $^ $(CLIBS)

combine_blobs.o: combine_blobs.c
	$(CC) $(CFLAGS) $(INCLUDES) -c $< 

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
