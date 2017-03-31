# ------------------------------------------------------------------------------
# FMS/FRE Project: Makefile to Build split_ncvars
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2017
# This program is distributed under the terms of the GNU General Public
# License. See the file COPYING contained in this directory
# ------------------------------------------------------------------------------
# Designed and written by Seth Underwood
# ------------------------------------------------------------------------------

# Bring in setting for the build
include ../../build.mk

# Default Directories
SRCDIR := .
PREFIX := .
PREFIX_EXEC = $(PREFIX)/$(SITE)/bin
PREFIX_SHARE = $(PREFIX)/share
PREFIX_DOCS = $(PREFIX_SHARE)
PREFIX_SHARE_EXEC = $(PREFIX_SHARE)/bin

VPATH = $(SRCDIR):.

TARGETS_SHARE = split_ncvars.csh split_ncvars.py split_ncvars.pl

all:

install: all $(TARGET_SHARE)
	install -m 755 -d $(PREFIX)
	install -m 755 -d $(PREFIX_SHARE)
	install -m 755 -d $(PREFIX_SHARE_EXEC)
	install -m 755 $(TARGETS_SHARE:%=$(SRCDIR)/%) $(PREFIX_SHARE_EXEC)

docs:

install-docs: docs

clean:
