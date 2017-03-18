# ------------------------------------------------------------------------------
# FRE netCDF Tools Project: Master Makefile
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2017
# Written by Seth Underwood
# ------------------------------------------------------------------------------

# This Makefile is used to build the FRE-NCTools packages.  This is designed to
# be used in conjuction with the make package script included in the
# fre-nctools package.

# Build site
SITE :=
export SITE

# Default directory locations
SRCDIR := .
PREFIX := .

toolsSRC := check_mask fregrid make_coupler_mosaic make_hgrid make_regional_mosaic
toolsSRC += make_quick_mosaic make_solo_mosaic make_topog make_vgrid ncexists
toolsSRC += remap_land river_regrid runoff_regrid transfer_to_mosaic_grid 
toolsSRC += mppncscatter make_land_domain

postpSRC := land_utils list_ncvars combine_blobs mppnccombine ncx plevel timavg
postpSRC += split_ncvars iceberg_comb combine_restarts

all: $(postpSRC)

docs:

install: $(postpSRC:%=%-install)

install-docs:

# Packages in the postprocessing directory
combine_blobs:
	make -C postprocessing/combine_blobs SRCDIR=$(SRCDIR)/postprocessing/combine_blobs

combine_restarts:
	make -C postprocessing/combine_restarts SRCDIR=$(SRCDIR)/postprocessing/combine_restarts

iceberg_comb:
	make -C postprocessing/iceberg_comb SRCDIR=$(SRCDIR)/postprocessing/iceberg_comb

land_utils:
	make -C postprocessing/land_utils SRCDIR=$(SRCDIR)/postprocessing/land_utils

list_ncvars:
	make -C postprocessing/list_ncvars SRCDIR=$(SRCDIR)/postprocessing/list_ncvars

mppnccombine:
	make -C postprocessing/mppnccombine SRCDIR=$(SRCDIR)/postprocessing/mppnccombine

ncx:
	make -C postprocessing/ncx SRCDIR=$(SRCDIR)/postprocessing/ncx

plevel:
	make -C postprocessing/plevel SRCDIR=$(SRCDIR)/postprocessing/plevel

split_ncvars:
	make -C postprocessing/split_ncvars SRCDIR=$(SRCDIR)/postprocessing/split_ncvars

timavg:
	make -C postprocessing/timavg SRCDIR=$(SRCDIR)/postprocessing/timavg

combine_blobs-install:
	make -C postprocessing/combine_blobs PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/combine_blobs install

combine_restarts-install:
	make -C postprocessing/combine_restarts PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/combine_restarts install

iceberg_comb-install:
	make -C postprocessing/iceberg_comb PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/iceberg_comb install

land_utils-install:
	make -C postprocessing/land_utils PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/land_utils install

list_ncvars-install:
	make -C postprocessing/list_ncvars PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/list_ncvars install

mppnccombine-install:
	make -C postprocessing/mppnccombine PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/mppnccombine install

ncx-install:
	make -C postprocessing/ncx PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/ncx install

plevel-install:
	make -C postprocessing/plevel PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/plevel install

split_ncvars-install:
	make -C postprocessing/split_ncvars PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/split_ncvars install

timavg-install:
	make -C postprocessing/timavg PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/timavg install

combine_blobs-docs:
	make -C postprocessing/combine_blobs SRCDIR=$(SRCDIR)/postprocessing/combine_blobs docs

combine_restarts-docs:
	make -C postprocessing/combine_restarts SRCDIR=$(SRCDIR)/postprocessing/combine_restarts docs

iceberg_comb-docs:
	make -C postprocessing/iceberg_comb SRCDIR=$(SRCDIR)/postprocessing/iceberg_comb docs

land_utils-docs:
	make -C postprocessing/land_utils SRCDIR=$(SRCDIR)/postprocessing/land_utils docs

list_ncvars-docs:
	make -C postprocessing/list_ncvars SRCDIR=$(SRCDIR)/postprocessing/list_ncvars docs

mppnccombine-docs:
	make -C postprocessing/mppnccombine SRCDIR=$(SRCDIR)/postprocessing/mppnccombine docs

ncx-docs:
	make -C postprocessing/ncx SRCDIR=$(SRCDIR)/postprocessing/ncx docs

plevel-docs:
	make -C postprocessing/plevel SRCDIR=$(SRCDIR)/postprocessing/plevel docs

split_ncvars-docs:
	make -C postprocessing/split_ncvars SRCDIR=$(SRCDIR)/postprocessing/split_ncvars docs

timavg-docs:
	make -C postprocessing/timavg SRCDIR=$(SRCDIR)/postprocessing/timavg docs

combine_blobs-install-docs:
	make -C postprocessing/combine_blobs PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/combine_blobs install-docs

combine_restarts-install-docs:
	make -C postprocessing/combine_restarts PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/combine_restarts install-docs

iceberg_comb-install-docs:
	make -C postprocessing/iceberg_comb PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/iceberg_comb install-docs

land_utils-install-docs:
	make -C postprocessing/land_utils PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/land_utils install-docs

list_ncvars-install-docs:
	make -C postprocessing/list_ncvars PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/list_ncvars install-docs

mppnccombine-install-docs:
	make -C postprocessing/mppnccombine PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/mppnccombine install-docs

ncx-install-docs:
	make -C postprocessing/ncx PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/ncx install-docs

plevel-install-docs:
	make -C postprocessing/plevel PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/plevel install-docs

split_ncvars-install-docs:
	make -C postprocessing/split_ncvars PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/split_ncvars install-docs

timavg-install-docs:
	make -C postprocessing/timavg PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/timavg install-docs

