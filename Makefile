# ------------------------------------------------------------------------------
# FRE netCDF Tools Project: Master Makefile
# ------------------------------------------------------------------------------
# Copyright (C) NOAA Geophysical Fluid Dynamics Laboratory, 2017
# Written by Seth Underwood
# ------------------------------------------------------------------------------

# This Makefile is used to build the FRE-NCTools packages.  This is designed to
# be used in conjuction with the make package script included in the
# fre-nctools package.

# Multiple make jobs
MAKEFLAGS += --jobs=8

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

postpSRC := land_utils list_ncvars combine_blobs mppnccombine plevel timavg
postpSRC += split_ncvars iceberg_comb combine_restarts

all: $(postpSRC) $(toolsSRC)

docs: $(postpSRC:%=%-docs) $(toolsSRC:%=%-docs)

install: $(postpSRC:%=%-install) $(toolsSRC:%=%-install)

install-docs: $(postpSRC:%=%-install-docs) $(toolsSRC:%=%-install-docs)

# ------------------------------------------------------------------------------
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

plevel-install-docs:
	make -C postprocessing/plevel PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/plevel install-docs

split_ncvars-install-docs:
	make -C postprocessing/split_ncvars PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/split_ncvars install-docs

timavg-install-docs:
	make -C postprocessing/timavg PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/postprocessing/timavg install-docs

# ------------------------------------------------------------------------------
# Packages in the postprocessing directory
check_mask:
	make -C tools/check_mask SRCDIR=$(SRCDIR)/tools/check_mask

fregrid:
	make -C tools/fregrid SRCDIR=$(SRCDIR)/tools/fregrid

make_coupler_mosaic:
	make -C tools/make_coupler_mosaic SRCDIR=$(SRCDIR)/tools/make_coupler_mosaic

make_hgrid:
	make -C tools/make_hgrid SRCDIR=$(SRCDIR)/tools/make_hgrid

make_land_domain:
	make -C tools/make_land_domain SRCDIR=$(SRCDIR)/tools/make_land_domain

make_quick_mosaic:
	make -C tools/make_quick_mosaic SRCDIR=$(SRCDIR)/tools/make_quick_mosaic

make_regional_mosaic:
	make -C tools/make_regional_mosaic SRCDIR=$(SRCDIR)/tools/make_regional_mosaic

make_solo_mosaic:
	make -C tools/make_solo_mosaic SRCDIR=$(SRCDIR)/tools/make_solo_mosaic

make_topog:
	make -C tools/make_topog SRCDIR=$(SRCDIR)/tools/make_topog

make_vgrid:
	make -C tools/make_vgrid SRCDIR=$(SRCDIR)/tools/make_vgrid

mppncscatter:
	make -C tools/mppncscatter SRCDIR=$(SRCDIR)/tools/mppncscatter

ncexists:
	make -C tools/ncexists SRCDIR=$(SRCDIR)/tools/ncexists

remap_land:
	make -C tools/remap_land SRCDIR=$(SRCDIR)/tools/remap_land

river_regrid:
	make -C tools/river_regrid SRCDIR=$(SRCDIR)/tools/river_regrid

runoff_regrid:
	make -C tools/runoff_regrid SRCDIR=$(SRCDIR)/tools/runoff_regrid

transfer_to_mosaic_grid:
	make -C tools/transfer_to_mosaic_grid SRCDIR=$(SRCDIR)/tools/transfer_to_mosaic_grid

check_mask-install:
	make -C tools/check_mask PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/check_mask install

fregrid-install:
	make -C tools/fregrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/fregrid install

make_coupler_mosaic-install:
	make -C tools/make_coupler_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_coupler_mosaic install

make_hgrid-install:
	make -C tools/make_hgrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_hgrid install

make_land_domain-install:
	make -C tools/make_land_domain PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_land_domain install

make_quick_mosaic-install:
	make -C tools/make_quick_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_quick_mosaic install

make_regional_mosaic-install:
	make -C tools/make_regional_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_regional_mosaic install

make_solo_mosaic-install:
	make -C tools/make_solo_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_solo_mosaic install

make_topog-install:
	make -C tools/make_topog PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_topog install

make_vgrid-install:
	make -C tools/make_vgrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_vgrid install

mppncscatter-install:
	make -C tools/mppncscatter PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/mppncscatter install

ncexists-install:
	make -C tools/ncexists PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/ncexists install

remap_land-install:
	make -C tools/remap_land PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/remap_land install

river_regrid-install:
	make -C tools/river_regrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/river_regrid install

runoff_regrid-install:
	make -C tools/runoff_regrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/runoff_regrid install

transfer_to_mosaic_grid-install:
	make -C tools/transfer_to_mosaic_grid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/transfer_to_mosaic_grid install

check_mask-docs:
	make -C tools/check_mask SRCDIR=$(SRCDIR)/tools/check_mask docs

fregrid-docs:
	make -C tools/fregrid SRCDIR=$(SRCDIR)/tools/fregrid docs

make_coupler_mosaic-docs:
	make -C tools/make_coupler_mosaic SRCDIR=$(SRCDIR)/tools/make_coupler_mosaic docs

make_hgrid-docs:
	make -C tools/make_hgrid SRCDIR=$(SRCDIR)/tools/make_hgrid docs

make_land_domain-docs:
	make -C tools/make_land_domain SRCDIR=$(SRCDIR)/tools/make_land_domain docs

make_quick_mosaic-docs:
	make -C tools/make_quick_mosaic SRCDIR=$(SRCDIR)/tools/make_quick_mosaic docs

make_regional_mosaic-docs:
	make -C tools/make_regional_mosaic SRCDIR=$(SRCDIR)/tools/make_regional_mosaic docs

make_solo_mosaic-docs:
	make -C tools/make_solo_mosaic SRCDIR=$(SRCDIR)/tools/make_solo_mosaic docs

make_topog-docs:
	make -C tools/make_topog SRCDIR=$(SRCDIR)/tools/make_topog docs

make_vgrid-docs:
	make -C tools/make_vgrid SRCDIR=$(SRCDIR)/tools/make_vgrid docs

mppncscatter-docs:
	make -C tools/mppncscatter SRCDIR=$(SRCDIR)/tools/mppncscatter docs

ncexists-docs:
	make -C tools/ncexists SRCDIR=$(SRCDIR)/tools/ncexists docs

remap_land-docs:
	make -C tools/remap_land SRCDIR=$(SRCDIR)/tools/remap_land docs

river_regrid-docs:
	make -C tools/river_regrid SRCDIR=$(SRCDIR)/tools/river_regrid docs

runoff_regrid-docs:
	make -C tools/runoff_regrid SRCDIR=$(SRCDIR)/tools/runoff_regrid docs

transfer_to_mosaic_grid-docs:
	make -C tools/transfer_to_mosaic_grid SRCDIR=$(SRCDIR)/tools/transfer_to_mosaic_grid docs

check_mask-install-docs:
	make -C tools/check_mask PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/check_mask install-docs

fregrid-install-docs:
	make -C tools/fregrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/fregrid install-docs

make_coupler_mosaic-install-docs:
	make -C tools/make_coupler_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_coupler_mosaic install-docs

make_hgrid-install-docs:
	make -C tools/make_hgrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_hgrid install-docs

make_land_domain-install-docs:
	make -C tools/make_land_domain PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_land_domain install-docs

make_quick_mosaic-install-docs:
	make -C tools/make_quick_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_quick_mosaic install-docs

make_regional_mosaic-install-docs:
	make -C tools/make_regional_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_regional_mosaic install-docs

make_solo_mosaic-install-docs:
	make -C tools/make_solo_mosaic PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_solo_mosaic install-docs

make_topog-install-docs:
	make -C tools/make_topog PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_topog install-docs

make_vgrid-install-docs:
	make -C tools/make_vgrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/make_vgrid install-docs

mppncscatter-install-docs:
	make -C tools/mppncscatter PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/mppncscatter install-docs

ncexists-install-docs:
	make -C tools/ncexists PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/ncexists install-docs

remap_land-install-docs:
	make -C tools/remap_land PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/remap_land install-docs

river_regrid-install-docs:
	make -C tools/river_regrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/river_regrid install-docs

runoff_regrid-install-docs:
	make -C tools/runoff_regrid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/runoff_regrid install-docs

transfer_to_mosaic_grid-install-docs:
	make -C tools/transfer_to_mosaic_grid PREFIX=$(PREFIX) SRCDIR=$(SRCDIR)/tools/transfer_to_mosaic_grid install-docs
