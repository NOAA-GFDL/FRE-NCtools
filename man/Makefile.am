#***********************************************************************
#                   GNU Lesser General Public License
#
# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).
#
# FRE-NCTools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# FRE-NCTools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
EXTRA_DIST = \
  check_mask.txt \
  combine-ncc.txt \
  combine_restarts.txt \
  decompress-ncc.txt \
  fregrid.txt \
  iceberg_comb.sh.txt \
  is-compressed.txt \
  list_ncvars.txt \
  list_ncvars.sh.txt \
  make_hgrid.txt \
  make_vgrid.txt \
  make_topog.txt \
  make_coupler_mosaic.txt \
  make_quick_mosaic.txt \
  make_regional_mosaic.txt \
  make_remap_file.txt \
  make_solo_mosaic.txt \
  mppnccombine.txt \
  mppncscatter.txt \
  ncexists.txt \
  plevel_interp.txt \
  plevel.sh.txt \
  river_regrid.txt \
  remap_land.txt \
  scatter-ncc.txt \
  split_ncvars.pl.txt \
  timavg.txt \
  time_average.txt \
  asciidoc.conf

# Documentation settings
A2X = a2x
ASCIIDOC_CFG = $(srcdir)/asciidoc.conf

man_MANS =

if WITH_ASCIIDOC
man_MANS += \
  check_mask.1 \
  combine-ncc.1 \
  combine_restarts.1 \
  decompress-ncc.1 \
  fregrid.1 \
  iceberg_comb.sh.1 \
  is-compressed.1 \
  list_ncvars.1 \
  list_ncvars.csh.1 \
  list_ncvars.sh.1 \
  make_hgrid.1 \
  make_vgrid.1 \
  make_topog.1 \
  make_coupler_mosaic.1 \
  make_quick_mosaic.1 \
  make_regional_mosaic.1 \
  make_remap_file.1 \
  make_solo_mosaic.1 \
  mppnccombine.1 \
  mppncscatter.1 \
  ncexists.1 \
  plevel_interp.1 \
  plevel.sh.1 \
  river_regrid.1 \
  remap_land.1 \
  scatter-ncc.1 \
  split_ncvars.pl.1 \
  timavg.1 \
  timavg.csh.1 \
  time_average.1
if WITH_MPI
man_MANS += \
  make_coupler_mosaic_parallel.1 \
  make_topog_parallel.1 \
  remap_land_parallel.1
endif
endif

.txt.1:
	$(A2X) -v -L\
          -D $(builddir) \
          -d manpage \
          -f manpage \
          -afre_version=$(PACKAGE_VERSION) \
          -apackage_bugreport=$(PACKAGE_BUGREPORT) \
          -acopyright_year=$(COPYRIGHT_YEAR) \
          --asciidoc-opts="-f $(ASCIIDOC_CFG)" \
          $<

make_alias_1 = @echo ".so man1/$<" > $@

list_ncvars.csh.1: list_ncvars.sh.1
	$(make_alias_1)

make_coupler_mosaic_parallel.1: make_coupler_mosaic.1
	$(make_alias_1)

make_topog_parallel.1: make_topog.1
	$(make_alias_1)

remap_land_parallel.1: remap_land.1
	$(make_alias_1)

timavg.csh.1: timavg.1
	$(make_alias_1)

clean-local:
	$(RM) *.1 *.xml
