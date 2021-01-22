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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
bin_PROGRAMS = combine-ncc decompress-ncc is-compressed scatter-ncc
dist_bin_SCRIPTS = combine-lumip-vars.py
# Helper library for the land utilties
noinst_LIBRARIES = libnfu.a

AM_FCFLAGS =  $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
AM_CFLAGS = $(NETCDF_CFLAGS)

libnfu_a_SOURCES = nfu.F90 nfu_compress.F90

combine_ncc_SOURCES = combine-ncc.F90
combine_ncc_LDADD = libnfu.a $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

decompress_ncc_SOURCES = decompress-ncc.F90
decompress_ncc_LDADD = libnfu.a $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

scatter_ncc_SOURCES = scatter-ncc.F90
scatter_ncc_LDADD = libnfu.a $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

is_compressed_SOURCES = is-compressed.c
is_compressed_LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

nfu_mod.$(FC_MODEXT): nfu.$(OBJEXT)
nfu_compress_mod.$(FC_MODEXT): nfu_compress.$(OBJEXT)

nfu_compress.$(OBJEXT): nfu_mod.$(FC_MODEXT)
combine-ncc.$(OBJEXT): nfu_mod.$(FC_MODEXT) nfu_compress_mod.$(FC_MODEXT)
decompress-ncc.$(OBJEXT): nfu_mod.$(FC_MODEXT) nfu_compress_mod.$(FC_MODEXT)
scatter-ncc.$(OBJEXT): nfu_mod.$(FC_MODEXT) nfu_compress_mod.$(FC_MODEXT)

clean-local:
	$(RM) *.$(FC_MODEXT)
