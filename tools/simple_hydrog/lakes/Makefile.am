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
bin_PROGRAMS = cr_lake_files

# Location of FMS library.
# TODO: allow builder to override and use a different libfms.
FMSROOT = $(top_builddir)/tools/simple_hydrog/libfmslite

AM_CFLAGS = $(NETCDF_CFLAGS)
AM_FCFLAGS = -I$(FMSROOT) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(NETCDF_FCLIBS)

cr_lake_files_SOURCES = cr_lake_files.f90
cr_lake_files_LDADD = $(FMSROOT)/libfms.a $(LDADD)

CLEANFILES = *.mod
