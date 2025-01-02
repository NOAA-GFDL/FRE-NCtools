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
bin_PROGRAMS = remap_land
if WITH_MPI
  bin_PROGRAMS += remap_land_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/lib/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)  $(RPATH_FLAGS)

remap_land_SOURCES = remap_land.c
remap_land_LDADD = $(top_builddir)/lib/libfrencutils/libfrencutils.a $(LDADD)

remap_land_parallel_SOURCES = $(remap_land_SOURCES)
remap_land_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
remap_land_parallel_LDADD = $(top_builddir)/lib/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
