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
bin_PROGRAMS = runoff_regrid
if WITH_MPI
  bin_PROGRAMS += runoff_regrid_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

runoff_regrid_SOURCES = runoff_regrid.c
runoff_regrid_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)

runoff_regrid_parallel_SOURCES = $(runoff_regrid_SOURCES)
runoff_regrid_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
runoff_regrid_parallel_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
