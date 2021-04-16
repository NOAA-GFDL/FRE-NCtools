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
bin_PROGRAMS = make_topog
if WITH_MPI
  bin_PROGRAMS += make_topog_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_topog_SOURCES = make_topog.c \
                     topog.c \
                     topog.h
make_topog_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)

make_topog_parallel_SOURCES = $(make_topog_SOURCES)
make_topog_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
make_topog_parallel_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
