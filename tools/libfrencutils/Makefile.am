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
noinst_LIBRARIES = libfrencutils.a
if WITH_MPI
  noinst_LIBRARIES += libfrencutils_mpi.a
endif

AM_CFLAGS = $(NETCDF_CFLAGS)

libfrencutils_a_SOURCES = affinity.c \
                          constant.h \
                          create_xgrid.c \
                          create_xgrid.h \
                          globals.h \
                          gradient_c2l.c \
                          gradient_c2l.h \
                          interp.c \
                          interp.h \
                          mosaic_util.c \
                          mosaic_util.h \
                          mpp_domain.c \
                          mpp_domain.h \
                          mpp_io.c \
                          mpp_io.h \
                          mpp.c \
                          mpp.h \
                          read_mosaic.c \
                          read_mosaic.h \
                          tool_util.c \
                          tool_util.h

libfrencutils_mpi_a_SOURCES = $(libfrencutils_a_SOURCES)
libfrencutils_mpi_a_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
