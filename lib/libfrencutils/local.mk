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

libfrencutils_libfrencutils_a_SOURCES = \
    libfrencutils/affinity.c \
    libfrencutils/constant.h \
    libfrencutils/create_xgrid.c \
    libfrencutils/create_xgrid.h \
    libfrencutils/globals.h \
    libfrencutils/gradient_c2l.c \
    libfrencutils/gradient_c2l.h \
    libfrencutils/interp.c \
    libfrencutils/interp.h \
    libfrencutils/mosaic_util.c \
    libfrencutils/mosaic_util.h \
    libfrencutils/mpp_domain.c \
    libfrencutils/mpp_domain.h \
    libfrencutils/mpp_io.c \
    libfrencutils/mpp_io.h \
    libfrencutils/mpp.c \
    libfrencutils/mpp.h \
    libfrencutils/read_mosaic.c \
    libfrencutils/read_mosaic.h \
    libfrencutils/tool_util.c \
    libfrencutils/tool_util.h

libfrencutils_libfrencutils_mpi_a_SOURCES = \
    $(libfrencutils_libfrencutils_a_SOURCES)
libfrencutils_libfrencutils_mpi_a_CFLAGS = \
    -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
