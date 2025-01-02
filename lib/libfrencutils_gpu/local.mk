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

libfrencutils_gpu_libfrencutils_gpu_a_SOURCES = \
    libfrencutils_gpu\create_xgrid_gpu.c \
    libfrencutils_gpu\create_xgrid_gpu.h \
    libfrencutils_gpu\create_xgrid_utils_gpu.c \
    libfrencutils_gpu\create_xgrid_utils_gpu.h \
    libfrencutils_gpu\general_utils_gpu.c \
    libfrencutils_gpu\general_utils_gpu.h

libfrencutils_gpu_libfrencutils_gpu_a_CFLAGS = \
    $(AM_CFLAGS) $(OPENACC_CFLAGS) -I$(top_srcdir)/lib/libfrencutils
libfrencutils_gpu_libfrencutils_gpu_a_LIBADD = \
    $(top_builddir)/lib/libfrencutils/libfrencutils.a
