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
bin_PROGRAMS = fregrid_gpu

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils $(NETCDF_CFLAGS) \
            -I$(top_srcdir)/tools/libfrencutils_gpu \
            -I$(top_srcdir)/tools/fregrid \
            $(OPENACC_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS)

fregrid_gpu_SOURCES = conserve_interp_gpu.c \
                      conserve_interp_gpu.h \
                      interp_utils_gpu.c \
                      interp_utils_gpu.h \
                      fregrid_utils_gpu.c \
                      fregrid_utils_gpu.h \
                      fregrid_gpu.c

fregrid_gpu_LDADD = $(top_builddir)/tools/fregrid/fregrid_util.o \
                    $(top_builddir)/tools/fregrid/conserve_interp.o \
                    $(top_builddir)/tools/fregrid/bilinear_interp.o \
                    $(top_builddir)/tools/libfrencutils/libfrencutils.a \
                    $(top_builddir)/tools/libfrencutils_gpu/libfrencutils_gpu.a \
                    $(LDADD)
