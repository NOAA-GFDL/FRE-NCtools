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

check_PROGRAMS = test_get_upbound_nxcells_2dx2d

AM_CFLAGS =  $(NETCDF_CFLAGS) \
            -I$(top_srcdir)/tools/fregrid_gpu \
            -I$(top_srcdir)/tools/libfrencutils \
            -I$(top_srcdir)/tools/libfrencutils_gpu -acc

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS) \
        $(top_builddir)/tools/fregrid_gpu/conserve_interp_gpu.o \
        $(top_builddir)/tools/fregrid_gpu/interp_utils_gpu.o \
        $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(top_builddir)/tools/libfrencutils_gpu/libfrencutils_gpu.a

test_get_upbound_nxcells_2dx2d_SOURCES = test_get_upbound_nxcells_2dx2d.c

TESTS = test_get_upbound_nxcells_2dx2d
