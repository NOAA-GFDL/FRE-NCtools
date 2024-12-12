#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

check_PROGRAMS = test_get_interp_order1

AM_CFLAGS =  $(NETCDF_CFLAGS) \
            -I$(top_srcdir)/tools/fregrid_gpu \
            -I$(top_srcdir)/tools/libfrencutils \
            -I$(top_srcdir)/tools/libfrencutils_gpu -acc

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS) \
        $(top_builddir)/tools/fregrid_gpu/interp_utils_gpu.o \
        $(top_builddir)/tools/fregrid_gpu/conserve_interp_gpu.o \
        $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(top_builddir)/tools/libfrencutils_gpu/libfrencutils_gpu.a

test_get_interp_order1_SOURCES = test_get_interp_order1.c

TESTS = test_get_interp_order1
