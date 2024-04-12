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

check_PROGRAMS = test_read_remap_file

AM_CFLAGS = $(NETCDF_CFLAGS) \
            -I$(top_srcdir)/tools/libfrencutils  \
            -I$(top_srcdir)/tools/libfrencutils_acc \
            -I$(top_srcdir)/tools/fregrid_acc -acc

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS) \
        $(top_builddir)/tools/fregrid_acc/conserve_interp_acc.o \
        $(top_builddir)/tools/fregrid_acc/fregrid_utils_acc.o \
        $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(top_builddir)/tools/libfrencutils_acc/libfrencutils_acc.a

test_read_remap_file_SOURCES = test_read_remap_file.c

TESTS = test_read_remap_file_conserve.sh

EXTRA_DIST = test_read_remap_file_conserve.sh \
             test_make_remap_file.py

cp_test_make_remap_file : test_make_remap_file_conserve.py
	cp $(top_srcdir)/t_acc/test_read_remap_file/test_make_remap_file_conserve.py \
    $(top_builddir)/t_acc/test_read_remap_file/test_make_remap_file_conserve.py


TEST_EXTENSIONS = .sh

CLEANFILES = *.nc *txt
