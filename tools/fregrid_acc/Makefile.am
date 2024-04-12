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
bin_PROGRAMS = fregrid_acc

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils $(NETCDF_CFLAGS) \
            -I$(top_srcdir)/tools/libfrencutils_acc \
            -I$(top_srcdir)/tools/fregrid -acc
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS)

fregrid_acc_SOURCES = conserve_interp_acc.c \
                      conserve_interp_acc.h \
                      fregrid_utils_acc.c \
                      fregrid_utils_acc.h \
                      fregrid_acc.c

fregrid_acc_LDADD = $(top_builddir)/tools/fregrid/fregrid_util.o \
                    $(top_builddir)/tools/fregrid/conserve_interp.o \
                    $(top_builddir)/tools/fregrid/bilinear_interp.o \
                    $(top_builddir)/tools/libfrencutils/libfrencutils.a \
                    $(top_builddir)/tools/libfrencutils_acc/libfrencutils_acc.a \
                    $(LDADD)
