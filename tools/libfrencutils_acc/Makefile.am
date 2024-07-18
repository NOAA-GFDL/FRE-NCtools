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
noinst_LIBRARIES = libfrencutils_acc.a

AM_CFLAGS = $(NETCDF_CFLAGS) $(OPENACC_CFLAGS) -I$(top_srcdir)/tools/libfrencutils

libfrencutils_acc_a_SOURCES = create_xgrid_acc.c \
                              create_xgrid_acc.h

libfrencutils_acc_a_LIBADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a
