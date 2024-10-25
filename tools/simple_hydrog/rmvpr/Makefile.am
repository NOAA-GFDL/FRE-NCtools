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
bin_PROGRAMS = rmv_parallel_rivers

AM_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS) -I$(top_builddir)/lib/libnfu

LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS) $(top_builddir)/lib/libnfu/libnfu.a $(RPATH_FLAGS)

rmv_parallel_rivers_SOURCES = rmv_parallel_rivers.F90
