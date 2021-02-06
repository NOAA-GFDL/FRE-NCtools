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
bin_PROGRAMS = TAVG.exe TAVG.r4.exe
dist_bin_SCRIPTS = timavg.csh

AM_FCFLAGS =  $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

TAVG_exe_SOURCES = time_average.f90
TAVG_exe_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG) $(AM_FCFLAGS)

TAVG_r4_exe_SOURCES = time_average.f90
TAVG_r4_exe_FCFLAGS = $(FC_DEFAULT_REAL_KIND4_FLAG) $(AM_FCFLAGS)

clean-local:
	$(RM) *.mod
