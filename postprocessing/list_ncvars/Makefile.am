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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
bin_PROGRAMS = list_ncvars.exe
dist_bin_SCRIPTS = list_ncvars.csh

list_ncvars_exe_SOURCES = list_ncvars.f90
list_ncvars_exe_FCFLAGS = $(NETCDF_FCFLAGS) $(FC_DEFAULT_REAL_KIND8_FLAG) $(AM_FCFLAGS)
list_ncvars_exe_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

clean-local:
	$(RM) *.mod
