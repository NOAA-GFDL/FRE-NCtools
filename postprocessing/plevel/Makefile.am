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
# Disable parallel make for these targets.  This is b/c the two executables use
# the same source, and in turn create module files with the same name -- which
# if run in parallel can be overwritten by the other compiles.
# This will stay in place until someone finds a way to do this.
MAKEFLAGS = -j 1

bin_PROGRAMS = PLEV.exe PLEV.r4.exe
dist_bin_SCRIPTS = plevel.sh

AM_FCFLAGS =  $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS) $(RPATH_FLAGS)

PLEV_exe_SOURCES = plev_constants.F90 \
                   moisture_convert.F90 \
                   pressure_interp.F90 \
                   pinterp_utilities.F90 \
                   run_pressure_interp.F90
PLEV_exe_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG) $(AM_FCFLAGS)

PLEV_r4_exe_SOURCES = plev_constants.F90 \
                      moisture_convert.F90 \
                      pressure_interp.F90 \
                      pinterp_utilities.F90 \
                      run_pressure_interp.F90
PLEV_r4_exe_FCFLAGS = $(FC_DEFAULT_REAL_KIND4_FLAG) $(AM_FCFLAGS)

plev_constants_mod.$(FC_MODEXT): plev_constants.$(OBJEXT)
moisture_convert_mod.$(FC_MODEXT): moisture_convert.$(OBJEXT)
pressure_interp_mod.$(FC_MODEXT): pressure_interp.$(OBJEXT)
pinterp_utilities_mod.$(FC_MODEXT): pinterp_utilities.$(OBJEXT)

moisture_convert.$(OBJEXT): plev_constants_mod.$(FC_MODEXT)
pressure_interp.$(OBJEXT): plev_constants_mod.$(FC_MODEXT) moisture_convert_mod.$(FC_MODEXT)


clean-local:
	$(RM) *.mod
