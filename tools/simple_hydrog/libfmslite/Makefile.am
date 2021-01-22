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
noinst_LIBRARIES = libfms.a

AM_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG)

libfms_a_SOURCES = constants.F90 \
  horiz_interp.F90 \
  horiz_interp_bilinear.F90 \
  horiz_interp_type.F90 \
  mpp.F90 \
  mpp_error_a_a.h \
  mpp_error_a_s.h \
  mpp_error_s_a.h \
  mpp_error_s_s.h

constants_mod.mod: constants.$(OBJEXT)
mpp_mod.mod: mpp.$(OBJEXT)
horiz_interp_type_mod.mod: mpp_mod.mod horiz_interp_type.$(OBJEXT)
horiz_interp_bilinear_mod.mod: constants_mod.mod mpp_mod.mod horiz_interp_type_mod.mod horiz_interp_bilinear.$(OBJEXT)
horiz_interp_mod.mod: constants_mod.mod mpp_mod.mod horiz_interp_type_mod.mod horiz_interp_bilinear_mod.mod horiz_intrep.$(OBJEXT)

mpp.$(OBJEXT): mpp.F90 mpp_error_a_a.h mpp_error_a_s.h mpp_error_s_a.h mpp_error_s_s.h
horiz_interp_type.$(OBJEXT): horiz_interp_type.F90 mpp_mod.mod
horiz_interp_bilinear.$(OBJEXT): horiz_interp_bilinear.F90 constants_mod.mod mpp_mod.mod horiz_interp_type_mod.mod
horiz_interp.$(OBJEXT): horiz_interp.F90 constants_mod.mod mpp_mod.mod horiz_interp_type_mod.mod horiz_interp_bilinear_mod.mod

CLEANFILES = *.mod
