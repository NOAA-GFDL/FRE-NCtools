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

noinst_LIBRARIES = \
    libnfu/libnfu.a


AM_FCFLAGS =  $(FC_MODOUT)$(@D) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(RPATH_FLAGS)

# Items to handle Fortran module files
SUFFIXES = _mod.$(FC_MODEXT) .$(FC_MODEXT)

.F90_mod.$(FC_MODEXT):
	@$(MKDIR_P) $(@D)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@D)/$(<F:.F90=.o) $<

.F90.$(FC_MODEXT):
	@$(MKDIR_P) $(@D)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@D)/$(<F:.F90=.o) $<

clean-local:
	$(RM) *.$(FC_MODEXT)

include $(top_srcdir)/lib/libnfu/local.mk