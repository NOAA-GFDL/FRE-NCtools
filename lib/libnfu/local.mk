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

libnfu_libnfu_a_SOURCES = \
    libnfu/nfu.F90 \
    libnfu/nfu_compress.F90

# nfu_mod.mod is a side effect of compiling nfu.F90 into nfu.o — both are
# produced by the same gfortran invocation.  Without this explicit rule,
# make -j would use the .F90.mod suffix rule to build nfu_mod.mod, which
# compiles nfu.F90 a second time in parallel with the .F90.o job.  The two
# simultaneous gfortran processes race to rename nfu_mod.mod0 → nfu_mod.mod
# and the loser dies with "Cannot rename … No such file or directory".
# This explicit no-op rule overrides the suffix rule and enforces that
# nfu_mod.mod is only produced once, as a side effect of nfu.o.
libnfu/nfu_mod.$(FC_MODEXT): libnfu/nfu.$(OBJEXT)
	@:

libnfu/nfu_compress.$(OBJEXT): libnfu/nfu_mod.$(FC_MODEXT)

CLEANFILES += libnfu/nfu_mod.$(FC_MODEXT) libnfu/nfu_compress_mod.$(FC_MODEXT)
