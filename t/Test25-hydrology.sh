#!/usr/bin/env bats

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
load input_util

@test "Wrapper complete hydrology test" {

  # Not all systems have /bin/csh, skip if it doesn't exist
  test ! -e /bin/csh && skip 'System does not have /bin/csh'
  # Skip if Test25-input is missing
  # TODO: Get a way to download Test25-input if missing
  test ! -e $top_srcdir/t/Test25-input/grid_spec.nc && skip 'Input directory does not exist on this system'

  run $top_srcdir/tools/simple_hydrog/share/make_simple_hydrog.csh
  [ "$status" -eq 1 ]

  mkdir tmpdir
  export TMPDIR=tmpdir

  #Run wrapper hydrology script
  run $top_srcdir/tools/simple_hydrog/share/make_simple_hydrog.csh -f 0. -t 1.e-5 -m $top_srcdir/t/Test25-input/grid_spec.nc
  [ "$status" -eq 0 ]

}
