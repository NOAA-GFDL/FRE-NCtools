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

# Test that mppnccombine reproduces a reference copy of a combined file.
# The reference file is currently from a Bronx-16 (GFDL system) mppnccombine.

load test_utils

@test "reference mppnccombine combines comparison to bronx-16 stored copy" {

  generate_all_from_ncl_num mppnccombine Test02-input

  #Combine the files into 1
   mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ -e mppnccombine_output.nc ]
  ncdump -h mppnccombine_output.nc

  [ -e $top_srcdir/t/Test02-reference/mppnccombine_output.nc ]

  nccmp -V

  nccmp -d mppnccombine_output.nc  $top_srcdir/t/Test02-reference/mppnccombine_output.nc
}
