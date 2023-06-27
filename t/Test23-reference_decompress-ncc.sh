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
load test_utils

@test "decompress input netcdf files" {

  cp $top_srcdir/t/Test18-input/decompress-ncc.atmos_daily.nc.copy .

  #Decompress compressed netcdf file(s) into 1
   decompress-ncc \
      decompress-ncc.atmos_daily.nc.copy \
      decompress-ncc_output.nc
  run_and_check [ -e decompress-ncc_output.nc ]
  ncdump -h decompress-ncc_output.nc

  run_and_check [ -e $top_srcdir/t/Test18-reference/decompress-ncc_output.nc ]

  run_and_check nccmp -V

  run_and_check nccmp -d decompress-ncc_output.nc  $top_srcdir/t/Test18-reference/decompress-ncc_output.nc
}
