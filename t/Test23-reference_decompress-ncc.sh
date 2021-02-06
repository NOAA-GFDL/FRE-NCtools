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


@test "decompress input netcdf files" {
  if [ ! -d "Test23" ] 
  then
  mkdir Test23
  fi

  cd Test23

  cp $top_srcdir/t/Test18-input/decompress-ncc.atmos_daily.nc.copy .  

  #Decompress compressed netcdf file(s) into 1
  run command decompress-ncc \
      decompress-ncc.atmos_daily.nc.copy \
      decompress-ncc_output.nc 
  [ "$status" -eq 0 ]
  [ -e decompress-ncc_output.nc ]
  run ncdump -h decompress-ncc_output.nc
  [ "$status" -eq 0 ]

  [ -e $top_srcdir/t/Test18-reference/decompress-ncc_output.nc ]

  run nccmp -V
  [ "$status" -eq 0 ]

  run nccmp -d decompress-ncc_output.nc  $top_srcdir/t/Test18-reference/decompress-ncc_output.nc
  [ "$status" -eq 0 ]


  #Clean up 
  cd ..
  rm -rf Test23
}
