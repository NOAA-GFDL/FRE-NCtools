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


@test "reference combine-ncc combines comparison to bronx-16 stored copy" {
  if [ ! -d "Test22" ] 
  then
  mkdir Test22
  fi

  cd Test22

  #ncks atmos daily *.nc file
 
  cp $top_srcdir/t/Test17-input/combine-ncc.atmos_daily.nc.copy .

  #Combine netcdf copy file 
  run command combine-ncc \
      combine-ncc.atmos_daily.nc.copy \
      combine-ncc_output.nc 
  [ "$status" -eq 0 ]
  [ -e combine-ncc_output.nc ]
  run ncdump -h combine-ncc_output.nc
  [ "$status" -eq 0 ]

  [ -e $top_srcdir/t/Test17-reference/combine-ncc_output.nc ]

  run nccmp -V
  [ "$status" -eq 0 ]

  run nccmp -d combine-ncc_output.nc  $top_srcdir/t/Test17-reference/combine-ncc_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test22
}
