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

@test "mppnccombine combines" {
  if [ ! -d "Test02" ] 
  then
  mkdir Test02
  fi

  cd Test02

  #Generate the netcdf files from the .ncl 
  for f in $top_srcdir/t/Test02-input/*.ncl.????
  do
    ncgen -o mppnccombine.nc.${f##*.} $f
  done

  #Combine the files into 1
  run command mppnccombine \
      mppnccombine_output.nc \
      mppnccombine.nc.????
  [ "$status" -eq 0 ]
  [ -e mppnccombine_output.nc ]
  run ncdump -h mppnccombine_output.nc
  [ "$status" -eq 0 ]

  #Clean up 
  cd ..
  rm -rf Test02
}
