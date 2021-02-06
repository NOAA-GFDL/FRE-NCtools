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

# first use mppncscatter to split the file
#             and then use mppnccombine to combine the file. The output should reproduce
#             original file

@test "Test mppncscatter and  mppnccombine" {
  if [ ! -d "Test12" ] 
  then
  		mkdir Test12
  fi

  cd Test12
  mkdir input
  cd input
  ncgen -o fv_core.res.tile1.nc $top_srcdir/t/Test12-input/fv_core.res.tile1.ncl
  ncgen -o fv_core.res.tile2.nc $top_srcdir/t/Test12-input/fv_core.res.tile2.ncl
  ncgen -o fv_core.res.tile3.nc $top_srcdir/t/Test12-input/fv_core.res.tile3.ncl
  ncgen -o fv_core.res.tile4.nc $top_srcdir/t/Test12-input/fv_core.res.tile4.ncl
  ncgen -o fv_core.res.tile5.nc $top_srcdir/t/Test12-input/fv_core.res.tile5.ncl
  ncgen -o fv_core.res.tile6.nc $top_srcdir/t/Test12-input/fv_core.res.tile6.ncl

  cd ..

#TO DO: For loop through tiles?
#Tile1 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile1.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile1.nc fv_core.res.tile1.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile1.nc input/fv_core.res.tile1.nc
#  [ "$status" -eq 0 ]

#Tile2 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile2.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile2.nc fv_core.res.tile2.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile2.nc input/fv_core.res.tile2.nc
#  [ "$status" -eq 0 ]

#Tile3 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile3.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile3.nc fv_core.res.tile3.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile3.nc input/fv_core.res.tile3.nc
#  [ "$status" -eq 0 ]

#Tile4 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile4.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile4.nc fv_core.res.tile4.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile4.nc input/fv_core.res.tile4.nc
#  [ "$status" -eq 0 ]

#Tile5 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile5.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile5.nc fv_core.res.tile5.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile5.nc input/fv_core.res.tile5.nc
#  [ "$status" -eq 0 ]

#Tile6 
#Split the file: 
  run command mppncscatter -i 2 -j 3 -x 2 -y 12 input/fv_core.res.tile6.nc
  [ "$status" -eq 0 ]

#Combine the file: 
  run command mppnccombine -64 fv_core.res.tile6.nc fv_core.res.tile6.nc.???? 
  [ "$status" -eq 0 ]

#Compare the two:
#  run command nccmp -md fv_core.res.tile6.nc input/fv_core.res.tile6.nc
#  [ "$status" -eq 0 ]

# TO DO: This works on GAEA, but can't run in on travis because the files are too large, 
#		(add ocean_temp_salt.res.nc and ice_model.res.nc) to put on git

#Try with ocean_temp_salt.res.nc
#  run command mppncscatter  -i 7 -j 7 -x 21 -y 14  input/ocean_temp_salt.res.nc
#  [ "$status" -eq 0 ]

#  run command mppnccombine -64 ocean_temp_salt.res.nc ocean_temp_salt.res.nc.???? 
#  [ "$status" -eq 0 ]

#  run command nccmp -md ocean_temp_salt.res.nc input/ocean_temp_salt.res.nc
#  [ "$status" -eq 0 ]

#Try with ice_model.res.nc
#  run command mppncscatter  -i 7 -j 7 -x 21 -y 14  input/ice_model.res.nc
#  [ "$status" -eq 0 ]

#  run command mppnccombine -64 ice_model.res.nc ice_model.res.nc.???? 
#  [ "$status" -eq 0 ]

#  run command nccmp -md ice_model.res.nc input/ice_model.res.nc
#  [ "$status" -eq 0 ]

  cd ..
  rm -rf Test12

}
