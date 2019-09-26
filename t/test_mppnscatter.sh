#!/usr/bin/env bats
# mppncscatter first use mppncscatter to split the file
#             and then use mppnccombine to combine the file. The output should reproduce
#             original file

@test "Test mppncscatter and  mppnccombine" {
  mkdir work_dir_11
  cd work_dir_11
  mkdir input
  cd input
  cp $top_srcdir/t/mppncscatter/* .
  cd ..

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

# To do: This works on GAEA, but can't run in on travis because the files are too large 
#		(add ocean_temp_salt.res.nc and ice_model.res.nc

  cd ..
  rm -rf work_dir_11

}
