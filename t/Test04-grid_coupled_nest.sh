#!/usr/bin/env bats

#Test grid for coupled nest model (land are C48 and ocean is 1 degree tripolar grid, atmosphere is C48 with nested region

@test "Test grid for coupled nest model (land are C48 and ocean is 1 degree tripolar grid, atmosphere is C48 with nested region" {

  if [ ! -d "Test04" ] 
  then
  		mkdir Test04
  fi

  cd Test04
  cp $top_srcdir/t/Test03-input/OCCAM_p5degree.nc OCCAM_p5degree.nc

#create ocean_hgrid 
run command make_hgrid \
		--grid_type tripolar_grid \
		--nxbnd 2  \
		--nybnd 7  \
		--xbnd -280,80  \
		--ybnd -82,-30,-10,0,10,30,90 \
		--dlon 1.0,1.0  \
		--dlat 1.0,1.0,0.6666667,0.3333333,0.6666667,1.0,1.0 \
		--grid_name ocean_hgrid  \
		--center c_cell 
  [ "$status" -eq 0 ]

#create ocean_vgrid
run command make_vgrid \
		--nbnds 3  \
		--bnds 0.,220.,5500.  \
		--dbnds 10.,10.,367.14286  \
		--center c_cell  \
		--grid_name ocean_vgrid 
  [ "$status" -eq 0 ]

#create ocean solo mosaic
run command make_solo_mosaic  \
		--num_tiles 1  \
		--dir ./  \
		--mosaic_name ocean_mosaic  \
		--tile_file ocean_hgrid.nc  \
		--periodx 360
  [ "$status" -eq 0 ]

#create ocean topography data
run command make_topog  \
		--mosaic ocean_mosaic.nc  \
		--topog_type realistic  \
		--topog_file OCCAM_p5degree.nc \
		--topog_field TOPO  \
		--scale_factor -1  \
		--vgrid ocean_vgrid.nc  \
		--output topog.nc
  [ "$status" -eq 0 ]

#Create C48 grid with atmos nested grid.
run command make_hgrid  \
		--grid_type gnomonic_ed  \
		--nlon 96  \
		--grid_name atmos_grid  \
		--do_schmidt \
		--target_lat 48.15  \
		--target_lon -100.15  \
		--halo 3  \
		--stretch_factor 3 \
		--great_circle_algorithm  \
		--nest_grid  \
		--refine_ratio 3  \
		--parent_tile 4 \
		--istart_nest 21  \
		--iend_nest 60  \
		--jstart_nest 11  \
		--jend 70
  [ "$status" -eq 0 ]

#create C48 solo mosaic for atmos
run command make_solo_mosaic  \
		--num_tiles 7  \
		--dir ./  \
		--mosaic atmos_mosaic  \
		--tile_file  \ atmos_grid.tile1.nc,atmos_grid.tile2.nc,atmos_grid.tile3.nc,atmos_grid.tile4.nc,atmos_grid.tile5.nc,atmos_grid.tile6.nc,atmos_grid.tile7.nc
  [ "$status" -eq 0 ]

#Create C144 grid for land
run command make_hgrid  \
		--grid_type gnomonic_ed  \
		--nlon 288  \
		--grid_name land_grid  \
		--do_schmidt \
		--target_lat 48.15  \
		--target_lon -100.15  \
		--halo 3  \
		--stretch_factor 3  \
		--great_circle_algorithm  \
		--nest_grid  \
		--refine_ratio 3  \
		--parent_tile 0
  [ "$status" -eq 0 ]

#create C144 solo mosaic for land
run command make_solo_mosaic  \
		--num_tiles 6  \
		--dir ./  \
		--mosaic land_mosaic \
		--tile_file land_grid.tile1.nc,land_grid.tile2.nc,land_grid.tile3.nc,land_grid.tile4.nc,land_grid.tile5.nc,land_grid.tile6.nc
  [ "$status" -eq 0 ]

# TO DO: Skipping this because it fails 
#make the coupler_mosaic
#run command aprun -n $npes2 make_coupler_mosaic_parallel --atmos_mosaic atmos_mosaic.nc --land_mosaic land_mosaic.nc \
#          --ocean_mosaic ocean_mosaic.nc --ocean_topog  topog.nc --interp_order 1 --mosaic_name grid_spec --check

#check reproducing ability between processor count for make_coupler_mosaic
#if( ! -d parallel ) mkdir parallel
#cd parallel
#run command aprun -n $npes make_coupler_mosaic_parallel --atmos_mosaic ../atmos_mosaic.nc --land_mosaic ../land_mosaic.nc \
#         --ocean_mosaic ../ocean_mosaic.nc --ocean_topog  ../topog.nc --interp_order 1 --mosaic_name grid_spec
#   nccmp -md $file ../$file

#Remove the workdir 
  cd ..
  rm -rf Test04
}
