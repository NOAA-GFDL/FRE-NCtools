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

@test "Test make_hgrid do_transform vs. FV3 generated grid spec" {

fv3_grid_name="fv3_48_grid_spec"   #name prefix for the internal fv3 grid spec files.
nct_ff_grid_name="C48_ff_grid_spec"   #name prefix for NCTools generated files from fv3 files
nct_grid_name="C48_grid_spec"   #name prefix for "the common analytic" make_hgrid generated files.
fv3_grids_filelist=""


 if [ ! -d "Test33" ]
  then
      mkdir Test33
  fi
 cd Test33

#0) Verify the existance of the FV3 grid spec files
for i in 1 2 3 4 5 6
do
    fv3_file=$fv3_grid_name".tile"$i".nc"
   run_and_check [ -e $top_srcdir/t/Test33-reference/$fv3_file ]
done

cp $top_srcdir/t/Test33-reference/*.nc .

# I) !Make the NCTools grid files.
#Ia) Prepare a list of FV3 files to pass to make_hgrid. These files
# can be created from SHiELD/FV3 runs with the do_cube, then (if neccesary)
# to modify the FV3 files with "ncks --fl_fmt=netcdf4_classic"
fv3_grids_filelist=$(get_csv_filename_list .  $fv3_grid_name"*.nc")

echo "about to call make_hgrid with file list: "
echo $fv3_grids_filelist
echo ""
#Ib) Call make_hgrid to generate the NCTools grids from FV3 files.
## TODO:  Note argument of --great_circle_algorithm option, which should not
## be necessary and code may be corrected in future.
make_hgrid --grid_type from_file \
	   --great_circle_algorithm \
	   --my_grid $fv3_grids_filelist \
	   --grid_name $nct_ff_grid_name

#II) "Analytically" create a equal distance gnomonic cubic grid using
#     the do_cube_transform option
make_hgrid \
    --grid_type gnomonic_ed \
    --nlon 96 \
    --do_cube_transform \
    --stretch_factor 3 \
    --target_lat 17.5 \
    --target_lon 172.5 \
    --grid_name $nct_grid_name


# III)  Compare the six tile files generated "analytically" to the corresponding ones
##    generated from FV3 grids. Note that we are comparing files with origins from
## two different systems, one Fortran based and another C based, and exact matching
# will generally not be possible. The two tolerances chosen below were determined
# by running on hardware, with FV3 origin files generated on Intel Skylake (C4) and
# the NCTools analytically generated files created on various AMD and Intel hardware
# (e.g. Intel analysis nodes, AMD T5 nodes, various home computers).
for i in 1 2 3 4 5 6
do
    fv3_file=$nct_ff_grid_name".tile"$i".nc"
    nct_file=$nct_grid_name".tile"$i".nc"
    run_and_check nccmp -d --variable=x,y,dx,dy --Tolerance=1.0e-9 $fv3_file $nct_file
    run_and_check nccmp -d --variable=area --Tolerance=1.0e-6 $fv3_file $nct_file
    # TODO: angle_dx and angle_dy may be done in future.
done

}
