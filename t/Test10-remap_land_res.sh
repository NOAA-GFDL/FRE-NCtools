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

# Test remap_land: remap land restart files.

# Prepare takes a directory name, and generates the required netCDF input
# files for a given test.  prepare_input_data expects a single argument <dir_name>
# which is the directory name in $top_srcdir/t/Test10-input/<dir_name>
# that contains the input data

load test_utils

prepare_input_data ()
{
  local inputdir=$top_srcdir/t/Test10-input
  for dir in $( ls -1 $inputdir )
  do
    local datadir=$inputdir/$dir
    mkdir -p $dir
    # Ensure the directory is empty
    rm -rf $dir/*
    # Generate the restart files
    for file in $( ls -1 $datadir/*.ncl )
    do
      ncgen -o $dir/$( basename ${file} .ncl ).nc  ${file}
    done
  done
}

@test "Test remap_land can remap land restart files" {

  prepare_input_data

   remap_land \
		--file_type land  \
		--src_mosaic C48_mosaic/C48_mosaic.nc \
		--dst_mosaic C192_mosaic/C192_mosaic.nc \
		--src_restart src_restart/land.res \
		--dst_restart land.res \
		--dst_cold_restart dst_cold_restart/land.res \
		--remap_file remap_file_C48_to_C192 --print_memory
}
