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

# Test procedure for mppncscatter is to first split the files.  Once scattered,
# run mppnccombine to re-combine the files.  The final, recombined file should
# match the unscattered file.

# The mppnccombine and mppncscatter commands should probably be tested in
# the same file, since here we assume mppnccombine is running correctly.
load test_utils

prepare_input_data ()
{
  # Directories for input and final output (re-combined) files
  mkdir input
  mkdir output
  test_file=fv_core.res.tile1.nc
  ncgen -o input/$test_file $top_srcdir/t/Test12-input/${test_file}l
}

@test "Test mppncscatter" {

  prepare_input_data

  # Scatter the file
   mppncscatter -i 2 -j 3 -x 2 -y 12 input/$test_file

  # Combine the file:
   mppnccombine -64 output/$test_file ${test_file}.????

  # Compare output
   run_and_check nccmp -w format -md output/$test_file input/$test_file

}
