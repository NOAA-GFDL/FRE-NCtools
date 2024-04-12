#!/bin/bash
#***********************************************************************
#*                   GNU Lesser General Public License
#*
#* This file is part of the GFDL Flexible Modeling System (FMS).
#*
#* FMS is free software: you can redistribute it and/or modify it under
#* the terms of the GNU Lesser General Public License as published by
#* the Free Software Foundation, either version 3 of the License, or (at
#* your option) any later version.
#*
#* FMS is distributed in the hope that it will be useful, but WITHOUT
#* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#* FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#* for more details.
#*
#* You should have received a copy of the GNU Lesser General Public
#* License along with FMS.  If not, see <http://www.gnu.org/licenses/>.
#***********************************************************************

#copy python script over
make cp_test_make_remap_file

echo -e "\ntest reading remap file for conserve_order1\n"
./test_make_remap_file_conserve.py 1
./test_read_remap_file 1

echo -e "\ntest reading remap file for conserve_order2\n"
./test_make_remap_file_conserve.py 2
./test_read_remap_file 2
