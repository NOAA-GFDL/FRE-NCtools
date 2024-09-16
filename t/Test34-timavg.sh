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

@test "Test timavg" {
cat <<_EOF > input.nml
 &input
    file_names(1) =  '20120101.ice_shelf.nc' ,
    file_name_out =  'new.nc' ,
    use_end_time  =  .false. ,
    verbose  =  .false. ,
    add_cell_methods =  .false. ,
    skip_tavg_errors =  .false. ,
    suppress_warnings =  .false.,
 &end
_EOF
    python3 $top_srcdir/t/test_time_avg.py
}
