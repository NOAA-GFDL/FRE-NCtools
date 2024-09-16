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
if [ -z "$skip_mpi" ]; then
cat <<_EOF > instantaneous.nml
 &input
    file_names(1) =  'timavg_instantaneous_in.nc' ,
    file_name_out =  'timavg_instantaneous_out.nc' ,
 &end
_EOF

cat <<_EOF > standard.nml
 &input
    file_names(1) =  'timavg_standard_in.nc' ,
    file_name_out =  'timavg_standard_out.nc' ,
 &end
_EOF

    python3 $top_srcdir/t/test_time_avg.py
fi
}
