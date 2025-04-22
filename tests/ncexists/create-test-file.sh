#!/bin/sh

# Copyright (C) 2024 Geophysical Fluid Dynamics Laboratory

# This file is part of the GFDL FRE NetCDF tools package (FRE-NCTools).

# FRE-NCtools is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# FRE-NCtools is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with FRE-NCTools.  If not, see
# <http://www.gnu.org/licenses/>.

# Create the test netCDF file.
ncgen -o test.nc << EOF
netcdf test {
dimensions:
  grid_xt = 2 ;
  grid_yt = 2 ;
  time = UNLIMITED ; // (1 currently)
variables:
  double grid_xt(grid_xt) ;
    grid_xt:long_name = "T-cell longitude" ;
    grid_xt:units = "degrees_E" ;
    grid_xt:axis = "X" ;
  double grid_yt(grid_yt) ;
    grid_yt:long_name = "T-cell latitude" ;
    grid_yt:units = "degrees_N" ;
    grid_yt:axis = "Y" ;
  double time(time) ;
    time:long_name = "time" ;
    time:units = "days since 0001-01-01 00:00:00" ;
    time:axis = "T" ;
    time:calendar_type = "NOLEAP" ;
    time:calendar = "noleap" ;
  float var(grid_yt, grid_xt) ;
    var:long_name = "Decomposed variable" ;
    var:units = "m" ;
    var:missing_value = 1.e+20f ;
    var:_FillValue = 1.e+20f ;
    var:dummy_attribute = "just a dummy attribute" ;
// global attributed:
  :filename = "test.nc" ;
  :NumFilesInSet = 2 ;
  :dummy_global_att = "a global dummy attribute" ;
data:
  grid_xt = 0.0, 1.0 ;
  grid_yt = 0.0, 1.0 ;
  time = 0 ;
  var =
    0, 1, 2, 3 ;
}
EOF

