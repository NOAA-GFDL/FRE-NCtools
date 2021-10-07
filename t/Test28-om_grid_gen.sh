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

# Test grid for multiple same level and telescoping nests 

@test "Test the ocean_model_grid_generator subproject python app" {

    if [ ! -d "Test28" ] 
  then
    mkdir Test28
  fi
  
  cd Test28
  
  oggappd=$top_srcdir/tools/ocean_model_grid_generator
  
  cp -p $oggappd/ocean_grid_generator.py .
  cp -p $oggappd/numpypi/numpypi_series.py .
  cp -p $oggappd/numpypi/ignore_this.py .

  chmod ugo+x *.py

  ./ocean_grid_generator.py -f ocean_hgrid_res4.0.nc -r 0.25 --even_j --no_changing_meta

  [ -e ocean_hgrid_res4.0.nc ]

  ./ocean_grid_generator.py -f ocean_hgrid_res1.0.nc -r 1.0  --south_cutoff_row 2 --no_changing_meta

  [ -e ocean_hgrid_res1.0.nc ]

  ./ocean_grid_generator.py -f ocean_hgrid_res0.5.nc -r 2  --no_changing_meta

  [ -e ocean_hgrid_res0.5.nc ]

  # The app may exit with error "Equator is not going to be a u-point" given the test args:
  #run python3 ./ocean_grid_generator.py -f ocean_hgrid_res0.5_equenh.nc -r 2 --south_cutoff_row 130
  #[ -e ocean_hgrid_res0.5_equenh.nc ]

  head -3 $top_srcdir/t/Test28-input/hash.md5 > hash.quick

  md5sum -c hash.quick

  cd ..
  rm -rf Test28
}
