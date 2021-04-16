#!/bin/sh
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
# License along with FRE-NCTools (LICENSE.md).  If not, see
# <http://www.gnu.org/licenses/>.
#***********************************************************************
#
# Copyright (c) 2021 - Seth Underwood (@underwoo)
#
# This script configures the environment using Environment modules
# for building FRE-NCtools.  This script can be run with the `eval`
# command to modify the environment.  Syntax is similar to the
# syntax used in modulefiles.

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DO NOT CHANGE THIS LINE
. $( dirname $( dirname $(readlink -f $0) ) )/env_functions.sh
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module load intel/2019.5
module load netcdf
module load hdf5
module load nccmp

# Set CONFIG_SITE to the correct config.site file for the system
setenv CONFIG_SITE $( dirname $(readlink -f $0) )/config.site
