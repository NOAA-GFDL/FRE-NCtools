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

intel_version=2019.5
ncc_version=4.7.4
hdf_version=1.10.6

# GCC is needed for icc to use newer C11 constructs.  However, on systems
# that use lmod, we cannot load two compiler family modules concurrently.
# This should add the GNU path information to the *PATH variables
module prepend-path PATH /apps/gcc-8/gcc-8.3.0//bin
module prepend-path LD_LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/gmp-6.1.2/lib
module prepend-path LD_LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/mpc-1.1.0/lib
module prepend-path LD_LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/mpfr-4.0.2/lib
module prepend-path LD_LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/lib
module prepend-path LD_LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/lib64
module prepend-path LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/lib
module prepend-path LIBRARY_PATH /apps/gcc-8/gcc-8.3.0/lib64
module prepend-path CPATH /apps/gcc-8/gcc-8.3.0/include

module load intel/${intel_version}
module load netcdf/${ncc_version}
module load hdf5/${hdf_version}

# Set CONFIG_SITE to the correct config.site file for the system
setenv CONFIG_SITE $( dirname $(readlink -f $0) )/config.site

# Include the netcdf-c/netcdf-fortran library paths during linking
setenv LD_RUN_PATH \$LD_LIBRARY_PATH
