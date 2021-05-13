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

# Variables to control versions used
intel_version=18.0.5.274
gnu_version=6.5.0
ncc_version=4.7.0
mpi_version=2018.4.274

# GCC is needed for icc to use newer C11 constructs.  However, on systems
# that use lmod, we cannot load two compiler family modules concurrently.
# This should add the GNU path information to the *PATH variables
module prepend-path PATH "/apps/gnu/gcc-9.2.0/bin"
module prepend-path MANPATH "/apps/gnu/gcc-9.2.0/share/man"
module prepend-path LD_LIBRARY_PATH "/apps/gnu/gcc-9.2.0/lib"
module prepend-path LIBRARY_PATH "/apps/gnu/gcc-9.2.0/lib"
module prepend-path LD_LIBRARY_PATH "/apps/gnu/gcc-9.2.0/lib64"
module prepend-path LIBRARY_PATH "/apps/gnu/gcc-9.2.0/lib64"
module prepend-path CPATH "/apps/gnu/gcc-9.2.0/include"
module prepend-path CMAKE_PREFIX_PATH "/apps/gnu/gcc-9.2.0/"

# Load the Intel compilers
module load intel/${intel_version}

# bats and nccmp are needed for tests
module prepend-path PATH /home/Seth.Underwood/opt/bats/0.4.0/bin
module load nccmp/1.8.5

# Load the Intel modules required for building
module load netcdf/$ncc_version
module load impi/$mpi_version

# Set CONFIG_SITE to the correct config.site file for the system
setenv CONFIG_SITE $( dirname $(readlink -f $0) )/config.site

# Include the netcdf-c/netcdf-fortran library paths during linking
setenv LD_RUN_PATH \$LD_LIBRARY_PATH
