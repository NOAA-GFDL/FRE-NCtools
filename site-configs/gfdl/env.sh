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
env_version=v0.15
intel_version=18.0.5
gcc_version=9.2.0
ncc_version=4.7.3
ncf_version=4.5.2
mpi_version=3.3.2

# Ensure the base spack modules are first in MODULEPATH
modpath_prepend /app/spack/${env_version}/modulefiles/linux-rhel7-x86_64
# GCC is needed for icc to use newer C11 constructs
module load gcc/$gcc_version
# bats is needed for tests
module load bats/0.4.0
# Need newer autoconf/automake than what pan has at the system level
module load autoconf/2.69
module load automake/1.16.2

# Load the Intel compilers
module load intel_compilers/${intel_version}

# Ensure the Intel spack modules are first in MODULEPATH
modpath_prepend /app/spack/${env_version}/modulefiles-intel-${intel_version}/linux-rhel7-x86_64
# Load the Intel modules required for building
module load netcdf-c/$ncc_version
module load netcdf-fortran/$ncf_version
module load mpich/$mpi_version

# Set CONFIG_SITE to the correct config.site file for the system
setenv CONFIG_SITE=$( dirname $(readlink -f $0) )/config.site
