# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source /usr/local/Modules/default/init/sh

module load intel_compilers/19.0.5

# gcc is needed for icc to use newer C11 constructs
# we want the spack-installed gcc/4.8.5 but configure fails
### checking whether /app/spack/v0.15/linux-rhel6-x86_64/gcc-4.4.7/gcc/4.8.5-uz2ynmmyxq4tsrf2tuyewx7dow4iszsx/bin/gfortran accepts -g... yes
### checking whether /app/spack/v0.15/linux-rhel6-x86_64/gcc-4.4.7/gcc/4.8.5-uz2ynmmyxq4tsrf2tuyewx7dow4iszsx/bin/gfortran understands -c and -o together... no
### checking for Fortran flag to compile .f90 files... unknown
### configure: error: Fortran could not compile .f90 files
#module load gcc/4.8.5
module load gcc/5.3.0

# Need a newer autoconf/automake than what is curerntly on the system.
module load autoconf automake

# only needed for testing
module load nccmp bats

# use intel-compiled libraries
module use /app/spack/v0.15/modulefiles-intel-19.0.5/linux-rhel6-x86_64
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
# this simplifies the setup but only works if you source this from the repo root dir
export CONFIG_SITE=`pwd`/site-configs/gfdl/config.site
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
