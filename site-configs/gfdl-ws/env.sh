# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source /usr/local/Modules/default/init/sh

module load intel_compilers/19.0.5

# gcc is needed for icc to use newer C11 constructs
# we want the spack-installed gcc/9.2.0 but configure fails:
### checking whether /app/spack/v0.15/linux-rhel7-x86_64/gcc-4.8.5/gcc/9.2.0-jlfyxwjebljjliixtaduabmcnh7zeg3l/bin/gfortran understands -c and -o together... no
### checking for Fortran flag to compile .f90 files... unknown
### configure: error: Fortran could not compile .f90 files
#module load gcc/9.2.0
module load gcc/6.2.0

# nccmp and bats needed only for testing
module load nccmp bats

# use intel-compiled libraries
module use /app/spack/v0.15/modulefiles-intel-19.0.5/linux-rhel7-x86_64
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
# this simplifies the setup but only works if you source this from the repo root dir
export CONFIG_SITE=`pwd`/site-configs/gfdl-ws/config.site
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
