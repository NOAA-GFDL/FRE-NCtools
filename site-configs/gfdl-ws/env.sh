# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source /usr/local/Modules/default/init/sh
module use /app/spack/default/modulefiles
module load intel_compilers/18.0.5
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2
# gcc is needed for icc to use newer C11 constructs
module load gcc/6.2.0
# nccmp and bats needed only for testing
module load nccmp bats

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
# this simplifies the setup but only works if you source this from the repo root dir
export CONFIG_SITE=`pwd`/site-configs/gfdl-ws/config.site
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
