# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
ENV_VERSION=v0.15
INTEL_VERSION=19.0.5

source /usr/local/Modules/default/init/sh

module load intel_compilers/${INTEL_VERSION}

# gcc is needed for icc to use newer C11 constructs
module load gcc/9.2.0

# nccmp and bats needed only for testing
module load nccmp bats

# use intel-compiled libraries
module use /app/spack/${ENV_VERSION}/modulefiles-intel-${INTEL_VERSION}/linux-rhel7-x86_64
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variables
# **********************************************************************
# this simplifies the setup but only works if you source this from the repo root dir
export CONFIG_SITE=`pwd`/site-configs/gfdl-ws/config.site
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
