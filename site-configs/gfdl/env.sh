# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
ENV_VERSION=v0.15
INTEL_VERSION=19.0.5

source /usr/local/Modules/default/init/sh
module use /app/spack/default/modulefiles

module load intel_compilers/${INTEL_VERSION}

# gcc>=4.9 is needed for icc to use newer C11 constructs
module load gcc/5.3.0

# Need a newer autoconf/automake than what is curerntly on the system.
module load autoconf automake

# needed for testing
module load nccmp bats

# use intel-compiled libraries
module use /app/spack/${ENV_VERSION}/modulefiles-intel-${INTEL_VERSION}/linux-rhel6-x86_64
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variables
# **********************************************************************
# this simplifies the setup but only works if you source this from the repo root dir
export CONFIG_SITE=`pwd`/site-configs/gfdl/config.site
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
