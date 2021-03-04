# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
ENV_VERSION=v0.15
INTEL_VERSION=18.0.5

source /usr/local/Modules/default/init/sh
module use /app/spack/${ENV_VERSION}/modulefiles/linux-rhel7-x86_64
# GCC is needed for icc to use newer C11 constructs
module load gcc/9.2.0
# bats is needed for tests
module load bats/0.4.0

module load intel_compilers/${INTEL_VERSION}
module use /app/spack/${ENV_VERSION}/modulefiles-intel-${INTEL_VERSION}/linux-rhel7-x86_64
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
# The GCC module sets CC and FC to gcc.  Since we want to use icc/ifort
# we will change CC and FC (and a few others)
FC=ifort
CC=icc
F77=ifort
CXX=icpc
export FC CC F77 CXX

