# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
modpath_prepend ()
{
  local myprepend=$1
  if $( module is-used $myprepend )
  then
    module remove-path MODULEPATH $myprepend
    module prepend-path MODULEPATH $myprepend
  else
    module use $myprepend
  fi
}

env_version=v0.15
intel_version=18.0.5

# Ensure the module environment is initialized
source /usr/local/Modules/default/init/sh

# Ensure the base spack modules are first in MODULEPATH
modpath_prepend /app/spack/${env_version}/modulefiles/linux-rhel7-x86_64
# GCC is needed for icc to use newer C11 constructs
module load gcc/9.2.0
# bats is needed for tests
module load bats/0.4.0

# Load the Intel compilers
module load intel_compilers/${intel_version}

# Ensure the Intel spack modules are first in MODULEPATH
modpath_prepend /app/spack/${env_version}/modulefiles-intel-${intel_version}/linux-rhel7-x86_64
# Load the Intel modules required for building
module load netcdf-c/4.7.3
module load netcdf-fortran/4.5.2
module load mpich/3.3.2

# **********************************************************************
# Set environment variables
# **********************************************************************
# None needed at this time

# **********************************************************************
# Clear temporary variables
# **********************************************************************
unset env_version
unset intel_version
unset -f modpath_prepend
