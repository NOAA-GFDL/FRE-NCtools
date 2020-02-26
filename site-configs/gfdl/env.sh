# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source /usr/local/Modules/default/init/sh
module load intel_compilers/18.0.5 mpich2/1.2.1p1 netcdf/4.2

# Need a newer autoconf/automake than what is curerntly on the system.
# This is needed until there is a module for these.
PATH=/app/spack/linux-rhel6-x86_64/gcc-4.4.7/autoconf/2.69-oz56rghtag37aapzrgmh2tzhyjek7o6c/bin:/app/spack/linux-rhel6-x86_64/gcc-4.4.7/automake/1.16.1-shm646foq3qrojtsqhc6uat77auqexx7/bin:${PATH}
# Add bats to PATH
PATH=${PATH}:/home/sdu/opt/bats/bin
export PATH

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
