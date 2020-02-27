# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source $MODULESHOME/init/sh
module load intel/2019.5
module load netcdf
module load hdf5

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
# Setting LD_RUN_PATH seems to make static executables that don't require 
# libraries to run. This is generally desired but doesn't work on orion so far.
# FRE loads intel/2019.5 (in order to load nco, which is in the lmod hierarchy)
# so the nctools works, but this is unoptimal
#LD_RUN_PATH=$LD_LIBRARY_PATH
#export LD_RUN_PATH
