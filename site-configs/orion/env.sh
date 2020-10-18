# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source $MODULESHOME/init/sh
module load intel/2019.5
module load netcdf
module load hdf5
module load nccmp

# **********************************************************************
# Set environment variables
# **********************************************************************
export LD_RUN_PATH=$LD_LIBRARY_PATH
