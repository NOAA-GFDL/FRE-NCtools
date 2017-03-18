# **********************************************************************
# Setup and Load the Modules
# **********************************************************************    
source $MODULESHOME/init/sh
module purge
module load intel/15.1.133
module load impi

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************    
LIBRARY_PATH=${LIBRARY_PATH}:/apps/intel/composer_xe_2015.1.133/compiler/lib/intel64
LD_RUN_PATH=${LD_LIBRARY_PATH}:/apps/hdf5/1.8.14-intel/lib:/apps/netcdf/4.3.0-intel/lib
export LIBRARY_PATH
export LD_RUN_PATH

# **********************************************************************
# Aliases
# **********************************************************************    

# **********************************************************************
# Other build configuration settings
# **********************************************************************    
