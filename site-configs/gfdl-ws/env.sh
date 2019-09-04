# **********************************************************************
# Setup and Load the Modules
# **********************************************************************    
source /usr/local/Modules/default/init/sh
module use -a /home/fms/local/modulefiles
module load netcdf/4.2
# netcdf/4.2 loads intel-11, so unload it now
module unload intel_compilers/11.1.073
module load intel_compilers/15.0.0
module use /home/sdu/publicmodules
module load mpich2/1.5b1

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************    
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH

