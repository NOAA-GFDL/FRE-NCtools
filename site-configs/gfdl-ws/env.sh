# **********************************************************************
# Setup and Load the Modules
# **********************************************************************
source /usr/local/Modules/default/init/sh
module use -a /home/fms/local/modulefiles
# netcdf/4.2 loads intel-11, so unload it now
module unload netcdf/4.2
module unload hdf5/1.8.8
module unload intel_compilers/11.1.073
module load intel_compilers/18.0.5
module load netcdf/4.2
module load nccmp
module use /home/sdu/publicmodules
module load mpich2/1.5b1
module load gcc/6.2.0

# Add bats to PATH
PATH=${PATH}:/home/sdu/opt/bats/bin
export PATH

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************
LD_RUN_PATH=$LD_LIBRARY_PATH
export LD_RUN_PATH
