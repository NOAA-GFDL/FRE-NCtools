# **********************************************************************
# Setup and Load the Modules
# **********************************************************************    
source /usr/local/Modules/default/init/sh
module use -a /home/fms/local/modulefiles
module load intel_compilers/15.0.1 mpich2/1.2.1p1 netcdf/4.2

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************    
FRE_SYSTEM_SITE=gfdl
LD_RUN_PATH=$LD_LIBRARY_PATH
export FRE_SYSTEM_SITE
export LD_RUN_PATH

# **********************************************************************
# Aliases
# **********************************************************************    
alias make="make HDF5_HOME=/usr/local/hdf5-1.8.8_optimized NETCDF_HOME=/usr/local/netcdf-4.2_optimized SITE=pan -f fre-nctools.mk"

# **********************************************************************
# Other build configuration settings
# **********************************************************************    
