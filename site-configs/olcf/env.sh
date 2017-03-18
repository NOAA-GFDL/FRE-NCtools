# **********************************************************************
# Setup and Load the Modules
# **********************************************************************    
source $MODULESHOME/init/sh
module rm PrgEnv-pgi PrgEnv-intel PrgEnv-gnu PrgEnv-cray
module load PrgEnv-intel
module swap intel intel/15.0.2.164
module load cray-netcdf/4.3.3.1

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************    
FRE_SYSTEM_SITE=olcf
export FRE_SYSTEM_SITE

# **********************************************************************
# Aliases
# **********************************************************************    
alias make="make HDF5_HOME=/opt/cray/hdf5/1.8.14/INTEL/140 NETCDF_HOME=/opt/cray/netcdf/4.3.3.1/INTEL/140 NC_BLKSZ=64K SITE=olcf -f fre-nctools.mk"

# **********************************************************************
# Other build configuration settings
# **********************************************************************    
