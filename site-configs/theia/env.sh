# **********************************************************************
# Setup and Load the Modules
# **********************************************************************    
source $MODULESHOME/init/sh
module use -a /home/fms/local/modulefiles
module purge
module load intel/15.1.133
module load impi
module load fre/bronx-9

# **********************************************************************
# Set environment variablesSetup and Load the Modules
# **********************************************************************    
FRE_SYSTEM_SITE=theia
LIBRARY_PATH=${LIBRARY_PATH}:/apps/intel/composer_xe_2015.1.133/compiler/lib/intel64
LD_RUN_PATH=${LD_LIBRARY_PATH}:/apps/hdf5/1.8.14-intel/lib:/apps/netcdf/4.3.0-intel/lib
export FRE_SYSTEM_SITE
export LIBRARY_PATH
export LD_RUN_PATH

# **********************************************************************
# Aliases
# **********************************************************************    
alias make="make HDF5_HOME=/apps/hdf5/1.8.14-intel NETCDF_HOME=/apps/netcdf/4.3.0-intel SITE=theia LDFLAGS='-Wl,-rpath=/apps/intel/composer_xe_2015.1.133/compiler/lib/intel64,-rpath=/apps/hdf5/1.8.14-intel/lib,-rpath=/apps/netcdf/4.3.0-intel/lib' -f fre-nctools.mk"

# **********************************************************************
# Other build configuration settings
# **********************************************************************    
