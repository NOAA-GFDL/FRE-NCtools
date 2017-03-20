# Make macros for the theia site

MPICC := mpiicc

CFLAGS_SITE :=
FFLAGS_SITE :=

CLIBS_SITE :=
FLIBS_SITE :=

LDFLAGS='-Wl,-rpath=/apps/intel/composer_xe_2015.1.133/compiler/lib/intel64,-rpath=/apps/hdf5/1.8.14-intel/lib,-rpath=/apps/netcdf/4.3.0-intel/lib'

NETCDF_HOME := /apps/netcdf/4.3.0-intel 
HDF5_HOME := /apps/hdf5/1.8.14-intel

STATIC :=
