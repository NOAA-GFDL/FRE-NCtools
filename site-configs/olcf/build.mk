# Make macros for the olcf site

MPICC := cc

CFLAGS_SITE := -msse2
FFLAGS_SITE := -msse2

CLIBS_SITE :=
FLIBS_SITE :=

HDF5_HOME := /opt/cray/hdf5/1.8.14/INTEL/140
NETCDF_HOME := /opt/cray/netcdf/4.3.3.1/INTEL/140

STATIC := -static
