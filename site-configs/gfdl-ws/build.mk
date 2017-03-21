# Make macros for the GFDL workstations

CFLAGS_SITE :=
FFLAGS_SITE :=

CLIBS_SITE :=
FLIBS_SITE :=

NETCDF_HOME := /usr/local/x64/netcdf-4.2_optimized
HDF5_HOME := /usr/local/x64/hdf5-1.8.8_optimized

STATIC :=

# NOPARALLEL controls if some tools also build the parallel (MPI) version of the
# executable.  If NOPARALLEL is non-blank, then the parallel version will not be
# built.
NOPARALLEL := t
