# Make macros for the ncrc4 site

MPICC := cc

CFLAGS_SITE := -msse2
FFLAGS_SITE := -msse2

CLIBS_SITE :=
FLIBS_SITE :=

NETCDF_HOME := $NETCDF_DIR
HDF5_HOME := $HDF5_DIR

STATIC := -static
