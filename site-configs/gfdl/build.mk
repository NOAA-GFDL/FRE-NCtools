# Make macros for the GFDL PP/AN system

MPICC := icc

CFLAGS_SITE := -I/app/mpich2-1.2.1p1/lib -DMPICH_IGNORE_CXX_SEEK
FFLAGS_SITE :=

CLIBS_SITE := -L/app/mpich2-1.2.1p1/lib -lmpich -lopa  -lpthread -lrt -lcurl
FLIBS_SITE := -lcurl

NETCDF_HOME := /usr/local/x64/netcdf-4.2_optimized
HDF5_HOME := /usr/local/x64/hdf5-1.8.8_optimized

STATIC :=
