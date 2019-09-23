bin_PROGRAMS = runoff_regrid
if WITH_MPI
  bin_PROGRAMS += runoff_regrid_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

runoff_regrid_SOURCES = runoff_regrid.c
runoff_regrid_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)

runoff_regrid_parallel_SOURCES = $(runoff_regrid_SOURCES)
runoff_regrid_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
runoff_regrid_parallel_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
