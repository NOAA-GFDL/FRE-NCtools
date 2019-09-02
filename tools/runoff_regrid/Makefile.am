bin_PROGRAMS = runoff_regrid runoff_regrid_parallel

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

runoff_regrid_SOURCES = runoff_regrid.c

runoff_regrid_parallel_SOURCES = $(runoff_regrid_SOURCES)
runoff_regrid_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
runoff_regrid_parallel_LDADD = $(LDADD) $(MPI_CLDFLAGS)
