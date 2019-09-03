bin_PROGRAMS = make_coupler_mosaic
if WITH_MPI
  bin_PROGRAMS += make_coupler_mosaic_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_coupler_mosaic_SOURCES = make_coupler_mosaic.c
make_coupler_mosaic_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)

make_coupler_mosaic_parallel_SOURCES = $(make_coupler_mosaic_SOURCES)
make_coupler_mosaic_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
make_coupler_mosaic_parallel_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
