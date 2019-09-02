bin_PROGRAMS = remap_land remap_land_parallel

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

remap_land_SOURCES = remap_land.c

remap_land_parallel_SOURCES = $(remap_land_SOURCES)
remap_land_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
remap_land_parallel_LDADD = $(LDADD) $(MPI_CLDFLAGS)
