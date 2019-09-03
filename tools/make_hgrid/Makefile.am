bin_PROGRAMS = make_hgrid
if WITH_MPI
  bin_PROGRAMS += make_hgrid_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_hgrid_SOURCES = create_conformal_cubic_grid.c \
                     create_gnomonic_cubic_grid.c \
                     create_grid_from_file.c \
                     create_hgrid.h \
                     create_lonlat_grid.c \
                     make_hgrid.c

make_hgrid_parallel_SOURCES = $(make_hgrid_SOURCES)
make_hgrid_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
make_hgrid_parallel_LDADD = $(LDADD) $(MPI_CLDFLAGS)
