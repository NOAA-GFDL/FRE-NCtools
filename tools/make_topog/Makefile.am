bin_PROGRAMS = make_topog
if WITH_MPI
  bin_PROGRAMS += make_topog_parallel
endif

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_topog_SOURCES = make_topog.c \
                     topog.c \
                     topog.h
make_topog_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)

make_topog_parallel_SOURCES = $(make_topog_SOURCES)
make_topog_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
make_topog_parallel_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils_mpi.a $(LDADD) $(MPI_CLDFLAGS)
