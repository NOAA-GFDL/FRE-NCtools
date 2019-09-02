bin_PROGRAMS = make_topog make_topog_parallel

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_topog_SOURCES = make_topog.c \
                     topog.c \
                     topog.h

make_topog_parallel_SOURCES = $(make_topog_SOURCES)
make_topog_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
make_topog_parallel_LDADD = $(LDADD) $(MPI_CLDFLAGS)
