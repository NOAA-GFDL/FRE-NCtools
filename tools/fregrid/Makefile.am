bin_PROGRAMS = fregrid fregrid_parallel

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

fregrid_SOURCES = bilinear_interp.c \
                  bilinear_interp.h \
                  conserve_interp.c \
                  conserve_interp.h \
                  fregrid_util.c \
                  fregrid_util.h \
                  fregrid.c \
                  globals.h

fregrid_parallel_SOURCES = $(fregrid_SOURCES)
fregrid_parallel_CFLAGS = -Duse_libMPI $(MPI_CFLAGS) $(AM_CFLAGS)
fregrid_parallel_LDADD = $(LDADD) $(MPI_CLDFLAGS)
