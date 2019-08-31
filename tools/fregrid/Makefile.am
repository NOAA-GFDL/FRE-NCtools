bin_PROGRAMS = fregrid

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
