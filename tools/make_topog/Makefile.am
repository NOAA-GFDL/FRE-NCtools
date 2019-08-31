bin_PROGRAMS = make_topog

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_topog_SOURCES = make_topog.c \
                     topog.c \
                     topog.h
