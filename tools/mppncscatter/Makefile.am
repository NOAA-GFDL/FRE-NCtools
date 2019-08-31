bin_PROGRAMS = mppncscatter

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

mppncscatter_SOURCES = main.c \
                       mppncscatter.c \
                       opt.c \
                       scatterdim.c \
                       strlist.c \
                       xmalloc.c
