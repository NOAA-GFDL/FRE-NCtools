bin_PROGRAMS = mppncscatter

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

mppncscatter_SOURCES = common.h \
                       main.c \
                       mppncscatter.c \
                       mppncscatter.h \
                       opt.c \
                       opt.h \
                       scatterdim.c \
                       scatterdim.h \
                       strlist.c \
                       strlist.h \
                       util.h \
                       xmalloc.c \
                       xmalloc.h
