bin_PROGRAMS = make_solo_mosaic

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_solo_mosaic_SOURCES = get_contact.c \
                           make_solo_mosaic.c
