bin_PROGRAMS = transfer_to_mosaic_grid

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

transfer_to_mosaic_grid_SOURCES = transfer_to_mosaic.c
