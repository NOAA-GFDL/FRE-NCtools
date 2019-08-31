bin_PROGRAMS = make_hgrid

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

make_hgrid_SOURCES = create_conformal_cubic_grid.c \
                     create_gnomonic_cubic_grid.c \
                     create_grid_from_file.c \
                     create_lonlat_grid.c \
                     make_hgrid.c
