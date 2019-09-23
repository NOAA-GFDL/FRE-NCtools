bin_PROGRAMS = ncexists

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
				$(NETCDF_CFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

ncexists_SOURCES = ncexists.c
