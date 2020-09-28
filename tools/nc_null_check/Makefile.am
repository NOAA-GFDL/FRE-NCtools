bin_PROGRAMS = nc_null_check

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
				$(NETCDF_FFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_LIBS)

nc_null_check_SOURCES = nc_null_check.F90
