bin_PROGRAMS = nc_null_check

AM_FCFLAGS = -I$(top_srcdir)/tools/libfrencutils \
				$(NETCDF_FCFLAGS)
LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a \
        $(NETCDF_LDFLAGS) $(NETCDF_FCLIBS)

nc_null_check_SOURCES = nc_null_check.F90
