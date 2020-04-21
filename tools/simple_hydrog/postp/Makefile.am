bin_PROGRAMS = cp_river_vars

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)

AM_FCFLAGS =  -I$(top_srcdir)/tools/libfrencutils \
            $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) -lnetcdff -lnetcdf -lhdf5_hl  -lhdf5 -lz -ldl

cp_river_vars_SOURCES = cp_river_vars.f90  

cp_river_vars_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)
 
