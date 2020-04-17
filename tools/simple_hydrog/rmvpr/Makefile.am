bin_PROGRAMS = rmv_parallel_rivers

AM_CFLAGS = -I$(top_srcdir)/tools/libfrencutils \
            $(NETCDF_CFLAGS)

AM_FCFLAGS =  -I$(top_srcdir)/tools/libfrencutils \
            $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)

LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) -lnetcdff -lnetcdf

rmv_parallel_rivers_SOURCES = rmv_parallel_rivers.f90  

rmv_parallel_rivers_LDADD = $(top_builddir)/tools/libfrencutils/libfrencutils.a $(LDADD)
 
