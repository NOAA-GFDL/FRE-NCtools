bin_PROGRAMS = cr_lake_files

# Location of FMS library.
# TODO: allow builder to override and use a different libfms.
FMSROOT = $(top_builddir)/tools/simple_hydrog/libfmslite

AM_CFLAGS = $(NETCDF_CFLAGS)
AM_FCFLAGS = -I$(FMSROOT) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(NETCDF_FCLIBS)

cr_lake_files_SOURCES = cr_lake_files.f90
cr_lake_files_LDADD = $(FMSROOT)/libfms.a $(LDADD)

CLEANFILES = *.mod
