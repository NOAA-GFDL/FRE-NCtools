bin_PROGRAMS = cr_lake_files

# Location of FMS library.
# TODO: allow builder to override and use a different libfms.
FMSROOT = $(top_builddir)/tools/simple_hydrog/libfmslite

LDADD = $(NETCDF_LDFLAGS) -lnetcdff $(NETCDF_LIBS) -lhdf5_hl -lhdf5 -lz -ldl
AM_FCFLAGS = -I$(FMSROOT) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)

cr_lake_files_SOURCES = cr_lake_files.f90
cr_lake_files_LDADD = $(FMSROOT)/libfms.a $(LDADD)

CLEANFILES = *.mod
