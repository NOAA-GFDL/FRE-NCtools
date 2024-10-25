BUILT_SOURCES =

combine_programs = \
  combine-ncc \
  is-compressed \
  mppnccombine \
  scatter-ncc

bin_PROGRAMS = \
  decompress-ncc \
  $(combine_programs)

noinst_LIBRARIES = libver.a

libver_a_SOURCES = \
  version.h \
  version.c

noinst_HEADERS = \
  version.h \
  version.inc

c_netcdf_ldadds = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
f_netcdf_ldadds = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

AM_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
AM_CFLAGS = $(NETCDF_CFLAGS)
LDADD = libver.a $(NETCDF_LDFLAGS) $(RPATH_FLAGS)

BUILT_SOURCES += version.inc
version.inc: Makefile
	$(AM_V_GEN)rm -f $@
	$(AM_V_at)printf 'interface\n' > $@t
	$(AM_V_at)printf '  subroutine print_version (command_name) bind(c)\n' >> $@t
	$(AM_V_at)printf '    use iso_c_binding\n' >> $@t
	$(AM_V_at)printf '    character(kind=c_char), intent(in) :: command_name(*)\n' >> $@t
	$(AM_V_at)printf '  end subroutine print_version\n' >> $@t
	$(AM_V_at)printf 'end interface\n\n' >>$@t
	$(AC_V_at)printf 'INTEGER, PARAMETER :: EXIT_SUCCESS = 0\n' >> $@t
	$(AC_V_at)printf 'INTEGER, PARAMETER :: EXIT_FAILURE = 1\n' >> $@t
	$(AM_V_at)chmod a-w $@t
	$(AM_V_at)mv $@t $@

BUILT_SOURCES += version.h
version.h: Makefile
	$(AM_V_GEN)rm -f $@
	$(AM_V_at)printf 'extern char const *Version;\n' > $@t
	$(AM_V_at)printf 'extern void print_version(const char *command_name);\n' >> $@t
	$(AC_V_at)chmod a-w $@t
	$(AM_V_at)mv $@t $@

DISTCLEANFILES = $(BUILT_SOURCES)
MAINTAINERCLEANFILES = $(BUILT_SOURCES)

# ********************************************************* #
# land_util programs                                        #
# ********************************************************* #
decompress_ncc_SOURCES = land_utils/decompress-ncc.F90
decompress_ncc_FCFLAGS = -I../lib/libnfu $(AM_FCFLAGS)
decompress_ncc_LDADD = $(f_netcdf_ldadds) $(top_builddir)/lib/libnfu/libnfu.a $(LDADD)

combine_ncc_SOURCES = land_utils/combine-ncc.F90
combine_ncc_FCFLAGS = -I../lib/libnfu $(AM_FCFLAGS)
combine_ncc_LDADD =  $(f_netcdf_ldadds) $(top_builddir)/lib/libnfu/libnfu.a $(LDADD)

is_compressed_SOURCES = land_utils/is-compressed.c
is_compressed_LDADD = $(c_netcdf_ldadds) $(LDADD)

scatter_ncc_SOURCES = land_utils/scatter-ncc.F90
scatter_ncc_FCFLAGS = -I../lib/libnfu $(AM_FCFLAGS)
scatter_ncc_LDADD = $(f_netcdf_ldadds) $(top_builddir)/lib/libnfu/libnfu.a $(LDADD)

# ********************************************************* #
# mppnccombine                                              #
# ********************************************************* #
mppnccombine_SOURCES = mpp-nccombine/mppnccombine.c
mppnccombine_LDADD = $(c_netcdf_ldadds) -lm $(LDADD)