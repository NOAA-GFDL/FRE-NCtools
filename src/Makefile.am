BUILT_SOURCES =

combine_programs = \
  combine-ncc \
  is-compressed \
  mppnccombine \
  scatter-ncc

bin_PROGRAMS = \
  decompress-ncc \
  $(combine_programs)

bin_SCRIPTS = \
  combine_restarts \
  iceberg_comb.sh \
  list_ncvars.csh \
  list_ncvars.sh \
  plevel.sh \
  timavg \
  timavg.csh

if HAVE_PLMOD_LIST_MOREUTILS
  bin_SCRIPTS += split_ncvars.pl
endif

pkglibexec_PROGRAMS = \
  list_ncvars \
  print_version \
  plevel-interp \
  time_average

noinst_LIBRARIES = \
  libver.a \
  libnctf.a

libver_a_SOURCES = \
  version.c
nodist_libver_a_SOURCES = \
  version.h

libnctf_a_SOURCES = \
  nct_nc_utils.F90
nodist_libnctf_a_SOURCES = \
  nct_nc_utils.$(FC_MODEXT)
BUILT_SOURCES += nct_nc_utils.$(FC_MODEXT)

c_netcdf_ldadds = $(NETCDF_LDFLAGS) $(NETCDF_LIBS)
f_netcdf_ldadds = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS) libnctf.a

AM_FCFLAGS = $(FC_MODOUT)$(MODDIR)$(@D) $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
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

EXTRA_DIST = \
  combine-restarts/combine_restarts.in \
  iceberg-comb/iceberg_comb.sh.in \
  list-ncvars/list_ncvars.csh.in \
  list-ncvars/list_ncvars.csh.in \
  plevel/plevel.sh.in \
  split_ncvars/split_ncvars.pl.in \
  timavg/timavg.csh.in \
  timeavg/timavg.sh.in

CLEANFILES = *.$(FC_MODEXT)
clean-local:
	$(RM) $(bin_SCRIPTS) $(BUILT_SOURCES)

# Perform several substitutions
do_subst = $(SED) -e 's,[@]prefix[@],$(prefix),g' \
               -e 's,[@]pkglibexecdir[@],$(pkglibexecdir),g' \
               -e 's,[@]NCDUMP[@],$(NCDUMP),g' \
               -e 's,[@]NCRCAT[@],$(NCRCAT),g' \
               -e 's,[@]NCATTED[@],$(NCATTED),g' \
               -e 's,[@]NCKS[@],$(NCKS),g' \
               -e 's,[@]SED[@],$(SED),g' \
               -e 's,[@]EGREP[@],$(EGREP),g'

# Items to handle Fortran module files
SUFFIXES = _mod.$(FC_MODEXT) .$(FC_MODEXT)

.F90_mod.$(FC_MODEXT):
	@$(MKDIR_P) $(@D)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@D)/$(<F:.F90=.o) $<

.F90.$(FC_MODEXT):
	@$(MKDIR_P) $(@D)
	$(AM_V_PPFC)$(PPFCCOMPILE) -c -o $(@D)/$(<F:.F90=.o) $<

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

# ********************************************************* #
# list_ncvars                                               #
# ********************************************************* #
list_ncvars_SOURCES = list-ncvars/list_ncvars.f90
list_ncvars_LDADD = $(f_netcdf_ldadds)

list_ncvars.sh: list-ncvars/list_ncvars.sh.in Makefile
	$(do_subst) < $(srcdir)/list-ncvars/list_ncvars.sh.in > list_ncvars.sh
	chmod +x list_ncvars.sh

list_ncvars.csh: list-ncvars/list_ncvars.csh.in Makefile
	$(do_subst) < $(srcdir)/list-ncvars/list_ncvars.csh.in > list_ncvars.csh
	chmod +x list_ncvars.csh

# ********************************************************* #
# Combine Scripts                                           #
# ********************************************************* #
combine_restarts: combine-restarts/combine_restarts.in
	$(do_subst) < $(srcdir)/combine-restarts/combine_restarts.in > combine_restarts
	chmod +x combine_restarts

iceberg_comb.sh: iceberg-comb/iceberg_comb.sh.in Makefile
	$(do_subst) < $(srcdir)/iceberg-comb/iceberg_comb.sh.in > iceberg_comb.sh
	chmod +x iceberg_comb.sh

# ********************************************************* #
# plevel                                                    #
# ********************************************************* #
plevel.sh: plevel/plevel.sh.in Makefile
	$(do_subst) < $(srcdir)/plevel/plevel.sh.in > plevel.sh
	chmod +x plevel.sh

plevel_interp_SOURCES = plevel/plev_constants.F90 \
                        plevel/moisture_convert.F90 \
                        plevel/pressure_interp.F90 \
                        plevel/pinterp_utilities.F90 \
                        plevel/run_pressure_interp.F90
nodist_plevel_interp_SOURCES = \
                        plevel/plev_constants_mod.$(FC_MODEXT) \
                        plevel/moisture_convert_mod.$(FC_MODEXT) \
                        plevel/pressure_interp_mod.$(FC_MODEXT) \
                        plevel/pinterp_utilities_mod.$(FC_MODEXT)

plevel_interp_LDADD = $(f_netcdf_ldadds)
BUILT_SOURCES += plevel/plev_constants_mod.$(FC_MODEXT) \
                 plevel/moisture_convert_mod.$(FC_MODEXT) \
                 plevel/pressure_interp_mod.$(FC_MODEXT) \
                 plevel/pinterp_utilities_mod.$(FC_MODEXT)
plevel/moisture_convert_mod.$(FC_MODEXT): plevel/plev_constants_mod.$(FC_MODEXT)
plevel/pressure_interp_mod.$(FC_MODEXT):  plevel/plev_constants_mod.$(FC_MODEXT) plevel/moisture_convert_mod.$(FC_MODEXT)
plevel/plev_constants.o: plevel/plev_constants_mod.$(FC_MODEXT)
plevel/moisture_convert.o: plevel/moisture_convert_mod.$(FC_MODEXT)
plevel/pressure_interp.o: plevel/pressure_interp_mod.$(FC_MODEXT)
plevel/pinterp_utilities.o: plevel/pinterp_utilities_mod.$(FC_MODEXT)
plevel/run_pressure_interp.o: plevel/plev_constants_mod.$(FC_MODEXT) plevel/pressure_interp_mod.$(FC_MODEXT) plevel/pinterp_utilities_mod.$(FC_MODEXT) nct_nc_utils.$(FC_MODEXT)

# ********************************************************* #
# time_average                                              #
# ********************************************************* #
time_average_SOURCES = timeavg/time_average.f90
time_average_LDADD = $(f_netcdf_ldadds)
timeavg/time_average.o: nct_nc_utils.$(FC_MODEXT)
timavg: timeavg/timavg.sh.in Makefile
	$(do_subst) < $(srcdir)/timeavg/timavg.sh.in > timavg
	chmod +x timavg
timavg.csh: timeavg/timavg.csh.in Makefile
	$(do_subst) < $(srcdir)/timeavg/timavg.csh.in > timavg.csh
	chmod +x timavg.csh


# ********************************************************* #
# split_ncvars.pl                                           #
# ********************************************************* #
split_ncvars.pl: split_ncvars/split_ncvars.pl.in Makefile
	$(do_subst) < $(srcdir)/split_ncvars/split_ncvars.pl.in > split_ncvars.pl
	chmod +x split_ncvars.pl
