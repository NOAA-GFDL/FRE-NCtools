bin_PROGRAMS = list_ncvars.exe
dist_bin_SCRIPTS = list_ncvars.csh

list_ncvars_exe_SOURCES = list_ncvars.f90
list_ncvars_exe_FCFLAGS = $(NETCDF_FCFLAGS) $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)
list_ncvars_exe_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

clean-local:
	$(RM) *.mod
