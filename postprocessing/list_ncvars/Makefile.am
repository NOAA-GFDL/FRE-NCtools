bin_PROGRAMS = list_ncvars.exe
bin_SCRIPTS = list_ncvars.csh

list_ncvars_exe_SOURCES = list_ncvars.f90
list_ncvars_exe_FCFLAGS = $(BASE_FCFLAGS) $(R8_FCFLAGS) $(AM_FCFLAGS)
list_ncvars_exe_LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

clean-local:
	$(RM) *.mod
