bin_PROGRAMS = PLEV.exe PLEV.r4.exe
bin_SCRIPTS = plevel.sh

AM_FCFLAGS =  $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

PLEV_exe_SOURCES = plev_constants.F90 moisture_convert.F90 pressure_interp.F90 pinterp_utilities.F90 run_pressure_interp.F90
PLEV_exe_FCFLAGS = $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)

PLEV_r4_exe_SOURCES = plev_constants.F90 moisture_convert.F90 pressure_interp.F90 pinterp_utilities.F90 run_pressure_interp.F90
PLEV_r4_exe_FCFLAGS = $(FC_REAL_KIND4_FLAG) $(AM_FCFLAGS)

clean-local:
	$(RM) *.mod
