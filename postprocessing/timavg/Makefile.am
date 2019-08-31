bin_PROGRAMS = TAVG.exe TAVG.r4.exe
bin_SCRIPTS = timavg.csh

AM_FCFLAGS =  $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_FCLDFLAGS) $(NETCDF_FCLIBS)

TAVG_exe_SOURCES = time_average.f90
TAVG_exe_FCFLAGS = $(FC_REAL_KIND8_FLAG) $(AM_FCFLAGS)

TAVG_r4_exe_SOURCES = time_average.f90
TAVG_r4_exe_FCFLAGS = $(FC_REAL_KIND4_FLAG) $(AM_FCFLAGS)

clean-local:
	$(RM) *.mod
