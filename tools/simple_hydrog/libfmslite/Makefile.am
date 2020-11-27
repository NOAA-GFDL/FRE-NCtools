noinst_LIBRARIES = libfms.a

AM_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG)

libfms_a_SOURCES = constants.F90 \
  horiz_interp.F90 \
  horiz_interp_bilinear.F90 \
  horiz_interp_type.F90 \
  mpp.F90 \
  mpp_error_a_a.h \
  mpp_error_a_s.h \
  mpp_error_s_a.h \
  mpp_error_s_s.h

mpp.$(OBJEXT): mpp_error_a_a.h mpp_error_a_s.h mpp_error_s_a.h mpp_error_s_s.h
horiz_interp_bilinear.$(OBJEXT): constants.$(OBJEXT) mpp.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp.$(OBJEXT): constants.$(OBJEXT) mpp.$(OBJEXT) horiz_interp_type.$(OBJEXT) horiz_interp_bilinear.$(OBJEXT)

CLEANFILES = *.mod
