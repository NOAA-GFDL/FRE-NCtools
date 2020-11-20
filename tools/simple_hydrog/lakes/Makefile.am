bin_PROGRAMS = cr_lake_files

AM_CFLAGS = $(NETCDF_CFLAGS)
AM_FCFLAGS = $(FC_DEFAULT_REAL_KIND8_FLAG) $(NETCDF_FCFLAGS)
LDADD = $(NETCDF_LDFLAGS) $(NETCDF_LIBS) $(NETCDF_FCLIBS)

cr_lake_files_SOURCES = affinity.c \
                      memuse.c \
                      gradient_c2l.c \
                      create_xgrid.c \
                      interp.c \
                      threadloc.c \
                      nsclock.c \
                      mosaic_util.c \
                      read_mosaic.c \
                      platform.F90 \
                      mpp_parameter.F90 \
                      mpp_data.F90 \
                      mpp.F90 \
                      mpp_memutils.F90 \
                      mpp_efp.F90 \
                      mpp_pset.F90 \
                      mpp_domains.F90 \
                      mpp_io.F90 \
                      fms_io.F90 \
                      memutils.F90 \
                      constants.F90 \
                      fms.F90 \
                      horiz_interp_type.F90 \
                      horiz_interp_conserve.F90 \
                      horiz_interp_spherical.F90 \
                      horiz_interp_bicubic.F90 \
                      horiz_interp_bilinear.F90 \
                      horiz_interp.F90 \
                      cr_lake_files.f90 \
                      create_xgrid.h \
                      file_version.h \
                      fms_platform.h \
                      gradient_c2l.h \
                      interp.h \
                      mosaic_util.h \
                      mpp_group_update.h \
                      read_mosaic.h \
                      fms_io_unstructured_field_exist.inc \
                      fms_io_unstructured_file_unit.inc \
                      fms_io_unstructured_get_field_size.inc \
                      fms_io_unstructured_get_file_name.inc \
                      fms_io_unstructured_get_file_unit.inc \
                      fms_io_unstructured_read.inc \
                      fms_io_unstructured_register_restart_axis.inc \
                      fms_io_unstructured_register_restart_field.inc \
                      fms_io_unstructured_save_restart.inc \
                      fms_io_unstructured_setup_one_field.inc \
                      group_update_pack.inc \
                      group_update_unpack.inc \
                      mpp_alltoall_mpi.h \
                      mpp_alltoall_nocomm.h \
                      mpp_alltoall_sma.h \
                      mpp_chksum.h \
                      mpp_chksum_int.h \
                      mpp_chksum_scalar.h \
                      mpp_comm.inc \
                      mpp_comm_mpi.inc \
                      mpp_comm_nocomm.inc \
                      mpp_comm_sma.inc \
                      mpp_data_mpi.inc \
                      mpp_data_nocomm.inc \
                      mpp_data_sma.inc \
                      mpp_define_nest_domains.inc \
                      mpp_do_check.h \
                      mpp_do_checkV.h \
                      mpp_do_get_boundary.h \
                      mpp_do_global_field.h \
                      mpp_do_redistribute.h \
                      mpp_do_update.h \
                      mpp_do_updateV.h \
                      mpp_do_updateV_ad.h \
                      mpp_do_updateV_nonblock.h \
                      mpp_do_update_ad.h \
                      mpp_do_update_nest.h \
                      mpp_do_update_nonblock.h \
                      mpp_domains_comm.inc \
                      mpp_domains_define.inc \
                      mpp_domains_misc.inc \
                      mpp_domains_reduce.inc \
                      mpp_domains_util.inc \
                      mpp_error_a_a.h \
                      mpp_error_a_s.h \
                      mpp_error_s_a.h \
                      mpp_error_s_s.h \
                      mpp_gather.h \
                      mpp_get_boundary.h \
                      mpp_global_field.h \
                      mpp_global_field_ug.h \
                      mpp_global_reduce.h \
                      mpp_global_sum.h \
                      mpp_global_sum_ad.h \
                      mpp_global_sum_tl.h \
                      mpp_io_connect.inc \
                      mpp_io_misc.inc \
                      mpp_io_read.inc \
                      mpp_io_unstructured_read.inc \
                      mpp_io_unstructured_write.inc \
                      mpp_io_util.inc \
                      mpp_io_write.inc \
                      mpp_read_2Ddecomp.h \
                      mpp_read_compressed.h \
                      mpp_read_distributed_ascii.h \
                      mpp_read_distributed_ascii.inc \
                      mpp_reduce_mpi.h \
                      mpp_reduce_nocomm.h \
                      mpp_reduce_sma.h \
                      mpp_scatter.h \
                      mpp_sum.inc \
                      mpp_sum_mpi.h \
                      mpp_sum_nocomm.h \
                      mpp_sum_sma.h \
                      mpp_transmit.inc \
                      mpp_transmit_mpi.h \
                      mpp_transmit_nocomm.h \
                      mpp_transmit_sma.h \
                      mpp_unstruct_domain.inc \
                      mpp_unstruct_pass_data.h \
                      mpp_update_domains2D.h \
                      mpp_update_domains2D_ad.h \
                      mpp_update_domains2D_nonblock.h \
                      mpp_update_nest_domains.h \
                      mpp_util.inc \
                      mpp_util_mpi.inc \
                      mpp_util_nocomm.inc \
                      mpp_util_sma.inc \
                      mpp_write.h \
                      mpp_write_2Ddecomp.h \
                      mpp_write_compressed.h \
                      mpp_write_unlimited_axis.h \
                      read_data_2d.inc \
                      read_data_3d.inc \
                      read_data_4d.inc \
                      system_clock.h \
                      write_data.inc

mpp_data.$(OBJEXT): mpp_data.F90 mpp_parameter.$(OBJEXT)
mpp.$(OBJEXT): mpp.F90 mpp_data.$(OBJEXT) mpp_parameter.$(OBJEXT)
mpp_memutils.$(OBJEXT): mpp_memutils.F90 mpp.$(OBJEXT)
mpp_efp.$(OBJEXT): mpp_efp.F90 mpp.$(OBJEXT)
mpp_pset.$(OBJEXT): mpp_pset.F90 mpp.$(OBJEXT)
mpp_domains.$(OBJEXT): mpp_domains.F90 mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) mpp_data.$(OBJEXT) mpp_memutils.$(OBJEXT) mpp_pset.$(OBJEXT) mpp_efp.$(OBJEXT)
mpp_io.$(OBJEXT): mpp_io.F90 mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) mpp_domains.$(OBJEXT)
fms_io.$(OBJEXT): fms_io.F90 mpp_io.$(OBJEXT) mpp_domains.$(OBJEXT) mpp.$(OBJEXT) mpp_parameter.$(OBJEXT) platform.$(OBJEXT)
memutils.$(OBJEXT): memutils.F90 mpp_io.$(OBJEXT) mpp.$(OBJEXT)
fms.$(OBJEXT): fms.F90 mpp.$(OBJEXT) mpp_domains.$(OBJEXT) mpp_io.$(OBJEXT) fms_io.$(OBJEXT) memutils.$(OBJEXT) constants.$(OBJEXT)
horiz_interp_type.$(OBJEXT): horiz_interp_type.F90 mpp.$(OBJEXT)
horiz_interp_conserve.$(OBJEXT): horiz_interp_conserve.F90 mpp.$(OBJEXT) fms.$(OBJEXT) fms_io.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_spherical.$(OBJEXT): horiz_interp_spherical.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_bicubic.$(OBJEXT): horiz_interp_bicubic.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp_bilinear.$(OBJEXT): horiz_interp_bilinear.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT)
horiz_interp.$(OBJEXT): horiz_interp.F90 mpp.$(OBJEXT) fms.$(OBJEXT) constants.$(OBJEXT) horiz_interp_type.$(OBJEXT) horiz_interp_conserve.$(OBJEXT) horiz_interp_bilinear.$(OBJEXT) horiz_interp_bicubic.$(OBJEXT) horiz_interp_spherical.$(OBJEXT)
cr_lake_files.$(OBJEXT): cr_lake_files.f90 horiz_interp.$(OBJEXT)

CLEANFILES = *.mod
