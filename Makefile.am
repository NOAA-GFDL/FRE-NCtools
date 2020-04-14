ACLOCAL_AMFLAGS = -I m4 
LDFLAGS="-L /usr/lib -lm"
SUBDIRS = postprocessing/combine_blobs \
          postprocessing/combine_restarts \
          postprocessing/iceberg_comb \
          postprocessing/land_utils \
          postprocessing/list_ncvars \
          postprocessing/mppnccombine \
          postprocessing/plevel \
          postprocessing/split_ncvars \
          postprocessing/timavg \
          tools/libfrencutils \
          tools/check_mask \
          tools/cubic_utils \
          tools/fregrid \
          tools/make_coupler_mosaic \
          tools/make_hgrid \
          tools/make_land_domain \
          tools/make_quick_mosaic \
          tools/make_regional_mosaic \
          tools/make_solo_mosaic \
          tools/make_topog \
          tools/make_vgrid \
          tools/mppncscatter \
          tools/ncexists \
          tools/remap_land \
          tools/river_regrid \
          tools/runoff_regrid \
          tools/transfer_to_mosaic_grid \
          tools/simple_hydrog/postp \
          tools/simple_hydrog/lakes \
          tools/simple_hydrog/rmvpr 

if WITH_CHECK_PROGS
SUBDIRS += t
endif
