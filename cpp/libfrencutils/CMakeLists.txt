
add_library(libFrencu  create_xgrid.C gradient_c2l.C interp.C mosaic_util.C
        mpp.C mpp_domain.C mpp_io.C read_mosaic.C tool_util.C config.h)

add_sycl_to_target(TARGET libFrencu SOURCES create_xgrid.C ) #may need other .C files later
add_dependencies(libFrencu libSearch)
target_link_libraries(libFrencu libSearch)
