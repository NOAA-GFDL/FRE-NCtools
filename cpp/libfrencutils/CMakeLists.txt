#add_library(libFrenc INTERFACE)
#target_include_directories(libFrenc INTERFACE ${CMAKE_SOURCE_DIR}/frenc)

add_library(libFrencu STATIC create_xgrid.C gradient_c2l.C interp.C mosaic_util.C
        mpp.C mpp_domain.C mpp_io.C
        read_mosaic.C tool_util.C config.h)
#target_include_directories(libFrenc STATIC  ${CMAKE_SOURCE_DIR}/frenc)
