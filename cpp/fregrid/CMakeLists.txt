
include_directories(${CMAKE_SOURCE_DIR}/../include)
include_directories(${CMAKE_SOURCE_DIR}/libfrencutils)
include_directories(${CMAKE_SOURCE_DIR}/libSearch)


link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(searchtest.x searchtest.C)
add_executable(fregrid.x fregrid.C conserve_interp.C fregrid_util.C bilinear_interp.C)

#target_link_libraries(fregrid libCommon)
target_link_libraries(searchtest.x libSearch)
target_link_libraries(fregrid.x  libFrencu libSearch netcdf)

