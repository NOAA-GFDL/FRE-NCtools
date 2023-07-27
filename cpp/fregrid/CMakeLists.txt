
include_directories(${CMAKE_SOURCE_DIR}/../include)
include_directories(${CMAKE_SOURCE_DIR}/libfrencutils)
include_directories(${CMAKE_SOURCE_DIR}/libSearch)
include_directories(/home/mzuniga/llvm/build/install/include)
include_directories(/home/mzuniga/llvm/build/install/include/sycl)
include_directories(/home/mzuniga/llvm/build/install/include/sycl/CL)



link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(fregrid.x fregrid.C conserve_interp.C fregrid_util.C bilinear_interp.C)
#add_sycl_to_target(TARGET fregrid.x SOURCES  fregrid.C)

add_executable(searchtest.x searchtest.C)
#add_sycl_to_target(TARGET searchtest.x SOURCES  searchtest.C) #in the future, may require other .C files

target_link_libraries(searchtest.x libSearch tbb)
target_link_libraries(fregrid.x  libFrencu libSearch netcdf tbb)

