

##if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
# including the format file if compiler is not gcc.
# nvc++ does not yet have this feature.
# Intel/CodePlay DPC++ has ?
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
##endif()

include_directories(${CMAKE_SOURCE_DIR}/libfrencutils)
include_directories(${CMAKE_SOURCE_DIR}/libSearch)
include_directories(/home/mzuniga/llvm/build/install/include)
include_directories(/home/mzuniga/llvm/build/install/include/sycl)
include_directories(/home/mzuniga/llvm/build/install/include/sycl/CL)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(fregrid.x fregrid.C conserve_interp.C fregrid_util.C bilinear_interp.C)
#        ../test/main.C
#        ../test/integration_tests.C)
#add_sycl_to_target(TARGET fregrid.x SOURCES  fregrid.C)

add_executable(searchtest.x searchtest.C)
#add_sycl_to_target(TARGET searchtest.x SOURCES  searchtest.C) #in the future, may require other .C files

#add tbb ?
target_link_libraries(searchtest.x libSearch tbb)
target_link_libraries(fregrid.x  libFrencu libSearch netcdf tbb)

