include_directories(${CMAKE_SOURCE_DIR}/../include)
include_directories(${CMAKE_SOURCE_DIR}/libfrencutils)
include_directories(${CMAKE_SOURCE_DIR}/libSearch)

include_directories(/home/mzuniga/llvm/build/install/include)
include_directories(/home/mzuniga/llvm/build/install/include/sycl)
include_directories(/home/mzuniga/llvm/build/install/include/sycl/CL)


#currently included in upper directory file:
#find_package(IntelSYCL REQUIRED)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(sycltest.x sycltest.C)
add_executable(synchro.x synchro.C)
#add_sycl_to_target(TARGET sycltest.x SOURCES  sycltest.C)
#add_sycl_to_target(TARGET synchro.x SOURCES  synchro.C)


#target_link_libraries(fregrid libCommon)
target_link_libraries(sycltest.x libSearch)


