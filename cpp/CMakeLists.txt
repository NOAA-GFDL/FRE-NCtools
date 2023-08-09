cmake_minimum_required(VERSION 3.21)
project(fregrid)

project(fregrid LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

##TODO: Upgrade to C++23 as soon as enough features are
# available in enough compilers.
set(CMAKE_CXX_STANDARD 20)

# These include directories below are only temporary as
# several C++23 lib features are being brought in for use
# with C++20 compilers. Some of them are purely in the "/include"
# directory, but the project source directory is also set as an
# include directory because  of the way the format lib was created.
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
#Perhaps will need to make it compiler dependent:
##if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
##endif()

#find_package(IntelSYCL REQUIRED)
option(USE_SYCL "Use SYCL compiling search and other components" OFF)
option(USE_NEXTAFTER "Use nextafter expanded bounding intervals" OFF)

##Location of some libs is external to this project.
#include_directories(/home/mzuniga/NetCdCPP/install/include)
#link_directories(/home/mzuniga/NetCdCPP/install/lib)
link_directories(/usr/lib/x86_64-linux-gnu)

##$nc-config --cflags
##-I/usr/include -I/usr/include/hdf5/serial
##$nc-config --libs
##-L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lnetcdf

#find_package(NetCDF REQUIRED COMPONENTS CPP)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")


add_subdirectory(common)
add_subdirectory(search)
add_subdirectory(libfrencutils)
add_subdirectory(fregrid)
if(USE_SYCL)
    add_subdirectory(sycl_tests)
endif(USE_SYCL)

#add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)

add_dependencies(searchtest.x libSearch)
add_dependencies(fregrid.x libSearch libFrencu)
