cmake_minimum_required(VERSION 3.21)
project(fregrid)

project(fregrid LANGUAGES CXX)

#TODO: C++20 or C++23?
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_VERBOSE_MAKEFILE ON)

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

#find_package(Boost REQUIRED)

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
