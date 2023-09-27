cmake_minimum_required(VERSION 3.21)
project(fregrid)

project(fregrid LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE ON)

##TODO: Upgrade to C++23 as soon as enough features are
# available in enough compilers that are being used.
set(CMAKE_CXX_STANDARD 20)

add_subdirectory(dependencies)

#Testing:
#NOTE: The initial framework are Catch2 + ApprovalTest. This should provide easy
# unit testing, acceptance testing, and micro benchmarking. However, it may be
# convenient to directly use and add cmakes Catch via "include(Catch)". E.g.
# can nccmp be directly called?"
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif(BUILD_TESTS)



# The " {CMAKE_SOURCE_DIR}/include" directories below is only temporary as
# several C++23 lib features are being brought in for use with C++20 compilers.
# Some of them are purely in the "/include" directory, but the project source
# directory is also set as an  include directory because  of the way the format
# lib was created.  The other possibility may be to add them (temporarily) in
# the dependencies CMaleList.txt.
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
#Perhaps will need to make it compiler dependent:
##if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
##endif()

##with-mpi option not yet available in the C++ codebase
##option(with-mpi "Use legacy MPI and legacy create_xgrid in fregrid build" OFF)
option(with-gpu "Compile for GPU. Only use with nvc++" OFF )
## These experimental options are only being used by developers:
##TODO: remove SYCL use from the C++ codebase
option(with-sycl "Compile with  SYCL *** WARNING: EXPERIMENTAL *** " OFF)
#find_package(IntelSYCL REQUIRED)

if(with-gpu)
    add_definitions(-DWITH_GPU)
endif()

if(with-sycl)
    add_definitions(-DWITH_SYCL)
endif()

##Location of some libs is external to this project.
#include_directories(/home/mzuniga/NetCdCPP/install/include)
#link_directories(/home/mzuniga/NetCdCPP/install/lib)
link_directories(/usr/lib/x86_64-linux-gnu)

#find_package(NetCDF REQUIRED COMPONENTS CPP)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../bin")

add_subdirectory(common)
add_subdirectory(search)
add_subdirectory(libfrencutils)
add_subdirectory(fregrid)

if(with-sycl)
    add_subdirectory(sycl_tests)
endif(with-sycl)

add_dependencies(searchtest.x libSearch)
add_dependencies(fregrid.x libSearch libFrencu)
