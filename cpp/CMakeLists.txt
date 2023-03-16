cmake_minimum_required(VERSION 3.24)
project(fregrid)

set(CMAKE_CXX_STANDARD 23)

option(USE_NEXTAFTER "Use nextafter expanded bounding intervals" OFF)

add_subdirectory(fregrid)
add_subdirectory(common)
add_subdirectory(search)

#add_subdirectory(vendor/edlib EXCLUDE_FROM_ALL)

add_dependencies(fregrid libSearch)
#add_dependencies(libSearch libCommon)
