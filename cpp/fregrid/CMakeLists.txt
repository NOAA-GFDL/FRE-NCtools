
include_directories(${CMAKE_SOURCE_DIR}/../include)

link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

add_executable(searchtest searchtest.C)

#target_link_libraries(fregrid libCommon)
target_link_libraries(searchtest libSearch)
