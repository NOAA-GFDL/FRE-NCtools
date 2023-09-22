
add_executable(bbox_tests bbox_tests.C)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(include)

#find_package(fmt REQUIRED)

target_link_libraries(bbox_tests libFrencu Catch2WithMain tbb)
include_directories(${CMAKE_SOURCE_DIR}/libfrencutils)



