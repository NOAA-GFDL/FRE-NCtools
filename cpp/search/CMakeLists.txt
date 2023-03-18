add_library(libSearch INTERFACE Polygon.h Point3D.h BruteBoxQuery.h)
target_include_directories(libSearch INTERFACE ${CMAKE_SOURCE_DIR}/search)