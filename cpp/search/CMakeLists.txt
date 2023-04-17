add_library(libSearch INTERFACE Polygon.h Point3D.h BruteBoxQuery.h DITreeNode.h DistanceInterval.h DITree.h
        Partition.h Comparators.h BoxedObj.h BBox3D.h)
target_include_directories(libSearch INTERFACE ${CMAKE_SOURCE_DIR}/search)

