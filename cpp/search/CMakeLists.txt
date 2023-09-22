#NOTE: the target sources was added as CLion sometimes would not display header files
#  in the Project widget.
set(MY_SEARCH_TSS BBox3D.h BoxedObj.h  DITree.h  Point3D.h Comparators.h DITreeNode.h  Polygon.h
        BruteBoxQuery.h  DistanceInterval.h  Partition.h   TreePerfStats.h)
add_library(libSearch INTERFACE )
target_sources(libSearch INTERFACE ${MY_SEARCH_TSS})
target_include_directories(libSearch INTERFACE .)



