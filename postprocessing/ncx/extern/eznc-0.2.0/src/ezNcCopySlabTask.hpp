/*
20111029 rsz Created.
*/
#ifndef EZNCCOPYSLABTASK_H
#define EZNCCOPYSLABTASK_H

#include <ezNc.hpp>
#include <ezNcSlab.hpp>

namespace ez {
template <typename T>
class ezNcCopySlabTask {
public:
  // Performs the task. Returns 0 if successful.
  int run(ezNcBuffers& buffer);
  
  // Source slab.
  ezNcSlab<T> slabSrc;
  // Destination slab.
  ezNcSlab<T> slabDst;
};
//#####################################################################
template <typename T>
int ezNcCopySlabTask<T>::run(ezNcBuffers& buffer) {
  if ( (slabSrc.var == 0) || (slabDst.var == 0) ) return 1;
  
  int status;
  status = slabSrc.var->nc->read(slabSrc.var->varid,
                                &buffer,
                                &(slabSrc.start[0]), 
                                &(slabSrc.count[0]));
  if (HandleNcStatus(status)) {
    fprintf(stderr, "ERROR: Task read failed for var \"%s\".\n", slabSrc.var->name.c_str());
    std::vector<size_t> dimsizes;
    slabSrc.var->nc->getDimSizes(slabSrc.var->varid, dimsizes);
    int i;
    int n = slabSrc.start.size();
    fprintf(stderr, "shape =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)dimsizes[i]);
    fprintf(stderr, "\nstart =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)slabSrc.start[i]);
    fprintf(stderr, "\ncount =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)slabSrc.count[i]);
    fprintf(stderr, "\n");
    return status;
  }

  status = slabDst.var->nc->write(slabDst.var->varid,
                                &buffer, 
                                &(slabDst.start[0]), 
                                &(slabDst.count[0]));
  if (HandleNcStatus(status)) {
    fprintf(stderr, "ERROR: Task write failed for var \"%s\".\n", slabDst.var->name.c_str());
    std::vector<size_t> dimsizes;
    slabDst.var->nc->getDimSizes(slabDst.var->varid, dimsizes);
    int i;
    int n = slabDst.start.size();
    fprintf(stderr, "shape =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)dimsizes[i]);
    fprintf(stderr, "\nstart =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)slabDst.start[i]);
    fprintf(stderr, "\ncount =");
    for(i=0; i < n; ++i) fprintf(stderr, " %d", (int)slabDst.count[i]);
    fprintf(stderr, "\n");
    return status;
  }
  return status;
};
//#####################################################################
};
#endif // EZNCCOPYSLABTASK_H