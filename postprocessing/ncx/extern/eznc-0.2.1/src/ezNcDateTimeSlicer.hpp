/*
Adds datetime slicing ability to standard slicer.

Usage:
  ezNcDateTimeSlicer<int> slice;
  slice.begin = "1900-01-01";
  slice.end = "1999-12-31";
  ezNcDateTime dt;
  // ...
  slice.build(dt);
  // Now use in slab iterator or use to get var/dim sizes, shapes, indices.

20120103 rsz Created.
*/
#ifndef EZNCDATETIMESLICER_H
#define EZNCDATETIMESLICER_H

#include "ezNc.hpp"
#include "ezNcSlicer.hpp"
#include "ezNcDateTime.hpp"

namespace ez {

template <typename T>
class ezNcDateTimeSlicer : public ezNcSlicer<T> {
public:
  // Ctor.
  ezNcDateTimeSlicer() : datetime(0) {};
  
  // Computes the date time index range, as well as any other dim slices.
  inline void build(ezNcDateTime*);
  
  // These datetime strings are used in "build".
  std::string begin, end;
  // This computes the indices from the begin,end strings.
  ezNcDateTime* datetime;
};

//##########################################################
template <typename T>
void ezNcDateTimeSlicer<T>::build(ezNcDateTime* dt) {
  datetime = dt;
  if (!datetime) return;  
  if (!datetime->var) return;
  if (!datetime->var->nc) return;
  
  // Get dim name for record axis.
  ezNcDim* dim = 0;
  if (datetime->var->dimids.size()) {
     dim = datetime->var->nc->getDim(datetime->var->dimids[0]);
  }
  
  if (!dim || dim->name.empty()) {
    printf("%s:%d WARNING: Failed to get name of record dimension.\n", __FILE__, __LINE__);
    return;
  }
  
  // Get the index range.
  int ibegin, iend;
  datetime->getRange(begin, end, ibegin, iend);
  //printf("%s:%d ibegin=%d, iend=%d\n", __FILE__, __LINE__, ibegin, iend);

  if ((ibegin == -1) && (iend == -1)) {
    // Invalid range, so don't slice the dim.
    ezNcSlicer<T>::indexSlices.erase(dim->name);
    //printf("%s:%d\n", __FILE__, __LINE__);
  } else {
    // Set the index slice for the rec dim.
    ez::ezSlice<T> & slice = ezNcSlicer<T>::indexSlices[dim->name];
    slice.set(ibegin, iend, 1, 0);
  }
  
  // Build the remaining slices.
  ezNcSlicer<T>::build(datetime->var->nc);
}
//##########################################################
};
#endif // EZNCDATETIMESLICER_H