/*
Iterates sequentially through possible slab "start" indices for a single var. The index sequences per dimension don't have to be linear, thus allowing noncontiguous, nonsequential slicing.

Static and record dimension incrementing can be exclusive to allow for flexible record-interleaved operations.

Template type should be signed, such as int or long long, to support negative strides and decreasing sequences.

20111029 rsz Created.
*/
#ifndef EZNCSLABITERATOR_H
#define EZNCSLABITERATOR_H

#include <ezNcSlicer.hpp>
#include <ezNcSlab.hpp>
#include <ezOdometer.hpp>

namespace ez {

template<typename T>
class ezNcSlabIterator {
public:
  // Ctor.
  inline ezNcSlabIterator();
  // Dtor.
  inline ~ezNcSlabIterator();
  // Clears all memory.
  inline void clear();
  // Fill an allocated slab instance with current state.
  inline void getSlab(ezNcSlab<size_t>*);
  // Try to increment using all digits.
  inline void inc() { ++odo; };
  // Increment only the rec dim counter, if any.
  inline void incRec();
  // Increment only the static dim counters, if any.
  inline void incStatic();
  // Return true if iterator reached end of its incrementing range.
  inline bool isDone() { return odo.isDone(); };
  // Check if record dim counter reached limit (if it uses rec dim).
  inline bool isDoneRec();
  // Check if all static dim counters reached limit (if it has any).
  inline bool isDoneStatic();
  // Normalizes the start counters so they all go from 0 to dim size minus 1.
  inline void normalize() { odo.normalize(); };
  // Set back counter so can iterate again.
  inline void reset();
  // Set back rec counter only to iterate again over rec dim.
  inline void resetRec();
  // Set back only static dim counters for re-iterating again.
  inline void resetStatic();
  // Sets up the odometer with proper values based on how var dims are sliced according to ezNcSlicer instance. sliceRec=true forces rec dim to have count=1 and be incrementable (desired for doing sequential I/O for many rec vars in a file). sliceRec=false does not override rec slice definition in the ezNcSlicer.
  inline void setVar(ezNcVar*, ezNcSlicer<T>*, bool sliceRec=false);
  
  // Shape of var won't change when set, and must be used to create slabs.
  std::vector<T> count;
  // Does the counter logic. odo.digits will be the slab 'start' index array. Must be signed integer to support negative strides.
  ezCustomOdometer<T> odo;
  // What var this is generating slabs for.
  ezNcVar* var;
};
//#####################################################################
template<typename T>
ezNcSlabIterator<T>::ezNcSlabIterator() {
  var = 0;
};
//#####################################################################
template<typename T>
ezNcSlabIterator<T>::~ezNcSlabIterator() {
  clear();
};
//#####################################################################
template<typename T>
void ezNcSlabIterator<T>::clear() {
  count.clear();
  odo.clear();
  var = 0;
};
//#####################################################################
template<typename T>
bool ezNcSlabIterator<T>::isDoneRec() {
  if ((var==0) || (var->hasrec==false)) return true;

  int oldend = odo.end;
  // rec is always first dim.
  odo.end = 0;
  bool res = odo.isDone();
  odo.end = oldend;
  
  return res;
}
//#####################################################################
template<typename T>
bool ezNcSlabIterator<T>::isDoneStatic() {
  int oldbegin = odo.begin;
  
  if (var && var->hasrec)
    // Skip first rec dim digit.
    odo.begin = 1;
  
  bool res = odo.isDone();
  odo.begin = oldbegin;      
  
  return res;
}
//#####################################################################
template <typename T>
void ezNcSlabIterator<T>::getSlab(ezNcSlab<size_t>* slab) {
  if (slab == 0) return;
  if (var == 0) return;
  
  slab->setVar(var);
  int i=0;
  int n = odo.size();
  
  for(; i < n; ++i) {
    slab->start[i] = (size_t)odo.digits[i];
    slab->count[i] = (size_t)count[i];
  }
};
//#####################################################################
template<typename T>
void ezNcSlabIterator<T>::incRec() { 
  if ((var==0) || (var->hasrec==false)) return;
  
  // If var has rec dim (index=0), then just inc first digit.
  int prevbegin = odo.begin;
  int prevend = odo.end;
  odo.begin = 0;
  odo.end = 0;
  ++odo;
  odo.begin = prevbegin;
  odo.end = prevend;
}
//#####################################################################
template<typename T>
void ezNcSlabIterator<T>::incStatic() { 
  int i = 0;
  if (var && var->hasrec)
    // Skip first digit which is rec dim.
    i = 1;
    
  // Increment all other dims.
  int prevbegin = odo.begin;
  odo.begin = i;
  ++odo;
  odo.begin = prevbegin;
}
//#####################################################################
template<typename T> 
void ezNcSlabIterator<T>::reset() {
  odo.reset();
};
//#####################################################################
template<typename T> 
void ezNcSlabIterator<T>::resetRec() {
  if ((var==0) || (var->hasrec==false)) return;
  
  int prevbegin = odo.begin;
  int prevend = odo.end;
  odo.begin = 0;
  odo.end = 0;
  odo.reset();
  odo.begin = prevbegin;
  odo.end = prevend;
};
//#####################################################################
template<typename T> 
void ezNcSlabIterator<T>::resetStatic() {
  int i = 0;
  if (var && var->hasrec)
    // Skip first digit which is rec dim.
    i = 1;
    
  // Increment all other dims.
  int prevbegin = odo.begin;
  odo.begin = i;
  odo.reset();
  odo.begin = prevbegin;
};
//#####################################################################
template<typename T>
void ezNcSlabIterator<T>::setVar(ezNcVar* _var, ezNcSlicer<T>* slicer, bool sliceRec) {
  if (_var == 0) return;
  if (slicer == 0) return;

  var = _var; 
  reset(); 

  int ndims = var->getNumDims();
  odo.init(ndims);
  count.resize(ndims);
  int i = 0;

  if (slicer->isVarSliced(var->name)) {
    if (var->nc == 0) return;
    ezNcDim* dim;
    std::vector<T> list;
    T min, max, stride;
    int dimid;
    
    // Set the indices per dim that will be sliced out of full var.
    for(; i < ndims; ++i) {
      dim = var->nc->getDim(var->dimids[i]);
      if (!dim) continue;
      switch(slicer->slicetypes[dim->name]) {
        case ez::ezNcSlicer<T>::INDEXSLICE:
          slicer->getDimIndices(dim->name, list);
          odo.setCustomValues(i, list); // start = some index.
          count[i] = 1; // count = single slice along dim.
          break;
        case ez::ezNcSlicer<T>::RANGESLICE:
          slicer->getDimMinMaxStride(dim->name, min, max, stride);
          ///printf("%s:%d min=%d, max=%d, stride=%d\n", __FILE__, __LINE__, (int)min, (int)max, (int)stride);
          if (stride == 1) {
            // Ignore stride because guaranteed to be 1.
            // Set up single contiguous slice that doesn't need to iterate.
            odo.setCustomValue(i, min);
            count[i] = max - min + 1; // Span from min to max inclusive.
          } else {
            // Stride must be negative, and Netcdf library can't do contiguous reversed read/write, so actually create singleton slices.
            odo.setrange(i, min, max, stride);
            count[i] = 1;
          }
          break;
        case ez::ezNcSlicer<T>::STRIDESLICE:
          slicer->getDimMinMaxStride(dim->name, min, max, stride);
          odo.setrange(i, min, max, stride);
          // Non contiguous, so can only slice one index at a time.
          count[i] = 1;
          break;
        default: //case ez::ezNcSlicer<T>::NOSLICE:
          odo.setCustomValue(i, 0); // start = 0.
          count[i] = slicer->getDimSize(dim->name); // count = n.
          break;
      }
    }
  } else {
    // Var not sliced. Set starts as all zeros and count as full shape.
    for(; i < ndims; ++i) {
      odo.setCustomValue(i, 0);
      count[i] = slicer->getDimSize(var->dimids[i]);
    }
  }

  // Slice rec only if not sliced already or singleton, i.e. count==1 means already sliced.
  if (sliceRec && var->hasrec && odo.size() && (count[0]>1)) {
    // Rec dim must be outermost index 0.
    // If count[0]>1, then slice it by making linear sequence range.
    // min = first rec digit.
    // max = first rec digit plus count minus one.
    // stride must be one, since slicing NOSLICE or RANGESLICE.
    odo.setrange(0, odo.digits[0], odo.digits[0]+count[0]-1, 1);
    count[0] = 1;
  }
};
//#####################################################################
};
#endif // EZNCSLABITERATOR_H