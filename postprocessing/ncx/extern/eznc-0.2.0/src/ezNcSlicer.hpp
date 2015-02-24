/*
A file-global slice maker such that all dim slices are set once and then shared by all variables. 4 slice modes are possible: indexed min/max/stride, value min/max/stride, indices, values.

Slices that are based on values (strided or list-based) are not valid until an active ezNc object is created and set so that dimension variables, aka coordinate variables, can be loaded and searched to define index-based slices.

Usage:
sl = ez::ezNcSlicer<size_t>();
sl.offset = 0; // or -1 for fortran.
sl.setDim...
sl.build(nc);
sl.getDim... // use details to build ezNcSlabIterator.

20111205 rsz Created.
*/
#ifndef EZNCSLICER_H
#define EZNCSLICER_H

#include <ezSlice.hpp>
#include "ezNc.hpp"

namespace ez {

template <typename T>
class ezNcSlicer {
typedef std::map<std::string, ez::ezSlice<T> > IndexSliceType;
typedef std::map<std::string, ez::ezSlice<double> > ValueSlicesType;
typedef std::map<std::string, std::vector<T> > IndexListType;
typedef std::map<std::string, std::vector<double> > ValueListType;

public:
  // Maps dim name to an index slice min/max/stride spec.
  IndexSliceType indexSlices;
  // Maps dim name to an value slice min/max/stride spec which must be processed to find index slice limits.
  ValueSlicesType valueSlices;
  // Maps dim name to an index list.
  IndexListType indexLists;
  // Maps dim name to list of values which must be processed to find indices.
  ValueListType valueLists;
  // What type of slice does dim use ultimately.
  std::map<std::string, char> slicetypes;
  // If using fortran 1-based indexing, set this to -1.
  int offset;
  // Operations are relative to this file.
  ez::ezNc* nc;
  
  // After slice optimization with "build", slices have 4 canonical categories.
  // NOSLICE means an entire dim is used as-is. No data from slicer is needed.
  // INDEXSLICE means that non-contiguous indices along a dim were chosen. Use getDimIndices to get list of indices, which are stored in "indexLists".
  // RANGESLICE means a stride=1 dim range from a min to max index. Use getDimMinMaxStride, and assume stride will be 1. Stored in indexSlices.
  // STRIDESLICE means that negative or greater than 1 stride from a min to max. Use getDimMinMaxStride.
  enum SLICETYPE { NOSLICE=0, INDEXSLICE, RANGESLICE, STRIDESLICE };
  
  // Builds index lists using values specific to nc instance and optimizes slice definitions if they're contiguous ranges, for example stride is 1 or indices aren't scattered or unordered. Helps make I/O contiguous. This may alter "indexSlices" and "indexLists" as they're used as inplace, for example a dim slice may be set as -9:-6:2, but after "build", it will become 1,4,2 so we don't repeatedly call ezSlice::get for each dim query.
  void build(ezNc*);
  // Frees all memory.
  void clear();
  // Ctor.
  ezNcSlicer() : nc(0), offset(0) {};
  // Dtor.
  ~ezNcSlicer() { clear(); }
  // Get final size of dim after slicing.
  T getDimSize(int id);
  T getDimSize(std::string& name);
  // Get list of indices for dim if it was index or value list sliced (slicetype==INDEXSLICE).
  void getDimIndices(int id, std::vector<T>& list);
  void getDimIndices(std::string&, std::vector<T>& list);
  // Get min/max/stride index slice info for dim (slicetype==RANGESLICE or STRIDESLICE).
  void getDimMinMaxStride(int id, T& mn, T& mx, T& stride);
  void getDimMinMaxStride(std::string&, T& mn, T& mx, T& stride);
  // Get list of dim sizes for a var.
  void getVarShape(int varid, std::vector<T>&);
  // Get list of dim sizes for a var.
  void getVarShape(std::string& varname, std::vector<T>&);
  // If is sliced convenience function.
  bool isDimSliced(std::string& name);
  // If is sliced convenience function.
  bool isDimSliced(int dimid);
  // If is sliced convenience function.
  bool isVarSliced(std::string& name);
  // If is sliced convenience function.
  bool isVarSliced(int varid);
  // Set the indices with strings, which can include elements from begin to end if the list has other strings such as from commmand-line parsing.
  void setDimIndices(std::string& name, std::vector<std::string>& list, int begin, int end);
  // Set the index values with strings, which can include elements from begin to end. These will then be searched in variable with the same name as dim in "build" function.
  void setDimValues(std::string& name, std::vector<std::string>& list, int begin, int end);
  // Set the min/max/stride index as strings.
  void setDimStride(std::string& name, std::string& mn, std::string& mx, std::string& stride);
  // Set the min/max/stride values as strings, which must then be searched to find their indices in the dimensions data array, if the file has one. Assumes the variable and dimension names are the same ala CF conventions.
  void setDimStrideValues(std::string& name, std::string& mn, std::string& mx, std::string& stride);
};
//#####################################################################
template <typename T>
void ezNcSlicer<T>::build(ezNc* _nc) {
  nc = _nc;
  
  // Loop over each dim in nc and check if was sliced by user, and if not set noslice.
  int ndims = nc->getNumDims();
  int i = 0;
  ezNcDim* dim;
  ezNcVar* var;
  ezNcBuffers buf;
  T min,max,stride;
  
  for(; i < ndims; ++i) {
    dim = nc->getDim(i);
    if (!dim) continue;
    
    if (indexSlices.count(dim->name) && indexSlices[dim->name].stride) {
      ezSlice<T> & slice = indexSlices[dim->name];
      slice.get(min, max, stride, dim->size);
      ///printf("%s:%d get min=%d max=%d stride=%d\n", __FILE__, __LINE__, (int)min, (int)max, (int)stride);
      // Now overwrite with final results that will persist throughout life of nc.
      // Even with non-unit stride, min/max doesn't not need to be inclusive bounds because any loop iteration will still terminate correctly.
      slice.min = min;
      slice.max = max;
      // Check if slice is simple linear sequence (i.e. abs(stride) is 1).
      if ((stride == 1) || (stride == -1)) {
        slicetypes[dim->name] = RANGESLICE;
      } else {
        slicetypes[dim->name] = STRIDESLICE;
      }
    } else if (valueSlices.count(dim->name) && valueSlices[dim->name].stride) {
      // Check if a var exists with dim's name (CF convention).
      var = nc->varsNameMap[dim->name];
      if (!var) {
        fprintf(stderr, "ERROR: %s:%d Failed to find variable with name \"%s\".\n", __FILE__, __LINE__, dim->name.c_str());
        slicetypes[dim->name] = NOSLICE;
        continue;
      }
      // Read its values, and find indices >= vmn and <= vmx.
      ezSlice<double> & vslice = valueSlices[dim->name];
      ezNcVariant variant;
      variant.type = var->type;
      variant.set<double>(vslice.min);
      nc->init(buf, var->varid);
      buf.allocate();
      nc->read(var->varid, &buf);
      // Coord var data can go "up" or "down", so be robust here.
      if (buf.increasing(var->type)) {
        min = buf.findFirstGreaterEqual(variant);
        variant.set<double>(vslice.max);
        max = buf.findLastLessEqual(variant);
      } else {
        // Reverse so indices still are increasing to follow data order in memory.
        max = buf.findLastGreaterEqual(variant);
        variant.set<double>(vslice.max);
        min = buf.findFirstLessEqual(variant);
      }
      if (vslice.strideMissing) vslice.stride = 1;
      ///printf("vslice.min=%g\n", vslice.min);
      ///printf("vslice.max=%g\n", vslice.max);
      ///printf("set min=%d max=%d stride=%d\n", (int)min, (int)max, (int)stride);

      // Create new slice, set values.
      ezSlice<T> & slice = indexSlices[dim->name];
      slice.set(min, max, (T)vslice.stride, offset);
      // Let it compute offset indices.
      slice.get(min, max, stride, dim->size);
      ///printf("get min=%d max=%d stride=%d\n", (int)min, (int)max, (int)stride);
      // Then check if a range or strided slice.
      if ((stride == 1) || (stride == -1)) {
        slicetypes[dim->name] = RANGESLICE;
      } else {
        slicetypes[dim->name] = STRIDESLICE;
      }
      buf.reset();
    } else if (indexLists.count(dim->name) && indexLists[dim->name].size()) {
      // Check if the indices are a sequence up or down and stride==1.
      // Don't sort list to preserve user's desired order, if they want to re-order dimensions, but check if contiguous to clump i/o.
      ///printf("%s:%d\n", __FILE__, __LINE__);
      if (isSeq<T>(indexLists[dim->name])) {
        // If seq, store the min/max indices and stride as a new slice.
        indexSlices[dim->name].min = indexLists[dim->name][0];
        indexSlices[dim->name].max = indexLists[dim->name][indexLists[dim->name].size()-1];
        indexSlices[dim->name].stride = 1;
        indexLists.erase(dim->name);
        slicetypes[dim->name] = RANGESLICE;
      } else if (isRevSeq<T>(indexLists[dim->name])) {
        ///printf("%s:%d\n", __FILE__, __LINE__);
        // If seq, store the min/max indices and stride as a new slice.
        indexSlices[dim->name].min = indexLists[dim->name][0];
        indexSlices[dim->name].max = indexLists[dim->name][indexLists[dim->name].size()-1];
        indexSlices[dim->name].stride = -1;
        indexLists.erase(dim->name);
        ///printf("%s:%d min=%d, max=%d, stride=%d\n", __FILE__, __LINE__, (int)indexSlices[dim->name].min, (int)indexSlices[dim->name].max, (int)indexSlices[dim->name].stride);
        slicetypes[dim->name] = RANGESLICE;
      } else if (isStrideSeq<T>(indexLists[dim->name])) {
        indexSlices[dim->name].min = indexLists[dim->name][0];
        indexSlices[dim->name].max = indexLists[dim->name][indexLists[dim->name].size()-1];
        indexSlices[dim->name].stride = (indexSlices[dim->name].max - indexSlices[dim->name].min)/(indexLists[dim->name].size()-1);
        indexLists.erase(dim->name);
        slicetypes[dim->name] = STRIDESLICE;      
      } else {
        // Not sequence. Use list as is.
        slicetypes[dim->name] = INDEXSLICE;
      }
    } else if (valueLists.count(dim->name) && valueLists[dim->name].size()) {
      // Test for equality may be bogus for real numbers if user gave a string, but do our best to find values in coordvar data.
      // Check if a var exists with dim's name (CF convention).
      var = nc->varsNameMap[dim->name];
      if (!var) {
        fprintf(stderr, "ERROR: %s:%d Failed to find variable with name \"%s\".\n", __FILE__, __LINE__, dim->name.c_str());
        slicetypes[dim->name] = NOSLICE;
        continue;
      }
      nc->init(buf, var->varid);
      buf.allocate();
      nc->read(var->varid, &buf);
      buf.find(var->type, valueLists[dim->name], indexLists[dim->name]);
      // Check if the indices are a sequence up or down and stride==1.
      if (isSeq<T>(indexLists[dim->name])) {
        // If seq, store the min/max indices and stride as a new slice.
        indexSlices[dim->name].min = indexLists[dim->name][0];
        indexSlices[dim->name].max = indexLists[dim->name][indexLists[dim->name].size()-1];
        indexSlices[dim->name].stride = 1;
        indexLists[dim->name].clear();
        slicetypes[dim->name] = RANGESLICE;
      } else if (isStrideSeq<T>(indexLists[dim->name])) {
        indexSlices[dim->name].min = indexLists[dim->name][0];
        indexSlices[dim->name].max = indexLists[dim->name][indexLists[dim->name].size()-1];
        indexSlices[dim->name].stride = (indexSlices[dim->name].max - indexSlices[dim->name].min)/(indexLists[dim->name].size()-1);
        indexLists[dim->name].clear();
        slicetypes[dim->name] = STRIDESLICE;      
      } else {
        // Not sequence. Use order of values as it.
        slicetypes[dim->name] = INDEXSLICE;
      }
      valueLists[dim->name].clear();
    } else {
      slicetypes[dim->name] = NOSLICE;
    }
  }
};
//#####################################################################
template <typename T>
void ezNcSlicer<T>::clear() {
  indexSlices.clear();
  valueSlices.clear();
  indexLists.clear();
  valueLists.clear();
};
//#####################################################################
template <typename T>
T ezNcSlicer<T>::getDimSize(int id) {
  if(!nc) return 0;
  ezNcDim* dim = nc->getDim(id);
  if (!dim) return 0;
  return getDimSize(dim->name);
};
//#####################################################################
template <typename T>
T ezNcSlicer<T>::getDimSize(std::string& dimname) {
  ezNcDim *dim, *dim2;
  ezSlice<T> * slice;
  
  switch(slicetypes[dimname]) {
    case NOSLICE:
      if(!nc) return 0;
      dim = nc->dimsNameMap[dimname];
      if (!dim) {
        fprintf(stderr, "ERROR: %s:%d Dim name \"%s\" not found.\n", __FILE__, __LINE__, dimname.c_str());
        return 0;
      }
      return dim->size;
      break;
    case INDEXSLICE:
      return indexLists[dimname].size();
      break;
    case RANGESLICE:
      slice = &indexSlices[dimname];
      if (slice->max > slice->min)
        return slice->max - slice->min + 1;
      else
        return (slice->min - slice->max) + 1;
      break;
    case STRIDESLICE:
      if(!nc) return 0;
      dim2 = nc->dimsNameMap[dimname];
      if (!dim2) {
        fprintf(stderr, "ERROR: %s:%d Dim name \"%s\" not found.\n", __FILE__, __LINE__, dimname.c_str());
        return 0;
      }
      return indexSlices[dimname].size(dim2->size);
      break;
    default: return 0;
  }
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getDimIndices(int id, std::vector<T>& list) {
  if(!nc) return;
  ezNcDim* dim = nc->getDim(id);
  if (!dim) return;
  getDimIndices(dim->name, list);
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getDimIndices(std::string& dimname, std::vector<T>& list) {
  list = indexLists[dimname];
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getDimMinMaxStride(int id, T& mn, T& mx, T& stride) {
  if(!nc) return;
  ezNcDim* dim = nc->getDim(id);
  if (!dim) return;
  getDimMinMaxStride(dim->name, mn, mx, stride);
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getDimMinMaxStride(std::string& dimname, T& mn, T& mx, T& stride) {
  ez::ezSlice<T> & slice = indexSlices[dimname];
  mn = slice.min;
  mx = slice.max;
  stride = slice.stride;
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getVarShape(int varid, std::vector<T>& shape) {
  if(!nc) return;
  ezNcVar* var = nc->getVar(varid);
  if(!var) return;
  shape.clear();
  int n = var->getNumDims();
  shape.resize( n );
  
  for(int i=0; i < n; ++i) {
    shape[i] = getDimSize(var->dimids[i]);
  }
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::getVarShape(std::string& varname, std::vector<T>& shape) {
  if(!nc) return;
  ezNcVar* var = nc->varsNameMap[varname];
  if(!var) return;
  shape.clear();
  int n = var->getNumDims();
  shape.resize( n );
  
  for(int i=0; i < n; ++i) {
    shape[i] = getDimSize(var->dimids[i]);
  }
}
//#####################################################################
template <typename T>
bool ezNcSlicer<T>::isDimSliced(int dimid) {
  if(!nc) return false;
  ezNcDim* dim = nc->getDim(dimid);
  if (!dim) return false;
  return isDimSliced(dim->name);
}
//#####################################################################
template <typename T>
bool ezNcSlicer<T>::isDimSliced(std::string& dimname) {
  switch(slicetypes[dimname]) {
    case NOSLICE: return false;
    default: return true;
  }
}
//#####################################################################
template <typename T>
bool ezNcSlicer<T>::isVarSliced(int varid) {
  if(!nc) return false;
  ezNcVar* var = nc->getVar(varid);
  if(!var) return false;
  return isVarSliced(var->name);
}
//#####################################################################
template <typename T>
bool ezNcSlicer<T>::isVarSliced(std::string& varname) {
  if(!nc) return false;
  ezNcVar* var = nc->varsNameMap[varname];
  if(!var) return false;
  int n = var->getNumDims();
  
  for(int i=0; i < n; ++i) {
    if (isDimSliced(var->dimids[i]))
      return true;
  }  
  
  return false;  
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::setDimIndices(std::string& dimname, std::vector<std::string>& list, int begin, int end) {
  if ( (begin >= list.size()) || (begin < 0) ) return;
  int n = end - begin + 1;
  if ((begin+n) > list.size()) return;

  std::vector<T> & indices = indexLists[dimname];
  std::stringstream ss;
  int i = 0;
  indices.resize(n);

  for(; i < n; ++i) {
    ss.clear();
    ss.str(list[begin+i]);
    ss >> indices[i];
    indices[i] += offset;
    ///printf("%s:%d indices[%d] = %d\n", __FILE__, __LINE__, i, (int)indices[i]);
  }
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::setDimValues(std::string& dimname, std::vector<std::string>& list, int begin, int end) {
  if ( (begin >= list.size()) || (begin < 0) ) return;
  int n = end - begin + 1;
  if ((begin+n) > list.size()) return;

  std::vector<double> & vals = valueLists[dimname];
  std::stringstream ss;
  int i = 0;
  vals.resize(n);
  
  for(; i < n; ++i) {
    ss.clear();
    ss.str(list[begin+i]);
    ss >> vals[i];
  }
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::setDimStride(std::string& dimname, std::string& min, std::string& max, std::string& stride) {
  if (dimname.empty()) return;
  ez::ezSlice<T> & slice = indexSlices[dimname];
  ///printf("%s:%d min=%s, max=%s, stride=%s, offset=%d\n", __FILE__, __LINE__, min.c_str(), max.c_str(), stride.c_str(), offset);
  slice.set(min.c_str(), max.c_str(), stride.c_str(), offset);
}
//#####################################################################
template <typename T>
void ezNcSlicer<T>::setDimStrideValues(std::string& dimname, std::string& min, std::string& max, std::string& stride) {
  if (dimname.empty()) return;
  ez::ezSlice<double> & slice = valueSlices[dimname];
  // Since we're dealing with values that must be searched, rather than simpler indices that can be inferred from dim size, do auto string substitution for empty args.
  slice.set(
    min.empty() ? "-1.79769e+308" : min.c_str(), 
    max.empty() ? "1.79769e+308" : max.c_str(), 
    stride.c_str(),
    offset);
}
//#####################################################################
};
#endif // EZNCSLICER_H
