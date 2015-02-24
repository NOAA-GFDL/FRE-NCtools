/*
20111029 rsz Created.
*/
#ifndef EZNCSLAB_H
#define EZNCSLAB_H

#include <ezNcVar.hpp>
#include <ezNcUtil.hpp>

namespace ez {

template <typename T>
class ezNcSlab {
public:
  // Ctor.
  inline ezNcSlab();
  // Ctor initializes arrays.
  inline ezNcSlab(ezNcVar*);
  // Dtor.
  inline ~ezNcSlab();
  // Resets to initial state.
  inline void clear();
  // Returns true if count arrays are same length with identical values.
  inline bool countEqual(ezNcSlab*);
  // Compute byte size of slab.
  inline size_t getNumBytes();
  // Initialize arrays to shape of var.
  inline void setVar(ezNcVar*);
  // Create string for printing.
  inline void toString(std::string& string);
  
  // What variable to slab.
  ezNcVar* var;
  // Start offsets into var data.
  std::vector<T> start;
  // Number of elements per dimension.
  std::vector<T> count;
};
//#####################################################################
template <typename T>
ezNcSlab<T>::ezNcSlab() {
  var = 0;
};
//#####################################################################
template <typename T>
ezNcSlab<T>::~ezNcSlab() {
  clear();
};
//#####################################################################
template <typename T>
bool ezNcSlab<T>::countEqual(ezNcSlab* s) {
  if (s == 0) return false;
  if (s->count.size() != count.size()) return false;
  int i=0;
  for(; i < count.size(); ++i)
    if (count[i] != s->count[i])
      return false;
      
  return true;
};
//#####################################################################
template <typename T>
void ezNcSlab<T>::clear() {
  var = 0;
  start.clear();
  count.clear();
};
//#####################################################################
template <typename T>
size_t ezNcSlab<T>::getNumBytes() {
  if (!var || count.empty()) return 0;
  
  int i = 0;
  int n = count.size();
  size_t total = 1;
  
  for(; i < n; ++i) {
    total *= count[i];
  }
  
  return total * sizeOf(var->type);
}
//#####################################################################
template <typename T>
void ezNcSlab<T>::setVar(ezNcVar* v) {
  clear();
  var = v;
  if (!v) return;
  
  int ndims = var->getNumDims();
  start.resize(ndims);
  count.resize(ndims);

  for(int i=0; i < ndims; ++i) {
    start[i] = 0;
    count[i] = 0;
  }
};
//#####################################################################
template <typename T>
void ezNcSlab<T>::toString(std::string& string) {
  int i = 0;
  int n = start.size();
  char tmp[16];
  string.append("start =");
  for(; i < n; ++i) {
    sprintf(tmp, " %d", (int)start[i]);
    string.append(tmp);
  }
  string.append(" count =");
  for(i=0; i < n; ++i) {
    sprintf(tmp, " %d", (int)count[i]);
    string.append(tmp);
  }
}
};
#endif // EZNCSLAB_H
