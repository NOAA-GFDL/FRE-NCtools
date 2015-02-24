/*
20111029 rsz Created. memtest passed.
*/
#ifndef EZPOINTERVECTOR_H
#define EZPOINTERVECTOR_H

#include <vector>

namespace ez {

/*
Smart pointer array.
Becomes owner of allocated pointers put into its vector and auto deallocates
them upon destruction.
*/
template <typename T>
class ezPointerVector {
public:
  // Dtor.
  ~ezPointerVector() { reset(); };
  // Frees pointers.
  void reset();
  
  // The real vector.
  std::vector<T*> vector;
};
//#####################################################################
template <typename T>
void ezPointerVector<T>::reset() {
  int i=0, n = vector.size();
  for(; i < n; ++i)
    delete vector[i];
    
  vector.clear();
};
//#####################################################################
};
#endif // EZPOINTERVECTOR_H