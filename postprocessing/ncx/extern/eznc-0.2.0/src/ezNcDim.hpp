/*

*/
/*
CHANGELOG
v0.0.0 20110523 rsz Created.
*/
#ifndef EZNCDIM_H
#define EZNCDIM_H

#include <vector>
#include <string>

namespace ez {
// Dummy forward definition so ezNcDim can declare a pointer.
class ezNc;

//####################################################################
class ezNcDim {
public:
  // The file wrapper this dim belongs to. 0 if there is no file.
  ezNc *nc;
  // The id of this dim in its parent file. -1 if no file uses it.
  int dimid;
  // If this dim is a record dimension.
  bool isrec;
  // Size of dimension; number of elements along its axis, regardless of stride/sample.
  size_t size;
  // The dimension's name.
  std::string name;
  
  ezNcDim() : nc(0), dimid(-1), isrec(false), size(0) {};
};
//###########################################################################
#ifdef DEBUG
void print(ezNcDim & d) {
  std::stringstream ss;
  ss << "Dim name=" << d.name << ", dimid=" << d.dimid << ", isrec=" << d.isrec << ", size=" << d.size << "\n";
  std::cout << ss;
};
#endif
}
#endif // EZNCDIM_H