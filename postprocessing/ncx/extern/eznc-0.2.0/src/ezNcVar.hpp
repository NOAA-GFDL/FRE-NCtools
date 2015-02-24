/*

*/
/*
CHANGELOG
v0.0.0 20110523 rsz Created.
*/
#ifndef EZNCVAR_H
#define EZNCVAR_H

#include <vector>
#include <string>
#include "netcdf.h"
#include <ezStringUtil.hpp>

namespace ez {
class ezNc;

//###########################################################################
class ezNcVar {
public:
  // The file wrapper of the variable owner. 0 if no such file exists.
  ezNc *nc;
  // The name of the variable.
  std::string name;
  // The id of the variable in the file. -1 if no file uses this variable.
  int varid;
  // Number of elements per record (or all the elements of non-record using variable).
  unsigned long long recsize;
  // Number of bytes this variable uses per record (or entirely for non-record using variables).
  unsigned long long recbytes;
  // Number of elements for the entire variable (including record dimension if any).
  unsigned long long size;
  // Number of bytes for entire variable (including record dimension if any).
  unsigned long long bytes;
  // Datatype code.
  nc_type type;
  // If uses a record dimension.
  bool hasrec;
  // Variable's attnames in same order as stored in file.
  std::vector<std::string> attnames;
  // Dimid list for dims this var uses.
  std::vector<int> dimids;
  
  // HDF feature.
  int chunkstorage;
  // HDF feature. Array must have one chunksize for each dimension of the variable.
  std::vector<size_t> chunksizes;
  // HDF feature.
  size_t chunksize;
  // HDF feature.
  size_t chunknelems;
  // HDF feature.
  float chunkpreemption;
  // Var specific fill mode.
  int no_fill;
  // For HDF style fill. Traditionally, set fill mode on file and use _FillValue attribute.
  union {
    char cfill;
    unsigned char ucfill;
    short sfill;
    unsigned short usfill;
    int ifill;
    unsigned int uifill;
    long long lfill;
    unsigned long long ulfill;
    float ffill;
    double dfill;
  };
  
  // HDF compression feature.
  int shuffle;
  // HDF compression feature.
  int deflate;
  // HDF compression feature.
  int deflate_level;
  // HDF feature.
  int checksum;
  // HDF feature.
  int endian;
  
  ezNcVar();
  ~ezNcVar();
  int getNumAtts() { return attnames.size(); }
  int getNumDims() { return dimids.size(); }
  // Returns true if att name is in list.
  bool hasAtt(const char* name);
  bool hasAtt(std::string& name);
  void reset();
};
//###########################################################################
ezNcVar::ezNcVar() {
  reset();
};
//###########################################################################
ezNcVar::~ezNcVar() {
  reset();
};
//###########################################################################
bool ezNcVar::hasAtt(std::string& name) {
  return hasAtt(name.c_str());
}
//###########################################################################
bool ezNcVar::hasAtt(const char* name) {
  if (name == 0) return false;
  
  return ez::find_first(attnames, name, true) != -1;
}
//###########################################################################
void ezNcVar::reset() {
  nc = 0;
  varid = -1;
  recsize = 0;
  recbytes = 0;
  size = 0;
  bytes = 0;
  hasrec = false;
  attnames.clear();
  dimids.clear();
  chunkstorage = NC_CONTIGUOUS;
  chunksizes.clear();
  chunksize = 0;
  chunknelems = 0;
  chunkpreemption = 0;
  no_fill = true;
  ulfill = 0;
  shuffle = 0;
  deflate = 0;
  deflate_level = 0;
  checksum = 0;
  endian = NC_ENDIAN_NATIVE;
};
//###########################################################################
#ifdef DEBUG
void print(ezNcVar & v) {
  std::stringstream ss;
  ss << "Var name=" << v.name << ", varid=" << v.varid << ", hasrec=" << v.hasrec << ", size=" << v.size << ", sizePerRecord=" << v.sizePerRecord << ", type=" << v.type << "\n";
  std::cout << ss;
};
#endif
}
#endif // EZNCVAR_H