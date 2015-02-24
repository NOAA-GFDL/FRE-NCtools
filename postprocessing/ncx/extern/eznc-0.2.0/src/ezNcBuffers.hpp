/*
CHANGELOG
v0.0.0 20110523 rsz Created.
v0.1.0 20111013 rsz Expanded datatypes to cover HDF4 additional native datatypes.
*/
#ifndef EZNCBUFFERS_H
#define EZNCBUFFERS_H

#include <stdio.h>
#include <algorithm> // std::fill
#include <string.h>
#include <ezNcVariant.hpp>
#include <argsort.hpp>

namespace ez {
//###########################################################################
// Contains pointers to allocated arrays for all netcdf datatypes for memory management of every datatype needed at runtime (not just one datatype, hence no use of union, however use only the types you need for max flexibility).
class ezNcBuffers {
public:
  char * c;
  short * s;
  int * i;
  long long * l;
  float * f;
  double * d;
  unsigned char * uc;
  unsigned short * us;
  unsigned int * ui;
  unsigned long long * ul;
  
  // Size of the arrays;
  size_t nc, ns, ni, nl, nf, nd, nuc, nus, nui, nul;
  
  ezNcBuffers() : nc(0), ns(0), ni(0), nl(0), nf(0), nd(0), nuc(0), nus(0), nui(0), nul(0), c(0), s(0), i(0), l(0), f(0), d(0), uc(0), us(0), ui(0), ul(0) {};
  ~ezNcBuffers();
  // Creates arrays for any non-zero sizes.
  void allocate();
  // True if arrays for type are equal. If no type given, compares all arrays.
  bool equals(ezNcBuffers& other, nc_type type=0);
  // Fill all arrays with a value.
  void fill(double, nc_type type=0);
  // Find indices where input list values (unsorted) equal an element in buffer array.
  // T is the value type, I is the index type.
  template <typename T, typename I>
  void find(nc_type type, std::vector<T>& list, std::vector<I>& indices);
  // Find index of first value greater-than-equal to a value in a buffer with same type as variant.
  int findFirstGreaterEqual(ezNcVariant&);
  int findFirstLessEqual(ezNcVariant&);
  int findLastGreaterEqual(ezNcVariant&);
  int findLastLessEqual(ezNcVariant&);
  // Check if an array's values are increasing. If all=false, only first 2 values are checked, otherwise all values are compared.
  bool increasing(nc_type, bool all=true);
  // Try to request raising ceiling size for type.
  void reserve(nc_type type, size_t size);
  // Clears and frees arrays.
  void reset();
  // Computes var element size excluding rec dimension and sets appropriate member (nc,ns,ni,etc.) if max for type.
  template <typename T>
  void setMaxSizePerRec(nc_type nctype, std::vector<T>& dimlens, std::vector<int>& dimids, int recid=-1);
  // Set number of elements for a specific type.
  void setSize(nc_type type, size_t size);
  // Fill all allocated arrays with zeros.
  void zeros();
};
//###########################################################################
ezNcBuffers::~ezNcBuffers() {
  reset();
}
//###########################################################################
void ezNcBuffers::allocate() {
  #define ALLOC(S, P, T) {\
    if (S) {\
      if (P)\
        delete [] P;\
      P = new T[S];\
    }\
  }
  
  ALLOC(nc,c,char);
  ALLOC(nuc,uc,unsigned char);
  ALLOC(ns,s,short);
  ALLOC(nus,us,unsigned short);
  ALLOC(ni,i,int);
  ALLOC(nui,ui,unsigned int);
  ALLOC(nl,l,long long);
  ALLOC(nul,ul,unsigned long long);
  ALLOC(nf,f,float);
  ALLOC(nd,d,double);
}
//###########################################################################
bool ezNcBuffers::equals(ezNcBuffers& other, nc_type type) {
  #define PEQ(P) {\
    if (n##P != other.n##P) return false;\
    for(int j=0; j < n##P; ++j) {\
      if (P[j] != other.P[j]) return false;\
    }\
  }
  
  switch(type) { 
    case NC_CHAR: PEQ(c); break;
    case NC_SHORT: PEQ(s); break;
    case NC_INT: PEQ(i); break;
    case NC_INT64: PEQ(l); break;
    case NC_FLOAT: PEQ(f); break;
    case NC_DOUBLE: PEQ(d); break;
    case NC_UBYTE: case NC_BYTE: PEQ(uc); break;
    case NC_USHORT: PEQ(us); break;
    case NC_UINT: PEQ(ui); break;
    case NC_UINT64: PEQ(ul); break;
    case 0:
       PEQ(c);
       PEQ(s);
       PEQ(i);
       PEQ(l);
       PEQ(f);
       PEQ(d);
       PEQ(uc);
       PEQ(us);
       PEQ(ui);
       PEQ(ul);
       break;
    default: break;
  }
  
  return true;
}
//###########################################################################
void ezNcBuffers::fill(double v, nc_type type) {
  switch(type) {
    case NC_CHAR: 
      if (nc) std::fill(c,c+nc,(char)v); 
      break;
    case NC_UBYTE: case NC_BYTE: 
      if (nuc) std::fill(uc,uc+nuc,(unsigned char)v); 
      break;
    case NC_SHORT: 
      if (ns) std::fill(s,s+ns,(short)v); 
      break;
    case NC_USHORT: 
      if (nus) std::fill(us,us+nus,(unsigned short)v); 
      break;
    case NC_INT:
      if (ni) std::fill(i,i+ni,(int)v); 
      break;
    case NC_UINT:
      if (nui) std::fill(ui,ui+nui,(unsigned int)v); 
      break;
    case NC_INT64:
      if (nl) std::fill(l,l+nl,(long long)v); 
      break;
    case NC_UINT64:
      if (nul) std::fill(ul,ul+nul,(unsigned long long)v); 
      break;
    case NC_FLOAT:
      if (nf) std::fill(f,f+nf,v); 
      break;
    case NC_DOUBLE:
      if (nd) std::fill(d,d+nd,v); 
      break;
    case 0:
      if (nc) std::fill(c,c+nc,(char)v);
      if (nuc) std::fill(uc,uc+nuc,(unsigned char)v);
      if (ns) std::fill(s,s+ns,(short)v);
      if (nus) std::fill(us,us+nus,(unsigned short)v);
      if (ni) std::fill(i,i+ni,(int)v);
      if (nui) std::fill(ui,ui+nui,(unsigned int)v);
      if (nl) std::fill(l,l+nl,(long long)v);
      if (nul) std::fill(ul,ul+nul,(unsigned long long)v);
      if (nf) std::fill(f,f+nf,v);
      if (nd) std::fill(d,d+nd,v);
      break;
    default: break;
  }
};
//###########################################################################
template <typename T, typename I>
void ezNcBuffers::find(nc_type type, std::vector<T>& list, std::vector<I>& indices) {
  indices.clear();
  int n = list.size();
  #define SORTFIND(D,P) {\
    std::vector<size_t> idx;\
    std::vector<D> unsorted(P,P+n##P);\
    std::vector<D> sorted;\
    argsort::sort(unsorted, sorted, idx);\
    indices.resize(n);\
    std::vector<D>::iterator b = sorted.begin(), e = sorted.end();\
    for(int i=0; i < n; ++i) {\
      indices[i] = idx[ std::lower_bound(b,e,list[i]) - b ];\
    }\
  }
  
  switch(type) { 
    case NC_CHAR: SORTFIND(char,c); break;
    case NC_SHORT: SORTFIND(short,s); break;
    case NC_INT: SORTFIND(int,i); break;
    case NC_INT64: SORTFIND(long long, l); break;
    case NC_FLOAT: SORTFIND(float, f); break;
    case NC_DOUBLE: SORTFIND(double, d); break;
    case NC_UBYTE: case NC_BYTE: SORTFIND(unsigned char, uc); break;
    case NC_USHORT: SORTFIND(unsigned short, us); break;
    case NC_UINT: SORTFIND(unsigned int, ui); break;
    case NC_UINT64: SORTFIND(unsigned long long, ul); break;
    default: break;
  }
}
//###########################################################################
int ezNcBuffers::findFirstGreaterEqual(ezNcVariant& v) {
  #define FGELOOP(T) {\
    for(int j=0; j < n##T; ++j) if (T[j] >= v.T) return j;\
    return n##T-1;\
  }

  switch(v.type) { 
    case NC_CHAR: FGELOOP(c); break;
    case NC_SHORT: FGELOOP(s); break;
    case NC_INT: FGELOOP(i); break;
    case NC_INT64: FGELOOP(l); break;
    case NC_FLOAT: FGELOOP(f); break;
    case NC_DOUBLE: FGELOOP(d); break;
    case NC_UBYTE: case NC_BYTE: FGELOOP(uc); break;
    case NC_USHORT: FGELOOP(us); break;
    case NC_UINT: FGELOOP(ui); break;
    case NC_UINT64: FGELOOP(ul); break;
    default: break;
  }
  return 0;
}
//###########################################################################
int ezNcBuffers::findFirstLessEqual(ezNcVariant& v) {
  #define FLELOOP(T) {\
    for(int j=0; j < n##T; ++j) if (T[j] <= v.T) return j;\
    return n##T-1;\
  }

  switch(v.type) { 
    case NC_CHAR: FLELOOP(c); break;
    case NC_SHORT: FLELOOP(s); break;
    case NC_INT: FLELOOP(i); break;
    case NC_INT64: FLELOOP(l); break;
    case NC_FLOAT: FLELOOP(f); break;
    case NC_DOUBLE: FLELOOP(d); break;
    case NC_UBYTE: case NC_BYTE: FLELOOP(uc); break;
    case NC_USHORT: FLELOOP(us); break;
    case NC_UINT: FLELOOP(ui); break;
    case NC_UINT64: FLELOOP(ul); break;
    default: break;
  }
  return 0;
}
//###########################################################################
int ezNcBuffers::findLastGreaterEqual(ezNcVariant& v) {
  #define LGELOOP(T) {\
    for(int j=n##T; j --> 0 ;) if (T[j] >= v.T) return j;\
    return 0;\
  }

  switch(v.type) { 
    case NC_CHAR: LGELOOP(c); break;
    case NC_SHORT: LGELOOP(s); break;
    case NC_INT: LGELOOP(i); break;
    case NC_INT64: LGELOOP(l); break;
    case NC_FLOAT: LGELOOP(f); break;
    case NC_DOUBLE: LGELOOP(d); break;
    case NC_UBYTE: case NC_BYTE: LGELOOP(uc); break;
    case NC_USHORT: LGELOOP(us); break;
    case NC_UINT: LGELOOP(ui); break;
    case NC_UINT64: LGELOOP(ul); break;
    default: break;
  }
  return 0;
}
//###########################################################################
int ezNcBuffers::findLastLessEqual(ezNcVariant& v) {
  #define LLELOOP(T) {\
    for(int j=n##T; j --> 0 ;) if (T[j] <= v.T) return j;\
    return 0;\
  }

  switch(v.type) { 
    case NC_CHAR: LLELOOP(c); break;
    case NC_SHORT: LLELOOP(s); break;
    case NC_INT: LLELOOP(i); break;
    case NC_INT64: LLELOOP(l); break;
    case NC_FLOAT: LLELOOP(f); break;
    case NC_DOUBLE: LLELOOP(d); break;
    case NC_UBYTE: case NC_BYTE: LLELOOP(uc); break;
    case NC_USHORT: LLELOOP(us); break;
    case NC_UINT: LLELOOP(ui); break;
    case NC_UINT64: LLELOOP(ul); break;
    default: break;
  }
  return 0;
}
//###########################################################################
bool ezNcBuffers::increasing(nc_type type, bool all) {
  int j=0, n;
  #define UPLOOP(T) {\
  if (all)\
    n = n##T-1;\
  else\
    n = n##T > 1 ? 1 : 0;\
  \
  for(; j < n; ++j) if (T[j] >= T[j+1]) return false;\
  }
  
  switch(type) { 
    case NC_CHAR: UPLOOP(c); break;
    case NC_SHORT: UPLOOP(s); break;
    case NC_INT: UPLOOP(i); break;
    case NC_INT64: UPLOOP(l); break;
    case NC_FLOAT: UPLOOP(f); break;
    case NC_DOUBLE: UPLOOP(d); break;
    case NC_UBYTE: case NC_BYTE: UPLOOP(uc); break;
    case NC_USHORT: UPLOOP(us); break;
    case NC_UINT: UPLOOP(ui); break;
    case NC_UINT64: UPLOOP(ul); break;
    default: break;
  }
  return true;
}
//###########################################################################
void print(ezNcBuffers & b) {
  printf("ezNcBuffer\n");
  printf("  nc=%zu, c is NULL? %d\n", b.nc, b.c==0);
  printf("  nuc=%zu, uc is NULL? %d\n", b.nuc, b.uc==0);
  printf("  ns=%zu, s is NULL? %d\n", b.ns, b.s==0);
  printf("  nus=%zu, us is NULL? %d\n", b.nus, b.us==0);
  printf("  ni=%zu, i is NULL? %d\n", b.ni, b.i==0);
  printf("  nui=%zu, ui is NULL? %d\n", b.nui, b.ui==0);
  printf("  nl=%zu, l is NULL? %d\n", b.nl, b.l==0);
  printf("  nul=%zu, ul is NULL? %d\n", b.nul, b.ul==0);
  printf("  nf=%zu, f is NULL? %d\n", b.nf, b.f==0);
  printf("  nd=%zu, d is NULL? %d\n", b.nd, b.d==0);
  
  int i;
  if (b.nf) {
    printf("  f=");
    for(i=0; i < b.nf; ++i)
      printf("%g ", b.f[i]);
    printf("\n");
  }
  
  if (b.nd) {
    printf("  d=");
    for(i=0; i < b.nd; ++i)
      printf("%g ", b.d[i]);
    printf("\n");
  }
};
//################################################################
void ezNcBuffers::reserve(nc_type type, size_t size) {
  #define MAX(S) if (size > S) S = size; 
  switch(type) { 
    case NC_CHAR: MAX(nc); break;
    case NC_SHORT: MAX(ns); break;
    case NC_INT: MAX(ni); break;
    case NC_INT64: MAX(nl); break;
    case NC_FLOAT: MAX(nf); break;
    case NC_DOUBLE: MAX(nd); break;
    case NC_UBYTE: case NC_BYTE: MAX(nuc); break;
    case NC_USHORT: MAX(nus); break;
    case NC_UINT: MAX(nui); break;
    case NC_UINT64: MAX(nul); break;
    default: break;
  }
};
//###########################################################################
void ezNcBuffers::reset() {
  #define FREE(S, P) if (P) { delete [] P; P = 0; S = 0; }

  FREE(nc,c);
  FREE(nuc,uc);
  FREE(ns,s);
  FREE(nus,us);
  FREE(ni,i);
  FREE(nui,ui);
  FREE(nl,l);
  FREE(nul,ul);
  FREE(nf,f);
  FREE(nd,d);
};
//###########################################################################
template <typename T>
void ezNcBuffers::setMaxSizePerRec(nc_type nctype, std::vector<T> & dimlens, std::vector<int> & dimids, int recid) {
  T total = 1;
  int i=0;
  int ndims = dimids.size();
  
  for(; i < ndims; ++i) {
    if (recid == dimids[i]) continue;
    total *= dimlens[i];
  }
  
  if ( (ndims == 1) && (recid == dimids[0]) )
    total = 1;
  
  reserve(nctype, total);
}    
//################################################################
void ezNcBuffers::setSize(nc_type type, size_t size) {
  switch(type) { 
    case NC_CHAR: nc = size; break;
    case NC_SHORT: ns = size; break;
    case NC_INT: ni = size; break;
    case NC_INT64: nl = size; break;
    case NC_FLOAT: nf = size; break;
    case NC_DOUBLE: nd = size; break;
    case NC_UBYTE: case NC_BYTE: nuc = size; break;
    case NC_USHORT: nus = size; break;
    case NC_UINT: nui = size; break;
    case NC_UINT64: nul = size; break;
    default: break;
  }
};
//################################################################
void ezNcBuffers::zeros() {
  if (nc) memset(c,0,nc);
  if (nuc) memset(uc,0,nuc);
  if (ns) memset(s,0,ns*sizeof(short));
  if (nus) memset(us,0,nus*sizeof(unsigned short));
  if (ni) memset(i,0,ni*sizeof(int));
  if (nui) memset(ui,0,nui*sizeof(unsigned int));
  if (nl) memset(l,0,nl*sizeof(long long));
  if (nul) memset(ul,0,nul*sizeof(unsigned long long));
  if (nf) std::fill(f,f+nf,0.0);
  if (nd) std::fill(d,d+nd,0.0);
};
}

#endif // EZNCBUFFERS_H