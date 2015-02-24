/*
20111125 rsz Created.
*/
#ifndef EZNCCOMPARE_H
#define EZNCCOMPARE_H

#include <ezNc.hpp>
#include <ezNcSlicer.hpp>
#include <ezNcSlabIterator.hpp>

namespace ez {
// Simple comparison implementation for unit tests.
// Compares between two files. Two compare within same file, supply same filename twice to 'init'.
class ezNcCompare {
public:
  // File wrappers.
  ezNc nc1, nc2;
  
  // Resets to empty state.
  inline void clear() { nc1.clear(); nc2.clear(); }
  // Compares size of dims for all dims in files. Both must have same dims.
  inline bool dimEquals();
  // Compares size of dims.
  inline bool dimEquals(std::string& dimname1, std::string& dimname2);
  // Compares all metadata and data.
  inline bool equals();
  // If formats equal.
  inline bool formatEquals() { return nc1.format == nc2.format; }
  // If all global atts equal.
  inline bool globalAttEquals();
  // If global atts equal.
  inline bool globalAttEquals(std::string& name1, std::string& name2);
  // Opens the files for comparison operations. Returns 0 if files opened ok.
  inline int init(std::string& filename1, std::string& filename2);
  // Compares formats, dims, global atts, var atts, var types, var shapes.
  inline bool metaDataEquals();
  // If all atts equal between two vars.
  inline bool varAttEquals(std::string& name1, std::string& name2);
  // If var atts equal.
  inline bool varAttEquals(std::string& vname1, std::string& aname1, std::string& vname2, std::string& aname2);
  // True if all data elements for all vars are equal. Both files must have same number of vars with identical names.
  inline bool varDataEquals();
  // True if all data elements are equal.
  inline bool varDataEquals(std::string& varname1, std::string& varname2);
  // If all vars meta data equal.
  inline bool varMetaDataEquals();
  // If metadata for two vars equal.
  inline bool varMetaDataEquals(std::string& varname1, std::string& varname2);
};
//######################################################################
bool ezNcCompare::dimEquals() {
  std::vector<std::string> names1, names2;
  nc1.getAllDimNames(names1);
  nc2.getAllDimNames(names2);
  std::sort(names1.begin(), names1.end());
  std::sort(names2.begin(), names2.end());
  if (names1 != names2) return false;
  
  int i=0;
  int n = names1.size();
  for(; i < n; ++i)
    if(!dimEquals(names1[i], names2[i])) return false;
    
  return true;
}
//######################################################################
bool ezNcCompare::dimEquals(std::string& dimname1, std::string& dimname2) {
  ezNcDim *dim1, *dim2;
  dim1 = nc1.getDim(dimname1);
  dim2 = nc2.getDim(dimname2);
  if (!dim1 && !dim2) return true;
  return dim1->size == dim2->size;
}
//######################################################################
bool ezNcCompare::equals() {
  return metaDataEquals() && varDataEquals();
}
//######################################################################
bool ezNcCompare::globalAttEquals() {
  std::vector<std::string> names1, names2;
  nc1.getAllGlobalAttNames(names1);
  nc2.getAllGlobalAttNames(names2);
  std::sort(names1.begin(), names1.end());
  std::sort(names2.begin(), names2.end());
  if (names1 != names2) return false;
  
  int i=0;
  int n = names1.size();
  for(; i < n; ++i)
    if(!globalAttEquals(names1[i], names2[i])) return false;
    
  return true;
}
//######################################################################
bool ezNcCompare::globalAttEquals(std::string& name1, std::string& name2) {
  ezNcBuffers buf1, buf2;
  nc1.getAtt(0, name1.c_str(), buf1);
  nc2.getAtt(0, name2.c_str(), buf2);
  return buf1.equals(buf2);
}
//######################################################################
int ezNcCompare::init(std::string& filename1, std::string& filename2) {
  int status;

  nc1.clear();
  nc2.clear();
  nc1.filename = filename1;
  nc2.filename = filename2;

  if(status = nc1.openReadOnly()) return status;
  if(status = nc1.loadAllMetadata()) return status;
  if(status = nc2.openReadOnly()) return status;
  if(status = nc2.loadAllMetadata()) return status;
  
  return 0;
}
//######################################################################
bool ezNcCompare::metaDataEquals() {
  return formatEquals() && dimEquals() && globalAttEquals() && varMetaDataEquals();
}
//######################################################################
bool ezNcCompare::varAttEquals(std::string& name1, std::string& name2) {
  ezNcVar *var1, *var2;
  var1 = nc1.getVar(name1);
  var2 = nc2.getVar(name2);
  if (!var1 && !var2) return true;
    
  if (var1->attnames != var2->attnames) return false;
  
  int i = 0;
  int n = var1->attnames.size();
  for(; i < n; ++i)
    if (!varAttEquals(name1, var1->attnames[i], name2, var2->attnames[i])) return false;
    
  return true;
}
//######################################################################
bool ezNcCompare::varAttEquals(std::string& vname1, std::string& aname1, std::string& vname2, std::string& aname2) {
  ezNcVar *var1, *var2;
  var1 = nc1.getVar(vname1);
  var2 = nc2.getVar(vname2);
  if (!var1 && !var2) return true;

  ezNcBuffers buf1, buf2;
  nc1.getAtt(vname1.c_str(), aname1.c_str(), buf1);
  nc2.getAtt(vname2.c_str(), aname2.c_str(), buf2);
  return buf1.equals(buf2);  
}
//######################################################################
bool ezNcCompare::varDataEquals() {
  std::vector<std::string> names1, names2;
  nc1.getAllVarNames(names1);
  nc2.getAllVarNames(names2);
  std::sort(names1.begin(), names1.end());
  std::sort(names2.begin(), names2.end());
  if (names1 != names2) return false;
  
  int n = names1.size();
  for(int i=0; i < n; ++i)
   if(!varDataEquals(names1[i], names2[i])) return false;
   
  return true;
}
//######################################################################
bool ezNcCompare::varDataEquals(std::string& varname1, std::string& varname2) {
  ezNcVar *var1, *var2;
  var1 = nc1.getVar(varname1);
  var2 = nc2.getVar(varname2);
  if (!var1 && !var2) return true;
  if (var1->type != var2->type) return false;
  if (var1->getNumDims() != var2->getNumDims()) return false;
  if (var1->size != var2->size) return false;
  
  ezNcBuffers buf1, buf2;
  nc1.reserve(buf1, var1->varid);
  nc2.reserve(buf2, var2->varid);
  buf1.allocate();
  buf2.allocate();

  ezNcSlicer<size_t> slicer;
  slicer.build(&nc1);
  
  ezNcSlabIterator<size_t> it;
  ezNcSlab<size_t> slab;
  it.setVar(var1, &slicer, true);
  while ( !it.isDone() ) {
    it.getSlab(&slab);
    nc1.read(var1->varid, &buf1, &(slab.start[0]), &(slab.count[0]));
    nc2.read(var2->varid, &buf2, &(slab.start[0]), &(slab.count[0]));
    if (!buf1.equals(buf2)) return false;
    it.inc();
  }
  
  return true;
}
//######################################################################
bool ezNcCompare::varMetaDataEquals() {
  std::vector<std::string> names1, names2;
  nc1.getAllVarNames(names1);
  nc2.getAllVarNames(names2);
  std::sort(names1.begin(), names1.end());
  std::sort(names2.begin(), names2.end());
  if (names1 != names2) return false;
  
  int n = names1.size();
  for(int i=0; i < n; ++i)
   if(!varMetaDataEquals(names1[i], names2[i])) return false;
   
  return true;
}
//######################################################################
bool ezNcCompare::varMetaDataEquals(std::string& varname1, std::string& varname2) {
  ezNcVar *var1, *var2;
  var1 = nc1.getVar(varname1);
  var2 = nc2.getVar(varname2);
  if (!var1 && !var2) return true;
  if (var1->type != var2->type) return false;
  if (var1->getNumDims() != var2->getNumDims()) return false;
  if (var1->size != var2->size) return false;
  return varAttEquals(varname1, varname2);
}
};

#endif // EZNCCOMPARE_H