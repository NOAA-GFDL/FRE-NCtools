/*
CHANGELOG
v0.0.0 20110523 rsz Created.
v0.1.0 20111120 rsz Housekeeping.
*/
#ifndef EZNC_H
#define EZNC_H

#include "netcdf.h"
#include <ezNcDim.hpp>
#include <ezNcVar.hpp>
#include <ezNcBuffers.hpp>
#include <ezStringUtil.hpp>
#include <ezNcUtil.hpp>
#include <map>
#include <vector>

#ifdef __DEBUG__
#define DEBUGLINE() fprintf(stderr, "%s:%d\n", __FILE__, __LINE__);
#endif

namespace ez {
typedef std::map<std::string, ezNcVar*> StringToVarMap;
typedef std::map<std::string, ezNcDim*> StringToDimMap;
typedef std::vector<ezNcVar*> VarVector;
typedef std::vector<ezNcDim*> DimVector;

class ezNc {
public:
  // File handle. -1 if none.
  int ncid;
  // The file path and name set by user.
  std::string filename;
  // Flavor of netcdf.
  int format;
  // Id for record dim in this file. -1 if missing.
  int recid;
  // Number of records.
  size_t nrec;
  // Initial size of file if creating new file.
  size_t initialsz;
  // Request desired blocksize for I/O.
  size_t bufrsizehint;
  // Sets the pad at the end of the "header" section.
  size_t h_minfree;
  // Controls the alignment of the beginning of the data section for fixed size variables.
  size_t v_align;
  // Sets the pad at the end of the data section for fixed size variables.
  size_t v_minfree;
  // Controls the alignment of the beginning of the data section for variables which have an unlimited dimension (record variables).
  size_t r_align;
  // Enable prefill of variables. Default is off.
  bool fill;
  
  // HDF feature, to be set before file open or create.
  size_t chunksize;
  // HDF feature, to be set before file open or create.
  size_t chunknelems;
  // HDF feature, to be set before file open or create.
  float chunkpreemption;
  
  // Maps names to variables.
  StringToVarMap varsNameMap;
  // Maps names to dimensions.
  StringToDimMap dimsNameMap;
  // Pointers to the vars indexed by their varid.
  VarVector varsIdVector;
  // Pointers to the dims indexed by their dimid.
  DimVector dimsIdVector;

  // Add newly created object to lookup data structures.
  inline void add(ezNcDim* dim);
  // Add newly created object to lookup structures.
  inline void add(ezNcVar* var);
  // Clears all data structures.
  inline void clear();
  // Closes file. Return 0 on success.
  inline int close();
  // Sets size info into var.
  inline void computeSize(ezNcVar*);
  // Sets size info into var.
  inline void computeSize(int varid);
  // Sets size info into var.
  inline void computeSize(std::string& varname);
  // Copy att from another file into this one. If varname is NULL, then att is global.
  inline int copyAtt(ezNc& fromnc, const char* fromvarname, const char* tovarname, const char* attname);
  // Copy all atts for a var or all global atts from one file into this one.
  inline int copyAtts(ezNc & fromnc, const char* fromvarname, const char* tovarname);
  // Copy a dimension from another file. dimid is not expected to be identical. Returns dimid or -1 if error.
  inline int copyDimDef(ezNc& fromnc, const char* name);
  // Copy dim using dimid.
  inline int copyDimDef(ezNc & fromnc, int dimid);
  // Copies all the dim defs.
  inline int copyDimDefs(ezNc& fromnc);
  // Copy file characteristics, but not any dims, vars, or atts.
  inline void copySettings(ezNc & from);
  // Copy a var from another file. varid is not expected to be identical. Returns varid or -1 if error. Dims for new var should already exist.
  inline int copyVarDef(ezNc& fromnc, const char* name);
  // Copy using varid.
  inline int copyVarDef(ezNc & fromnc, int varid);
  // Copies all the var defs.
  inline int copyVarDefs(ezNc& fromnc);
  // Copies all the var defs for names listed.
  inline int copyVarDefs(ezNc& fromnc, std::vector<std::string>& names);
  // Creates 64bit netcdf3 style file.
  inline int create64();
  inline int create64Clobber();
  // Creates HDF style file.
  inline int createHDF();
  inline int createHDFClobber();
  // Creates new file. Returns 0 if successful. mode is by default to clobber and use classic format. Applies initialsz and bufrsize tuning parameters.
  inline int create(int mode=0);
  // Creates classic netcdf3 style file.
  inline int createClobber();
  // Define dim in file.
  inline int create(ezNcDim *);
  // Define var in file.
  inline int create(ezNcVar *);
  // Returns dimid or -1 if error.
  inline int createDim(const char* name, size_t size, bool isRecord=false);
  // Returns varid or -1 on error.
  inline int createVar(const char* name, nc_type type, int ndims, int* dimids, int storage=0, size_t* chunksizes=0, size_t size=0, size_t nelems=0, float preemption=0.0, int no_fill=1, void* fill_value=0, int shuffle=0, int deflate=0, int deflate_level=0, int checksum=0, int endian=NC_ENDIAN_NATIVE);
  // End definition mode.
  inline int enddef();
  // Cter.
  ezNc() { ncid = -1; clear(); }
  // Dter.
  ~ezNc() { clear(); }
  // Get sorted dimids.
  inline void getAllDimIds(std::vector<int> &);
  // Get all names of dims sorted by dimids.
  inline void getAllDimNames(std::vector<std::string> &);
  // Gets list of names sorted as stored in file.
  inline void getAllGlobalAttNames(std::vector<std::string> &);
  // Return size in bytes for all vars, static and record.
  inline size_t getAllVarBytes();
  // Get sorted ids of vars (both static and record based).
  inline void getAllVarIds(std::vector<int> &);
  // Get all names of vars (both static and record based) sorted by varids.
  inline void getAllVarNames(std::vector<std::string> &);
  // Get var or global att value.
  inline int getAtt(const char* varname, const char* attname, ezNcBuffers& values);
  // Gets string contents for att, if exists. If not, returns empty string.
  inline int getAttString(const char* varname, const char* attname, std::string& string);
  // Get all var names that have the same name as var's dim names (CF convention). Use varid and get varids.
  inline void getCoordVars(int varid, std::vector<int>& varids);
  // Use varname and return var names.
  inline void getCoordVars(std::string& varname, std::vector<std::string>& varnames);
  // Get pointer from id by safely checking if id is valid. Returns 0 if invalid.
  ezNcDim* getDim(unsigned int id) { if (id >= dimsIdVector.size()) return 0; else return dimsIdVector[id]; }
  ezNcDim* getDim(const char* name) { if (name==0) return 0; return dimsNameMap[name]; }
  ezNcDim* getDim(std::string& name) { return dimsNameMap[name]; }
  // Get dimids used by all vars listed in varid vector.
  inline void getDimIds(std::set<int> & varids, std::set<int> & dimids);
  inline void getDimIds(std::vector<int> & varids, std::set<int> & dimids);
  // Get list of dimids in same order as names list.
  inline void getDimIds(std::vector<std::string> & dimnames, std::vector<int> & dimids);
  // Get dims names in proper order for a var.
  inline void getDimNames(int varid, std::vector<std::string> & dims);
  // Return vector of unique dim names used by each var listed.
  inline void getDimNamesForVars(std::vector<std::string> & vars, std::set<std::string> & dims);
  // Get dim sizes for a var, aka shape.
  inline void getDimSizes(int varid, std::vector<size_t> &);
  // Get dim sizes for a var, aka shape.
  inline void getDimSizes(std::string & name, std::vector<size_t> &);
  // Returns max number of dimensions used by list of var names.
  inline int getMaxNumDims(std::vector<std::string> &);
  // List of dims in the file.
  inline int getNumDims() { return dimsIdVector.size(); };
  // Get total number of variables.
  inline int getNumVars() { return getNumStaticVars() + getNumRecVars(); };
  // Get total number of static variables.
  inline int getNumStaticVars();
  // Get total number of record variables.
  inline int getNumRecVars();
  // Return size in bytes for all record vars.
  inline size_t getRecVarBytes();
  // Get pointer from id by safely checking if id is valid. Returns 0 if invalid.
  inline ezNcVar* getVar(unsigned int id);
  inline ezNcVar* getVar(const char* name);
  inline ezNcVar* getVar(std::string& name);
  // Get var id from name. -1 if invalid.
  inline int getVarId(const char* name);
  inline int getVarId(const std::string& name);
  // Returns list of record var names, sorted by varid.
  inline void getRecVarNames(std::vector<std::string> &);
  // Returns list of record var ids sorted.
  inline void getRecVarIds(std::vector<int> &);
  // Return size in bytes for all static vars.
  inline size_t getStaticVarBytes();
  // Returns list of non-record var names, sorted by varid.
  inline void getStaticVarNames(std::vector<std::string> &);
  // Returns list of non-record var ids sorted.
  inline void getStaticVarIds(std::vector<int> &);
  // Get size of each dimension size assuming tzyx dimension order for typical Coards var.
  inline void getVarSizeTZYX(int id, size_t & nt, size_t & nz, size_t & ny, size_t & nx);
  // Get size of each dimension size assuming tzyx dimension order for typical Coards var.
  inline void getVarSizeTZYX(std::string& name, size_t & nt, size_t & nz, size_t & ny, size_t & nx);
  // Get vector of ezNcVar pointers that exactly match (or exclude) names or patterns, optionally with their coord vars (vars that have same name as their dims), all sorted by their varids.
  inline void getVars(std::vector<ezNcVar*>& out, std::vector<std::string>& patterns, bool regex, bool regex_nocase, bool coordvars, bool exclude);
  // Check nc return status and print message if an error.
  inline int handle_status(int);
  // True if file has record dimension.
  inline bool hasrec() const { return recid != -1; }
  // Return true if underlying format is HDF flavor (false if nc3 classic or 64bit).
  inline bool isHDF();
  // Load all the dim metadata from file. Returns false if failed.
  inline int loadDimsMetadata();
  // Loads globals, dim, and var metadata. Returns false if failed.
  inline int loadAllMetadata();
  // Load all the var metadata from file. Returns false if failed.
  inline int loadVarsMetadata();
  // Sets filename, format and ncid. Uses buffersizehint.
  inline int openReadOnly();
  // For updating metadata or appending variables or records. Uses buffersizehint.
  inline int openWritable();
  // Write attribute data. If name is null, then writes global attribute.
  inline int putAtt(const char* varname, const char* attname, nc_type type, std::vector<std::string> & values, int first, int last);
  // More direct and efficient write of att data.
  inline int putAtt(const char* varname, const char* attname, nc_type type, ezNcBuffers& values);
  // Read portion of variable.
  inline int read(const char* name, ezNcBuffers*, size_t * start, size_t * count);
  inline int read(int varid, ezNcBuffers*, size_t * start, size_t * count);
  // Read entire variable.
  inline int read(int varid, ezNcBuffers*);

  // Configure buffers for all vars in file.
  inline int reserve(ezNcBuffers &);
  // Configure buffers for var id for just single record. Doesn't reset or allocate buffer.
  inline int reserve(ezNcBuffers &, int varid);
  // Configure buffers for var with name. Doesn't reset or allocate buffer.
  inline int reserve(ezNcBuffers & buffers, const char* varname);
  // Configure buffers for list of var ids. Doesn't reset or allocate buffer.
  inline int reserve(ezNcBuffers &, std::vector<int> & );
  // Configure buffers for list of var names. Doesn't reset or allocate buffer.
  inline int reserve(ezNcBuffers &, std::vector<std::string> & names);
  // Configure buffers for var id for entire var shape (all dims including all records). Doesn't reset or allocate buffer.
  inline int reserveEntire(ezNcBuffers &, int varid);

  // Returns human readable string for format type.
  inline static void typeToString(nc_type, std::string&);
  // Write portion of variable.
  inline int write(const char* name, ezNcBuffers*, size_t * start, size_t * count);
  inline int write(int varid, ezNcBuffers*, size_t * start, size_t * count);
};
//###########################################################################
void ezNc::add(ezNcDim* dim) {
  if (dim == 0) return;
  
  dim->nc = this;
  dimsNameMap[dim->name] = dim;
  if (dimsIdVector.size() <= dim->dimid) dimsIdVector.resize(dim->dimid+1);
  dimsIdVector[dim->dimid] = dim;
}
//###########################################################################
void ezNc::add(ezNcVar* var) {
  if (var == 0) return;
  
  var->nc = this;
  varsNameMap[var->name] = var;
  if (varsIdVector.size() <= var->varid) varsIdVector.resize(var->varid+1);
  varsIdVector[var->varid] = var;
}
//###########################################################################
void ezNc::clear() {
  if (ncid != -1) 
    close();
  
  filename.clear();
  format = -1;
  recid = -1;
  nrec = 0;
  initialsz = 0;
  bufrsizehint = NC_SIZEHINT_DEFAULT;
  h_minfree = 0;
  v_align = 4;
  v_minfree = 0;
  r_align = 4;
  fill = false;
  chunksize = 0;
  chunknelems = 0;
  chunkpreemption = 0;
  
  varsNameMap.clear();
  dimsNameMap.clear();
  
  int n = varsIdVector.size();
  int i = 0;
  for(; i < n; ++i) {
    if (varsIdVector[i])
      delete varsIdVector[i];
  }
  varsIdVector.clear();

  n = dimsIdVector.size();
  for(i=0; i < n; ++i) {
    if (dimsIdVector[i])
      delete dimsIdVector[i];
  }
  dimsIdVector.clear();
};
//###########################################################################
int ezNc::close() {
  if (ncid != -1) {
    int status = nc_close(ncid);
    if(handle_status(status)) return status;
    ncid = -1;
  }
  return 0;
};
//###########################################################################
void ezNc::computeSize(int varid) {
  ezNcVar* v = getVar(varid);
  if (v == 0) return;
  computeSize(v);
}
//###########################################################################
void ezNc::computeSize(std::string& varname) {
  ezNcVar* v = getVar(varname);
  if (v == 0) return;
  computeSize(v);
}
//###########################################################################
void ezNc::computeSize(ezNcVar* var) {
  if (var==0) return;
  
  // Compute sizes for convenience.
  size_t total = 1;
  int i = 0;
  int ndims = var->getNumDims();
  ezNcDim* dim;
  
  for(; i < ndims; ++i) {
    if ((i==0) && (var->hasrec)) {
      continue;
    } else {
      dim = getDim(var->dimids[i]);
      if (dim)
        total *= dim->size;
    }
  }
  
  if (var->hasrec) {
    var->recsize = total;
    var->recbytes = total*sizeOf(var->type);
    var->size = total*nrec;
    var->bytes = var->recbytes*nrec;
  } else {
    var->recsize = var->size = total;
    var->recbytes = var->bytes = total*sizeOf(var->type);
  }
};
//###########################################################################
int ezNc::copyAtts(ezNc & fromnc, const char* fromvarname, const char* tovarname) {
  if (fromvarname) {
    ezNcVar* fromv = fromnc.varsNameMap[fromvarname];
    if (!fromv) return 1;
    int n = fromv->attnames.size();
    for(int i=0; i < n; ++i)
      copyAtt(fromnc, fromvarname, tovarname, fromv->attnames[i].c_str());
  } else {
    // Global atts.
    int ngatts;
    int status = nc_inq_natts(fromnc.ncid, &ngatts);
    if(handle_status(status)) return status;
    char attname[256];
    
    for(int i=0; i < ngatts; ++i) {
      status = nc_inq_attname(fromnc.ncid, NC_GLOBAL, i, attname);
      if(handle_status(status)) return status;
      status = copyAtt(fromnc, 0, tovarname, attname);
      if (status) return status;
    }
  }
  return 0;
}
//###########################################################################
int ezNc::copyAtt(ezNc & fromnc, const char* fromvarname, const char* tovarname, const char* attname) {
  int tovarid = NC_GLOBAL;
  int fromvarid = NC_GLOBAL;
  
  if (tovarname) {
    ezNcVar* v = varsNameMap[tovarname];
    if (!v) return 1;
    tovarid = v->varid;
  }
  
  if (fromvarname) {
    ezNcVar* fromv = fromnc.varsNameMap[fromvarname];
    if (!fromv) return 1;
    fromvarid = fromv->varid;
  }

  int status;
  status = nc_copy_att(fromnc.ncid, fromvarid, attname, ncid, tovarid);
  if(handle_status(status)) return status;
 
  return 0;
};
//###########################################################################
int ezNc::copyDimDef(ezNc & fromnc, const char* name) {
  ezNcDim *fromdim = fromnc.dimsNameMap[name];
  if (fromdim == 0) {
    fprintf(stderr, "ERROR:%s:%d Failed to find dim \"%s\".\n", __FILE__, __LINE__, name);
    return -1;
  }
  
  ezNcDim *dim = new ezNcDim;
  dim->name = name;
  dim->size = fromdim->size;
  dim->isrec = fromdim->isrec;

  return create(dim) == -1;
};
//###########################################################################
int ezNc::copyDimDef(ezNc & fromnc, int dimid) {
  ezNcDim *fromdim = fromnc.getDim(dimid);
  if (fromdim == 0) {
    fprintf(stderr, "ERROR:%s:%d Failed to find dimid \"%d\".\n", __FILE__, __LINE__, dimid);
    return -1;
  }
  
  return createDim(fromdim->name.c_str(), fromdim->size, fromdim->isrec)  == -1;
};
//###########################################################################
int ezNc::copyDimDefs(ezNc & from) {
  int n = from.dimsIdVector.size();
  for(int i=0; i < n; ++i) {
    copyDimDef(from, i);
  }
};
//###########################################################################
void ezNc::copySettings(ezNc & from) {
  ncid = from.ncid;
  filename = from.filename;
  format = from.format;
  recid = from.recid;
  nrec = from.nrec;
  initialsz = from.initialsz;
  bufrsizehint = from.bufrsizehint;
  h_minfree = from.h_minfree;
  v_align = from.v_align;
  v_minfree = from.v_minfree;
  r_align = from.r_align;
  fill = from.fill;
  chunksize = from.chunksize;
  chunknelems = from.chunknelems;
  chunkpreemption = from.chunkpreemption;
};
//###########################################################################
int ezNc::copyVarDef(ezNc & fromnc, const char* name) {
  ezNcVar *fromv = fromnc.varsNameMap[name];
  if (fromv == 0) {
    fprintf(stderr, "ERROR:%s:%d Failed to find var \"%s\".\n", __FILE__, __LINE__, name);
    return 1;
  }
  
  return copyVarDef(fromnc, fromv->varid);
}
//###########################################################################
int ezNc::copyVarDef(ezNc & fromnc, int varid) {
  ezNcVar *fromv = fromnc.getVar(varid);
  if (fromv == 0) {
    fprintf(stderr, "ERROR:%s:%d Failed to find varid \"%d\" in %s.\n", __FILE__, __LINE__, varid, fromnc.filename.c_str());
    return -1;
  }

  ezNcVar *v = new ezNcVar;
  v->name = fromv->name;
  v->nc = this;
  v->hasrec = fromv->hasrec;
  v->type = fromv->type;
  v->attnames = fromv->attnames;
  
  // Copy the dim dependencies.
  int ndims = fromv->getNumDims();
  int i=0;
  
  for(; i < ndims; ++i) {
    ezNcDim * dim1 = fromnc.getDim(fromv->dimids[i]);
    if(dim1 == 0) {
      fprintf(stderr, "ERROR:%s:%d Failed to find dim id %d.\n", __FILE__, __LINE__, fromv->dimids[i]);
      delete v;
      return -1;
    }
    
    ezNcDim * dim2 = dimsNameMap[dim1->name];
    if(dim2 == 0) {
      fprintf(stderr, "ERROR:%s:%d Failed to find dim \"%s\".\n", __FILE__, __LINE__, dim1->name.c_str());
      delete v;
      return -1;
    }
    v->dimids.push_back(dim2->dimid);
  }

  if (isHDF()) {
    v->chunkstorage = fromv->chunkstorage;
    v->chunksizes = fromv->chunksizes;
    v->no_fill = fromv->no_fill;
    v->lfill = fromv->lfill;
    v->shuffle = fromv->shuffle;
    v->deflate = fromv->deflate;
    v->deflate_level = fromv->deflate_level;
    v->checksum = fromv->checksum;
    v->endian = fromv->endian;
  }

  // This assigns a new varid, which is needed to copy atts, otherwise default (-1) will cause atts to be copied to global att list.
  if (create(v) == -1) {
    delete v;
    return -1;
  }

  // Copy the atts.
  int status;
  for(i=0; i < v->attnames.size(); ++i) {
    status = nc_copy_att(fromnc.ncid, fromv->varid, v->attnames[i].c_str(), ncid, v->varid);
    if (handle_status(status)) {
      delete v;
      return -1;
    }
  }
      
  return v->varid;
};
//###########################################################################
int ezNc::copyVarDefs(ezNc & from) {
  int n = from.varsIdVector.size();
  for(int i=0; i < n; ++i) {
    copyVarDef(from, i);
  }
};
//###########################################################################
int ezNc::copyVarDefs(ezNc& from, std::vector<std::string>& names) {
  int n = names.size();
  for(int i=0; i < n; ++i) {
    copyVarDef(from, names[i].c_str());
  }
}
//###########################################################################
int ezNc::createClobber() {
  format = NC_FORMAT_CLASSIC;
  return create(NC_CLOBBER);
};
//###########################################################################
int ezNc::create64() {
  format = NC_FORMAT_64BIT;
  return create(NC_64BIT_OFFSET | NC_NOCLOBBER);
};
//###########################################################################
int ezNc::create64Clobber() {
  format = NC_FORMAT_64BIT;
  return create(NC_CLOBBER|NC_64BIT_OFFSET);
};
//###########################################################################
int ezNc::createHDF() {
  format = NC_FORMAT_NETCDF4;
  return create(NC_NETCDF4 | NC_NOCLOBBER);
};
//###########################################################################
int ezNc::createHDFClobber() {
  format = NC_FORMAT_NETCDF4;
  return create(NC_CLOBBER|NC_NETCDF4);
};
//###########################################################################
int ezNc::create(int mode) {
  int status;
  
  if (chunksize > 0) {
    status = nc_set_chunk_cache(chunksize, chunknelems, chunkpreemption);
    if (handle_status(status)) return status;
  }

  status = nc__create(filename.c_str(), mode, initialsz, &bufrsizehint, &ncid);
  if (handle_status(status)) return status;
  
  int dummy;
  status = nc_set_fill(ncid, fill ? NC_FILL : NC_NOFILL, &dummy);
  return handle_status(status);
};
//###########################################################################
int ezNc::create(ezNcDim * dim) {
  if ((ncid == -1) || (dim == 0)) return -1;

  size_t len = NC_UNLIMITED;
  if (!dim->isrec)
    len = dim->size;
    
  int dimid;
  int status = nc_def_dim(ncid, dim->name.c_str(), len, &dimid);
  if (handle_status(status)) return -1;
  
  dim->dimid = dimid;
  add(dim);
  
  return dimid;
};
//###########################################################################
int ezNc::create(ezNcVar * var) {
  if ( (ncid == -1) || (var == 0) || (var->getNumDims()<1) ) return -1;

  var->nc = this;
  int ndims = var->getNumDims();
  int status = nc_def_var(ncid, var->name.c_str(), var->type, ndims, &(var->dimids[0]), &(var->varid));
  if (handle_status(status)) return -1;
  
  computeSize(var);

  if (isHDF()) {

    if ( (var->chunkstorage != 0) && !var->chunksizes.empty() ) {
      status = nc_def_var_chunking(ncid, var->varid, var->chunkstorage, &(var->chunksizes[0]));
      if (handle_status(status)) return -1;
    }

    if (var->chunksize != 0) {
      status = nc_set_var_chunk_cache(ncid, var->varid, var->chunksize, var->chunknelems, var->chunkpreemption);
      if (handle_status(status)) return -1;
    }

    if (!var->no_fill) {
      status = nc_def_var_fill(ncid, var->varid, (int)var->no_fill, (void*)&(var->cfill));
      if (handle_status(status)) return -1;
    }

    if (var->shuffle) {
      status = nc_def_var_deflate(ncid, var->varid, var->shuffle, var->deflate, var->deflate_level);
      if (handle_status(status)) return -1;
    }

    if (var->checksum) {
      status = nc_def_var_fletcher32(ncid, var->varid, var->checksum);
      if (handle_status(status)) return -1;
    }

    if (var->endian) {
      status = nc_def_var_endian(ncid, var->varid, var->endian);
      if (handle_status(status)) return -1;
    }
  }
  
  add(var);
  
  return var->varid;
};
//###########################################################################
int ezNc::createDim(const char* name, size_t size, bool isRecord) {
  if ((ncid == -1) || (name == 0)) return -1;
  int dimid;
  int status;

  status = nc_def_dim(ncid, name, isRecord ? NC_UNLIMITED : size, &dimid);
  if (handle_status(status)) return -1;

  ezNcDim * dim = new ezNcDim;
  dim->name = name;
  dim->dimid = dimid;
  dim->size = size;
  
  if (isRecord) {
    dim->isrec = 1;
    recid = dimid;
    nrec = size;
  }
  
  add(dim);
  
  return dimid;
};
//###########################################################################
int ezNc::createVar(const char* name, nc_type type, int ndims, int* dimids, int storage, size_t* chunksizes, size_t chunksize, size_t chunknelems, float chunkpreemption, int no_fill, void* fill_value, int shuffle, int deflate, int deflate_level, int checksum, int endian) {
  if ( (ncid == -1) || (name == 0) || (ndims == 0) || (dimids == 0) ) return -1;

  ezNcVar *var = new ezNcVar;
  var->name = name;
  var->type = type;
  int i=0;
  ezNcDim * dim;
  
  for(; i < ndims; ++i) {
    var->dimids.push_back(dimids[i]);
    dim = getDim(dimids[i]);
    if (dim == 0) {
      fprintf(stderr, "ERROR: Dimension with id \"%d\" not found for var \"%s\".\n", dimids[i], name);
      delete var;
      return -1;
    }
    if (! var->hasrec) var->hasrec = (dimids[i] == recid);
    
    if (chunksizes != 0) var->chunksizes.push_back(chunksizes[i]);
  }
  
  if (isHDF()) {
    var->chunkstorage = storage;
    var->chunksize = chunksize;
    var->chunknelems = chunknelems;
    var->chunkpreemption = chunkpreemption;
    var->no_fill = no_fill;

    if (fill_value != 0) {
      switch(type) {
        case NC_CHAR: var->cfill = *((char*)fill_value); break;
        case NC_BYTE: var->ucfill = *((unsigned char*)fill_value); break;
        case NC_SHORT: var->sfill = *((short*)fill_value); break;
        case NC_USHORT: var->usfill = *((unsigned short*)fill_value); break;
        case NC_INT: var->ifill = *((int*)fill_value); break;
        case NC_UINT: var->uifill = *((unsigned int*)fill_value); break;
        case NC_INT64: var->lfill = *((long long*)fill_value); break;
        case NC_UINT64: var->ulfill = *((unsigned long long*)fill_value); break;
        case NC_FLOAT: var->ffill = *((float*)fill_value); break;
        case NC_DOUBLE: var->dfill = *((double*)fill_value); break;
        default: break;
      }
    }
    
    var->shuffle = shuffle;
    var->deflate = deflate;
    var->deflate_level = deflate_level;
    var->checksum = checksum;
    var->endian = endian;
  }
  
  if (create(var) == -1) {
    delete var;
    return -1;
  }
  
  return var->varid;
};
//###########################################################################
int ezNc::enddef() {
  if (ncid < 0) return 0;
  
  int status = nc__enddef(ncid, h_minfree, v_align, v_minfree, r_align);
  handle_status(status);
  
  return status;
};
//###########################################################################
void ezNc::getAllDimIds(std::vector<int> & ids) {
  ids.clear();
  int n = dimsIdVector.size();
  ids.reserve(n);
  
  for(int i=0; i < n; ++i) {
    if (dimsIdVector[i] != 0)
      ids.push_back(i);
  }
};
//###########################################################################
void ezNc::getAllDimNames(std::vector<std::string> & names) {
  std::vector<int> ids;
  getAllDimIds(ids);
  names.clear();
  int i=0, n=ids.size();
  
  for(; i < n; ++i) {
    names.push_back(dimsIdVector[ids[i]]->name);
  }
}
//###########################################################################
void ezNc::getAllGlobalAttNames(std::vector<std::string> & names) {
  int ngatts;
  int status = nc_inq_natts(ncid, &ngatts);
  if (handle_status(status)) return;
  
  char name[256];
  for(int i=0; i < ngatts; ++i) {
    status = nc_inq_attname(ncid, NC_GLOBAL, i, name);
    handle_status(status);
    names.push_back(name);
  }
}
//###########################################################################
size_t ezNc::getAllVarBytes() {
  return getStaticVarBytes() + getRecVarBytes();
};
//###########################################################################
void ezNc::getAllVarIds(std::vector<int> & ids) {
  ids.clear();
  int n = varsIdVector.size();
  ids.reserve(n);
  
  for(int i=0; i < n; ++i) {
    if (varsIdVector[i] != 0)
      ids.push_back(i);
  }
};
//###########################################################################
void ezNc::getAllVarNames(std::vector<std::string> & names) {
  std::vector<int> ids;
  getAllVarIds(ids);
  names.clear();
  int i=0, n=ids.size();
  
  for(; i < n; ++i) {
    names.push_back(varsIdVector[ids[i]]->name);
  }
};
//###########################################################################
int ezNc::getAtt(const char* varname, const char* attname, ezNcBuffers& values) {
  nc_type type;
  int varid = NC_GLOBAL;
  size_t len;
  
  if (varname) {
    ezNcVar* v = varsNameMap[varname];
    if (!v) return 1;
    varid = v->varid;
  }
  
  int status = nc_inq_att(ncid, varid, attname, &type, &len);
  if (handle_status(status)) return status;
  values.reset();
  values.setSize(type, len);
  values.allocate();
  
  switch(type) {
    case NC_CHAR: status = nc_get_att_text(ncid, varid, attname, values.c); break; 
    case NC_BYTE: case NC_UBYTE: status = nc_get_att_uchar(ncid, varid, attname, values.uc); break;
    case NC_SHORT: status = nc_get_att_short(ncid, varid, attname, values.s); break;
    case NC_USHORT: status = nc_get_att_ushort(ncid, varid, attname, values.us); break;
    case NC_INT: status = nc_get_att_int(ncid, varid, attname, values.i); break;
    case NC_UINT: status = nc_get_att_uint(ncid, varid, attname, values.ui); break;
    case NC_INT64: status = nc_get_att_longlong(ncid, varid, attname, values.l); break;
    case NC_UINT64: status = nc_get_att_ulonglong(ncid, varid, attname, values.ul); break;
    case NC_FLOAT: status = nc_get_att_float(ncid, varid, attname, values.f); break;
    case NC_DOUBLE: status = nc_get_att_double(ncid, varid, attname, values.d); break;
    default: return 0;
  }    
  
  return status;
}
//###########################################################################
int ezNc::getAttString(const char* varname, const char* attname, std::string& string) {
  nc_type type;
  int varid = NC_GLOBAL;
  size_t len;
  string.clear();
  
  if (varname) {
    ezNcVar* v = varsNameMap[varname];
    if (!v) return 1;
    varid = v->varid;
  }

  int status = nc_inq_att(ncid, varid, attname, &type, &len);
  if (handle_status(status)) return status;
  
  if (type == NC_CHAR) {
    char *tmp = new char[len];
    status = nc_get_att_text(ncid, varid, attname, tmp);
    if (handle_status(status)) return status;
    string.assign(tmp, len);
    delete [] tmp;
  }
  
  return 0;
}
//###########################################################################
void ezNc::getCoordVars(int varid, std::vector<int>& varids) {
  ezNcVar* var = getVar(varid);
  if (!var) return;
  
  int n = var->getNumDims();
  int i = 0;
  ezNcDim *dim;
  ezNcVar *dimvar;
  varids.clear();
  varids.reserve(n);
  
  for(; i < n; ++i) {
    dim = getDim(var->dimids[i]);
    if (dim) {
      dimvar = getVar(dim->name);
      if (dimvar)
        varids.push_back(dimvar->varid);
    }
  }
}
//###########################################################################
void ezNc::getCoordVars(std::string& varname, std::vector<std::string>& varnames) {
  ezNcVar* var = getVar(varname);
  if (!var) return;

  int n = var->getNumDims();
  int i = 0;
  ezNcDim *dim;
  ezNcVar *dimvar;
  varnames.clear();
  varnames.reserve(n);

  for(; i < n; ++i) {
    dim = getDim(var->dimids[i]);
    if (dim) {
      varnames.push_back(dim->name);
    }
  }
}
//###########################################################################
void ezNc::getDimIds(std::set<int> & varids, std::set<int> & dimids) {
  std::set<int>::iterator it = varids.begin(), end = varids.end();
  ezNcVar *var;
  for(; it != end; ++it) {
    var = getVar(*it);
    if (!var) continue;
    std::copy(var->dimids.begin(), var->dimids.end(), std::inserter(dimids, dimids.end()));  
  }
}
//###########################################################################
void ezNc::getDimIds(std::vector<int> & varids, std::set<int> & dimids) {
  int i=0, n=varids.size();
  ezNcVar *var;
  for(; i < n; ++i) {
    var = getVar(varids[i]);
    if (!var) continue;
    std::copy(var->dimids.begin(), var->dimids.end(), std::inserter(dimids, dimids.end()));  
  }
}
//###########################################################################
void ezNc::getDimIds(std::vector<std::string> & dimnames, std::vector<int> & dimids) {
  int i=0, n=dimnames.size();
  ezNcDim * dim;
  dimids.resize(n);
  
  for(; i < n; ++i) {
    dim = dimsNameMap[dimnames[i]];
    if (dim)
      dimids[i] = dim->dimid;
    else
      dimids[i] = -1;
  }
};
//###########################################################################
void ezNc::getDimNames(int varid, std::vector<std::string> & dims) {
  ezNcVar * var = getVar(varid);
  if (var == 0) return;
  dims.clear();
  ezNcDim * dim;
  int n = var->getNumDims();
  for(int i=0; i < n; ++i) {
    dim = getDim(var->dimids[i]);
    if (dim)
      dims.push_back(dim->name);
  }
}
//###########################################################################
void ezNc::getDimNamesForVars(std::vector<std::string> & vars, std::set<std::string> & dims) {
  int n = vars.size();
  int i=0;
  int ndims,j;
  
  for(; i < n; ++i) {
    ezNcVar *v = varsNameMap[vars[i]];
    if (v==0) continue;
    
    ndims = v->getNumDims();
    for(j=0; j < ndims; ++j) {
      ezNcDim *d = getDim(v->dimids[j]);
      if (d==0) continue;
      dims.insert(d->name);
    }
  }
};
//###########################################################################
void ezNc::getDimSizes(int varid, std::vector<size_t>& dimsizes) {
  ezNcVar * var = getVar(varid);
  if (var == 0) return;
  dimsizes.clear();
  ezNcDim * dim;
  int n = var->dimids.size();
  dimsizes.resize(n);
  for(int i=0; i < n; ++i) {
    dim = getDim(var->dimids[i]);
    if (dim) {
      dimsizes[i] = dim->size;
    } else {
      dimsizes[i] = 0;
    }
  }
};
//###########################################################################
void ezNc::getDimSizes(std::string & name, std::vector<size_t>& dimsizes) {
  if (varsNameMap.count(name))
    getDimSizes(varsNameMap[name]->varid, dimsizes);
};
//###########################################################################
int ezNc::getMaxNumDims(std::vector<std::string> & names) {
  int maxndims=0;
  int nvars = names.size();
  int size;
  
  ezNcVar * var;
  for(int i=0; i < nvars; ++i) {
    var = varsNameMap[names[i]];
    if (var==0) continue;
    size = var->getNumDims();
    if (size > maxndims)
      maxndims = size;
  }
  
  return maxndims;
};
//###########################################################################
int ezNc::getNumStaticVars() {
  int n = varsIdVector.size();
  int ctr = 0;
  ezNcVar* var;
  
  for(int i=0; i < n; ++i) {
    var = getVar(i);
    if (var && !var->hasrec)
      ++ctr;
  }
  
  return ctr;
};
//###########################################################################
int ezNc::getNumRecVars() {
  int n = varsIdVector.size();
  int ctr = 0;
  ezNcVar* var;
  
  for(int i=0; i < n; ++i) {
    var = getVar(i);
    if (var && var->hasrec)
      ++ctr;
  }
  
  return ctr;
};
//###########################################################################
size_t ezNc::getRecVarBytes() {
  size_t total = 0;
  int n = varsIdVector.size();
  int ctr = 0;
  ezNcVar* var;
  
  for(int i=0; i < n; ++i) {
    var = getVar(i);
    if (var && var->hasrec)
      total += var->bytes;
  }

  return total;
}; 
//###########################################################################
void ezNc::getRecVarNames(std::vector<std::string> & names) {
  std::vector<int> ids;
  getRecVarIds(ids);
  names.clear();
  
  int i=0, n=ids.size();
  
  for(; i < n; ++i) {
    names.push_back(varsIdVector[ids[i]]->name);
  }
};
//###########################################################################
void ezNc::getRecVarIds(std::vector<int> & varids) {
  varids.clear();
  int n = varsIdVector.size();
  int i=0;
  ezNcVar* var;

  for(; i < n; ++i) {
    var = getVar(i);
    if (var && var->hasrec)
      varids.push_back(var->varid);
  }
  
  std::sort(varids.begin(), varids.end());
};
//###########################################################################
size_t ezNc::getStaticVarBytes() {
  size_t total = 0;
  int n = varsIdVector.size();
  int ctr = 0;
  ezNcVar* var;
  
  for(int i=0; i < n; ++i) {
    var = getVar(i);
    if (var && !var->hasrec)
      total += var->bytes;
  }

  return total;
};
//#######################################################################
inline ezNcVar* ezNc::getVar(unsigned int id) {
  if (id >= varsIdVector.size()) return 0; 
  else return varsIdVector[id]; 
}
//#######################################################################
inline ezNcVar* ezNc::getVar(const char* name) {
  if (name==0) return 0; 

  return varsNameMap[name]; 
}
//#######################################################################
inline ezNcVar* ezNc::getVar(std::string& name) {
  return varsNameMap[name];
}
//#######################################################################
// Get var id from name. -1 if invalid.
inline int ezNc::getVarId(const char* name) {
  if (name==0) return -1;
  ezNcVar* v = varsNameMap[name];
  if (v)
    return v->varid;
  else
    return -1; 
}
//#######################################################################
// Get var id from name. -1 if invalid.
inline int ezNc::getVarId(const std::string& name) {
  if (name.empty()) return -1;
  ezNcVar* v = varsNameMap[name];
  if (v)
    return v->varid;
  else
    return -1; 
}
//#######################################################################
void ezNc::getStaticVarNames(std::vector<std::string> & names) {
  std::vector<int> ids;
  getStaticVarIds(ids);
  names.clear();
  
  int i=0, n=ids.size();
  
  for(; i < n; ++i) {
    names.push_back(varsIdVector[ids[i]]->name);
  }
};
//###########################################################################
void ezNc::getStaticVarIds(std::vector<int> & varids) {
  varids.clear();
  int n = varsIdVector.size();
  int i=0;
  ezNcVar* var;

  for(; i < n; ++i) {
    var = getVar(i);
    if (var && !var->hasrec)
      varids.push_back(var->varid);
  }
  
  std::sort(varids.begin(), varids.end());
};
//###########################################################################
void ezNc::getVarSizeTZYX(int id, size_t & nt, size_t & nz, size_t & ny, size_t & nx) {
  nt = 0;
  nz = 0;
  ny = 0;
  nx = 0;
  ezNcVar* var = getVar(id);
  if (var==0)
    return;
  
  int ndims = var->getNumDims();
  ezNcDim* dim;
  
  if (var->hasrec) {
    dim = getDim(var->dimids[0]);
    if (dim) nt = dimsIdVector[dim->dimid]->size;
    
    switch(ndims) {
      case 4:
        dim = getDim(var->dimids[1]);
        if (dim) nz = dimsIdVector[dim->dimid]->size;
        dim = getDim(var->dimids[2]);
        if (dim) ny = dimsIdVector[dim->dimid]->size;
        dim = getDim(var->dimids[3]);
        if (dim) nx = dimsIdVector[dim->dimid]->size;
        break;
      case 3:
        dim = getDim(var->dimids[1]);
        if (dim) ny = dimsIdVector[dim->dimid]->size;
        dim = getDim(var->dimids[2]);
        if (dim) nx = dimsIdVector[dim->dimid]->size;
        break;
      case 2:
        dim = getDim(var->dimids[1]);
        if (dim) nx = dimsIdVector[dim->dimid]->size;
        break;
      }
  } else {
  switch(ndims) {
    case 3:
      dim = getDim(var->dimids[0]);
      if (dim) nz = dimsIdVector[dim->dimid]->size;
      dim = getDim(var->dimids[1]);
      if (dim) ny = dimsIdVector[dim->dimid]->size;
      dim = getDim(var->dimids[2]);
      if (dim) nx = dimsIdVector[dim->dimid]->size;
      break;
    case 2:
      dim = getDim(var->dimids[0]);
      if (dim) ny = dimsIdVector[dim->dimid]->size;
      dim = getDim(var->dimids[1]);
      if (dim) nx = dimsIdVector[dim->dimid]->size;
      break;
    case 1:
      dim = getDim(var->dimids[0]);
      if (dim) nx = dimsIdVector[dim->dimid]->size;
      break;
    }
  }
};
//#########################################################
void ezNc::getVars(std::vector<ezNcVar*>& out, std::vector<std::string>& patterns, bool regex, bool regex_nocase, bool coordvars, bool exclude) {
  std::vector<std::string> allNames;
  getAllVarNames(allNames);
  int j,m,i;
  int n = patterns.size();
  // Unsorted list of regex matches of varnames.
  std::vector<std::string> matches;
  std::set<int> matchIds;

  if (regex) {
    for(i=0; i < n; ++i) {
      matches.clear();
      ez::regex(patterns[i], allNames, matches, 1/* Extended regex */, regex_nocase);
      
      m = matches.size();
      for(j=0; j < m; ++j)
        matchIds.insert(getVarId(matches[j]));
    }
  } else {
    for(i=0; i < n; ++i) {
      if (ez::find_first(allNames, patterns[i]) >= 0) {
        matchIds.insert(getVarId(patterns[i]));
      }
    }
  }
  
  if (exclude) {
    // Remove from entire list all those that were matched.
    std::vector<int> allIds;
    std::set<int> notExcludedIds;
    getAllVarIds(allIds);
    std::set_difference(allIds.begin(), allIds.end(), matchIds.begin(), matchIds.end(), std::inserter(notExcludedIds, notExcludedIds.begin()));
    matchIds = notExcludedIds;
  }
  
  //------------------------------------------------------
  // Do coordvar lookups to get set of vars to also include or exclude.
  std::vector<int> allCoordVarIds;
  std::set<int> coordVarIds;
  std::set<int>::iterator it = matchIds.begin(), end = matchIds.end();
  for(; it != end; ++it) {

    getCoordVars(*it, allCoordVarIds);
    n = allCoordVarIds.size();
    // Check that coordvar was not in pattern list, otherwise we would accidentally remove it if coordvars are not desired.
    for(i=0; i < n; ++i) {

    // If coordvar is in match list, then skip it.
      if (matchIds.count( allCoordVarIds[i]) ) continue;
      // Otherwise, add to unique set of coordvar ids that will later be kept or used for exclusion.
      coordVarIds.insert( allCoordVarIds[i] );
    }
  }
  
  it = coordVarIds.begin();
  end = coordVarIds.end();

  if (coordvars) {
    // Include coord vars.
    for(; it != end; ++it) {
      matchIds.insert(*it);
    }
  } else {
    // Exclude the coord varids we found.
    for(; it != end; ++it) {
      matchIds.erase(*it);
    }
  }

  // Prepare for output.
  it = matchIds.begin();
  end = matchIds.end();
  for(; it != end; ++it) 
    out.push_back(getVar(*it));
}
//###########################################################################
void ezNc::getVarSizeTZYX(std::string& name, size_t & nt, size_t & nz, size_t & ny, size_t & nx) {
  ezNcVar* var = varsNameMap[name];
  if(var==0) return;
  
  getVarSizeTZYX(var->varid, nt, nz, ny, nx);
};
//###########################################################################
int ezNc::handle_status(int status) {
  if (status != NC_NOERR) {
    fprintf(stderr, "%s\n", nc_strerror(status));
  }
  return status;
};
//###########################################################################
bool ezNc::isHDF() {
  switch(format) {
    case NC_FORMAT_NETCDF4: 
    case NC_FORMAT_NETCDF4_CLASSIC:
      return true;
    default: 
      return false;
  }
};
//###########################################################################
int ezNc::loadAllMetadata() {
  int status = loadDimsMetadata();
  if (handle_status(status)) return status;

  status = loadVarsMetadata();
  handle_status(status);

  return status;
};
//###########################################################################
int ezNc::loadDimsMetadata() {
  int status, ndims;
  
  status = nc_inq_unlimdim(ncid, &recid);
  if (handle_status(status)) return status;

  if (recid > 0) {
    status = nc_inq_dimlen(ncid, recid, &nrec);
    if (handle_status(status)) return status;
  } else {
    nrec = 0;
  }
  
  status = nc_inq_ndims(ncid, &ndims);
  if (handle_status(status)) return false;
  
  char name[256];
  size_t len;
  dimsIdVector.resize(ndims);
  
  for(int i=0; i < ndims; ++i) {
    status = nc_inq_dim(ncid, i, name, &len);    
    if (handle_status(status)) return status;

    ezNcDim *dim = new ezNcDim;
    dim->nc = this;
    dim->dimid = i;
    dim->name = name;
    dim->isrec = (i == recid);
    dim->size = len;      
    dimsNameMap[dim->name] = dim;
    dimsIdVector[i] = dim;
  }
  
  return 0;
};
//###########################################################################
int ezNc::loadVarsMetadata() {
  int status, nvars, ndims;
  size_t len, dimlen;
  status = nc_inq_nvars(ncid, &nvars);
  if (handle_status(status)) return status;
  
  char name[256];
  nc_type type;
  int nvardims, dimids[10], natts, j, bytesize;
  ezNcVar *var;
  varsIdVector.resize(nvars);
  
  for(int i=0; i < nvars; ++i) {
    status = nc_inq_var(ncid, i, name, &type, &nvardims, dimids, &natts);
    if (handle_status(status)) return status;

    var = new ezNcVar;
    var->nc = this;
    var->varid = i;
    var->name = name;
    var->type = type;
    varsNameMap[var->name] = var;
    varsIdVector[i] = var;
    
    for(j=0; j < nvardims; ++j) {
      var->dimids.push_back(dimids[j]);
      if (dimids[j] == recid)
        var->hasrec = true;
    }

    computeSize(var);

    for(j=0; j < natts; ++j) {
      status = nc_inq_attname(ncid, i, j, name);
      if (handle_status(status)) return status;

      var->attnames.push_back(name);
    }
  }

  if (isHDF()) {
    for(int i=0; i < nvars; ++i) {
      var = varsIdVector[i];
      var->chunksizes.resize(var->getNumDims());
      status = nc_inq_var_chunking(ncid, i, &var->chunkstorage, &(var->chunksizes[0]));
      if (handle_status(status)) return status;

      status = nc_get_var_chunk_cache(ncid, i, &var->chunksize, &var->chunknelems, &var->chunkpreemption);
      if (handle_status(status)) return status;
      
      status = nc_inq_var_fill(ncid, i, &var->no_fill, (void*)(&var->ifill) );
      if (handle_status(status)) return status;

      status = nc_inq_var_deflate(ncid, i, &var->shuffle, &var->deflate, &var->deflate_level);
      if (handle_status(status)) return status;
      
      status = nc_inq_var_fletcher32(ncid, i, &var->checksum);
      if (handle_status(status)) return status;
      
      status = nc_inq_var_endian(ncid, i, &var->endian);
      if (handle_status(status)) return status;
    }
  }
  
  return 0;  
};
//###########################################################################
int ezNc::openReadOnly() {
  int status = nc__open(filename.c_str(), NC_NOWRITE, &bufrsizehint, &ncid);
  if (handle_status(status)) return status;
  
  status = nc_inq_unlimdim(ncid, &recid);
  if (handle_status(status)) return status;

  status = nc_inq_format(ncid, &format);
  if (handle_status(status)) return status;
  
  return status;
};
//###########################################################################
int ezNc::openWritable() {
  int status = nc__open(filename.c_str(), NC_WRITE, &bufrsizehint, &ncid);
  if (handle_status(status)) return status;
  
  status = nc_inq_unlimdim(ncid, &recid);
  if (handle_status(status)) return status;

  status = nc_inq_format(ncid, &format);
  if (handle_status(status)) return status;
  
  status = nc_redef(ncid);
  if (handle_status(status)) return status;

  return status;
};
//###########################################################################
#ifdef DEBUG
void ezNc::print(ezNc & nc) {
  std::stringstream ss;
  ss << "filename=" << nc.filename << "\n";
  ss << "ncid=" << nc.ncid << "\n";
  ss << "format=" << nc.format << "\n";
  ss << "recid=" << nc.recid << "\n";
  ss << "initialsz=" << nc.initialsz << "\n";
  ss << "bufrsizehint=" << nc.bufrsizehint << "\n";
  ss << "h_minfree=" << nc.h_minfree << "\n";
  ss << "v_align=" << nc.v_align << "\n";
  ss << "v_minfree=" << nc.v_minfree << "\n";
  ss << "r_align=" << nc.r_align << "\n";
  std::cout << ss;
  
  int i;
  for(i=0; i < dims.size(); ++i)
    print(dims[i]);
    
  for(i=0; i < vars.size(); ++i)
    print(vars[i]);
};
#endif
//###########################################################################
int ezNc::putAtt(const char* varname, const char* attname, nc_type type, std::vector<std::string> & values, int first, int last) {
  ezNcBuffers buf;
  int n = last-first+1;
  if (n<1) return 0;
  
  if (type != NC_CHAR) {
    buf.setSize(type, n);
    buf.allocate();
  }
  
  int status;
      
  #define STR2BUF(P,T)  {\
    StringsToValues<T>(values, first, last, buf.P);\
  }
  
  switch(type) {
    case NC_BYTE: case NC_UBYTE: STR2BUF(uc,unsigned char); break;
    case NC_CHAR: 
      buf.setSize( type, values[first].size() );
      buf.allocate();
      strncpy(buf.c, values[first].c_str(), buf.nc);
      break;
    case NC_SHORT: STR2BUF(s,short); break;
    case NC_USHORT: STR2BUF(us,unsigned short); break;
    case NC_INT: STR2BUF(i,int); break;
    case NC_UINT: STR2BUF(ui,unsigned int); break;
    case NC_INT64: STR2BUF(l,long long); break;
    case NC_UINT64: STR2BUF(ul,unsigned long long); break;
    case NC_FLOAT: STR2BUF(f,float); break;
    case NC_DOUBLE: STR2BUF(d,double); break;
    default: break;
  }
  return putAtt(varname, attname, type, buf);
};
//###########################################################################
int ezNc::putAtt(const char* varname, const char* attname, nc_type type, ezNcBuffers& values) {
  int varid = NC_GLOBAL; // Default to global att if no var name.
  ezNcVar* var;
  int status;
  
  if ( (varname != 0) && (var = varsNameMap[varname]) )
    varid = var->varid;

  #define PUTBUFATT(T,P)  {\
    status = nc_put_att_##T(ncid, varid, attname, type, values.n##P, values.P);\
    handle_status(status); return status; \
  }
  
  switch(type) {
    case NC_BYTE: case NC_UBYTE: PUTBUFATT(uchar,uc); break;
    case NC_CHAR: 
      status = nc_put_att_text(ncid, varid, attname, values.nc, values.c);
      handle_status(status);
      return status;
    case NC_SHORT: PUTBUFATT(short,s); break;
    case NC_USHORT: PUTBUFATT(ushort,us); break;
    case NC_INT: PUTBUFATT(int,i); break;
    case NC_UINT: PUTBUFATT(uint,ui); break;
    case NC_INT64: PUTBUFATT(longlong,l); break;
    case NC_UINT64: PUTBUFATT(ulonglong,ul); break;
    case NC_FLOAT: PUTBUFATT(float,f); break;
    case NC_DOUBLE: PUTBUFATT(double,d); break;
    default: break;
  }

  return 0;
}
//###########################################################################
int ezNc::read(const char* name, ezNcBuffers * buf, size_t * start, size_t * count) {
  if (!name) return 1;
  ezNcVar * var = varsNameMap[name];
  if (!var) return 1;
  
  return read(var->varid, buf, start, count);
};
//###########################################################################
int ezNc::read(int varid, ezNcBuffers * buf, size_t * start, size_t * count) {
  #define GET(N,P) {\
    status = nc_get_vara_##N(ncid, varid, start, count, buf->P);\
    if(handle_status(status)) return status;\
  }
  if ( (varid < 0) || (varid >= varsIdVector.size()) ) 
    return 1;
    
  int status;
  switch(varsIdVector[varid]->type) {
    case NC_BYTE: case NC_UBYTE: GET(uchar,uc); break;
    case NC_CHAR: GET(text,c); break;
    case NC_SHORT: GET(short,s); break;
    case NC_USHORT: GET(ushort,us); break;
    case NC_INT: GET(int,i); break;
    case NC_UINT: GET(uint,ui); break;
    case NC_INT64: GET(longlong,l); break;
    case NC_UINT64: GET(ulonglong,ul); break;
    case NC_FLOAT: GET(float,f); break;
    case NC_DOUBLE: GET(double,d); break;
    default: return 1;
  }
  
  return 0;
};
//###########################################################################
int ezNc::read(int varid, ezNcBuffers * buf) {
  #define GETENTIRE(N,P) {\
    status = nc_get_var_##N(ncid, varid, buf->P);\
    if(handle_status(status)) return status;\
  }
  if ( (varid < 0) || (varid >= varsIdVector.size()) ) 
    return 1;
    
  int status;
  switch(varsIdVector[varid]->type) {
    case NC_BYTE: case NC_UBYTE: GETENTIRE(uchar,uc); break;
    case NC_CHAR: GETENTIRE(text,c); break;
    case NC_SHORT: GETENTIRE(short,s); break;
    case NC_USHORT: GETENTIRE(ushort,us); break;
    case NC_INT: GETENTIRE(int,i); break;
    case NC_UINT: GETENTIRE(uint,ui); break;
    case NC_INT64: GETENTIRE(longlong,l); break;
    case NC_UINT64: GETENTIRE(ulonglong,ul); break;
    case NC_FLOAT: GETENTIRE(float,f); break;
    case NC_DOUBLE: GETENTIRE(double,d); break;
    default: return 1;
  }
  
  return 0;
};
//###########################################################################
int ezNc::reserve(ezNcBuffers & buffers) {
  std::vector<int> varids;
  getAllVarIds(varids);
  return reserve(buffers, varids);
}
//###########################################################################
int ezNc::reserve(ezNcBuffers & buffers, const char* varname) {
  if (!varname) return 1;
  ezNcVar* var = varsNameMap[varname];
  if (!var) return 1;
  return reserve(buffers, var->varid);
}
//###########################################################################
int ezNc::reserve(ezNcBuffers & buffers, int varid) {
  std::vector<size_t> dimsizes;
  ezNcVar* var;
  
  var = getVar(varid);
  if (!var) return 1;
  getDimSizes(varid, dimsizes);
  buffers.setMaxSizePerRec<size_t>(var->type, dimsizes, var->dimids, recid);
    
  return 0;
};
//###########################################################################
int ezNc::reserveEntire(ezNcBuffers & buffers, int varid) {
  std::vector<size_t> dimsizes;
  ezNcVar* var;
  
  var = getVar(varid);
  if (!var) return 1;
  getDimSizes(varid, dimsizes);

  // Compute size of entire variable across all dims.
  int n = dimsizes.size();
  size_t size = n ? 1 : 0;
  for(int i=0; i < n; ++i) { 
    size *= dimsizes[i];
  }
  buffers.reserve(var->type, size);
    
  return 0;
};
//###########################################################################
int ezNc::reserve(ezNcBuffers & buffers, std::vector<int> & varids) {
  int i=0;
  std::vector<size_t> dimsizes;
  int n = varids.size();
  ezNcVar* var;
  int status = 0;
  
  for(; i < n; ++i) {
    var = getVar(varids[i]);
    if (!var) { status = 1; continue; }
    
    getDimSizes(varids[i], dimsizes);
    buffers.setMaxSizePerRec(var->type, dimsizes, var->dimids, recid);
  }
  return status;
};
//###########################################################################
int ezNc::reserve(ezNcBuffers & buffers, std::vector<std::string> & names) {
  std::vector<int> varids;
  int n = names.size();
  int i = 0;
  ezNcVar* var;
  
  for(; i < n; ++i)
    if (var = varsNameMap[names[i]])
      varids.push_back( var->varid );
    
  return reserve(buffers, varids);
};
//###########################################################################
void ezNc::typeToString(nc_type t, std::string& s) {
  s.clear();
  switch(t) {
    case NC_CHAR: s = "char"; break; 
    case NC_BYTE: case NC_UBYTE: s = "unsigned char"; break;
    case NC_SHORT: s = "short"; break;
    case NC_USHORT: s = "unsigned short"; break;
    case NC_INT: s = "int"; break;
    case NC_UINT: s = "unsigned int"; break;
    case NC_INT64: s = "long long"; break;
    case NC_UINT64: s = "unsigned long long"; break;
    case NC_FLOAT: s = "float"; break;
    case NC_DOUBLE: s = "double"; break;
    default: break;
  }
};
//###########################################################################
int ezNc::write(const char* name, ezNcBuffers * buf, size_t * start, size_t * count) {
  ezNcVar * var = varsNameMap[name];
  if (!var) return 1;
  return write(var->varid, buf, start, count);
};  
//###########################################################################
int ezNc::write(int varid, ezNcBuffers * buf, size_t * start, size_t * count) {
  #define PUT(N,P) {\
    status = nc_put_vara_##N(ncid, varid, start, count, buf->P);\
    if(handle_status(status)) return status;\
  }
  
  int status;
  switch(varsIdVector[varid]->type) {
    case NC_BYTE: case NC_UBYTE: PUT(uchar,uc); break;
    case NC_CHAR: PUT(text,c); break;
    case NC_SHORT: PUT(short,s); break;
    case NC_USHORT: PUT(ushort,us); break;
    case NC_INT: PUT(int,i); break;
    case NC_UINT: PUT(uint,ui); break;
    case NC_INT64: PUT(longlong,l); break;
    case NC_UINT64: PUT(ulonglong,ul); break;
    case NC_FLOAT: PUT(float,f); break;
    case NC_DOUBLE: PUT(double,d); break;
    default: return 1;
  }
  
  return 0;
};

}
#endif // EZNC_H
