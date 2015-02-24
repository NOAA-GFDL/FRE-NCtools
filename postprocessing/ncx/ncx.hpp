/*
ncx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ncx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ncx.  If not, see <http://www.gnu.org/licenses/>.

Copyright 2011,2012 Remik Ziemlinski

20111029 rsz Created.
20111212 rsz Export coordvars by default.
20111212 rsz Fix record coordvar data export. Add ncx version global att string.
20111219 rsz Fix 0D recvar data export.
*/
#ifndef NCX_H
#define NCX_H

#define NCX_VERSION "1.1.1"

#include <ezOptionParser.hpp>
#include <ezProgressBar.hpp>
#include <ezRateProgressBar.hpp>
#include <ezNc.hpp>
#include <ezNcSlicer.hpp>
#include <ezPointerVector.hpp>
#include <ezNcCopySlabTask.hpp>
#include <ezNcSlabIterator.hpp>
#include <time.h> // To make date-time history global attribute.

class NcX {
public:
  // Ctor.
  NcX() : debug(0), mem_virt(0), mem_rss(0), verbose(1) {};
  
  // See if we hit peak memory and store results.
  void GetPeakMemory();
  // Process input options to designate input source vars to export.
  int ChooseSrcVars();
  // Closes all files. Returns 0 on success.
  int CloseFiles();
  // Create all data transfer tasks between files. Returns 0 on success.
  int CreateTasks();
  // Make history string for global attribute in output file.
  void MakeHistory();
  // Returns 0 if successfully opens input and output files.
  int OpenFiles();
  // Sets up the metadata in the output file. Returns 0 on success.
  int ProcessMetadata();
  // Print summary of peak memory usage.
  void ReportMemory();
  // Call this after SetOptions. Does all the real work.
  int Run();
  // Runs the i/o tasks between files. Returns 0 on success.
  int RunTasks();
  // Sets up filenames and optionally for append mode, and if not append and a file exists, then ensure overwrite is enabled. Returns 0 on success, 1 on failure.
  int SetFileNames();
  // Sets options from command-line input. Returns 0 if successful.
  int SetOptions(int argc, const char *argv[]);
  // Do slicing of dims requested by user to determine output shapes.
  void SliceSrcDims();
  // Prints usage message to stdout.
  void Usage();
  // Get version message.
  void Version(std::string&);
  
  // Command line args. Keep so we can make "history" global attribute of the command that was invoked.
  int argc;
  const char **argv;
  // Arena memory (pool) for data transfers between files.
  ez::ezNcBuffers buffers;
  // Print detailed messages to debug.
  bool debug;
  // Rich progress indicator.
  ez::ezRateProgressBar<float> eta;
  // Peak memory used. Reported with -m option.
  double mem_virt, mem_rss;
  // Command-line options.
  ez::ezOptionParser opt;
  // Simple progress indicator.
  ez::ezProgressBar progress;
  // File with data that will be appended into, but is really just another input source.
  ez::ezNc ncAppendSrc;
  // One and only output file sink.
  ez::ezNc ncDst;
  // Required input file.
  ez::ezNc ncSrc;
  // Does all the hard slicing logic for input source file.
  ez::ezNcSlicer<int> ncSrcSlicer;
  // Final (sorted) set of dimids of ncSrc to export.
  std::set<int> srcDimIdsToExport;
  // Final (sorted) set of varids of ncSrc to export.
  std::set<int> srcVarIdsToExport;
  // Data read/write tasks in execution order.
  ez::ezPointerVector< ez::ezNcCopySlabTask<size_t> > tasks;
  // If verbose info or warning messages. All errors will always print.
  bool verbose;
};
//#####################################################################
void NcX::GetPeakMemory() {
  if (!opt.isSet("-m")) return;
  
  double vm, rss;
  ez::process_mem_usage(vm, rss);
  if (vm > mem_virt) mem_virt = vm;
  if (rss > mem_rss) mem_rss = rss;
}
//#####################################################################
int NcX::ChooseSrcVars() {
  if(debug) {
    printf("%s:%d ChooseSrcVars\n", __FILE__, __LINE__); fflush(stdout);
  }

  // Find which dims and vars we need to copy from input file. Dims/vars from append source are dealt with in CreateTasks since they're mandatory if doing append.
  std::vector<std::string> allSrcVars;
  ncSrc.getAllVarNames(allSrcVars);
  
  // Command line var list (which could be a regex pattern list).
  std::vector<std::string> varPatterns;
  opt.get("-v")->getStrings(varPatterns);
  
  // Unsorted list of regex matches of varnames.
  std::vector<std::string> varsToExport;
  
  bool regex = opt.isSet("-R");
  int i, j;
  
  int n = varPatterns.size();
  if (regex) {
    for(i=0; i < n; ++i) {
      ez::regex(varPatterns[i], allSrcVars, varsToExport, 1/* Extended regex? */, 0 /* case insensitive? */);
    }
  } else {
    for(i=0; i < n; ++i) {
      j = ez::find_first(allSrcVars, varPatterns[i]);
      if (j >= 0)
        varsToExport.push_back(varPatterns[i]);
      else {
        // Varname not found. Abort like "ncks".
        fprintf(stderr, "ERROR: User specified variable \"%s\" is not in input file \"%s\". Aborting.\n", varPatterns[i].c_str(), ncSrc.filename.c_str());
        return 1;
      }
    }
  }
    
  std::vector<std::string> tmpStrings;
  std::set<int>::iterator it, end;
  
  if (opt.isSet("-x")) {
    // Exclude var names found.
    ez::not_in_second(allSrcVars, varsToExport, tmpStrings);
    varsToExport = tmpStrings;
    tmpStrings.clear();
  }
  
  n = varsToExport.size();
  // Convert names to varids for faster lookups.
  int varid;
  for(i=0; i < n; ++i) {
    varid = ncSrc.getVarId(varsToExport[i]);
    if (varid != -1) srcVarIdsToExport.insert(varid);
  }

  //------------------------------------------------------
  // Do coordvar lookups to get set of vars to also include or exclude.
  std::set<int> coordSet;
  std::vector<int> coordVector;
  it = srcVarIdsToExport.begin();
  end = srcVarIdsToExport.end();
  for(; it != end; ++it) {
    ncSrc.getCoordVars(*it, coordVector);
    // Check that a found coordvar was not in user's -v list, otherwise we would accidentally remove it when -C option is used.
    n = coordVector.size();
    for(i=0; i < n; ++i) {
      // If coordvar is in user's list, then skip it.
      if (srcVarIdsToExport.count( coordVector[i]) ) continue;
      // Grow unique set of coordvar ids that will then be included or excluded.
      coordSet.insert( coordVector[i] );
    }
  }

  it = coordSet.begin();
  end = coordSet.end();

  if (opt.isSet("-C")) { 
    // Exclude from srcVarIdsToExport the coord varids we found.
    for(; it != end; ++it) {
      srcVarIdsToExport.erase(*it);
    }
  } else {
    // Include coord vars.
    for(; it != end; ++it) {
      srcVarIdsToExport.insert(*it);
    }
  }
  
  // Build set of dim ids used by the exportable vars.
  ncSrc.getDimIds(srcVarIdsToExport, srcDimIdsToExport);
  
  if(debug) {
    printf("%s:%d Exporting these %d var names (ids) from input file:", __FILE__, __LINE__, (int)srcVarIdsToExport.size());
    it = srcVarIdsToExport.begin();
    end = srcVarIdsToExport.end();
    for(; it != end; ++it) {
      ez::ezNcVar *pvar = ncSrc.getVar(*it);
      printf(" %s (%d)", (pvar ? pvar->name.c_str() : "null"), *it);
    }
    printf("\n");
    fflush(stdout);
  }
  
  return 0;
}
//#####################################################################
int NcX::CloseFiles() {
  if (debug) printf("%s:%d CloseFiles\n", __FILE__, __LINE__);
  int status = 0;
  
  if (status = ncSrc.close()) return status;
  if (status = ncDst.close()) return status;
  
  bool append = !ncAppendSrc.filename.empty();
  
  if (append) {
    if (status = ncAppendSrc.close()) return status;
    // Remove old file that has name of what we're appending into.
    if (remove(ncAppendSrc.filename.c_str())) {
      fprintf(stderr, "ERROR: Failed to remove old file \"%s\" that is being appended into. Results were saved to \"%s\". Aborting.\n", ncAppendSrc.filename.c_str(), ncDst.filename.c_str());
      return 1;
    }
    // Rename tmp output file to final output name used by append source.
    if (rename(ncDst.filename.c_str(), ncAppendSrc.filename.c_str())) {
      fprintf(stderr, "ERROR: Failed to rename temporary file \"%s\" to appended output file \"%s\".\n", ncDst.filename.c_str(), ncAppendSrc.filename.c_str());
      return 1;	
    }
  }
  
  if (opt.isSet("-p") && verbose) { ++progress; printf("\n"); }
  if (opt.isSet("-e") && verbose) eta += 1.0; // Signal finally done.
  if (opt.isSet("-m")) GetPeakMemory();

  if (debug) printf("%s:%d CloseFiles done\n", __FILE__, __LINE__);
  return status;
}
//#####################################################################
int NcX::CreateTasks() {
  if(debug) printf("%s:%d CreateTasks\n", __FILE__, __LINE__);
  bool append = !ncAppendSrc.filename.empty();
  ez::ezNcSlabIterator<int> slitSrc, slitDst;
  std::set<int>::iterator idit, idend;
  int i=0, n=0;
  int nSrcRecVars=0, nAppendSrcRecVars=0;
  ez::ezNcVar *varSrc=0, *varDst=0;
  ez::ezNcCopySlabTask<size_t> *task=0;
  std::vector<int> vInt;
  std::string tmpstring;
  
  // Guestimate number of tasks before making the tasks.
  if (append) {
    // Vars in append source aren't sliced, so each static var is a single task, and each record per recvar is a task.
    nAppendSrcRecVars = ncAppendSrc.getNumRecVars();
    i = ncAppendSrc.getNumStaticVars() + nAppendSrcRecVars * ncAppendSrc.nrec;     
  }  
  i += ncSrc.getNumStaticVars() + ncSrc.getNumRecVars() * ncSrc.nrec;
  tasks.vector.reserve(i);

  if (append) {
    ncAppendSrc.getStaticVarIds(vInt);
    n = vInt.size();
    if(debug) printf("%s:%d append has %d static vars\n", __FILE__, __LINE__, n);

    for(i=0; i < n; ++i) {
      task = new ez::ezNcCopySlabTask<size_t>;
      tasks.vector.push_back(task);
      varSrc = ncAppendSrc.getVar(vInt[i]);
      if(debug) printf("%s:%d append file static var %s\n", __FILE__, __LINE__, varSrc->name.c_str());

      task->slabSrc.setVar( varSrc );
      varDst = ncDst.getVar(varSrc->name);
      if (!varDst) {
        fprintf(stderr, "%s:%d ERROR: Failed to get output file var \"%s\"\n", __FILE__, __LINE__, varSrc->name.c_str());
        return 1;
      }
      task->slabDst.setVar( varDst );
      
      // Set the slabs' counts to copy entire static var. We can't slice vars in a file we're appending into. Start vectors will be already all zeros.
      ncAppendSrc.getDimSizes(varSrc->name, task->slabSrc.count);
      task->slabDst.count = task->slabSrc.count;
      
      if(debug) {
        tmpstring.clear();
        task->slabSrc.toString(tmpstring);
        printf("%s:%d %s slabIn: %s\n", __FILE__, __LINE__, task->slabSrc.var->name.c_str(), tmpstring.c_str());
        tmpstring.clear();
        task->slabDst.toString(tmpstring);
        printf("%s:%d %s slabOut: %s\n", __FILE__, __LINE__, task->slabDst.var->name.c_str(), tmpstring.c_str());
      }
    }
  }
  
  // Input source static vars.
  idit = srcVarIdsToExport.begin();
  idend = srcVarIdsToExport.end();
  nSrcRecVars = srcVarIdsToExport.size();
  for(; idit != idend; ++idit) {
    i = *idit;
    varSrc = ncSrc.getVar(i);
    if (!varSrc || varSrc->hasrec) {
      --nSrcRecVars;
      continue;
    }
    
    if(debug) printf("%s:%d input file static var %s\n", __FILE__, __LINE__, varSrc->name.c_str());

    varDst = ncDst.getVar(varSrc->name);
    if (!varDst) {
      fprintf(stderr, "%s:%d ERROR: Failed to get output file var \"%s\"\n", __FILE__, __LINE__, varSrc->name.c_str());
      return 1;
    }
    slitSrc.setVar(varSrc, &ncSrcSlicer);
    slitDst = slitSrc;
    slitDst.var = varDst;
    slitDst.normalize();
    while(! slitSrc.isDone() ) {
      task = new ez::ezNcCopySlabTask<size_t>;
      tasks.vector.push_back(task);
      slitSrc.getSlab(&task->slabSrc);
      slitDst.getSlab(&task->slabDst);
      slitSrc.inc();
      slitDst.inc();
      if(debug) {
        tmpstring.clear();
        task->slabSrc.toString(tmpstring);
        printf("%s:%d %s slabIn: %s\n", __FILE__, __LINE__, task->slabSrc.var->name.c_str(), tmpstring.c_str());
        tmpstring.clear();
        task->slabDst.toString(tmpstring);
        printf("%s:%d %s slabOut: %s\n", __FILE__, __LINE__, task->slabDst.var->name.c_str(), tmpstring.c_str());
      }
    }
  }

  slitSrc.clear();
  slitDst.clear();
  
  // All record vars interleaved, append source first, then input source.
  // Create iterators for each rec var, then create tasks using iterators.
  ez::ezPointerVector< ez::ezNcSlabIterator<int> > srcIters,dstIters;
  ez::ezNcSlabIterator<int> *pSrcSlit, *pDstSlit;
  srcIters.vector.reserve(nAppendSrcRecVars + nSrcRecVars);
  dstIters.vector.reserve(nAppendSrcRecVars + nSrcRecVars);

  if (append) {
    ez::ezNcSlicer<int> *ncAppendSrcSlicer = new ez::ezNcSlicer<int>;
    ncAppendSrcSlicer->build(&ncAppendSrc);
    vInt.clear();
    ncAppendSrc.getRecVarIds(vInt);
    n = vInt.size();
    for(i=0; i < n; ++i) {
      pSrcSlit = new ez::ezNcSlabIterator<int>;
      srcIters.vector.push_back(pSrcSlit);
      varSrc = ncAppendSrc.getVar(vInt[i]);
      if(debug) printf("%s:%d append file rec var %s\n", __FILE__, __LINE__, varSrc->name.c_str());

      pSrcSlit->setVar(varSrc, ncAppendSrcSlicer, true);
      // Make analog for dst output file.
      pDstSlit = new ez::ezNcSlabIterator<int>;
      dstIters.vector.push_back(pDstSlit);
      // Make copy of iteration map, then set dst var.
      *pDstSlit = *pSrcSlit;
      pDstSlit->var = ncDst.getVar(varSrc->name);
      pDstSlit->normalize();
    }
    if (ncAppendSrcSlicer) delete ncAppendSrcSlicer;
  }
  
  idit = srcVarIdsToExport.begin();
  idend = srcVarIdsToExport.end(); 
  for(; idit != idend; ++idit) {
    i = *idit;
    varSrc = ncSrc.getVar(i);
    if (!varSrc || !varSrc->hasrec) continue;
    if(debug) printf("%s:%d input file rec var %s\n", __FILE__, __LINE__, varSrc->name.c_str());

    pSrcSlit = new ez::ezNcSlabIterator<int>;
    srcIters.vector.push_back(pSrcSlit);
    pSrcSlit->setVar(varSrc, &ncSrcSlicer, true);
    // Make analog for dst output file.
    pDstSlit = new ez::ezNcSlabIterator<int>;
    dstIters.vector.push_back(pDstSlit);
    // Make copy of iteration map, then set dst var.
    *pDstSlit = *pSrcSlit;
    pDstSlit->var = ncDst.getVar(varSrc->name);
    pDstSlit->normalize();
  }
  
  n = srcIters.vector.size();
  if (n < 1) return 0;
  
  #define COPYSLABTASK() { /**/ \
    task = new ez::ezNcCopySlabTask<size_t>; \
    tasks.vector.push_back(task); \
    srcIters.vector[i]->getSlab(&task->slabSrc); \
    dstIters.vector[i]->getSlab(&task->slabDst); \
    srcIters.vector[i]->incStatic(); \
    dstIters.vector[i]->incStatic(); \
    if (debug) { \
      tmpstring = task->slabSrc.var->name; \
      tmpstring.append(" slabIn: "); \
      task->slabSrc.toString(tmpstring); \
      tmpstring.append(" slabOut: "); \
      task->slabDst.toString(tmpstring); \
      printf("%s:%d %s\n", __FILE__, __LINE__, tmpstring.c_str()); \
    } \
  }
  
  // Create the rec tasks by simply iterating. If one is done, all are done.
  while(!srcIters.vector[0]->isDone()) {
    //printf("%s:%d\n", __FILE__, __LINE__);
    for(i=0; i < n; ++i) {
      //printf("%s:%d\n", __FILE__, __LINE__);
      srcIters.vector[i]->resetStatic();
      dstIters.vector[i]->resetStatic();
      
      // For all static dims of var.
      while(!srcIters.vector[i]->isDoneStatic()) {
        COPYSLABTASK();
      }
      
      // If there are no static dims, var can still be 0D recvar.
      if (srcIters.vector[i]->var->hasrec && (srcIters.vector[i]->var->getNumDims() == 1))
        COPYSLABTASK();
      
      srcIters.vector[i]->incRec();
      dstIters.vector[i]->incRec();
    }
  }
  //printf("%s:%d\n", __FILE__, __LINE__);
  GetPeakMemory();
  if(debug) printf("%s:%d CreateTasks done\n", __FILE__, __LINE__);
  
  return 0;
}
//#####################################################################
void NcX::MakeHistory() {
  ez::ezNcBuffers attbuf;
  std::string history = NCX_VERSION;

  // Add ncx version global attribute.
  attbuf.nc = history.size();
  attbuf.allocate();
  memcpy(attbuf.c, history.c_str(), attbuf.nc*sizeof(char));
  ncDst.putAtt(0, "NCX", NC_CHAR, attbuf);  
  attbuf.reset();
  history.clear();
  
  if (opt.isSet("-h")) return;
  
  // Append to history global att, if it exists.
  if (!ncAppendSrc.filename.empty()) {
    std::vector<std::string> strings;
    ncAppendSrc.getAllGlobalAttNames(strings);
      
    history = "history";
    if (ez::find_first(strings, history) >= 0) {
      ncAppendSrc.getAtt(0, "history", attbuf);
      if (attbuf.nc) {
        history.assign(attbuf.c, attbuf.nc);
        history.append("\n");
      } else
        history = "";
    } else
      history = "";
  }
    
  attbuf.reset();
  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  history.append( asctime(timeinfo) );
  // Remove trailing newline auto inserted by asctime.
  history.erase( history.rfind('\n') );
  
  for(int i=0; i < argc; ++i) {
    history.append(" ");
    history.append(argv[i]);
  }
  
  // Write att.
  attbuf.nc = history.size();
  attbuf.allocate();
  memcpy(attbuf.c, history.c_str(), attbuf.nc*sizeof(char));
  ncDst.putAtt(0, "history", NC_CHAR, attbuf);  
}
//#####################################################################
int NcX::OpenFiles() {
  int status = 0;
  if(debug) {
    printf("%s:%d OpenFiles\n", __FILE__, __LINE__); fflush(stdout);
  }

  // Open input file as read only.
  if (status = ncSrc.openReadOnly()) {
    fprintf(stderr, "ERROR: Failed to open input file \"%s\". Aborting.\n", ncSrc.filename.c_str());
    return status;
  }
  
  // If append file, open it as read only.
  bool append = !ncAppendSrc.filename.empty();
  if (append) {
    if (status = ncAppendSrc.openReadOnly()) {
      fprintf(stderr, "ERROR: Failed to open append file \"%s\". Aborting.\n", ncAppendSrc.filename.c_str());
      return status;
    }    
  }
  
  //---------------------------------------------
  // Create output file. Assume clobber ok at this point. Overwrite already checked in 'SetFileNames'.
  
  // Set blocksize.
  int format;
  opt.get("-b")->getInt(format);
  ncDst.bufrsizehint = (size_t)format;
  
  // Set header pad.
  opt.get("-P")->getInt(format);
  ncDst.h_minfree = (size_t)format;
  
  // If no -3/-6/-H given, use append file format or use input file format.
  if (opt.isSet("-3")) {
    format = NC_FORMAT_CLASSIC;
  } else if (opt.isSet("-6")) {
    format = NC_FORMAT_64BIT;
  } else if (opt.isSet("-H")) {
    format = NC_FORMAT_NETCDF4;
  } else if (append) {
    format = ncAppendSrc.format;
  } else {
    format = ncSrc.format;
  }
  
  if(debug) {
    std::string formatstr;
    ez::formatToString(format, formatstr);
    printf("%s:%d Setting output format to %s.\n", __FILE__, __LINE__, formatstr.c_str()); fflush(stdout);
  }
  
  switch(format) {
    case NC_FORMAT_64BIT:
      if(status = ncDst.create64Clobber()) {
        fprintf(stderr, "ERROR: Failed to create 64bit output file \"%s\". Aborting.\n", ncDst.filename.c_str());
        return status;
      }
      break;
    case NC_FORMAT_NETCDF4: case NC_FORMAT_NETCDF4_CLASSIC:
      if(status = ncDst.createHDFClobber()) {
        fprintf(stderr, "ERROR: Failed to create netcdf4 HDF output file \"%s\". Aborting.\n", ncDst.filename.c_str());
        return status;
      }
      break;
    default: 
      if(status = ncDst.createClobber()) {
        fprintf(stderr, "ERROR: Failed to create classic output file \"%s\". Aborting.\n", ncDst.filename.c_str());
        return status;
      }
      break;
  }
  
  return status;
}
//#####################################################################
int NcX::ProcessMetadata() {
  int status = 0;
  if(debug) {
    printf("%s:%d ProcessMetadata\n", __FILE__, __LINE__); fflush(stdout);
  }

  // If we're doing append.
  bool append = !ncAppendSrc.filename.empty();
  
  //--------------------------------------
  // Load all the metadata in each input file.
  if (status = ncSrc.loadAllMetadata()) return status;
  if (append && (status = ncAppendSrc.loadAllMetadata())) return status;
  
  if (status = ChooseSrcVars()) return status;
  SliceSrcDims();
  
  std::set<int>::iterator it, end;
  ez::ezNcDim *dim1, *dim2;
  int id, i;
  size_t size;
  ez::ezNcVar *var1, *var2;

  if (append) {
    if(debug) {
      printf("%s:%d copy dim defs from ncAppendSrc\n", __FILE__, __LINE__);
      fflush(stdout);
    }
    // Copy dims of append file to output file.
    if (status = ncDst.copyDimDefs(ncAppendSrc)) return status;
    
    if(debug) {
      printf("%s:%d copy var defs from ncAppendSrc\n", __FILE__, __LINE__);
      fflush(stdout);
    }
    // Copy vars of append file to output file,
    // but skip any that will be overrided by source.
    std::vector<std::string> appendSrcVarNames, userVarNames, diffVarNames;
    ncAppendSrc.getAllVarNames(appendSrcVarNames);
    
    // Get the user selected (and depenedent) var names.
    it = srcVarIdsToExport.begin();
    end = srcVarIdsToExport.end();
    for(; it != end; ++it) {
      id = *it;
      var1 = ncSrc.getVar(id);
      if (var1) userVarNames.push_back(var1->name);
    }
    
    ez::not_in_second(appendSrcVarNames, userVarNames, diffVarNames);
    if (status = ncDst.copyVarDefs(ncAppendSrc, diffVarNames)) return status;
    
    if(debug) {
      printf("%s:%d copy global atts from ncAppendSrc\n", __FILE__, __LINE__);
      fflush(stdout);
    }
    // Copy global atts.
    if (status = ncDst.copyAtts(ncAppendSrc,0,0)) return status;
  }
  
  // Define dims for output.
  if(debug) {
    printf("%s:%d define dims\n", __FILE__, __LINE__);
    fflush(stdout);
  }

  // Define dims from input src (and not already defined from append src).
  it = srcDimIdsToExport.begin();
  end = srcDimIdsToExport.end();
  for(; it != end; ++it) {
    id = *it;
    dim1 = ncSrc.getDim(id);
    if(!dim1) continue;
    size = (size_t)ncSrcSlicer.getDimSize(id);

    // Check if dim already exists due to append.
    dim2 = ncDst.getDim(dim1->name);
    if (dim2) {
      if (dim2->size == size)
        continue;
      else {
        fprintf(stderr, "ERROR: Dimension \"%s\" has size of %d in input file while output file has it as size %d. Aborting.\n", dim1->name.c_str(), (int)size, (int)dim2->size);
        return 1;
       }
    }
    
    id = ncDst.createDim(dim1->name.c_str(), size, dim1->isrec);
    if (id < 0) {
      fprintf(stderr, "ERROR: Failed to define dimension \"%s\" in output file.\n", dim1->name.c_str());
      return 1;
    }
  }
  
  // Define vars for output.
  if(debug) {
    printf("%s:%d define vars\n", __FILE__, __LINE__);
    fflush(stdout);
  }
  
  it = srcVarIdsToExport.begin();
  end = srcVarIdsToExport.end();
  std::vector<int> dimids;
  std::vector<std::string> strings;
  
  for(; it != end; ++it) {
    id = *it;
    var1 = ncSrc.getVar(id);
    if (!var1) continue;
    var2 = ncDst.getVar(var1->name);
    if (var2) {
      // This shouldn't happen. Append mode wouldn't copy old var into dst.
      continue;
    }
    
    ncSrc.getDimNames(id, strings);
    ncDst.getDimIds(strings, dimids);
    if(debug) {
      printf("%s:%d Defining output var \"%s\" with type %d with dimids = [", __FILE__, __LINE__, var1->name.c_str(), var1->type);
      for(i=0; i < dimids.size(); ++i) printf(" %d", (int)dimids[i]);
      printf(" ] and dimlens = [");
      for(i=0; i < dimids.size(); ++i) printf(" %d", (int)ncDst.getDim(dimids[i])->size);
      printf(" ]\n");
    }
    
    id = ncDst.createVar(var1->name.c_str(), var1->type, dimids.size(), &(dimids[0]));
    if (id < 0) {
      fprintf(stderr, "ERROR: Failed to define variable \"%s\" in output file.\n", var1->name.c_str());
      continue;
    }
    
    // Copy atts.
    if (status = ncDst.copyAtts(ncSrc, var1->name.c_str(), var1->name.c_str())) return status;
  }
  
  MakeHistory();
  
  status = ncDst.enddef();
 
  return status;
}
//#####################################################################
void NcX::ReportMemory() {
  if (!opt.isSet("-m")) return;
  if (!verbose) return;
    
  int mem = (int)(mem_virt/1000);
  std::ostringstream stm;
  stm << mem;
  std::string str = stm.str();
  for(int i=str.size()-3; i>0; i-=3) str.insert(i, 1, ',');
  printf("Peak Virtual Memory: %s MB\n", str.c_str());
  
  mem = (int)(mem_rss/1000);
  stm.str(""); stm.clear();
  stm << mem;
  str = stm.str();
  for(int i=str.size()-3; i>0; i-=3) str.insert(i, 1, ',');
  printf("Peak Resident Set Size: %s MB\n", str.c_str());
}
//#####################################################################
int NcX::Run() {
  int status;
  if(status = SetFileNames()) return status;
  if(status = OpenFiles()) return status;
  if(status = ProcessMetadata()) return status;
  if(status = CreateTasks()) return status;
  if(status = RunTasks()) return status;
  if(status = CloseFiles()) return status;
  ReportMemory();
  
  return 0;
}
//#####################################################################
int NcX::RunTasks() {
  if(debug) printf("%s:%d RunTasks\n", __FILE__, __LINE__);
  int status = 0;
  ez::ezNcCopySlabTask<size_t> *task=0;
  int i = 0;
  int n = tasks.vector.size();
  bool doProgress = opt.isSet("-p") && verbose;
  bool doEta = opt.isSet("-e") && verbose;
  
  // Make buffers to fit all exporting variables' atomic slabs.
  ez::ezNcBuffers buf;
  
  if (!ncAppendSrc.filename.empty())
    ncAppendSrc.reserve(buf);
    
  std::vector<int> ncSrcVarIds(srcVarIdsToExport.begin(), srcVarIdsToExport.end());
  ncSrc.reserve(buf, ncSrcVarIds);
  ncSrcVarIds.clear();
  
  buf.allocate();
  
  if (doProgress) {
    // Number of tasks plus one for closing the files.
    progress.n = n + 1;
    progress.start();
  }
  
  std::vector<float> taskSizes;
  if (doEta) {
    // Compute task sizes.
    taskSizes.reserve(n);
    float totalMB = 0;
    float size;
    for(i=0; i < n; ++i) {
      size = tasks.vector[i]->slabSrc.getNumBytes() / 1.0e6;
      taskSizes.push_back( size );
      totalMB += size;
    }
    
    eta.n = totalMB + 1.0; // Add one fake MB as marker when files close.
    eta.units = "MB";
    eta.start();
  }
    
  for(i=0; i < n; ++i) {
    task = tasks.vector[i];
    if (task) {
      if (status = task->run(buf)) 
        return status;
    }
    
    if (doProgress) ++progress;
    if (doEta) eta += taskSizes[i];
  }
  
  GetPeakMemory();
  if(debug) printf("%s:%d RunTasks done\n", __FILE__, __LINE__);
  
  return status;
}
//#####################################################################
int NcX::SetFileNames() {
  if(debug) {
    printf("%s:%d SetFileNames\n", __FILE__, __LINE__); fflush(stdout);
  }
  
  bool append = opt.isSet("-A");
  bool overwrite = opt.isSet("-O");
  
  ncSrc.filename = (opt.lastArgs[0] ? *opt.lastArgs[0] : "");
  ncDst.filename = (opt.lastArgs[1] ? *opt.lastArgs[1] : "");

  if (!ez::FileExists(ncSrc.filename.c_str())) {
    std::cerr << "ERROR: Input file \"" << ncSrc.filename << "\" doesn't exist. Aborting.\n\n";
    return 1;
  }
  
  bool dstExists = ez::FileExists(ncDst.filename.c_str());
  if (append) {
    if (!dstExists) {
      if (verbose || debug)
        std::cout << "WARNING: Append file \"" << ncDst.filename << "\" doesn't exist. Will create new file.\n\n";
    } else {
      // Append mode enabled.
      ncAppendSrc.filename = ncDst.filename;
      // Make tmp filename for new append-into file.
      // Generate temporary file so we don't overwrite existing output file until we copy appended contents. Gives illusion of appending.
      // Get any possible directory path prefix from existing file to re-use for temporary file's path prefix.
      std::string dir;
      ez::DirName(ncAppendSrc.filename, dir);
      ez::TmpFilename(dir, ".nc");      
      ncDst.filename = dir;
      if (debug) {
        std::cout << "Temporary filename for appending: " << ncDst.filename.c_str() << std::endl;
      }
      if (ncDst.filename.empty()) {
        std::cerr << "ERROR: Failed to create temporary filename for output in append mode.\n\n";
        return 1;
      }
    }
  } else if (dstExists && !overwrite) {
    std::cerr << "ERROR: Output file \"" << ncDst.filename << "\" already exists and overwrite not enabled with -O. Aborting.\n\n";
    return 1;
  }
  
  return 0;
}
//#####################################################################
int NcX::SetOptions(int argc, const char *argv[]) {
  // Save these just to make history global attribute.
  this->argc = argc;
  this->argv = argv;
  
  opt.overview = "ncx -- NetCDF extractor and appender.";
  opt.syntax = "ncx [OPTIONS] in.nc out.nc";
  opt.example = "# Extract any variable that has ps or temp in name and show peak memory, throughput and time to process.\nncx -m --eta --regex -v ps,temp in.nc out.nc\n\n# Subset and reverse dim order using fortran indices.\nncx -F -i level,9,8,7,6 -v rh in.nc out.nc\n\n# Select only specific pressure levels, exclude var precip and overwrite output.\nncx --val level,800,900,1000 -x -v precip -O in.nc out.nc\n\n# Append all vars with names that start with 's' and don't alter global history attribute.\nncx -A --regex -v '^s' -h in.nc out.nc\n\n";
  Version(opt.footer);

  // Used by header_pad.
  ez::ezOptionValidator *validateUINT32 = new ez::ezOptionValidator(ez::ezOptionValidator::UINT32);

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Force output file to netCDF3 classic (32-bit offset) format. Default will use append file's or input file's format, in that order.", // Help description.
    "-3",
    "--32",
    "--classic"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Force output file to netCDF3 (64-bit offset) format. Default will use append file's or input file's format, in that order.", // Help description.
    "-6",
    "--64",
    "--64bit"
  );
  
  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Append to existing output file.", // Help description.
    "-A",
    "--apn",
    "--append"
  );

  opt.add(
    "65536", // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Request custom I/O blocksize, in bytes. Default is 65536 (64KB).", // Help description.
    "-b",
    "--blksz",
    "--blocksize",
    validateUINT32
  );
  
  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Don't process associated coordinate variables. Coordinate variables are 1D variables that have a dimension name.", // Help description.
    "-C",
    "--nocoords"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Print detailed messages for debugging. Overrides -q.", // Help description.
    "-D",
    "--dbg",
    "--debug"
  );

  opt.add(
    "", // Default.
    0, // Required?
    -1, // Number of args expected.
    ',', // Delimiter if expecting multiple args.
    "Dimension name with optional slice limits and stride. Syntax:\ndim,min[,max[,stride]]\n\nIf min/max are integers, then they will be used as indices. If they are reals (have a decimal), then they will be used as values only if there exists a coordinate var with the same name as the dim. See also --ind and --val. Extended regular expression supported for dim name with -R.", // Help description.
    "-d",
    "--dim",
    "--dimension"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Show progress, I/O rate, elapsed and ETA time. See --progress for simpler version.", // Help description.
    "-e",
    "--eta"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Use 1-based fortran indexing.", // Help description.
    "-F",
    "--ftn",
    "--fortran"
  );

  opt.add(
    "0", // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Display usage instructions.\nThere is a choice of three different layouts for description alignment. Your choice can be any one of the following to suit your style:\n\n0 - align (default)\n1 - interleave\n2 - stagger", // Help description.
    "-u",
    "--usage",
    "--help" // Flag token.
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Force output file to netCDF4 (HDF) format. Default will use append file's or input file's format, in that order.", // Help description.
    "-H",
    "--HDF",
    "--hdf"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Do not append to \"history\" global attribute.", // Help description.
    "-h",
    "--hst",
    "--history"
  );

  opt.add(
    "0", // Default.
    0, // Required?
    1, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Allow extra header space, in bytes.", // Help description.
    "-P",
    "--hdr_pad",
    "--header_pad",
    validateUINT32
  );

  opt.add(
    "", // Default.
    0, // Required?
    -1, // Number of args expected.
    ',', // Delimiter if expecting multiple args.
    "List of dim indices to extract. See also --dim and --val. Syntax:\ndim,i1[,i2[...]]", // Help description.
    "-i",
    "--ind",
    "--index"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Report memory usage (virtual and resident set size).", // Help description.
    "-m",
    "--mem",
    "--memory"
  );
  
  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Overwrite existing output file.", // Help description.
    "-O",
    "--ovr",
    "--overwrite"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Show simple progress bar. See --eta for advanced version.", // Help description.
    "-p",
    "--progress"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Turn off all info and warning messages (errors will always print).", // Help description.
    "-q",
    "--quiet"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Enable regular expressions for dim and var names.", // Help description.
    "-R",
    "--regex"
  );

  opt.add(
    "", // Default.
    0, // Required?
    -1, // Number of args expected.
    ',', // Delimiter if expecting multiple args.
    "Variable(s) to process, separated by commas.\nExtended regular expressions supported with -R.", // Help description.
    "-v",     // Flag token.
    "--var",
    "--variable"
  );

  opt.add(
    "", // Default.
    0, // Required?
    -1, // Number of args expected.
    ',', // Delimiter if expecting multiple args.
    "List of dim values to extract from coordinate var. See also --dim and --ind. Syntax:\ndim,val1[,val2[...]]", // Help description.
    "-V",
    "--val",
    "--values"
  );

  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Display version info.", // Help description.
    "-ver",
    "--version"
  );
  
  opt.add(
    "", // Default.
    0, // Required?
    0, // Number of args expected.
    0, // Delimiter if expecting multiple args.
    "Extract all variables EXCEPT those specified with -v.", // Help description.
    "-x",
    "--xcl",
    "--exclude"
  );

  opt.parse(argc, argv);

  if (opt.isSet("-u")) {
    Usage();
    return 1;
  }

  if (opt.isSet("--version")) {
    std::string ver;
    Version(ver);
    std::cout << ver << std::endl;
    return 1;
  }

  if (opt.lastArgs.size() != 2) {
    std::cerr << "ERROR: Expected one input file name and one output file name. See --help.\n\n";
    return 1;
  } 
  
  if (!opt.isSet("-v")) {
    fprintf(stderr, "ERROR: No variables specified with -v. Aborting. See --help.\n");
    return 1;
  }

  std::vector<std::string> badOptions;
  int i;
  if(!opt.gotRequired(badOptions)) {
    for(i=0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Missing required option " << badOptions[i] << ". See --help.\n\n";

    return 1;
  }

  if(!opt.gotExpected(badOptions)) {
    for(i=0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Got unexpected number of arguments for option " << badOptions[i] << ". See --help.\n\n";
            
    return 1;
  }

  std::vector<std::string> badArgs;
  if(!opt.gotValid(badOptions, badArgs)) {
    for(i=0; i < badOptions.size(); ++i)
      std::cerr << "ERROR: Got invalid argument \"" << badArgs[i] << "\" for option " << badOptions[i] << ". See --help.\n\n";

    return 1;
  }
  
  debug = opt.isSet("-D");
  verbose = !opt.isSet("-q");
  
  if (debug) {
    std::string options;
    opt.prettyPrint(options);
    std::cout << options << std::endl;
  }
   
  return 0;
}
//#####################################################################
void NcX::SliceSrcDims() {
  if(debug) printf("%s:%d SliceSrcDims\n", __FILE__, __LINE__);
  
  std::vector< std::vector<std::string> > strings2d;
  std::vector<std::string> dimnames, allSrcDimNames;
  int i, j, nlists, ndims;
  
  bool regex = opt.isSet("-R");
  ncSrc.getAllDimNames(allSrcDimNames);

  // Are these fortran 1-based indices?
  if (opt.isSet("-F"))
    ncSrcSlicer.offset = -1;
  
  if (opt.isSet("-d")) {
    // Slice by index or value min/max/stride.
    std::string min, max, stride;
    opt.get("-d")->getMultiStrings(strings2d);
    // Each vector is for a single dim's -d tuple.
    nlists = strings2d.size();
    for(i=0; i < nlists; ++i) {
      switch(strings2d[i].size()) {
        case 4: stride = strings2d[i][3];
        case 3: max = strings2d[i][2];
        case 2: min = strings2d[i][1]; break;
        default: continue; // No min,max, or stride set, so skip.
      }
      if (debug) {
        printf("%s:%d i=%d, dimpattern=%s, min=%s, max=%s, stride=%s\n", __FILE__, __LINE__, i, strings2d[i][0].c_str(), min.c_str(), max.c_str(), stride.c_str());
      }
      
      bool useValues = ((min.find('.') != std::string::npos) ||
                        (max.find('.') != std::string::npos));
      dimnames.clear();
      if (regex) {
        // Do regex lookup of dim name(s).
        ez::regex(strings2d[i][0], allSrcDimNames, dimnames, 1, 0);
      } else {
        if ((j = ez::find_first(allSrcDimNames, strings2d[i][0])) >= 0)
          dimnames.push_back(strings2d[i][0]);
      }
      ndims = dimnames.size();
      for(j=0; j < ndims; ++j) {
        // Check if min or max has decimal point to denote values.
        if (useValues) {
          ncSrcSlicer.setDimStrideValues(dimnames[j], min, max, stride);
        } else {
          // Otherwise, they're just plain indices.
          ///printf("%s:%d min=%s, max=%s, stride=%s\n", __FILE__, __LINE__, min.c_str(), max.c_str(), stride.c_str());
          ncSrcSlicer.setDimStride(dimnames[j], min, max, stride);
        }
      }
    }
  }
  
  if (opt.isSet("-i")) {
    // Slice by index list.
    opt.get("-i")->getMultiStrings(strings2d);
    nlists = strings2d.size();
    for(i=0; i < nlists; ++i) {
      // Do regex lookup of dim name(s).
      dimnames.clear();
      if (regex) {
        ez::regex(strings2d[i][0], allSrcDimNames, dimnames, 1, 0);
      } else {
        if ((j = ez::find_first(allSrcDimNames, strings2d[i][0])) >= 0)
          dimnames.push_back(strings2d[i][0]);
      }
      ndims = dimnames.size();
      for(j=0; j < ndims; ++j) {
        // Have it parse the entire list, except for the dim name itself.
        ncSrcSlicer.setDimIndices(dimnames[j], strings2d[i], 1, strings2d[i].size()-1);
      }
    }
  }
  
  if (opt.isSet("-V")) {
    // Slice by value list.
    opt.get("-V")->getMultiStrings(strings2d);
    nlists = strings2d.size();
    for(i=0; i < nlists; ++i) {
      // Do regex lookup of dim name(s).
      dimnames.clear();
      if (regex) {
        ez::regex(strings2d[i][0], allSrcDimNames, dimnames, 1, 0);
      } else {
        if ((j = ez::find_first(allSrcDimNames, strings2d[i][0])) >= 0)
          dimnames.push_back(strings2d[i][0]);
      }
      ndims = dimnames.size();
      for(j=0; j < ndims; ++j) {
        // Have it parse the entire list, except for the dim name itself.
        ncSrcSlicer.setDimValues(dimnames[j], strings2d[i], 1, strings2d[i].size()-1);
      }
    }
  }
  
  ncSrcSlicer.build(&ncSrc);
  if (debug) {
    std::vector<int> alldimids;
    ncSrc.getAllDimIds(alldimids);
    ndims = alldimids.size();
    for(j=0; j < ndims; ++j) {
      printf("%s:%d dimname = %s, sliced? %d, size = %d\n", __FILE__, __LINE__, ncSrc.getDim(alldimids[j])->name.c_str(), ncSrcSlicer.isDimSliced(alldimids[j]), (int)ncSrcSlicer.getDimSize(alldimids[j]));
    }
  }

  if(debug) printf("%s:%d SliceSrcDims done\n", __FILE__, __LINE__);
}
//#####################################################################
void NcX::Usage() {
  std::string usage;
  int layout = 0;
  opt.get("--help")->getInt(layout);
  
  switch(layout) {
  case 0:
    opt.getUsage(usage,78,ez::ezOptionParser::ALIGN);
    break;
  case 1:
    opt.getUsage(usage,78,ez::ezOptionParser::INTERLEAVE);
    break;
  case 2:
    opt.getUsage(usage,78,ez::ezOptionParser::STAGGER);
    break;
  }
  
  std::cout << usage << std::endl;
};
//#####################################################################
void NcX::Version(std::string& s) {
  s.append("ncx ");
  s.append(NCX_VERSION);
  s.append(" Copyright (C) 2011,2012 Remik Ziemlinski\nThis program is free and without warranty.\nBuilt with NetCDF ");
  s.append(nc_inq_libvers());
  s.append("\n");
}
//#####################################################################
#endif // NCX_H
