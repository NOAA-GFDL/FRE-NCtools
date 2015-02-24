// 
/*
Utility class to hold and process a collection of Nc files.

Example:
  ez::ezNcFiles files;
  std::vector<std::string> names;
  std::vector<ez::ezNc*> sorted;
  files.setFileNames(names);
  files.load();
  files.getSortedByDatetime(sorted);

v0.2.1 20111231 rsz Created.
*/
#ifndef EZNCFILES_H
#define EZNCFILES_H

#include <ezNc.hpp>
#include <ezNcDateTime.hpp>
#include <map>
#include <vector>

namespace ez {

typedef std::map<std::string,ezNc*> NameToNcMap;
typedef std::map<ezNc*,ezNcDateTime*> NcToDateTimeMap;

class ezNcFiles {
public:
  // File name to file object.
  NameToNcMap files;
  // Files with time coordvar will have datetime info available.
  NcToDateTimeMap datetimes;
  
  // Dtor.
  ~ezNcFiles() { clear(); }
  
  // Closes all the file objects and clears memory.
  inline void clear();
  // Closes the file objects.
  inline void closeAll();
  // Returns vector of items sorted by datetime var data. Files without a datetime coordvar are not output. All files must have been already loaded with "loadAll".
  inline void getSortedByDatetime(std::vector<ezNc*>& sorted);
  // Sorts by file names.
  inline void getSortedByName(std::vector<ezNc*>& sorted);
  // Loads all the metadata info for all the files.
  inline void load();
  // Creates datetime objects for files with time coordvar, only after "load".
  inline void loadDateTime();
  // Set names and allocate ezNc objects, but do nothing else.
  inline void setFileNames(std::vector<std::string>& names);
};
//######################################################
void ezNcFiles::clear() {
  closeAll();

  NameToNcMap::iterator it = files.begin(), end = files.end();
  ezNc* nc = 0;

  for(; it != end; ++it) {
    nc = it->second;
    if (nc) delete nc;
  }
  
  files.clear();
  
  ezNcDateTime* dt = 0;
  NcToDateTimeMap::iterator dtit = datetimes.begin(), dtend = datetimes.end();
  for(; dtit != dtend; ++dtit) {
    nc = dtit->first;
    dt = dtit->second;
    if (nc && dt) delete dt;
  }  
  
  datetimes.clear();
}
//######################################################
void ezNcFiles::closeAll() {
  NameToNcMap::iterator it = files.begin(), end = files.end();
  for(; it != end; ++it)
    if (it->second) it->second->close();
}
//######################################################
void ezNcFiles::getSortedByDatetime(std::vector<ezNc*>& sorted) {
  sorted.clear();
  NameToNcMap::iterator it = files.begin(), end = files.end();
  ezNcDateTime* dt = 0;
  ezNc* nc = 0;

  if (datetimes.empty()) loadDateTime();

  std::vector<ezNcDateTime*> unsorted;
  for(; it != end; ++it) {
    nc = it->second;
    if (nc)
      unsorted.push_back( datetimes[nc] );
  }
  
  std::sort(unsorted.begin(), unsorted.end(), ezNcDateTime::comparePtr);
  
  int i = 0, n = unsorted.size();
  for(; i < n; ++i) {
    dt = unsorted[i];
    if (dt && dt->var)
      sorted.push_back( dt->var->nc );
  }
}
//######################################################
void ezNcFiles::getSortedByName(std::vector<ezNc*>& sorted) {
  sorted.clear();
  NameToNcMap::iterator it = files.begin(), end = files.end();
  // Map should already have file names sorted when iterating.
  for(; it != end; ++it)
    if (it->second) sorted.push_back(it->second);
}
//######################################################
void ezNcFiles::load() {
  NameToNcMap::iterator it = files.begin(), end = files.end();
  ezNc* nc = 0;
  
  for(; it != end; ++it) {
    nc = it->second;
    if (nc) {
      nc->openReadOnly();
      nc->loadAllMetadata();
    }
  }
}
//######################################################
void ezNcFiles::loadDateTime() {
  NameToNcMap::iterator it = files.begin(), end = files.end();
  ezNc* nc = 0;
  
  for(; it != end; ++it) {
    nc = it->second;
    if (nc) {
      ezNcDateTime* dt = new ezNcDateTime;
      datetimes[nc] = dt;
      dt->build(nc);
    }
  }
}
//######################################################
void ezNcFiles::setFileNames(std::vector<std::string>& names) {
  int i = 0;
  int n = names.size();
  ezNc* nc = 0;
  
  clear();
  
  for(; i < n; ++i) {
    nc = new ezNc;
    nc->filename = names[i];
    files[ names[i] ] = nc;
  }
}
//######################################################
};
#endif // EZNCFILES_H