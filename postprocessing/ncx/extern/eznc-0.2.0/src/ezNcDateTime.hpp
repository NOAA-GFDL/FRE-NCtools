/*
Utility class to handle time coordinate variables.

v0.2.0 20111203 rsz Created.
*/
#ifndef EZNCDATETIME_H
#define EZNCDATETIME_H

#include <ezNc.hpp>
#include <ezNcCF.hpp>
#include <ezDateTime.hpp>

namespace ez {

class ezNcDateTime {
public:
  // Ctor.
  ezNcDateTime() : unit(ezDateTime::YEARS), var(0) {};
  // Dtor.
  ~ezNcDateTime() { clear(); }
  
  // Tries to find time coord var in nc. Returns null if not found.
  // First tests vars using CF rules. If fails, checks 1D vars with a name like "time", "record", "rec", "t".
  ezNcVar* autodetect(ezNc* nc);
  // Loads data for var. If var is null, autodetects time var. Returns 0 on success.
  int build(ezNc* nc);
  // Resets memory and members.
  void clear();
  // True if first comes before the second. For sorting in ascending datetime order. Note, that only the front elements of the value arrays are compared, which is fast, but also assumes no overlap between datetime value arrays.
  static bool compare(ezNcDateTime& first, ezNcDateTime& second);
  // True if value falls withing range of values in vector.
  bool contains(double datetime);
  // Get inclusive index range for values that fall within datetime range begin,end.
  void getRange(double begin, double end, int& ibegin, int& iend);
  // Converts file raw values (accounting for units & calendar) to ezDateTime double.
  template <typename T> void 
  fromUnits(T* t, int nt, std::vector<double>& out);
  // Same as above, but from buffer array of certain type.
  void fromUnits(ezNcBuffers& buf, nc_type type, std::vector<double>& out);
  // Parses units string to define raw value unit and sets baseline ezDateTime.
  void setUnits(std::string& units);
  
  // Encoded values converted from raw file values to ezDateTime doubles.
  std::vector<double> values;
  // The time coord var.
  ezNcVar *var;
  // Encodes var's units attribute baseline date time and var's calendar.
  ezDateTime baseline;
  // What delta unit used by file's raw values. Set to an enum from ezDateTime, such as DAY, HOUR, or YEAR.
  char unit;
};
//##########################################################
ezNcVar* ezNcDateTime::autodetect(ezNc* nc) {
  if (nc == 0) return 0;
  
  // Check each var and find first that passes CF "isT" test.
  std::vector<int> varids;
  std::vector<ezNcVar*> vars;
  int i, n, nrecvars;
  ezNcVar* var;
  
  // Start with rec vars. Get only 1d vars.
  nc->getRecVarIds(varids);
  n = varids.size();
  for(i=0; i < n; ++i) {
    var = nc->getVar(varids[i]);
    if (var==0) continue;
    if (var->getNumDims() == 1)
      vars.push_back(var);
  }
  
  nrecvars = vars.size();
  for(i=0; i < nrecvars; ++i) {
    var = vars[i];
    if (ezNcCF::isT(vars[i])) return var;
  }
  
  // And try static vars.
  nc->getStaticVarIds(varids);
  n = varids.size();
  for(i=0; i < n; ++i) {
    var = nc->getVar(varids[i]);
    if (var==0) continue;
    if (var->getNumDims() == 1)
      vars.push_back(var);
  }

  // Already checked recvars, so skip those.
  n = vars.size();
  for(i=nrecvars; i < n; ++i) {
    var = vars[i];
    if (ezNcCF::isT(vars[i])) return var;
  }

  // If none found, ad hoc check 1D var with name "time", "rec", "record", "t" regardless if recvar or static.
  std::string lower;
  for(i=0; i < n; ++i) {
     var = vars[i];
     lower = var->name;
     ez::lower(lower);
     if (lower.compare("time")==0) return vars[i];
     else if (lower.compare("rec")==0) return vars[i];
     else if (lower.compare("record")==0) return vars[i];
     else if (lower.compare("t")==0) return vars[i];
  }
}
//##########################################################
int ezNcDateTime::build(ezNc* nc) {
  if (nc==0) return 1;
  
  // No time var was set, so find it ourselves using best guess.
  if (var==0) 
    var = autodetect(nc);
    
  // Anything found?
  if (var==0) return 1;
   
  // Set calendar type.
  std::string tmp;
  if (var->hasAtt("calendar")) {
    if (nc->getAttString(var->name.c_str(), "calendar", tmp) == 0) {
      if (tmp.size() >= 3) {
        baseline.calendar = ezDateTime::getCalendarType(tmp);
      }
    }    
  }
  
  int status;
  
  // Get the units string.
  tmp.clear();
  if (status = nc->getAttString(var->name.c_str(), "units", tmp)) 
    return status;

  setUnits(tmp);
  
  // Load all the data for the time var.
  ezNcBuffers buf;

  if (status = nc->initEntire(buf, var->varid)) 
    return status;

  if (status = nc->read(var->varid, &buf)) 
    return status;

  // Encode the data into our custom format.
  fromUnits(buf, var->type, values);
}
//##########################################################
void ezNcDateTime::clear() {
  baseline.calendar = 0;
  baseline.value = 0;
  values.clear();
  var = 0;
}
//##########################################################
template <typename T>
void ezNcDateTime::fromUnits(T* t, int nt, std::vector<double>& out) {
  int i = 0;
  ezTimeDelta delta;
  out.clear();
  out.reserve(nt);
  
  #define TADDLOOP(UNIT) { \
    for(i; i < nt; ++i) { \
      delta.UNIT = (double)t[i]; \
      out.push_back( baseline.add(delta) ); \
    } \
  }
  
  switch(unit) {
    case ezDateTime::DAYS: TADDLOOP(days); break;
    case ezDateTime::HOURS: TADDLOOP(hours); break;
    case ezDateTime::MINUTES: TADDLOOP(minutes); break;
    case ezDateTime::SECONDS: TADDLOOP(seconds); break;
    case ezDateTime::MONTHS: TADDLOOP(months); break;
    default: TADDLOOP(years); break;
  }
}
//##########################################################
void ezNcDateTime::fromUnits(ezNcBuffers& buf, nc_type type, std::vector<double>& out) {
  switch(type) {
    case NC_BYTE: case NC_UBYTE: fromUnits<unsigned char>(buf.uc, buf.nuc, values); break;
    case NC_CHAR: fromUnits<char>(buf.c, buf.nc, values); break;
    case NC_SHORT: fromUnits<short>(buf.s, buf.ns, values); break;
    case NC_USHORT: fromUnits<unsigned short>(buf.us, buf.nus, values); break;
    case NC_INT: fromUnits<int>(buf.i, buf.ni, values); break;
    case NC_UINT: fromUnits<unsigned int>(buf.ui, buf.nui, values); break;
    case NC_INT64: fromUnits<long long>(buf.l, buf.nl, values); break;
    case NC_UINT64: fromUnits<unsigned long long>(buf.ul, buf.nul, values); break;
    case NC_FLOAT: fromUnits<float>(buf.f, buf.nf, values); break;
    case NC_DOUBLE: fromUnits<double>(buf.d, buf.nd, values); break;
    default: break;
  }
}
/*//##########################################################
double ezNcDateTime::fromString(std::string& u) {
  double res = 0.0;
  if (u.empty()) return res;
  
  size_t i = u.find("since");

  return fromString(u, this->calendar);
}*/
//##########################################################
void ezNcDateTime::setUnits(std::string& units) {
  if (units.empty()) return;
  
  size_t i=0;
  int n = units.size();
  // Parse delta unit string before "since" keyword. 
  // Skip leading whitespace.
  while ((i < n) && (units[i] == ' ')) ++i;
  
  switch(units[i]) {
    case 'd': case 'D': unit = ezDateTime::DAYS; break;
    case 'h': case 'H': unit = ezDateTime::HOURS; break;
    case 'm': case 'M': 
      if ((i+1) < n) {
        switch(units[i+1]) {
          case 'i': case 'I': unit = ezDateTime::MINUTES; break;
          case 'o': case 'O': unit = ezDateTime::MONTHS; break;
        }
      } else
        unit = ezDateTime::MINUTES;
    case 's': case 'S': unit = ezDateTime::SECONDS; break;
    default: unit = ezDateTime::YEARS; break;
  }
  
  i = units.find("since", i);
  if (i == std::string::npos) return;

  // Convert date/time to converted baseline value.
  std::string substr = units.substr(i+5);
  baseline.setValue(substr);
}
//##########################################################
};
#endif // EZNCDATETIME_H