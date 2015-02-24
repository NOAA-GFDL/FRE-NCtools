/*
Utility class to handle time coordinate variables.
"Raw values" referred to below are those in the netcdf file.
"Values" are the double datatype encoded equivalent of "raw values" that allow unified comparison and math.

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
  ezNcDateTime() { clear(); };
  // Dtor.
  ~ezNcDateTime() { clear(); }
  
  // Tries to find time coord var in nc. Returns null if not found.
  // First tests vars using CF rules. If fails, checks 1D vars with a name like "time", "record", "rec", "t".
  inline ezNcVar* autodetect(ezNc* nc);
  // Loads data for set var. If var is null, autodetects time var in nc. Returns 0 on success.
  inline int build(ezNc* nc);
  // Resets memory and members.
  inline void clear();
  // True if first comes before the second. For sorting in ascending datetime order. Note, that only the front elements of the value arrays are compared, which is fast, but also assumes no overlap between datetime value arrays and arrays are sorted and increasing.
  inline static bool compare(ezNcDateTime& first, ezNcDateTime& second);
  inline static bool comparePtr(ezNcDateTime* first, ezNcDateTime* second);
  // Compares front time values as encoded in files.
  // static bool compareRawValues(ezNcDateTime& first, ezNcDateTime& second);
  // True if value falls withing range of values in vector.
  inline bool contains(double datetime);
  /// True if raw file value falls withing range of raw values.
  /// bool containsRaw(double datetime);
  // Returns true if every attribute and value are identical.
  inline bool equals(ezNcDateTime& other);
  // Get inclusive index range for values that fall within datetime range begin,end. Invalids indices get -1. If both ibegin/iend are -1, then given time does not overlap file datetime range.
  inline void getRange(double begin, double end, int& ibegin, int& iend);
  inline void getRange(std::string& begin, std::string& end, int& ibegin, int& iend);
  // Converts file raw values (accounting for units, timezone and calendar) to ezDateTime double.
  template <typename T> void 
  inline fromUnits(T* t, int nt, std::vector<double>& out);
  // Same as above, but from buffer array of certain type.
  inline void fromUnits(ezNcBuffers& buf, nc_type type, std::vector<double>& out);
  // Parses units string to define raw value unit and sets baseline ezDateTime and timezone.
  inline void setUnits(std::string& units);
  // Returns true if the utc values are equal.
  inline bool valuesEqual(ezNcDateTime& o) { return values == o.values; }
  
  // What calendar is time var using.
  char calendar;
  // Encoded values converted from raw file values to ezDateTime doubles.
  std::vector<double> values;
  // Values read straight from file without any conversion.
  ezNcBuffers rawValues;
  // The time coord var.
  ezNcVar *var;
  // Encodes var's units attribute baseline date time and var's calendar.
  double baselineValue;
  // What delta unit used by file's raw values. Set to an enum from ezDateTime, such as DAY, HOUR, or YEAR.
  char unit;
  // Units timezone offsets hours and minutes (both can be negative).
  char tzHours, tzMinutes;
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
     if ( (lower.compare("time")==0) ||
          (lower.compare("rec")==0) ||
          (lower.compare("record")==0) ||
          (lower.compare("t")==0) )
      return vars[i];
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
        calendar = ezDateTime::getCalendarType(tmp);
      } else
        tmp.clear();
    } else
      tmp.clear();
  }
  
  if (tmp.empty()) {
    printf("WARNING: No calendar attribute found for \"%s\" in \"%s\". Assuming no-leap years with 365 days.\n", var->name.c_str(), nc->filename.c_str());
  }
      
  int status;
  
  // Get the units string.
  tmp.clear();
  if (status = nc->getAttString(var->name.c_str(), "units", tmp)) 
    return status;

  if (tmp.empty()) {
    printf("WARNING: No time units attribute found for \"%s\" in \"%s\". Assuming years since 0000-01-01 00:00Z.\n", var->name.c_str(), nc->filename.c_str());
    tmp = "years since 0000-01-01 00:00";
  }

  setUnits(tmp);
  
  // Load all the data for the time var.
  rawValues.reset();

  if (status = nc->reserveEntire(rawValues, var->varid)) 
    return status;

  rawValues.allocate();
  if (status = nc->read(var->varid, &rawValues)) 
    return status;

  // Encode the data into our custom format.
  fromUnits(rawValues, var->type, values);
}
//##########################################################
void ezNcDateTime::clear() {
  baselineValue = 0;
  calendar = 0;
  tzHours = 0;
  tzMinutes = 0;
  unit = 0;
  values.clear();
  var = 0;
}
//##########################################################
bool ezNcDateTime::compare(ezNcDateTime& first, ezNcDateTime& second) {
  if (first.values.empty()) return true;
  if (second.values.empty()) return true;
  return first.values.front() < second.values.front();
}
//##########################################################
bool ezNcDateTime::comparePtr(ezNcDateTime* first, ezNcDateTime* second) {
  if (!first || !second) return false;

  if (first->values.empty()) return true;
  if (second->values.empty()) return true;
  return first->values.front() < second->values.front();
}
//##########################################################
bool ezNcDateTime::contains(double datetime) {
  int n = values.size();
  
  switch(n) {
    case 0: return 0;
    case 1: return values[0] == datetime;
    default: return (values[0] <= datetime) && (datetime <= values[n-1]);
  }
}
//##########################################################
/*
bool ezNcDateTime::containsRaw(double datetime) {
  if (!var) return false;
  int n = rawValues.getSize(var->type);
  
  switch(n) {
    case 0: return 0;
    case 1: 
      return values[0] == datetime;
    default: return (values[0] <= datetime) && (datetime <= values[n-1]);
  }
}*/
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
      out.push_back( ezDateTime::add(baselineValue, delta, calendar) ); \
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
//##########################################################
void ezNcDateTime::getRange(double begin, double end, int& ibegin, int& iend) {
  if (values.empty() || 
    (begin > values.back()) || 
    (end < values.front()) ) {
    ibegin = -1;
    iend = -1;
    return;
  }

  std::vector<double>::iterator it, b=values.begin(), e=values.end();
  
  if (begin < values.front()) {
    ibegin = 0;
  } else {
    it = std::lower_bound(b, e, begin);
    ibegin = int(it - b);
  }
  
  if (end > values.back()) {
    iend = values.size()-1;
  } else {
    it = std::upper_bound(b, e, end);
    iend = int(it - b - 1);
  }
}
//##########################################################
void ezNcDateTime::getRange(std::string& begin, std::string& end, int& ibegin, int& iend) {
  double dBegin, dEnd;
  char ignore;

  ezDateTime::fromString(dBegin, begin, calendar, ignore, ignore);
  ezDateTime::fromString(dEnd, end, calendar, ignore, ignore);
  
  getRange(dBegin, dEnd, ibegin, iend);
}
//##########################################################
bool ezNcDateTime::equals(ezNcDateTime& other) {
  return (calendar == other.calendar) &&
    (baselineValue == other.baselineValue) &&
    (unit == other.unit) &&
    (tzHours == other.tzHours) &&
    (tzMinutes == other.tzMinutes) &&
    (values == other.values);
}
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
      break;
    case 's': case 'S': unit = ezDateTime::SECONDS; break;
    default: unit = ezDateTime::YEARS; break;
  }

  i = units.find("since", i);
  if (i == std::string::npos) return;

  // Convert date/time to converted baseline value.
  std::string substr = units.substr(i+5);
  ezDateTime::fromString(baselineValue, substr, calendar, tzHours, tzMinutes);
}
//##########################################################
};
#endif // EZNCDATETIME_H