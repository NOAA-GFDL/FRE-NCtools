/*
Utility class following the NetCDF Climate and Forecast Metadata Conventions (CF) 1.0.

v0.2.0 20111202 rsz Created.
*/
#ifndef EZNCCF_H
#define EZNCCF_H

#include <ezNc.hpp>
#include <ezStringUtil.hpp>

namespace ez {

class ezNcCF {
public:  
  // Return true if attributes designate var as Longitude or X coord var.
  inline static bool isX(ezNcVar* var);
  // Return true if attributes designate var as Latitude or Y coord var.
  inline static bool isY(ezNcVar* var);
  // Return true if attributes designate var as Vertical or Z coord var.
  inline static bool isZ(ezNcVar* var);
  // Return true if attributes designate var as Time or T coord var.
  inline static bool isT(ezNcVar* var);
  // Return true if var has attribute "coordinates" for curvilinear grid.
  inline static bool isWarped(ezNcVar* var);
};
//#####################################################################
bool ezNcCF::isX(ezNcVar* var) {
  if ((var == 0) || (var->nc == 0 )) return false;
  
  const char* name = var->name.c_str();
  std::string string;
  
  if (var->hasAtt("standard_name")) {
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    if (string.compare("longitude") == 0) return true;
  }
  
  if (var->hasAtt("long_name")) {
    string.clear();
    var->nc->getAttString(name, "long_name", string);
    ez::lower(string);
    if (string.compare("longitude") == 0) return true;
  }
  
  if (var->hasAtt("axis")) {
    string.clear();
    var->nc->getAttString(name, "axis", string);
    if (string.compare("X") == 0) return true;
    if (string.compare("x") == 0) return true;
  }

  if (var->hasAtt("units")) {
    const char *units[] = {"degrees_east", "degree_east", "degree_e", "degrees_e", "degreee", "degreese", 0};
    string.clear();
    var->nc->getAttString(name, "units", string);
    ez::lower(string);
    for (const char **it = units; *it; ++it) {
      ///printf("string=%s, it=%s\n", string.c_str(), *it);
      if (string.compare(*it) == 0) return true;
    }
  }
  
  if (var->hasAtt("standard_name")) {
    string.clear();
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    if (string.compare("grid_longitude") == 0) return true;
  }

  return false;
};
//#####################################################################
bool ezNcCF::isY(ezNcVar* var) {
  if ((var == 0) || (var->nc == 0 )) return false;
  
  const char* name = var->name.c_str();
  std::string string;
  
  if (var->hasAtt("standard_name")) {
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    if (string.compare("latitude") == 0) return true;
  }

  if (var->hasAtt("long_name")) {
    string.clear();
    var->nc->getAttString(name, "long_name", string);
    ez::lower(string);
    if (string.compare("latitude") == 0) return true;
  }

  if (var->hasAtt("axis")) {
    string.clear();
    var->nc->getAttString(name, "axis", string);
    if (string.compare("Y") == 0) return true;
    if (string.compare("y") == 0) return true;
  }
  
  if (var->hasAtt("units")) {
    const char *units[] = {"degrees_north", "degree_north", "degree_n", "degrees_n", "degreen", "degreesn", 0};
    string.clear();
    var->nc->getAttString(name, "units", string);
    ez::lower(string);
    for (const char **it = units; *it; ++it) {
      if (string.compare(*it) == 0) return true;
    }
  }

  if (var->hasAtt("standard_name")) {
    string.clear();
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    if (string.compare("grid_latitude") == 0) return true;
  }

  return false;
}
//#####################################################################
bool ezNcCF::isZ(ezNcVar* var) {
  if ((var == 0) || (var->nc == 0 )) return false;
  
  const char* name = var->name.c_str();
  std::string string;

  if (var->hasAtt("standard_name")) {
    const char *standard_names[] = {"atmosphere_sigma_coordinate","atmosphere_hybrid_sigma_pressure_coordinate", "atmosphere_hybrid_height_coordinate", "atmosphere_sleve_coordinate", "ocean_sigma_coordinate", "ocean_s_coordinate", "ocean_sigma_z_coordinate", "ocean_double_sigma_coordinate", 0};
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    for (const char **it = standard_names; *it; ++it) {
      if (string.compare(*it) == 0) return true;
    }
  }

  if (var->hasAtt("long_name")) {
    string.clear();
    var->nc->getAttString(name, "long_name", string);
    if (string.compare("elevation") == 0) return true;
  }

  if (var->hasAtt("axis")) {
    string.clear();
    var->nc->getAttString(name, "axis", string);
    if (string.compare("Z") == 0) return true;
    if (string.compare("z") == 0) return true;
  }

  if (var->hasAtt("units")) {
    const char *units[] = {"bar", "millibar", "decibar", "atmosphere", "atm", "pascal", "pa", "hpa", "meter", "m", "metre", "kilometer", "km", "level", "layer", "sigma_level", "bars", "millibars", "decibars", "atmospheres", "pascals",  "meters", "metres", "kilometers", "levels", "layers", "sigma_levels", 0};
    string.clear();
    var->nc->getAttString(name, "units", string);
    ez::lower(string);
    for (const char **it = units; *it; ++it) {
      if (string.compare(*it) == 0) return true;
    }
  }

  if (var->hasAtt("positive")) {
    string.clear();
    var->nc->getAttString(name, "positive", string);
    if (string.compare("down") == 0) return true;
    if (string.compare("up") == 0) return true;
  }
  
  return false;
}
//#####################################################################
/*  According to CF, the name itself is arbitrary.
    Only the attributes are used. See CF 1.0, Section 4.4.
*/
bool ezNcCF::isT(ezNcVar* var) {
  if ((var == 0) || (var->nc == 0 )) return false;
  
  const char* name = var->name.c_str();
  std::string string;

  if (var->hasAtt("standard_name")) {
    var->nc->getAttString(name, "standard_name", string);
    ez::lower(string);
    if (string.compare("time") == 0) return true;
  }

  if (var->hasAtt("long_name")) {
    string.clear();
    var->nc->getAttString(name, "long_name", string);
    ez::lower(string);
    if (string.compare("time") == 0) return true;
  }
  
  if (var->hasAtt("axis")) {
    string.clear();
    var->nc->getAttString(name, "axis", string);
    if (string.compare("T") == 0) return true;
    if (string.compare("t") == 0) return true;
  }

  if (var->hasAtt("units")) {
    const char *units[] = {"day", "d", "hour", "hr", "h", "minute", "min", "second","sec", "s", "year", "month", "days", "hours", "minutes", "seconds", "years", "months", 0};
    string.clear();
    var->nc->getAttString(name, "units", string);
    ez::lower(string);
    std::string substr;
    for (const char **it = units; *it; ++it) {
      substr = *it;
      substr.append(" since");
      if (string.find(substr) != std::string::npos) return true;
    }
  }

  return var->hasAtt("month_lengths") || var->hasAtt("calendar") || var->hasAtt("leap_year") || var->hasAtt("leap_month");
}
//#####################################################################
bool ezNcCF::isWarped(ezNcVar* var) {
  if ((var == 0) || (var->nc == 0 )) return false; 
  return var->hasAtt("coordinates");
}
};
#endif // EZNCCF_H