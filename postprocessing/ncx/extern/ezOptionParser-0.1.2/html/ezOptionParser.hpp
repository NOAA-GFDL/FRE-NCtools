/*
This file is part of ezOptionParser.

ezOptionParser is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ezOptionParser is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with ezOptionParser.  If not, see <http://www.gnu.org/licenses/>.

Copyright (C) 2011 Remik Ziemlinski <first d0t surname att gmail>
*/
/*
CHANGELOG

v0.0.0 20110505 rsz Created.
*/
#ifndef EZ_OPTION_PARSER_H
#define EZ_OPTION_PARSER_H

#include <vector>
#include <list>
#include <map>
#include <string>
#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace ez {
/* ################################################################### */
/* 
Compare strings for opts, so short opt flags come before long format flags.
For example, -d < --dimension < --dmn, and also lower come before upper. The default STL std::string compare doesn't do that.
*/
bool CmpOptStringPtr(std::string * s1, std::string * s2) { 
  int c1,c2;
	const char *s=s1->c_str();
	for(c1=0; c1 < s1->size(); ++c1)
		if (isalpha(s[c1])) // locale sensitive.
			break;
			
	s=s2->c_str();
	for(c2=0; c2 < s2->size(); ++c2)
		if (isalpha(s[c2]))
			break;

	// Test which has more symbols before its name.
  if (c1 > c2)
    return false;
  else if (c1 < c2)
    return true;
    
	// Both have same number of symbols, so compare first letter.
	char char1 = s1->at(c1);
	char char2 = s2->at(c2);
  char lo1 = tolower(char1);
  char lo2 = tolower(char2);

  if (lo1 != lo2)
    return lo1 < lo2;

	// Their case doesn't match, so find which is lower.
  char up1 = isupper(char1);
  char up2 = isupper(char2);
    
  if (up1 && !up2)
    return false;
  else if (!up1 && up2)
    return true;
    
  return (s1->compare(*s2)<0); 
};
/* ################################################################### */
/* 
Makes a vector of strings from one string,
splitting at (and excluding) delimiter "token".
*/
void SplitDelim( const std::string& s, const char token, std::vector<std::string*> * result) {
  std::string::const_iterator i = s.begin();
  std::string::const_iterator j = s.begin();
  const std::string::const_iterator e = s.end();

  while(i!=e) {
    while(i!=e && *i++!=token);
    std::string *newstr = new std::string(j, i);
    if (newstr->at(newstr->size()-1) == token) newstr->erase(newstr->size()-1);
    result->push_back(newstr);
    j = i;
  }
};
/* ################################################################### */
// Variant that uses deep copies and references instead of pointers (less efficient).
void SplitDelim( const std::string& s, const char token, std::vector<std::string> & result) {
  std::string::const_iterator i = s.begin();
  std::string::const_iterator j = s.begin();
  const std::string::const_iterator e = s.end();

  while(i!=e) {
    while(i!=e && *i++!=token);
    std::string newstr(j, i);
    if (newstr.at(newstr.size()-1) == token) newstr.erase(newstr.size()-1);
    result.push_back(newstr);
    j = i;
  }
};
/* ################################################################### */
// Variant that uses list instead of vector for efficient insertion, etc.
void SplitDelim( const std::string& s, const char token, std::list<std::string*> & result) {
  std::string::const_iterator i = s.begin();
  std::string::const_iterator j = s.begin();
  const std::string::const_iterator e = s.end();

  while(i!=e) {
    while(i!=e && *i++!=token);
    std::string *newstr = new std::string(j, i);
    if (newstr->at(newstr->size()-1) == token) newstr->erase(newstr->size()-1);
    result.push_back(newstr);
    j = i;
  }
};
/* ################################################################### */
void StringsToInts(std::vector<std::string> & strings, std::vector<int> & out) {
  for(int i=0; i < strings.size(); ++i) {
    out.push_back(atoi(strings[i].c_str()));
  }
};
/* ################################################################### */
void StringsToInts(std::vector<std::string*> * strings, std::vector<int> * out) {
  for(int i=0; i < strings->size(); ++i) {
    out->push_back(atoi(strings->at(i)->c_str()));
  }
};
/* ################################################################### */
void StringsToLongs(std::vector<std::string> & strings, std::vector<long> & out) {
  for(int i=0; i < strings.size(); ++i) {
    out.push_back(atol(strings[i].c_str()));
  }
};
/* ################################################################### */
void StringsToLongs(std::vector<std::string*> * strings, std::vector<long> * out) {
  for(int i=0; i < strings->size(); ++i) {
    out->push_back(atol(strings->at(i)->c_str()));
  }
};
/* ################################################################### */
void StringsToULongs(std::vector<std::string> & strings, std::vector<unsigned long> & out) {
  for(int i=0; i < strings.size(); ++i) {
    out.push_back(strtoul(strings[i].c_str(),0,0));
  }
};
/* ################################################################### */
void StringsToULongs(std::vector<std::string*> * strings, std::vector<unsigned long> * out) {
  for(int i=0; i < strings->size(); ++i) {
    out->push_back(strtoul(strings->at(i)->c_str(),0,0));
  }
};
/* ################################################################### */
void StringsToFloats(std::vector<std::string> & strings, std::vector<float> & out) {
  for(int i=0; i < strings.size(); ++i) {
    out.push_back(atof(strings[i].c_str()));
  }
};
/* ################################################################### */
void StringsToFloats(std::vector<std::string*> * strings, std::vector<float> * out) {
  for(int i=0; i < strings->size(); ++i) {
    out->push_back(atof(strings->at(i)->c_str()));
  }
};
/* ################################################################### */
void StringsToDoubles(std::vector<std::string> & strings, std::vector<double> & out) {
  for(int i=0; i < strings.size(); ++i) {
    out.push_back(atof(strings[i].c_str()));
  }
};
/* ################################################################### */
void StringsToDoubles(std::vector<std::string*> * strings, std::vector<double> * out) {
  for(int i=0; i < strings->size(); ++i) {
    out->push_back(atof(strings->at(i)->c_str()));
  }
};
/* ################################################################### */
char** CommandLineToArgvA(char* CmdLine, int* _argc) {
  char** argv;
  char*  _argv;
  unsigned long   len;
  unsigned long   argc;
  char   a;
  unsigned long   i, j;

  bool  in_QM;
  bool  in_TEXT;
  bool  in_SPACE;

  len = strlen(CmdLine);
  i = ((len+2)/2)*sizeof(void*) + sizeof(void*);

  argv = (char**)malloc(i + (len+2)*sizeof(char));

  _argv = (char*)(((unsigned char*)argv)+i);

  argc = 0;
  argv[argc] = _argv;
  in_QM = false;
  in_TEXT = false;
  in_SPACE = true;
  i = 0;
  j = 0;

  while( a = CmdLine[i] ) {
    if(in_QM) {
      if( (a == '\"') ||
					(a == '\'')) // rsz. Added single quote.
			{
        in_QM = false;
      } else {
        _argv[j] = a;
        j++;
      }
    } else {
      switch(a) {
      case '\"': 
			case '\'': // rsz. Added single quote.
        in_QM = true;
        in_TEXT = true;
        if(in_SPACE) {
          argv[argc] = _argv+j;
          argc++;
        }
        in_SPACE = false;
        break;
      case ' ':
      case '\t':
      case '\n':
      case '\r':
        if(in_TEXT) {
          _argv[j] = '\0';
          j++;
        }
        in_TEXT = false;
        in_SPACE = true;
        break;
      default:
        in_TEXT = true;
        if(in_SPACE) {
          argv[argc] = _argv+j;
          argc++;
        }
        _argv[j] = a;
        j++;
        in_SPACE = false;
        break;
      }
    }
    i++;
  }
  _argv[j] = '\0';
  argv[argc] = NULL;

  (*_argc) = argc;
  return argv;
};
/* ################################################################### */
class OptionGroup {
public:
  OptionGroup() : delim(0), expectArgs(0), isSet(false), isRequired(false) { }

  ~OptionGroup() {
    int i, j;
    for(i=0; i < flags.size(); ++i)
      delete flags[i];
      
    flags.clear();
    clearArgs();
  };

	void clearArgs();
  void getInt(int&);
  void getLong(long&);
  void getULong(unsigned long&);
  void getFloat(float&);
  void getDouble(double&);
  void getString(std::string&);
  void getInts(std::vector<int>&);
  void getLongs(std::vector<long>&);
  void getULongs(std::vector<unsigned long>&);
  void getFloats(std::vector<float>&);
  void getDoubles(std::vector<double>&);
  void getStrings(std::vector<std::string>&);
  void getMultiInts(std::vector< std::vector<int> >&);
  void getMultiLongs(std::vector< std::vector<long> >&);
  void getMultiULongs(std::vector< std::vector<unsigned long> >&);
  void getMultiFloats(std::vector< std::vector<float> >&);
  void getMultiDoubles(std::vector< std::vector<double> >&);
  void getMultiStrings(std::vector< std::vector<std::string> >&);
  
  // defaults value regardless of being set by user.
  std::string defaults;
  // If expects arguments, this will delimit arg list.
  char delim;
  // If not 0, then number of delimited args. -1 for arbitrary number.
  int expectArgs;
  // Descriptive help message shown in usage instructions for option.
  std::string help;
  // 0 or 1.
  bool isRequired;
  // A list of flags that denote this option, i.e. -d, --dimension.
  std::vector< std::string* > flags;
  // If was set (or found).
  bool isSet;
  // Lists of arguments, per flag instance, after splitting by delimiter.
  std::vector< std::vector< std::string* > * > args;
};
/* ################################################################### */
void OptionGroup::clearArgs() {
	int i,j;
	for(i=0; i < args.size(); ++i) {
		for(j=0; j < args[i]->size(); ++j)
			delete args[i]->at(j);
			
		delete args[i];
	}

	args.clear();
	isSet = false;
};
/* ################################################################### */
void OptionGroup::getInt(int & out) {
  if (!isSet) {
    if (defaults.empty())
      out = 0;
    else
      out = atoi(defaults.c_str());
  } else {
    if (args.empty() || args[0]->empty())
      out = 0;
    else {
      out = atoi(args[0]->at(0)->c_str());
    }
  }
};
/* ################################################################### */
void OptionGroup::getLong(long & out) {
  if (!isSet) {
    if (defaults.empty())
      out = 0;
    else
      out = atoi(defaults.c_str());
  } else {
    if (args.empty() || args[0]->empty())
      out = 0;
    else {
      out = atol(args[0]->at(0)->c_str());
    }
  }
};
/* ################################################################### */
void OptionGroup::getULong(unsigned long & out) {
  if (!isSet) {
    if (defaults.empty())
      out = 0;
    else
      out = atoi(defaults.c_str());
  } else {
    if (args.empty() || args[0]->empty())
      out = 0;
    else {
      out = strtoul(args[0]->at(0)->c_str(),0,0);
    }
  }
};
/* ################################################################### */
void OptionGroup::getFloat(float & out) {
  if (!isSet) {
    if (defaults.empty())
      out = 0.0;
    else
      out = (float)atof(defaults.c_str());
  } else {
    if (args.empty() || args[0]->empty())
      out = 0.0;
    else {
      out = (float)atof(args[0]->at(0)->c_str());
    }
  }
};
/* ################################################################### */
void OptionGroup::getDouble(double & out) {
  if (!isSet) {
    if (defaults.empty())
      out = 0.0;
    else
      out = atof(defaults.c_str());
  } else {
    if (args.empty() || args[0]->empty())
      out = 0.0;
    else {
      out = atof(args[0]->at(0)->c_str());
    }
  }
};
/* ################################################################### */
void OptionGroup::getString(std::string & out) {
  if (!isSet) {
    out = defaults;
  } else {
    if (args.empty() || args[0]->empty())
      out = "";
    else {
      out = *args[0]->at(0);
    }
  }
};
/* ################################################################### */
void OptionGroup::getInts(std::vector<int> & out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      StringsToInts(strings, out);
    }
  } else {
    if (!(args.empty() || args[0]->empty()))
      StringsToInts(args[0], &out);
  }
};
/* ################################################################### */
void OptionGroup::getLongs(std::vector<long> & out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      StringsToLongs(strings, out);
    }
  } else {
    if (!(args.empty() || args[0]->empty()))
      StringsToLongs(args[0], &out);
  }
};
/* ################################################################### */
void OptionGroup::getULongs(std::vector<unsigned long> & out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      StringsToULongs(strings, out);
    }
  } else {
    if (!(args.empty() || args[0]->empty()))
      StringsToULongs(args[0], &out);
  }
};
/* ################################################################### */
void OptionGroup::getFloats(std::vector<float> & out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      StringsToFloats(strings, out);
    }
  } else {
    if (!(args.empty() || args[0]->empty()))
      StringsToFloats(args[0], &out);
  }
};
/* ################################################################### */
void OptionGroup::getDoubles(std::vector<double> & out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      StringsToDoubles(strings, out);
    }
  } else {
    if (!(args.empty() || args[0]->empty()))
      StringsToDoubles(args[0], &out);
  }
};
/* ################################################################### */
void OptionGroup::getMultiInts(std::vector< std::vector<int> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      StringsToInts(strings, out[0]);
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);
      for(int i=0; i < n; ++i) {
        StringsToInts(args[i], &out[i]);
      }
    }
  }
};
/* ################################################################### */
void OptionGroup::getMultiLongs(std::vector< std::vector<long> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      StringsToLongs(strings, out[0]);
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);
      for(int i=0; i < n; ++i) {
        StringsToLongs(args[i], &out[i]);
      }
    }
  }
};
/* ################################################################### */
void OptionGroup::getMultiULongs(std::vector< std::vector<unsigned long> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      StringsToULongs(strings, out[0]);
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);
      for(int i=0; i < n; ++i) {
        StringsToULongs(args[i], &out[i]);
      }
    }
  }
};
/* ################################################################### */
void OptionGroup::getMultiFloats(std::vector< std::vector<float> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      StringsToFloats(strings, out[0]);
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);
      for(int i=0; i < n; ++i) {
        StringsToFloats(args[i], &out[i]);
      }
    }
  }
};
/* ################################################################### */
void OptionGroup::getMultiDoubles(std::vector< std::vector<double> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      StringsToDoubles(strings, out[0]);
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);
      for(int i=0; i < n; ++i) {
        StringsToDoubles(args[i], &out[i]);
      }
    }
  }
};
/* ################################################################### */
void OptionGroup::getMultiStrings(std::vector< std::vector<std::string> >& out) {
  if (!isSet) {
    if (!defaults.empty()) {
      std::vector< std::string > strings;
      SplitDelim(defaults, delim, strings);
      if (out.size() < 1) out.resize(1);
      out[0] = strings;
    }
  } else {
    if (!args.empty()) {
      int n = args.size();
      if (out.size() < n) out.resize(n);

      for(int i=0; i < n; ++i) {
				for(int j=0; j < args[i]->size(); ++j) 
					out[i].push_back( *args[i]->at(j) );
      }
    }
  }
};
/* ################################################################### */
class ezOptionParser {
public:
  // How to layout usage descriptions with the option flags.
  enum Layout { ALIGN, INTERLEAVE, STAGGER };

	~ezOptionParser();
	
	void add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1);
	void add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2);
	void add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2, const char * flag3);
	void add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2, const char * flag3, const char * flag4);
	bool exportFile(const char * filename, bool all=false);
	OptionGroup * get(const char * name);
  void getUsage(std::string & usage, int width=80, Layout layout=ALIGN);
	void getUsageDescriptions(std::string & usage, int width=80, Layout layout=STAGGER);
  bool gotExpected(std::vector<std::string> & badOptions);
  bool gotRequired(std::vector<std::string> & badOptions);
	bool importFile(const char * filename, char comment='#');
	int isSet(const char * name);
	int isSet(std::string & name);
	void parse(int argc, const char * argv[]);
	void prettyPrint(std::string & out); 
	void reset();
	void resetArgs();
		
  // General description in human language on what the user's tool does.
  // It's the first section to get printed in the full usage message.
  std::string overview;
  // A synopsis of command and options usage to show expected order of input arguments.
  // It's the second section to get printed in the full usage message.
  std::string syntax;
  // Example (third) section in usage message.
  std::string example;
  // Final section printed in usage message. For contact, copyrights, version info.
  std::string footer;
	// Map from an option to an Id of its parent group.
	std::map< std::string, int > optionGroupIds;
	// Unordered collection of the option groups.
	std::vector< OptionGroup* > groups;
	// Store unexpected args in input.
	std::vector< std::string* > unknownArgs;
	// List of args that occur left-most before first option flag.
	std::vector< std::string* > firstArgs;
	// List of args that occur after last right-most option flag and its args.
	std::vector< std::string* > lastArgs;
};
/* ################################################################### */
ezOptionParser::~ezOptionParser() {
	reset();
}
/* ################################################################### */
void ezOptionParser::reset() {
	int i;
	for(i=0; i < groups.size(); ++i)
		delete groups[i];
	groups.clear();
	
	for(i=0; i < unknownArgs.size(); ++i)
		delete unknownArgs[i];
	unknownArgs.clear();
	
	for(i=0; i < firstArgs.size(); ++i)
		delete firstArgs[i];
	firstArgs.clear();

	for(i=0; i < lastArgs.size(); ++i)
		delete lastArgs[i];
	lastArgs.clear();
	
	optionGroupIds.clear();
};
/* ################################################################### */
void ezOptionParser::resetArgs() {
	int i;
	for(i=0; i < groups.size(); ++i)
		groups[i]->clearArgs();
	
	for(i=0; i < unknownArgs.size(); ++i)
		delete unknownArgs[i];
	unknownArgs.clear();
	
	for(i=0; i < firstArgs.size(); ++i)
		delete firstArgs[i];
	firstArgs.clear();

	for(i=0; i < lastArgs.size(); ++i)
		delete lastArgs[i];
	lastArgs.clear();	
};
/* ################################################################### */
void ezOptionParser::add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1) {
	int id = this->groups.size();
	OptionGroup * g = new OptionGroup;
  g->defaults = defaults;
  g->isRequired = required;
	g->expectArgs = expectArgs;
	g->delim = delim;
	g->isSet = 0;
	g->help = help;
	std::string *f1 = new std::string(flag1);
	g->flags.push_back( f1 );
	this->optionGroupIds[flag1] = id;

	this->groups.push_back(g);
};
/* ################################################################### */
void ezOptionParser::add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2) {
	int id = this->groups.size();
	OptionGroup * g = new OptionGroup;
  g->defaults = defaults;
  g->isRequired = required;
	g->expectArgs = expectArgs;
	g->delim = delim;
	g->isSet = 0;
	g->help = help;
	std::string *f1 = new std::string(flag1);
	g->flags.push_back( f1 );
	std::string *f2 = new std::string(flag2);
	g->flags.push_back( f2 );
	this->optionGroupIds[flag1] = id;
	this->optionGroupIds[flag2] = id;

	this->groups.push_back(g);
};
/* ################################################################### */
void ezOptionParser::add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2, const char * flag3) {
	int id = this->groups.size();
	OptionGroup * g = new OptionGroup;
  g->defaults = defaults;
  g->isRequired = required;
	g->expectArgs = expectArgs;
	g->delim = delim;
	g->isSet = 0;
	g->help = help;
	std::string *f1 = new std::string(flag1);
	g->flags.push_back( f1 );
	std::string *f2 = new std::string(flag2);
	g->flags.push_back( f2 );
	std::string *f3 = new std::string(flag3);
	g->flags.push_back( f3 );
	this->optionGroupIds[flag1] = id;
	this->optionGroupIds[flag2] = id;
	this->optionGroupIds[flag3] = id;

	this->groups.push_back(g);
};
/* ################################################################### */
void ezOptionParser::add(const char * defaults, bool required, int expectArgs, char delim, const char * help, const char * flag1, const char * flag2, const char * flag3, const char * flag4) {
	int id = this->groups.size();
	OptionGroup * g = new OptionGroup;
  g->defaults = defaults;
  g->isRequired = required;
	g->expectArgs = expectArgs;
	g->delim = delim;
	g->isSet = 0;
	g->help = help;
	std::string *f1 = new std::string(flag1);
	g->flags.push_back( f1 );
	std::string *f2 = new std::string(flag2);
	g->flags.push_back( f2 );
	std::string *f3 = new std::string(flag3);
	g->flags.push_back( f3 );
	std::string *f4 = new std::string(flag4);
	g->flags.push_back( f4 );
	this->optionGroupIds[flag1] = id;
	this->optionGroupIds[flag2] = id;
	this->optionGroupIds[flag3] = id;
	this->optionGroupIds[flag4] = id;

	this->groups.push_back(g);
};
/* ################################################################### */
bool ezOptionParser::exportFile(const char * filename, bool all) {
	int i;
	std::string out;
	bool quote;
	
	// Export the first args, except the program name, so start from 1.
	for(i=1; i < firstArgs.size(); ++i) {
		quote = ((firstArgs[i]->find_first_of(" \t") != std::string::npos) && (firstArgs[i]->find_first_of("\'\"") == std::string::npos));

		if (quote)
			out.append("\"");

		out.append(*firstArgs[i]);
		if (quote)
			out.append("\"");

		out.append(" ");
	}
	
	if (firstArgs.size() > 1)
		out.append("\n");
	
	std::vector<std::string* > stringPtrs(groups.size());
	int j,m;
	int n = groups.size();
	for(i=0; i < n; ++i) {
		stringPtrs[i] = groups[i]->flags[0];
	}

	OptionGroup *g;
	// Sort first flag of each group with other groups.
	std::sort(stringPtrs.begin(), stringPtrs.end(), CmpOptStringPtr);
	for(i=0; i < n; ++i) {
		g = get(stringPtrs[i]->c_str());
		if (g->isSet || all) {
			if (!g->isSet || g->args.empty()) {
				if (!g->defaults.empty()) {
					out.append(*stringPtrs[i]);
					out.append(" ");
					quote = ((g->defaults.find_first_of(" \t") != std::string::npos) && (g->defaults.find_first_of("\'\"") == std::string::npos));
					if (quote)
						out.append("\"");
						
					out.append(g->defaults);
					if (quote)
						out.append("\"");
						
					out.append("\n");
				}
			} else {
				int n = g->args.size();
				for(int j=0; j < n; ++j) {
					out.append(*stringPtrs[i]);
					out.append(" ");
					m = g->args[j]->size();
					
					for(int k=0; k < m; ++k) {
						quote = ( (*g->args[j]->at(k)).find_first_of(" \t") != std::string::npos );
						if (quote)
							out.append("\"");
							
						out.append(*g->args[j]->at(k));
						if (quote)
							out.append("\"");
							
						if ((g->delim) && ((k+1) != m))
							out.append(1,g->delim);
					}
					out.append("\n");
				}
			}
		}
	}
	
	// Export the last args.
	for(i=0; i < lastArgs.size(); ++i) {
		quote = ( lastArgs[i]->find_first_of(" \t") != std::string::npos );
		if (quote)
			out.append("\"");

		out.append(*lastArgs[i]);
		if (quote)
			out.append("\"");

		out.append(" ");
	}
	
	std::ofstream file(filename);
  if (!file.is_open())
		return false;

	file << out;
	file.close();
	
  return true;
};
/* ################################################################### */
// Does not overwrite current options.
// Returns true if file was read successfully.
// So if this is used before parsing CLI, then option values will reflect
// this file, but if used after parsing CLI, then values will contain
// both CLI values and file's values.
//
// Comment lines are allowed if prefixed with #.
// Strings should be quoted as usual.
bool ezOptionParser::importFile(const char * filename, char comment) {
	std::ifstream file (filename, std::ios::in | std::ios::ate);
  if (!file.is_open())
    return false;

  // Read entire file contents.
  std::ifstream::pos_type size = file.tellg();
  char * memblock = new char[(int)size+1]; // Add one for end of string.
  file.seekg (0, std::ios::beg);
  file.read (memblock, size);
	memblock[size] = '\0';
  file.close();
	
  // Find comment lines.
  std::list<std::string*> lines;
  std::string memblockstring(memblock);
  delete[] memblock;
  SplitDelim(memblockstring, '\n', lines);
  int i,j,n;
  std::list<std::string*>::iterator iter;
  std::vector<int> sq, dq; // Single and double quote indices.
	std::vector<int>::iterator lo; // For searching quote indices.
  size_t pos;
	const char *str;
	std::string *line;
  // Find all single and double quotes to correctly handle comment tokens.
  for(iter=lines.begin(); iter != lines.end(); ++iter) {
    line = *iter;
    str = line->c_str();
    n = line->size();
    sq.clear();
    dq.clear();
    if (n) {
      // If first char is comment, then erase line and continue.
      pos = line->find_first_not_of(" \t\r");
      if ((pos==std::string::npos) || (line->at(pos)==comment)) {
				line->erase();
        continue; 
      } else {
				// Erase whitespace prefix.
				line->erase(0,pos);
        n = line->size();
			}

      if (line->at(0)=='"')
        dq.push_back(0);
    
      if (line->at(0)=='\'')
        sq.push_back(0);
    } else { // Empty line.
      continue; 
    }

    for(i=1; i < n; ++i) {
      if ( (str[i]=='"') && (str[i-1]!='\\') )
        dq.push_back(i);
      else if ( (str[i]=='\'') && (str[i-1]!='\\') )
        sq.push_back(i);
    }
		// Scan for comments, and when found, check bounds of quotes.
		// Start with second char because already checked first char.
		for(i=1; i < n; ++i) {
			if ( (line->at(i)==comment) && (line->at(i-1)!='\\') ) {
				// If within open/close quote pair, then not real comment.
				if (sq.size()) {
					lo = std::lower_bound(sq.begin(), sq.end(), i);
					// All start of strings will be even indices, closing quotes is odd indices.
					j = (int)(lo-sq.begin());
					if ( (j % 2) == 0) { // Even implies comment char not in quote pair.
						// Erase from comment char to end of line.
						line->erase(i);
						break;
					}
				}	else if (dq.size()) {
					// Repeat tests for double quotes.
					lo = std::lower_bound(dq.begin(), dq.end(), i);
					j = (int)(lo-dq.begin());
					if ( (j % 2) == 0) {
						line->erase(i);
						break;
					}
				} else {				
					// Not in quotes.
					line->erase(i);
					break;
				}
			}
		}	
  }
			
	std::string cmd;
  // Convert list to string without newlines to simulate commandline.
  for(iter=lines.begin(); iter != lines.end(); ++iter) {
		if (! (*iter)->empty()) {
			cmd.append(**iter);
			cmd.append(" ");
		}
	}
  
  // Now parse as if from command line. 
  int argc=0;
  char** argv = CommandLineToArgvA((char*)cmd.c_str(), &argc);
  
  // Parse.
  parse(argc, (const char**)argv);
  if (argv) free(argv);
  for(iter=lines.begin(); iter != lines.end(); ++iter)
    delete *iter;
    
  return true;
};
/* ################################################################### */
int ezOptionParser::isSet(const char * name) {
	std::string sname(name);
	
	if (this->optionGroupIds.count(sname)) {
		return this->groups[this->optionGroupIds[sname]]->isSet;
	}
	
	return 0;
};
/* ################################################################### */
int ezOptionParser::isSet(std::string & name) {
	if (this->optionGroupIds.count(name)) {
		return this->groups[this->optionGroupIds[name]]->isSet;
	}
	
	return 0;
};
/* ################################################################### */
OptionGroup * ezOptionParser::get(const char * name) {
	if (optionGroupIds.count(name)) {
		return groups[optionGroupIds[name]];
	}
	
	return 0;
};
/* ################################################################### */
void ezOptionParser::getUsage(std::string & usage, int width, Layout layout) {
  usage.append(overview);
  usage.append("\n\n");
  usage.append("USAGE: ");
  usage.append(syntax);
  usage.append("\n\nOPTIONS:\n\n");
  getUsageDescriptions(usage, width, layout);
  
  if (!example.empty()) {
    usage.append("EXAMPLES:\n\n");
    usage.append(example);
  }
  
  if (!footer.empty()) {
    usage.append(footer);
  }
};
/* ################################################################### */
// Creates 2 column formatted help descriptions for each option flag.
void ezOptionParser::getUsageDescriptions(std::string & usage, int width, Layout layout) {
	// Sort each flag list amongst each group.
	int i;
	// Store index of flag groups before sort for easy lookup later.
	std::map<std::string*, int> stringPtrToIndexMap;
	std::vector<std::string* > stringPtrs(groups.size());
	
	for(i=0; i < groups.size(); ++i) {
		std::sort(groups[i]->flags.begin(), groups[i]->flags.end(), CmpOptStringPtr);
		stringPtrToIndexMap[groups[i]->flags[0]] = i;
		stringPtrs[i] = groups[i]->flags[0];
	}
	
	size_t j, k, n;
	std::string opts;
	std::vector<std::string> sortedOpts;
	// Sort first flag of each group with other groups.
	std::sort(stringPtrs.begin(), stringPtrs.end(), CmpOptStringPtr);
 	for(i=0; i < groups.size(); ++i) {
		//printf("DEBUG:%d: %d %d %s\n", __LINE__, i, stringPtrToIndexMap[stringPtrs[i]], stringPtrs[i]->c_str());
		k = stringPtrToIndexMap[stringPtrs[i]];
		opts.clear();
		for(j=0; j < groups[k]->flags.size()-1; ++j) {
			opts.append(*groups[k]->flags[j]);
			opts.append(", ");
			
			if (opts.size() > width)
				opts.append("\n");
		}
		// The last flag. No need to append comma anymore.
		opts.append( *groups[k]->flags[j] );
		
		if (groups[k]->expectArgs) {
			opts.append(" ARG");
			
			if (groups[k]->delim) {
				opts.append("1[");
				opts.append(1, groups[k]->delim);
				opts.append("ARGn]");
      }
		}
			
		sortedOpts.push_back(opts);
	}	
	
	// Each option group will use this to build multiline help description.
  std::list<std::string*> desc; 
  // Number of whitespaces from start of line to description (interleave layout) or
  // gap between flag names and description (align, stagger layouts).
	int gutter = 3;
	
  // Find longest opt flag string to set column start for help usage descriptions.
  int maxlen=0;
  if (layout == ALIGN) {
    for(i=0; i < groups.size(); ++i) {
      if (maxlen < sortedOpts[i].size())
        maxlen = sortedOpts[i].size();
    }	
  }
  
  // The amount of space remaining on a line for help text after flags.
  int helpwidth;
  std::list<std::string*>::iterator cIter, insertionIter;
  size_t pos;
  for(i=0; i < groups.size(); ++i) {
    k = stringPtrToIndexMap[stringPtrs[i]];
    
    if (layout == STAGGER)
      maxlen = sortedOpts[i].size();
      
    int pad = gutter + maxlen;
    helpwidth = width - pad;        
    
    // All the following split-fu could be optimized by just using substring (offset, length) tuples, but just to get it done, we'll do some not-too expensive string copying.
    SplitDelim(groups[k]->help, '\n', desc);
    // Split lines longer than allowable help width.
    for(insertionIter=desc.begin(), cIter=insertionIter++; 
        cIter != desc.end(); 
        cIter=insertionIter++) {
      if ((*cIter)->size() > helpwidth) {
        // Get pointer to next string to insert new strings before it.
        std::string *rem = *cIter;
        // Remove this line and add back in pieces.
        desc.erase(cIter);
        // Loop until remaining string is short enough.
        while (rem->size() > helpwidth) {
          // Find whitespace to split before helpwidth.
          if (rem->at(helpwidth) == ' ') {
            // If word ends exactly at helpwidth, then split after it.
            pos = helpwidth;
          } else {
            // Otherwise, split occurs midword, so find whitespace before this word.
            pos = rem->rfind(" ", helpwidth);
          }
          // Insert split string.
          desc.insert(insertionIter, new std::string(*rem, 0, pos));
          // Now skip any whitespace to start new line.
          pos = rem->find_first_not_of(' ', pos);
          rem->erase(0, pos);
        }
        
        if (rem->size())
          desc.insert(insertionIter, rem);
        else
          delete rem;
      }
    }

    usage.append(sortedOpts[i]);
    if (layout != INTERLEAVE)
      // Add whitespace between option names and description.
      usage.append(pad - sortedOpts[i].size(), ' ');
    else {
      usage.append("\n\n");
      usage.append(gutter, ' ');
    }
    
    // First line already padded above (before calling SplitDelim) after option flag names.
    cIter = desc.begin();
    usage.append(**cIter);
    usage.append("\n");
    // Now inject the pad for each line.
    for(++cIter; cIter != desc.end(); ++cIter) {			
      usage.append(pad, ' ');
      usage.append(**cIter);
      usage.append("\n");
    }

    usage.append("\n");
    
    if (desc.size()) {
      for(cIter=desc.begin(); cIter != desc.end(); ++cIter)		
        delete *cIter;
      
      desc.clear();
    }
  }	
};
/* ################################################################### */
bool ezOptionParser::gotExpected(std::vector<std::string> & badOptions) {
  int i,j;
  
  for(i=0; i < groups.size(); ++i) {
		OptionGroup *g = groups[i];
    // If was set, ensure number of args is correct.
    if (g->isSet) {
      if ((g->expectArgs != 0) && g->args.empty()) {
        badOptions.push_back(*g->flags[0]);
        continue;
      }
      
      for(j=0; j < g->args.size(); ++j) {
        if ((g->expectArgs != -1) && (g->expectArgs != g->args[j]->size()))
          badOptions.push_back(*g->flags[0]);
      }
		}
  }

  return badOptions.empty();
};
/* ################################################################### */
bool ezOptionParser::gotRequired(std::vector<std::string> & badOptions) {
  int i;
  
  for(i=0; i < groups.size(); ++i) {
		OptionGroup *g = groups[i];
    // Simple case when required but user never set it.
    if (g->isRequired && (!g->isSet)) {
      badOptions.push_back(*g->flags[0]);
      continue;
    }
  }

  return badOptions.empty();
};
/* ################################################################### */
void ezOptionParser::parse(int argc, const char * argv[]) {
	if (argc < 1) return;
	
	/*
	std::map<std::string,int>::iterator it;
	for ( it=optionGroupIds.begin() ; it != optionGroupIds.end(); it++ )
		std::cout << (*it).first << " => " << (*it).second << std::endl;
	*/
	
	int found=0, i, k, firstOptIndex=0, lastOptIndex=0;
	std::string s;
	OptionGroup *g;
	
	for(i=0; i < argc; ++i) {
		s = argv[i];

		if (optionGroupIds.count(s))
			break;
	}
	
	firstOptIndex = i;

	if (firstOptIndex == argc) {
    // No flags encountered, so set last args.
    this->firstArgs.push_back(new std::string(argv[0]));

    for(k=1; k < argc; ++k)
      this->lastArgs.push_back(new std::string(argv[k]));
      
    return;
  }

	// Store initial args before opts appear.
	for(k=0; k < i; ++k) {
		this->firstArgs.push_back(new std::string(argv[k]));
	}
	
	for(; i < argc; ++i) {
		s = argv[i];
		
		if (optionGroupIds.count(s)) {
			k = optionGroupIds[s];
			g = groups[k];
			g->isSet = 1;
			if (g->expectArgs) {
				// Read ahead to get args.
				++i;
				if (i >= argc) return;
				g->args.push_back(new std::vector<std::string*>);
				SplitDelim(argv[i], g->delim, g->args.back());
			}
			lastOptIndex = i;
		}
	}
	
	// Scan for unknown opts/arguments.
	for(i=firstOptIndex; i <= lastOptIndex; ++i) {
		s = argv[i];
		
		if (optionGroupIds.count(s)) {
			k = optionGroupIds[s];
			g = groups[k];
			g->isSet = 1;
			if (g->expectArgs) {
				// Read ahead for args.
				++i;
			}
		} else {
			unknownArgs.push_back(new std::string(argv[i]));
		}
	}
	
	if ( lastOptIndex >= (argc-1) ) return;
	
	// Store final args without flags.
	for(k=lastOptIndex + 1; k < argc; ++k) {
		this->lastArgs.push_back(new std::string(argv[k]));
	}
};
/* ################################################################### */
void ezOptionParser::prettyPrint(std::string & out) {
	char tmp[256];
	int i,j,k;
	
	out += "First Args:\n";
	for(i=0; i < firstArgs.size(); ++i) {
		sprintf(tmp, "%d: %s\n", i+1, firstArgs[i]->c_str());
		out += tmp;
	}

	// Sort the option flag names.
	int n = groups.size();
	std::vector<std::string* > stringPtrs(n);
	for(i=0; i < n; ++i) {
		stringPtrs[i] = groups[i]->flags[0];
	}

	// Sort first flag of each group with other groups.
	std::sort(stringPtrs.begin(), stringPtrs.end(), CmpOptStringPtr);

	out += "\nOptions:\n";
	OptionGroup *g;
	for(i=0; i < n; ++i) {						
		g = get(stringPtrs[i]->c_str());
    out += "\n";
		// The flag names:
		for(j=0; j < g->flags.size()-1; ++j) {
			sprintf(tmp, "%s, ", g->flags[j]->c_str());
			out += tmp;
		}
		sprintf(tmp, "%s:\n", g->flags.back()->c_str());
		out += tmp;

		if (g->isSet) {
			if (g->expectArgs) {
				if (g->args.empty()) {
					sprintf(tmp, "%s (default)\n", g->defaults.c_str());
					out += tmp;
				} else {
					for(k=0; k < g->args.size(); ++k) {						
						for(j=0; j < g->args[k]->size()-1; ++j) {
							sprintf(tmp, "%s%c", g->args[k]->at(j)->c_str(), g->delim);
							out += tmp;
						}
						sprintf(tmp, "%s\n", g->args[k]->back()->c_str(), g->delim);
						out += tmp;
					}
				}
			} else { // Set but no args expected.
				sprintf(tmp, "Set\n");
				out += tmp;
			}
		} else {
			sprintf(tmp, "Not set\n");
			out += tmp;
		}
	}
	
	out += "\nLast Args:\n";
	for(i=0; i < lastArgs.size(); ++i) {
		sprintf(tmp, "%d: %s\n", i+1, lastArgs[i]->c_str());
		out += tmp;
	}
	
	out += "\nUnknown Args:\n";
	for(i=0; i < unknownArgs.size(); ++i) {
		sprintf(tmp, "%d: %s\n", i+1, unknownArgs[i]->c_str());
		out += tmp;
	}
};
}
/* ################################################################### */
#endif /* EZ_OPTION_PARSER_H */