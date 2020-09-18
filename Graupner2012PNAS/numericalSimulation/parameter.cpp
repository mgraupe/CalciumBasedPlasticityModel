/*************************************************************************
file:         parameter.cpp
author:       Gernot Schaller
mail:         schaller@theory.phy.tu-dresden.de
version:      0.5
Time-stamp:   <03/11/25 16:40:52 schaller>
---------------------------------------------------------------------
description:  implementation of class Parameter (parameter.h)
This class implements all the basic features of working with
parameters stored in files etc.
**************************************************************************/  

#include "parameter.hpp"  // parameter management

#include <sstream>    // conversion string <-> <template>
#include <fstream>    // save/load data files
#include <iostream>   // error messages
using namespace std;

Parameter::Parameter(){
}

Parameter::Parameter(const string& dateiname){
  init(dateiname);
}

/* gets parameter filename and stores the read parameters within a
   map using a global string format. The object which is initialized
   using the parameter class should know which of the builtin 
   data types it is expecting and use the corresponding function.
   The parameter file must contain the parameters in the form
   'name' = 'value', empty lines and comments #... will be ignored.
   one line must not contain more than 1 parameter assignment */

void Parameter::init(const string& dateiname){
  string par_name="", par_value="", linestr="";
  unsigned int linectr=0;
  ifstream inputfile(dateiname.c_str());
  if (!inputfile.good()){
    cerr << "File '" << dateiname << "' not found!" << endl;
    cerr << "ERROR found in method Parameter::init(*), aborting ..." << endl;
    exit(100);
  }
  else while (inputfile.good()){
    // read the parameter file and store the data in the map
    ++linectr;
    par_name="";
    par_value="";
    linestr="";
    getline(inputfile, linestr);
    replace(linestr.begin(), linestr.end(), '\t', ' '); 
    unsigned int n=0;
    while (n<linestr.size() && linestr[n]!='#')
      ++n; // ignore comments at end of line
    /* the following routine copies the parameter name found in
       line.substr into par_name and the value into par_value */
    linestr=linestr.substr(0,n);
    if (linestr.find("=")!=linestr.npos){
      /* Then the sign '=' has been found in linestr. */
      extract_param(par_name, par_value, linestr);
      if ((par_name=="") || (par_value=="")){
	// error-handling
	cerr << "Invalid data format at line " << linectr 
	     << " in file '" << dateiname << "'!\n"
	     << "(parameter name='" << par_name 
	     << "' parameter value='" << par_value << "')"
	     << " ignoring ..." << endl;
      }
      else{
	// save parameters in map 
	par_map[upcase(par_name)]=par_value;
      }
    }
    else{
      n=0;
      while (n<linestr.size() && linestr[n]==' ')
	++n; // determine if line contains characters
      if (n<linestr.size()){
	cerr << "Invalid data format at line " << linectr 
	     << " in file '" << dateiname << "'!\n"
	     << "('" << linestr << "' does not contain '=')"
	     << " ignoring ..." << endl;
      }
    }
  } 
  inputfile.close();
}

Parameter::~Parameter(){
}

void Parameter::extract_param(string& name, string& value, 
			      const string& line) const{
  /* this routine looks for a '=' char in the string line and assigns the
     text before to name and the text after to value after stripping tabs and ' ' */
  unsigned int eqpos=line.find("=");
  if (eqpos!=line.npos){
    unsigned int startpos=eqpos-1;
    unsigned int endpos=eqpos-1;
    while (line[endpos]==' ' && endpos>0) --endpos;
    startpos=line.rfind(" ", endpos);
    if (startpos==line.npos) name=line.substr(0, endpos+1);
    else name=line.substr(startpos+1, endpos-startpos);
    startpos=eqpos+1;
    while (line[startpos]==' ' && startpos<line.size())
      startpos++;
    endpos=line.find(" ", startpos);
    if (endpos==line.npos)
      value=line.substr(startpos, line.size()-startpos);
    else
      value=line.substr(startpos, endpos-startpos);
    // the '=' char should not be assigned to value
  }
}

string Parameter::upcase(const string& small) const{
  string dummy=small;
  unsigned int endpos=strlen(small.c_str());
  for (unsigned int i=0; i<=endpos; i++)
    dummy[i]=toupper(small[i]);
  return dummy;
}

/* It would be nice to be able to use command line parameters to overwrite 
   the values found in the parameter file. Consequently, the key word 'name' 
   in the command line option "-'name'='value'" must match the one given in  
   parameter file. If not, another new parameter will be added in the map.
   Command line options not containing an assignment will be ignored.  */

void Parameter::set_cmdlineopts(const unsigned int& argc, char** argv){
  string line="", par_name="", par_value="";
  unsigned int pos1=0, pos2=0;
  for (unsigned int i=1; i<argc; i++){
    line=argv[i];
    pos2=line.size();
    pos1=pos2;
    while (pos1>0 && line[pos1]!='-')
      pos1--;
    if (line.find("=")!=line.npos){
      /* Then the sign '=' has been found in line. */
      extract_param(par_name, par_value, line.substr(pos1+1, pos2-pos1));
      if (par_name=="" || par_value==""){
	// error-handling
	cerr << "Invalid data format in " << i
	     << ". command line option: '" << par_name 
	     << "'='" << par_value << "'" 
	     << " ignoring ..." << endl;
      }
      else{
	/* save parameters in map */
	par_map[upcase(par_name)]=par_value;
      }
    }
  }
}

/* In some cases it might be useful to overwrite or add parameters while the
   program is running. This should overwrite both parameter file values and 
   command line options. */
void Parameter::set_param(const string& par_name, const string& par_value){
  if (par_name=="" || par_value==""){
    // error-handling
    cerr << "Invalid data format in attempted assignment: '" 
	 << par_name << "'='" << par_value << "'" << endl;
    cerr << "FATAL ERROR occured in method Parameter::set_param(*,*)" 
	 << " aborting ..." << endl;
    exit(100);
  }
  else{
    par_map[upcase(par_name)]=par_value;
  }
}

/* The following functions return the value saved as a string in the
   map in the desired - converted - form to the calling object. */
template<class T> T Parameter::get_param(const string& name) const{
  string NAME=upcase(name);
  T value;
  map<string,string>::const_iterator pos=par_map.find(NAME);
  if (pos!=par_map.end()){
    /* then the key 'NAME' has been found in the map and the parameter 
       'value' should be assigned the corresponding long int conversion */
    if (pos->second!="" && NAME!=""){
      stringstream str;
      str << pos->second;
      str >> value;
      return value;
    }
    else{
      cerr << "Parameter '" << NAME << "' has value '"
	   << pos->second << "' which cannot be converted!" << endl;
      cerr << "ERROR occured in method Parameter::get_param(*)"
	   << " aborting ..." << endl;
      exit(100);
    }
  }
  else{
    cerr << "Parameter '" << NAME << "' has not been defined!" << endl
	 << "ERROR occured in method Parameter::get_param<T>(*)" 
	 << " aborting ..." << endl;
    exit(100);
  }
}

// bool conversion needs a special function
template<> bool Parameter::get_param<bool>(const string& name) const{
  string NAME=upcase(name);
  map<string,string>::const_iterator pos=par_map.find(NAME);
  if (pos!=par_map.end()){
    /* then the key 'NAME' has been found in the map and the parameter 
       'value' should be assigned the corresponding long int conversion */
    if (NAME!="" && (pos->second=="TRUE" || pos->second=="1")){
      return true;
    }
    else if (NAME!="" && (pos->second=="FALSE" || pos->second=="0")){
      return false;
    }
    else{
      cerr << "Boolean parameter '" << NAME << "' has value '"
	   << pos->second << "' which cannot be converted!"
	   << "ERROR occured in method Parameter::get_param<bool>(*)"
	   << " aborting ..." << endl;
      exit(100);
    }
  }
  else{
    cerr << "Boolean parameter '" << NAME << "' has not been defined!"
	 << "ERROR occured in method Parameter::get_param<bool>(*)"
	 << " aborting ..." << endl;
    exit(100);
  }
}

bool Parameter::if_param(const string& name) const{
  string NAME=upcase(name);
  if (par_map.find(NAME)!=par_map.end())
    return true;
  else
    return false;
}

void Parameter::write_status(const string& dateiname, 
			     const string& remark="") const{
  ofstream file(dateiname.c_str());
  file << "# -------------- parameter backup file --------------" << endl
       << "# " << remark << endl
       << "# file contains " << par_map.size() << " parameters..." << endl
       << "# ---------------------------------------------------" << endl;
  for(map<string,string>::const_iterator pos=par_map.begin();
      pos!=par_map.end(); ++pos)
    file << pos->first << " = " << pos->second << endl;
  file.close();
}


/* explicit instanciation: The template get_param is creating object
   code for all these instanciations. */
template int           Parameter::get_param<int>(const string&) const;
template unsigned int  Parameter::get_param<unsigned int>(const string&) const;
template unsigned long Parameter::get_param<unsigned long>(const string&) const;
template long          Parameter::get_param<long>(const string&) const;
template double        Parameter::get_param<double>(const string&) const;
template float         Parameter::get_param<float>(const string&) const;
template string        Parameter::get_param<string>(const string&) const;
template char          Parameter::get_param<char>(const string&) const;
