/*************************************************************************
file:         parameter.h
author:       Gernot Schaller
mail:         schaller@theory.phy.tu-dresden.de
version:      0.5
Time-stamp:   <03/11/25 16:41:20 schaller>
---------------------------------------------------------------------- 
description:  implementation of class Parameter (parameter.h)
This class implements all the basic features of working with
parameters stored in files etc.
**************************************************************************/  

#ifndef __PARAMETER__
#define __PARAMETER__

#include <map>
#include <string>
#include <sstream>
#include <fstream>    // save/load data files
#include <iostream>   // error messages
#include <algorithm>
#include <cstring>


using namespace std;

class Parameter{
public:
  /* gets parameter filename and stores the read parameters within a
     map using a global string format. The object which is initialized
     using the parameter class should know which of the builtin 
     data types it is expecting and use the corresponding function.
     The parameter file must contain the parameters in the form
     'name' = 'value', empty lines and comments beginning with '#' will be 
     ignored. All non-comment-lines must contain an assignment sign ('=').
     Important: Note that internally the program uses capital letters for the
     keywords 'name', i.e. "volume" and "VoLUme" are identical keywords. */
  Parameter();
  Parameter(const string&);
  ~Parameter();

  void init(const string&);

  /* command line parameters should overwrite the values found in the
     parameter file or add new ones. Consequently, the key word given by 
     'name'='value' must match the one given in the parameter file */
  void set_cmdlineopts(const unsigned int&, char**);

  /* New parameters can also be set in the program. This should overwrite
     both parameter file values and command line options. However, this 
     depends on which of the set_... routines is called first. */
  void set_param(const string&, const string&);

  /* The following functions return the value saved under "string" in the
     map in the desired - converted - form to the calling object. */
  template <class T>
    T get_param(const string&) const;

  /* one might be interested in how many parameters have been saved in the
     map. */
  unsigned int count_param(void) const{
      return par_map.size();
  }

  // one can check if a parameter has already been stored in par_map
  bool if_param(const string&) const;

  // the parameters in par_map can be written in a separate file
  // of course without comments
  void write_status(const string& file, const string&) const;

private: 

  /* The map should not differentiate between upper- and lowercase letters
     in the keywords. One can therefore also save everything in uppercase
     letters. */
  string upcase(const string&) const;

  void extract_param(string&, string&, const string&) const;

  map<string,string> par_map;
  unsigned int number;
 
};

#endif

