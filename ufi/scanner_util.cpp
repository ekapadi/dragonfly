// $Source: /usr/data0/leipzig_work/tmat_cvs/src/scanner_util.cpp,v $

/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2010  Kort Travis                                         */
/*                                                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */

#if defined(__USE_PTHREAD)
#include <pthread.h>
#endif

#if defined(__USE_MPI)
#include <mpi.h>
#endif

#include <typeinfo>

#include <iostream>
#include <iomanip>

#include <complex>
#include <string>
#include <sstream>
#include <list>

// #include <functional>
// #include <algorithm>
// #include <numeric>
#include <unordered_map>

#include "cross_platform.h"
#include "scanner_util.h"

namespace scanner_util{

/*
 * command-line arguments as single string (excluding argv[0]):
 */
void argument_string(std::string& s, int argc, const char* argv[])
{
  s.clear();
  for(int n = 1; n < argc; ++n){    
    if (n != 1)
      s += ' ';
    s += argv[n];
  }
}

// load all available contents of a std::istream to a std::string:
// helper method (based on a utility from boost::regex examples)    
void load_from_stream(std::string& s, std::istream& is)
{
  s.erase();
  s.reserve(is.rdbuf()->in_avail()); // if this doesn't work, gradually increase capacity (below)
  char c;
  while(is.get(c))
  {
     if(s.capacity() == s.size()) 
        s.reserve(s.capacity() * 3); 
     s.append(1, c);
  }
}

// load contents from a stream to a string:
//   contents are delimited by start-delim, end-delim
//   leading whitespace is ignored
void load_group_from_stream(std::string& s, std::istream& is, char ldelim, char rdelim)
{ 
  using std::isspace;
  // allow a failed read without destroying stream contents
  std::streampos savePos(is.tellg());
  
  s.erase();
  #if 0  // ------- assume that a group is relatively small in comparison to the entire input file: ------------
  s.reserve(is.rdbuf()->in_avail()); // if this doesn't work, gradually increase capacity (below)
  #endif // -------
  char c;
  
  while(is.get(c) && isspace(c));
  int nesting_level(0);
  if (c == ldelim){
    do{
      if (c == ldelim)
        ++nesting_level;
      else
      if (c == rdelim)
        --nesting_level;
        
      if(s.capacity() == s.size()) 
         s.reserve(s.capacity() * 3); 
      s.append(1, c); 
      
      if (nesting_level == 0)
        break;
    } while(is.get(c));
  }  
  else{
   is.setstate(std::ios_base::failbit); // return status in the stream
   is.seekg(savePos, std::ios::beg);    // restore entry position
  }
}

} // namespace scanner_util
