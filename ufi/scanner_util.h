#if !defined(__scanner_util__h)
#define __scanner_util__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/scanner_util.h,v $

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


//! Utility methods for scanner implementatation.
/**
 * Utility methods to assist in scanner implementation.
 */
namespace scanner_util{

/*
 * command-line arguments as single string (excluding argv[0]):
 */
void argument_string(std::string& s, int argc, const char* argv[]);


// load all available contents of a std::istream to a std::string:
// helper method (based on a utility from boost::regex examples)    
void load_from_stream(std::string& s, std::istream& is);


// load contents from a stream to a string:
//   contents are delimited by start-delim, end-delim
//   leading whitespace is ignored
void load_group_from_stream(std::string& s, std::istream& is, char ldelim, char rdelim);


// commonly used regular expression strings:
const std::string 
  identifier_regex_str("[_$\\w]+[_$\\w\\d]*"),
  integer_regex_str("[+\\-]?\\d*"),
  real_regex_str("[+\\-]?\\d*(\\.\\d+)?[eEdD][+\\-]?\\d+"),
  complex_regex_str(real_regex_str + "\\s*\\+\\s*" + real_regex_str + "j"),
  number_regex_str(real_regex_str + "|" + complex_regex_str);
      
} // namespace scanner_util


#if 1 // ---------------------- moved from numericalFunctor.h: 08.08.2009 --------------------------------
// Specialization of hash<T> functor required in order to use hash_map with std::string as a key:
namespace _STL_EXT_NAMESPACE_
{
  #if defined(__PGI)
  #else
  using std::size_t;
  #endif

  inline size_t hash_std_string(const std::string& s)
  {
    const char* __s(s.c_str());

    unsigned long __h = 0; 
    for ( ; *__s; ++__s)
      __h = 5*__h + *__s;

    return size_t(__h);
  }

  template<> struct hash<std::string>
  {
    size_t operator()(const std::string& __s) const { return hash_std_string(__s); }
  };

  template<> struct hash<const std::string>
  {
    size_t operator()(const std::string& __s) const { return hash_std_string(__s); }
  };

} // namespace _STL_EXT_NAMESPACE_
#endif


#include <scanner_util_module.h>

#endif // __scanner_util__h
