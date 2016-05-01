
/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2016  Kort Travis                                         */
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

/**
 * @file: commUtil_stub.cpp
 *   *Temporary* reduction of "commUtil" dependency prior to *upgrade* of "commUtil" classes to use std::iostream-derived binary I/O.
 */

#include <assert.h>
#include <stdexcept>

#include <iostream>
#include <iomanip>
#include <string>
#include <complex>
#include <sstream>
#include <vector>

#include <typeinfo>
#include <type_traits>

#include "cross_platform.h"
#include "commUtil_stub.h"

using std::cout;
using std::endl;

namespace commUtil{


comm_error& comm_error::operator=(const comm_error& other)
{ 
  base_class::operator=(static_cast<const base_class&>(other));     
  return *this;
}

comm_error::~comm_error(void)
{ }

comm_error::comm_error(void)
  : base_class("unspecified comm_error")
{ }

comm_error::comm_error(const comm_error& other)
  : base_class(static_cast<const std::runtime_error&>(other))
{ }

comm_error::comm_error(const std::string& msg)
  : base_class(msg)
{ }


// cannot use _other_ read/write methods which may not have been declared yet:
template <>
bool readBinary(std::istream &in, std::string& s)
{
  bool status(true);
  size_t N(0);
  status = status && in.read(reinterpret_cast<char *>(&N), sizeof(size_t));
  s.resize(N);
  for(std::string::iterator itS = s.begin(), itSEnd = s.end();
      status && (itS != itSEnd);
		  ++itS)
	  status = status && in.read(&(*itS), sizeof(char));	

  return status;	  
}

template <>
bool writeBinary(std::ostream &out, const std::string& s)
{
  bool status(true);
  const size_t N(s.size());
  status = status && out.write(reinterpret_cast<const char *>(&N), sizeof(size_t));
  for(std::string::const_iterator itS = s.begin(), itSEnd = s.end();
      status && (itS != itSEnd);
		  ++itS)
	  status = status && out.write(&(*itS), sizeof(char));	

  return status;	  
}

template <>
size_t binarySize(const std::string& s)
{
  size_t val(0);

  val += sizeof(size_t);

  // assume 8-bit characters:
  val += s.size();

	return val;
}  


} // namespace commUtil
