
// $Source: /usr/data0/leipzig_work/tmat_cvs/src/scanner_util_module.cpp,v $

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


#include "cross_platform.h"

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
#if defined(__PGI)
  // TACC port: see discussion in "portability.h" at definition of _STL_EXT_NAMESPACE_ macro:
  // #include <slist>
  #include <hash_map>
#else
  // #include <ext/algorithm>
  // #include <ext/numeric>
  // #include <ext/slist>
  #include <ext/hash_map>
#endif

// use the _outer_ header, in order to obtain the definition of the hash function for std::string:
#include <scanner_util.h> 

namespace scanner_util{
namespace implementation_module{



} // namespace implementation_module
} // namespace scanner_util
