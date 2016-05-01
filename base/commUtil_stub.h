#if !defined(__commUtil_stub__h)
#define __commUtil_stub__h

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
 * @file: commUtil_stub.h
 *   *Temporary* reduction of "commUtil" dependency prior to *upgrade* of "commUtil" classes to use std::iostream-derived binary I/O.
 */
 
/**
 * @namespace commUtil
 * @ingroup parallel
 * @brief Unified low-level binary interface to disk-resident, memory-resident, and MPI inter-process data transfer.
 *
 * Using these methods to implement ``readBinary'', ''writeBinary'', and ``binarySize'' member or static methods for an end-user class
 * is sufficient to provide unified disk, memory, and/or inter-process I/O of class instances.
 */

namespace commUtil{


/**
 * @brief Base-class for binary-I/O exceptions.
 */
class comm_error: public std::runtime_error{
  private:

  public:
    typedef std::runtime_error base_class;
    
    comm_error& operator=(const comm_error& other);    
    
    virtual ~comm_error(void);
    
    comm_error(void);
    comm_error(const comm_error& other);
    comm_error(const std::string& msg);
}; // class comm_error 


// **********************************************************************
//
// _generic_ writeBinary, readBinary and binarySize:
//

template<class U>
bool writeBinary(std::ostream &out, const U& u);

template<class U>
bool readBinary(std::istream &in, U& u );

// in consideration of variable-precision number types, the following "binarySize" method requires an instance of class U:
template <class U>
size_t binarySize(const U& u);

// **********************************************************************

// "std::complex<double>" is *not* a POD-type, however with respect to I/O it can be treated as a POD type,
//    therefore, "readBinary", "writeBinary", and "binarySize" need to be specialized to deal with this type.

template <>
bool readBinary(std::istream &in, std::complex<double> &c);

template <>
bool writeBinary(std::ostream &out, const std::complex<double> &c);

template <>
size_t binarySize(const std::complex<double> &c);


// The following string methods are declared _here_ in order to allow this file to be included _first_, prior to other read/write declarations.
// These methods are required _now_ for the _transfer_ of error messages (e.g. between threads or processes):

template <>
bool readBinary(std::istream &in, std::string& s);

template <>
bool writeBinary(std::ostream &out, const std::string& s);

template <>
size_t binarySize(const std::string& s);


} // namespace commUtil

#include "commUtil_stub_template.h"

#endif // __commUtil_stub__h
