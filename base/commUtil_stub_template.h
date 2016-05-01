#if !defined(__commUtil_stub_template__h)
#define __commUtil_stub_template__h

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


namespace commUtil{
      
//
// _generic_ linalg writeBinary and readBinary.
//

template <class U>
inline bool writeBinary_POD_dispatch_(std::ostream &out, const U& u, std::false_type)
{ return u.writeBinary(out); }

template <class U>
inline bool writeBinary_POD_dispatch_(std::ostream &out, const U& u, std::true_type)
{ 
  bool status(true);
  status = status && out.write(reinterpret_cast<const char *>(&u), sizeof(u));	

  return status; 
}

template<class U>
inline bool writeBinary(std::ostream & out, const U& u)
{ return writeBinary_POD_dispatch_(out, u, typename std::is_pod<U>::type()); }

template <class U>
inline bool readBinary_POD_dispatch_(std::istream &in, U& u, std::false_type)
{ return u.readBinary(in); }

template <class U>
inline bool readBinary_POD_dispatch_(std::istream &in, U& u, std::true_type)
{ 
  bool status(true);
  status = status && in.read(reinterpret_cast<char *>(&u), sizeof(u));	

  return status; 
}

template<class U>
inline bool readBinary(std::istream &in, U& u)
{ return readBinary_POD_dispatch_(in, u, typename std::is_pod<U>::type()); }

template <class U>
inline size_t binarySize_POD_dispatch_(U& u, std::false_type)
{ return u.binarySize(); }

template <class U>
inline size_t binarySize_POD_dispatch_(U& u, std::true_type)
{ return sizeof(u); }

template <class U>
inline size_t binarySize(const U& u)
{ return binarySize_POD_dispatch_(u, typename std::is_pod<U>::type()); }


// "std::complex<double>" is *not* a POD-type, however with respect to I/O it can be treated as a POD type,
//    therefore, "readBinary", "writeBinary", and "binarySize" need to be specialized to deal with this type.

template <>
inline bool readBinary(std::istream &in, std::complex<double> &c)
{ return readBinary_POD_dispatch_(in, c, std::true_type()); }

template <>
inline bool writeBinary(std::ostream &out, const std::complex<double> &c)
{ return writeBinary_POD_dispatch_(out, c, std::true_type()); }

template <>
inline size_t binarySize(const std::complex<double> &c)
{ return binarySize_POD_dispatch_(c, std::true_type()); }

				
} // namespace commUtil

#endif // __commUtil_stub_template__h
