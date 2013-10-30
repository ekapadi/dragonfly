#if !defined(__particle_packing_inline__h)
#define __particle_packing_inline__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/particle_packing_inline.h,v $

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


namespace particle_packing{

// -------------- allow this module to be independent of the TMatrix namespace: ------------------------
#if !defined(__numericalConstants__h)
// --------------- explicit specializations: -----------------------------
#if defined(__USE_MERE)
template <>
inline mere::C integer<mere::C,long>(const long& n)
{
  return mere::C(n,0); 
}
template <>
inline const mere::C& zero<mere::C>(void)
{ return mere::C::zero(); }
template <>
inline const mere::C& epsilon<mere::C>(void)
{ return mere::C::epsilon(); }


template <>
inline mere::R integer<mere::R,long>(const long& n)
{
  return mere::R(n); 
}
template <>
inline const mere::R& zero<mere::R>(void)
{ return mere::R::zero(); }
template <>
inline const mere::R& epsilon<mere::R>(void)
{ return mere::R::epsilon(); }
#else
template <>
inline std::complex<double> integer<std::complex<double>,long>(const long& n)
{
  return std::complex<double>(static_cast<double>(n)); 
}
template <>
inline const std::complex<double> & zero<std::complex<double> >(void)
{ return std::complex<double>(0.0); }
template <>
inline const std::complex<double> & epsilon<std::complex<double> >(void)
{ return std::complex<double>(std::numeric_limits<std::complex<double> >::epsilon()); }

template <>
double integer<double,long>(const long& n)
{ return static_cast<double>(n); }
template <>
inline const double& zero<double>(void)
{ return 0.0; }
template <>
inline const double& epsilon<double>(void)
{ return std::numeric_limits<double>::epsilon(); }
#endif
// -----------------------------------------------------------------------
#endif
// ----------------------------- end: module independence section ------------------------------------------------

} // namespace particle_packing

#endif // __particle_packing_inline__h
