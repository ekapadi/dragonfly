#if !defined(__commUtil_template__h)
#define __commUtil_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/commUtil_template.h,v $

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


namespace commUtil{


// ----- global methods: ---------

// ---- these won't work correctly until after inclusion of "commUtil_template_forward.h": ----

#if defined(__USE_MPI) // ------------------------------------------- MPI ONLY SECTION ------------------------------------------------------------
template <class U>
inline void gather(const U& u, abstractCommHandle *h)
{  
 processCommHandle *h_(dynamic_cast<processCommHandle*>(h));
 if (NULL == h_)
   throw std::runtime_error("commUtil::gather: cast to derived processCommHandle* fails");
 h_->gather(u); 
}

template <class U>
inline void gather(const U& u, std::vector<U>& vU, abstractCommHandle *h)
{ 
 processCommHandle *h_(dynamic_cast<processCommHandle*>(h));
 if (NULL == h_)
   throw std::runtime_error("commUtil::gather: cast to derived processCommHandle* fails");
 h_->gather(u, vU); 
}

template <class U>
inline void all_gather(const U& u, std::vector<U>& vU, abstractCommHandle *h)
{ 
 processCommHandle *h_(dynamic_cast<processCommHandle*>(h));
 if (NULL == h_)
   throw std::runtime_error("commUtil::all_gather: cast to derived processCommHandle* fails");
 h_->all_gather(u, vU); 
}

template <class U>
inline void scatter(U& u, abstractCommHandle *h)
{ 
 processCommHandle *h_(dynamic_cast<processCommHandle*>(h));
 if (NULL == h_)
   throw std::runtime_error("commUtil::scatter: cast to derived processCommHandle* fails");
 h_->scatter(u); 
}

template <class U>
inline void scatter(const std::vector<U>& vU, U& u, abstractCommHandle *h)
{ 
 processCommHandle *h_(dynamic_cast<processCommHandle*>(h));
 if (NULL == h_)
   throw std::runtime_error("commUtil::scatter: cast to derived processCommHandle* fails");
 h_->scatter(vU, u); 
}

#endif // ---------------------------------------end:  MPI ONLY SECTION ------------------------------------------------------------

} // namespace commUtil

#endif // __commUtil_template__h
