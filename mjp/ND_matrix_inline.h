#if !defined(__ND_matrix_inline__h)
#define __ND_matrix_inline__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/ND_matrix_inline.h,v $

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


namespace linalg{



/**
 * @brief number of neighbor points participating in the interpolation
 */
template <class U, size_t NDIM>
inline size_t linear_interpolator<U,NDIM>::N_neighbor_(void)const
{ 
  // note: > max-dimensions test at constructor, not here.
  return 1UL<<NDIM; 
}


/**
 * @brief index-tuple n0, and co-ordinate-tuple x0 of lowest-index corner of hypercube containing point x_i:
 */
template <class U, size_t NDIM>
inline void linear_interpolator<U,NDIM>::base_index_(const ntuple<R,NDIM>& x_i, ntuple<size_t,NDIM>& n0, ntuple<R,NDIM>& x0)const
{
  // note: out-of-domain test at apply method, not here.
  for(size_t nd = 0; nd < NDIM; ++nd){
    n0[nd] = static_cast<size_t>(floor(((x_i[nd] - domain_.start()[nd])/scale_[nd])*integer<R>(shape_[nd]-1)));
    
    // allow x_i[nd] == domain_.end()[nd]:
    // (max for _base_, however, is shape[nd] - 2)
    if (n0[nd] >= shape_[nd]-1)
      n0[nd] = shape_[nd] - 2; 
  
    x0[nd] = domain_.start()[nd] + ratio<R,size_t>(n0[nd], shape_[nd]-1) * scale_[nd];
  }
}


/**
 * @brief interpolation-kind includes periodic boundary-conditions.
 */
template <class U, size_t NDIM>
inline bool linear_interpolator<U,NDIM>::periodic_BC_(void)const
{ return (kind_ == NEAREST_PERIODIC) || (kind_ == LINEAR_PERIODIC); }
 
/**
 * @brief index-tuple n of neighbor from base index-tuple n0 and packed-offset index nn:
 * (neighbors are relative to lowest-index corner of hypercube containing x_i:
 *    offset bits of "nn", by dimension as:
 *      0: no offset 
 *      1: +1 offset)
 * @return  true if neighbor exists (note: if periodic B.C., neighbor always exists).
 */
template <class U, size_t NDIM>
inline bool linear_interpolator<U,NDIM>::neighbor_(const ntuple<size_t,NDIM>& n0, size_t nn, ntuple<size_t,NDIM>& n)const
{
  bool wrap(false);
    
  for(size_t nd = 0; nd < NDIM; ++nd){
    if (!(nn & (1<<nd)))
      n[nd] = n0[nd];
    else{
      n[nd] = n0[nd] + 1;
      if (n[nd] >= shape_[nd]){
        n[nd] = 0;
        wrap = true;
      }
    }  
  }
  
  return !wrap || periodic_BC_();
}

/**
 * @brief linear basis functions
 *   @param[in] dx  N-dimensional offset from lowest-index corner of hypercube (i.e. containing x_i)
 *   @param[in] nn  index of basis-function: corresponds to packed-bit neighbor offset-index  
 */
template <class U, size_t NDIM>
inline typename linear_interpolator<U,NDIM>::S linear_interpolator<U,NDIM>::B_(const ntuple<R,NDIM>& dx, size_t nn)const
{
  // notes: out-of-domain test at apply method, not here.
  //   R: magnitude-type
  //   S: scalar-type (e.g. S may be complex => convert as late as possible from known real values)
  
  S val(one<S>());

  for(size_t nd = 0; nd < NDIM; ++nd){
    if (nn & (1<<nd))
      val *= dx[nd]*integer<R>(shape_[nd]-1)/scale_[nd];
    else
      val *= one<R>() - dx[nd]*integer<R>(shape_[nd]-1)/scale_[nd];
  }
  
  return val;
}


template <class U, size_t NDIM>
inline typename linear_interpolator<U,NDIM>::value_type linear_interpolator<U,NDIM>::operator()(const ntuple<R,NDIM>& x_i)const throw(std::runtime_error)
{ 
  value_type y_i;
  apply(x_i, y_i);
  return y_i;
}


} // namespace linalg

#endif // __ND_matrix_inline__h
