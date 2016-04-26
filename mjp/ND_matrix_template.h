#if !defined(__ND_matrix_template__h)
#define __ND_matrix_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/ND_matrix_template.h,v $

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


// ---------------------- moved from particle_packing_template.h: -----------------------------

template <class N, size_t NDIM>   
template <class N1>     
inline N1 row_major_index<N,NDIM>::product_(const ntuple<N1,NDIM>& t, size_t offset)
{
 #if 0
 N1 P(static_cast<N1>(1));
 for(size_t n = offset; n < NDIM; ++n)
   P *= t[n];   
 return P;
 #else
 return product(t, offset, NDIM);
 #endif
}
  
template <class N, size_t NDIM>        
inline row_major_index<N,NDIM> row_major_index<N,NDIM>::operator+(const row_major_index& other)const
{ 
  row_major_index val(*this);
  val += other;
  return val;
}

template <class N, size_t NDIM>        
inline row_major_index<N,NDIM> row_major_index<N,NDIM>::operator-(const row_major_index& other)const
{ 
  row_major_index val(*this);
  val -= other;
  return val;
}

template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>& row_major_index<N,NDIM>::operator+=(const row_major_index& other)
{ 
  base_class::operator+=(other); 
  return *this;
}

template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>& row_major_index<N,NDIM>::operator-=(const row_major_index& other)
{ 
  base_class::operator-=(other); 
  return *this;
}

// apply periodic B.C. as constrained by shape:
template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>& row_major_index<N,NDIM>::mod_assign(const shape_type& shape)
{
  // mod(a, b) \defas a - floor(a/b)*b
  //   where floor(x): smallest integer n s.t. n <= x
  
  for(size_t n = 0; n < NDIM; ++n){
    typename base_class::value_type 
      &x(base_class::operator[](n)),
      ns(static_cast<value_type>(shape[n])),
      nq(static_cast<value_type>(floor(static_cast<double>(x)/static_cast<double>(ns))));
    x -= nq*ns;
  }
  return *this;
}

template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>& row_major_index<N,NDIM>::operator=(const row_major_index& other)
{ 
  base_class::operator=(other); 
  return *this;
}


#if !defined(__INTEL_COMPILER) && !defined(__PGI)
template <class N, size_t NDIM>        
template <class N1>
inline row_major_index<N,NDIM>& row_major_index<N,NDIM>::operator=(const ntuple<N1,NDIM>& other)
{ 
  return base_class::operator=(other); 
  return *this;
}
#endif

/*
 * row_major_index<N,NDIM>::inverse(shape, shift=0)
 *   linear index from index-tuple constrained by shape
 */
 
template <class N, size_t NDIM>        
inline size_t row_major_index<N,NDIM>::inverse(const shape_type& shape, value_type shift)const
{ 
  size_t n_linear(0);
  #if 0
  for(size_t n = 0; n < NDIM; ++n){
    size_t n_i(static_cast<size_t>(base_class::operator[](n) - shift));
    assert(n_i < shape[n]);
    n_linear += n_i * product_(shape, n+1);
  }  
  #else
  // incremental product:
  size_t stride(1);
  for(size_t nd = NDIM-1; 
      nd != static_cast<size_t>(-1); 
      stride *= shape[nd], --nd){
    size_t n_i(static_cast<size_t>(base_class::operator[](nd) - shift));
    assert(n_i < shape[nd]);
    n_linear += n_i * stride;
  }
  #endif
  return n_linear;  
}

/*
 * test wrap-around condition for index "x2" with respect to index "x1";
 *   initialize dimensions corresponding to wrapping with 0 or +-1, corresponding to wrap direction;
 *   return true if any dimensions wrap.
 */
template <class N, size_t NDIM>        
bool row_major_index<N,NDIM>::wrap(const row_major_index& x1, const row_major_index& x2, const shape_type& shape,
                                           row_major_index& wrap_dim)
{
  bool test(false);

  #if 0 // --------------------- previous version: 2.7e8 counts -------------------------
  for(size_t n = 0; n < NDIM; ++n){
    value_type &xw(wrap_dim[n]);
    xw = 0;
    if ((x1[n] == static_cast<value_type>(shape[n])-1) && (x2[n] == 0)){
      test = true;
      xw = 1;
    }
    else
    if ((x1[n] == 0) && (x2[n] == static_cast<value_type>(shape[n])-1)){
      test = true;
      xw = -1;
    }
  }  
  #endif
  
  #if 1 // ------------------------- new version 0: 2.5e8 counts  ---------------------------
        // (note: speed-up this version, only about 10%;
        //    all other versions were _much_ worse)
  for(size_t n = 0; n < NDIM; ++n){
    value_type &xw(wrap_dim[n]);
    xw = 0;
    if (x1[n] == 0){
      if (x2[n] == static_cast<value_type>(shape[n])-1){
        test = true;
        xw = -1;      
      }
    }
    else
    if (x2[n] == 0){
      if (x1[n] == static_cast<value_type>(shape[n])-1){
        test = true;
        xw = 1;      
      }
    }
  }  
  #endif // ------------------------ end: new version 0 ----------------------

  return test;  
}                                           

template <class N, size_t NDIM>        
inline bool row_major_index<N,NDIM>::writeBinary(commUtil::abstractCommHandle* fp)const
{ return base_class::writeBinary(fp); }

template <class N, size_t NDIM>        
inline bool row_major_index<N,NDIM>::readBinary(commUtil::abstractCommHandle* fp)
{ return base_class::readBinary(fp); }

template <class N, size_t NDIM>        
inline size_t row_major_index<N,NDIM>::binarySize(void)const
{ return base_class::binarySize(); }

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
/*
 * return the row_major_index corresponding to the position "p" within the hyper-cube "domain_cube"
 *   given the array "shape" and index-base "shift"
 */
template <class N, size_t NDIM>        
template <class R>
row_major_index<N,NDIM> row_major_index<N,NDIM>::ND_bin(
  const ntuple<R,NDIM>& p, const ntuple_interval<R,NDIM>& domain_cube,
  const shape_type& shape, value_type shift)
{
  row_major_index result(0,false);
  
  // "left_closure" \leftarrow  x \in [start, end)
  if (!(domain_cube.left_closure(p)))
    throw std::runtime_error("row_major_index<N,NDIM>::ND_bin: ntuple not in domain hypercube");
    
  for(size_t n = 0; n < NDIM; ++n){
    R L((domain_cube.end()[n] - domain_cube.start()[n])/integer<R>(shape[n])), 
      nx((p[n] - domain_cube.start()[n])/L);
    conv(result[n], floor(nx));
    result[n] += shift;  
  }
  
  return result;  
}  
#endif

#if 0  // completely screws-up the compiler for some reason (in unrelated code sections):
/*
 * return the sub-hypercube corresponding to the row_major_index within the hyper-cube "domain_cube"
 *   given the array "shape" and index-base "shift"
 */
template <class N, size_t NDIM>        
ntuple_interval<R,NDIM> row_major_index<N,NDIM>::ND_edges(
  const row_major_index& indices, const ntuple_interval<R,NDIM>& domain_cube,
  const shape_type& shape, value_type shift, const R& epsilon_)
{
  ntuple<R,NDIM> start_, end_;

  if (!indices.in_domain(shape, shift))
    throw std::runtime_error("row_major_index<N,NDIM>::ND_edges: indices out-of-range for specified shape and base shift");

  for(size_t n = 0; n < NDIM; ++n){
    R L_bin((domain_cube.end()[n] - domain_cube.start()[n])/integer<R>(shape[n])); 
    start_[n] = integer<R>(indices[n] - shift) * L_bin;
    end_[n] = integer<R>(indices[n] - shift + 1) * L_bin;
  }  

  return ntuple_interval<R,NDIM>(start_, end_, epsilon_);  
}  
#endif


/*
 * row_major index in-domain of given shape and base-shift
 */
template <class N, size_t NDIM>        
inline bool row_major_index<N,NDIM>::in_domain(const shape_type& shape, value_type shift)const
{
  bool test(true);

  for(size_t n = 0; n < NDIM; ++n){
    const value_type x(base_class::operator[](n) - shift);
    if ((static_cast<typename shape_type::value_type>(x) >= shape[n]) || (x < 0)){
      test = false;
      break;
    }
  }
  
  return test;
} 

template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>::row_major_index(value_type shift, bool clear)
  :base_class(clear)
{
  if (clear && (shift != 0)) 
    base_class::operator+=(shift); 
}

/*
 * row_major_index(<raw linear index>, <shape>, shift=0)
 *   index-tuple from linear index constrained by shape
 *   shift => starting offset for index intervals
 */
template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>::row_major_index(size_t n_linear, const shape_type& shape, value_type shift)
  : base_class(false)
{  

  size_t n_rem(n_linear), N_elts(0), x(0);
  for(size_t n = 0; n < NDIM; ++n){
    // number of elements per increment at this ndim:
    N_elts = product_(shape, n+1);
    x = n_rem/N_elts; // truncation
    base_class::operator[](n) = static_cast<value_type>(x) + shift;
    n_rem -= x * N_elts;
  }
}


template <class N, size_t NDIM>        
inline row_major_index<N,NDIM>::row_major_index(const row_major_index& other)
  :base_class(other)
{ }

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
template <class N, size_t NDIM>        
template <class N1>
inline row_major_index<N,NDIM>::row_major_index(const ntuple<N1,NDIM>& other)
  :base_class(other)  
{ }
#endif

// ---------------------- end: moved from particle_packing_template.h -----------------------------


/**
 * @brief Dereference an N-dimensional matrix.
 *   @tparam U a linear container class with fixed-allocation value-type (e.g. scalar, ntuple, or tensor)
 *   @tparam NDIM  number of dimensions of this matrix
 *   @param[in]  u  container
 *   @param[in]  shape  dimension limits
 *   @param[in]  n  ntuple index of element to reference
 */    
template <class U, size_t NDIM>
typename gmm::linalg_traits<U>::reference
  ND_deref(U& u, const ntuple<size_t,NDIM>& shape, const ntuple<size_t,NDIM>& n)
{
  const row_major_index<size_t,NDIM>& n_(static_cast<const row_major_index<size_t,NDIM>&>(n));
  return u[n_.inverse(shape)]; 
}

/**
 * @brief Dereference an N-dimensional matrix.
 *   @tparam U a linear container class with fixed-allocation value-type (e.g. scalar, ntuple, or tensor)
 *   @tparam NDIM  number of dimensions of this matrix
 *   @param[in]  u  container
 *   @param[in]  shape  dimension limits
 *   @param[in]  n  ntuple index of element to reference
 */    
template <class U, size_t NDIM>
const typename gmm::linalg_traits<U>::reference
  ND_deref(const U& u, const ntuple<size_t,NDIM>& shape, const ntuple<size_t,NDIM>& n)
{
  const row_major_index<size_t,NDIM>& n_(static_cast<const row_major_index<size_t,NDIM>&>(n));
  return const_cast<U&>(u)[n_.inverse(shape)]; // force by-reference return.
}

/** 
 * @brief Weights for central-difference derivative calculation.
 *   param[in]  Np  number-of-points in calculation
 *   param[in]  deriv_order  derivative order
 */
template <class T>
std::vector<T> central_diff_weights(size_t Np, size_t deriv_order)
{
  // based on python implementation from scipy-0.7.1: "common.py: central_diff_weights":
  using linalg::factorial;
  using number::pow_n;
  
  if ((Np < deriv_order+1) || (Np % 2 == 0))
    throw std::runtime_error("central_diff_weights: Number of points must be odd, and at least the derivative order + 1.");

  gmm::dense_matrix<T> X(Np, Np);
  for(size_t nr = 0; nr < Np; ++nr)
  for(size_t nc = 0; nc < Np; ++nc)
    X(nr,nc) = pow_n<T>(integer<T>(static_cast<long>(nr) - static_cast<long>(Np/2)), static_cast<long>(nc));

  gmm::lu_inverse(X);  

  std::vector<T> w(Np, zero<T>());
  gmm::copy(gmm::mat_row(X, deriv_order), w);

  long prefactor = factorial<long>(static_cast<size_t>(deriv_order));
  gmm::scale(w, integer<T>(prefactor));

  return w;   
}

/** 
 * @brief Central-difference numerical derivative calculation.
 *   @param[in]  src data points
 *   @param[in]  dest output points
 *   @param[in]  shape  dimension limits for src and dest  
 *   @param[in]  dim dimension over which to take derivative
 *   @param[in]  Np  number-of-points included in calculation of each derivative point
 *   @param[in]  div_order  derivative order
 *   @param[in]  boundary_condition  enum: 0: periodic_BC 
 *     - Np >= div_order+1, and must be odd
 *     - At present, only implemented B.C. is periodic
 *     .
 */
template <class U1, class U2, size_t NDIM>
void central_diff(
  U1& dest, const ntuple<size_t, NDIM>& shape, 
  const U2& src, size_t dim, size_t Np, size_t deriv_order, int BC_enum)
{
  // assume U1 and U2 are _linear_ container-types 
  //   (this is what the present N-dimensional array implementation expects, in addition to the "shape" ntuple).
  // assume U1::value_type and U2::value_type are "compatible" 
  //   (e.g. support direct conversion and mixed arithmetic operations)
  // data-reduction of source => implement all calculations using scalar_type of source.
  typedef typename gmm::linalg_traits<U1>::value_type T1;
  typedef typename gmm::linalg_traits<U2>::value_type T2;
  typedef typename tensor_traits<T2>::scalar_type S;
  
  if (dim >= NDIM)
    throw std::runtime_error("central_diff: value of \"dim\" parameter out-of-range; \"dim\" must be < NDIM");
  
  if (0 != BC_enum)
    throw std::runtime_error("central_diff: value of boundary-condition enum not recognized, supported values are: 0: periodic B.C. (only)");
    
  if (Np > shape[dim])
    throw std::runtime_error("central_diff: source has less than number-of-points for derivative dimension");

  std::vector<S> weight(central_diff_weights<S>(Np, deriv_order));
  
  // resize dest only if destination is _not_ a reference type:
  // see note at "simple_object_base::is" regarding type_info comparison by-name, to enable cross-module use with current GNU C++ ABI:
  if (0 == strcmp(typeid(typename gmm::linalg_traits<U1>::is_reference).name(), typeid(gmm::linalg_false).name()))
    gmm::linalg_traits<U1>::resize(dest, gmm::vect_size(src)); 
  else
  if (gmm::vect_size(dest) != gmm::vect_size(src))
    throw std::runtime_error("central_diff: for reference types, destination must be pre-allocated to size of source");
         
  gmm::clear(dest); // initialize to zero

  /* 
   * Implementation notes: 
   *   - this implementation assumes that  benefits for BLAS1 utilization (i.e. by-row processing) are significant, 
   *     but that speedup for processing with 2D array blocks is not particularly so.  If this is _not_ the case, then
   *     optimization should be applied for 2D block processing when NDIM >= 2.
   *   - processing is separated into the common middle section, and the boundary sections;
   *       this will eventually allow different boundary-condition types to be implemented.
   *   .
   */
  const size_t Ncol(shape[NDIM-1]); 
  if (NDIM-1 == dim){
    for(size_t nraw=0, nraw_max = src.size()-1; nraw <= nraw_max; nraw += Ncol){
      T1 *dest_row_begin(&*(dest.begin() + nraw)), *dest_row_end(dest_row_begin + Ncol);      
      const T2 *src_row_begin(&*(src.begin() + nraw)), *src_row_end(src_row_begin + Ncol); 
      for(size_t nw = 0; nw < Np; ++nw){
        // process central section:
        //   y += a*x: BLAS1: axpy(n, a, x, incx, y, incy):
        BLAS::axpy<T2,T1,S>(Ncol-Np, weight[nw], src_row_begin+nw, 1, dest_row_begin+Np/2, 1);
        
        // process boundary sections: Np/2 .... Np/2
        //  source section: [-Np+nw:],wrap,[:nw]  ::=  s1, s2
        //    n-points: 
        //      s1: Np-nw
        //      s2: nw
        //  dest.  section: [-Np/2:],wrap,[:Np/2]  ::=  d1, d2
        //   n-points:
        //      d1: Np/2
        //      d2: Np/2
        
        if (nw <= Np/2){
          BLAS::axpy<T2,T1,S>(Np/2, weight[nw], src_row_end-Np+nw, 1, dest_row_end-Np/2, 1);
          BLAS::axpy<T2,T1,S>(Np/2-nw, weight[nw], src_row_end-Np/2+nw, 1, dest_row_begin, 1);
          BLAS::axpy<T2,T1,S>(nw, weight[nw], src_row_begin, 1, dest_row_begin+Np/2-nw, 1);
        }
        else{
          BLAS::axpy<T2,T1,S>(Np-nw, weight[nw], src_row_end-Np+nw, 1, dest_row_end-Np/2, 1);
          // Np/2 - (Np - nw) --> dest_row_end - Np/2 + (Np - nw): 
          BLAS::axpy<T2,T1,S>(nw-Np/2, weight[nw], src_row_begin, 1, dest_row_end+Np/2-nw, 1);                
          BLAS::axpy<T2,T1,S>(Np/2, weight[nw], src_row_begin+nw, 1, dest_row_begin, 1);        
        }
      }    
    }
  }
  else{
    const size_t 
      Nslice(product(shape, 0,dim+1)),
      slice_stride(product(shape, dim+1,NDIM));

    // process central section:
    for(size_t n_slice = Np/2; n_slice < Nslice-Np/2; ++n_slice){
      T1 *dest_slice_begin(&*(dest.begin() + n_slice*slice_stride)),
         *dest_slice_end(dest_slice_begin + slice_stride);      
      for(size_t nw = 0; nw < Np; ++nw){
        const T2 
          *src_slice_begin(&*(src.begin() + (n_slice - Np/2 + nw)*slice_stride)),
          *it_src(src_slice_begin);
        for(T1 *it_dest = dest_slice_begin, *it_destEnd = dest_slice_end;
            it_dest != it_destEnd;
            it_dest += Ncol, it_src += Ncol)
          BLAS::axpy<T2,T1>(Ncol, weight[nw], it_src, 1, it_dest, 1);          
      }
    }
    
    // process boundary sections:
    for(size_t n_slice = 0; n_slice < Np/2; ++n_slice){
      T1 *dest_slice_begin(&*(dest.begin() + n_slice*slice_stride)),
         *dest_slice_end(dest_slice_begin + slice_stride);      
      for(size_t nw = 0; nw < Np; ++nw){
        // periodic B.C.:
        size_t n_wrap = static_cast<size_t>(mod(static_cast<long>(nw) + static_cast<long>(n_slice) - static_cast<long>(Np/2), static_cast<long>(Nslice)));
        const T2 
          *src_slice_begin(&*(src.begin() + n_wrap*slice_stride)),
          *it_src(src_slice_begin);
        for(T1 *it_dest = dest_slice_begin, *it_destEnd = dest_slice_end;
            it_dest != it_destEnd;
            it_dest += Ncol, it_src += Ncol)
          BLAS::axpy<T2,T1>(Ncol, weight[nw], it_src, 1, it_dest, 1);          
      }
    }

    for(size_t n_slice = Nslice-Np/2; n_slice < Nslice; ++n_slice){
      T1 *dest_slice_begin(&*(dest.begin() + n_slice*slice_stride)),
         *dest_slice_end(dest_slice_begin + slice_stride);      
      for(size_t nw = 0; nw < Np; ++nw){
        // periodic B.C.:
        size_t n_wrap = static_cast<size_t>(mod(static_cast<long>(nw) + static_cast<long>(n_slice) - static_cast<long>(Np/2), static_cast<long>(Nslice)));
        const T2 
          *src_slice_begin(&*(src.begin() + n_wrap*slice_stride)),
          *it_src(src_slice_begin);
        for(T1 *it_dest = dest_slice_begin, *it_destEnd = dest_slice_end;
            it_dest != it_destEnd;
            it_dest += Ncol, it_src += Ncol)
          BLAS::axpy<T2,T1>(Ncol, weight[nw], it_src, 1, it_dest, 1);          
      }
    }    
  }
}




template <class U, size_t NDIM>
void linear_interpolator<U,NDIM>::apply(const ntuple<R,NDIM>& x_i, value_type& y_i)const
{
  if (!domain_.closure(x_i))
    throw std::runtime_error("linear_interpolator<U,NDIM>::apply: coordinates out of instantiated domain");
    
  // index and co-ordinates of lowest-index corner of hypercube containing x_i:
  ntuple<size_t,NDIM> n0;
  ntuple<R,NDIM> x0;
  base_index_(x_i, n0, x0);          

  // offset from x0 corner:        
  ntuple<R,NDIM> dx(x_i);
  dx -= x0;
 
  
  // absolute neighbor index:
  ntuple<size_t,NDIM> n;

  
  switch (kind_){
    
    case NEAREST:
    case NEAREST_PERIODIC:
    {
      // absolute neighbor index:
      ntuple<size_t,NDIM> n;
      neighbor_(n0, 0, n);
      
      y_i = ND_deref<U,NDIM>(data_, shape_, n);         
    }
    break;
    
    case LINEAR:
    case LINEAR_PERIODIC:
    {
      // offset from x0 corner:        
      ntuple<R,NDIM> dx(x_i);
      dx -= x0;
      
      // value_type shall support assigment to scalar_type (i.e. meaning initialize all components to scalar).      
      y_i = zero<S>();
      
      // absolute neighbor index:
      ntuple<size_t,NDIM> n;
      for(size_t nn = 0, nnEnd = N_neighbor_();
          nn != nnEnd;
          ++nn)
        if (neighbor_(n0, nn, n)){        
          const S &yn(ND_deref<U,NDIM>(data_, shape_, n));
          y_i += B_(dx,nn) * yn;
        }         
    }
    break;
    
    default:
      throw std::runtime_error("linear_interpolator<U,NDIM>::apply: unrecognized interp_kind enum");
  }
}  


#if !defined(__INTEL_COMPILER) && !defined(__PGI)
/**
* @brief interpolation over N-dimensional matrix of co-ordinates.
 *   @tparam  U1  a linear container with value_type compatible with  ntuple<R,NDIM>
 *   @tparam  U2  a linear container with value_type compatible with linear_interpolator<U, NDIM>::value_type
 *   @param[in] x_i coordinate values
 *   @param[out] y_i  interpolated data values
 */
template <class U, size_t NDIM>
template <class U1, class U2>
void linear_interpolator<U,NDIM>::apply(const U1& x_i, U2& y_i)const
{
  // resize dest only if destination is _not_ a reference type:
  // see note at "simple_object_base::is" regarding type_info comparison by-name, to enable cross-module use with current GNU C++ ABI:
  if (0 == strcmp(typeid(typename gmm::linalg_traits<U2>::is_reference).name(), typeid(gmm::linalg_false).name()))
    gmm::linalg_traits<U2>::resize(y_i, gmm::vect_size(x_i)); 
  else
  if (gmm::vect_size(y_i) != gmm::vect_size(x_i))
    throw std::runtime_error("linear_interpolator<U,NDIM>::apply: for reference types, destination must be pre-allocated to size of source");
 
  ntuple<R,NDIM> x_i_; // temporaries
  value_type y_i_; 
  typename gmm::linalg_traits<U1>::const_iterator it_x(gmm::vect_const_begin(x_i));
  for(typename gmm::linalg_traits<U2>::iterator it_y = gmm::vect_begin(y_i), it_yEnd = gmm::vect_end(y_i);
      it_y != it_yEnd;
      ++it_y, ++it_x){
    conv(x_i_, *it_x);
    apply(x_i_, y_i_);
    conv(*it_y, y_i_);  
  }     
}
#endif


template <class U, size_t NDIM>
linear_interpolator<U,NDIM>::linear_interpolator
  (const ntuple_interval<R,NDIM>& domain, const U& y, const ntuple<size_t,NDIM>& shape, interp_kind e)
  : domain_(domain), scale_(domain.scale()), data_(y), shape_(shape), kind_(e)
{
  if (NDIM+1 > sizeof(size_t)*8)
    throw std::runtime_error("linear_interpolator<U,NDIM>::linear_interpolator: NDIM exceeds capacity of packed bit-width");
}


template <class U, size_t NDIM>
linear_interpolator<U,NDIM>::linear_interpolator
  (const ntuple_interval<R,NDIM>& domain, U& y, const ntuple<size_t,NDIM>& shape, interp_kind e, bool transfer_ownership)
  : domain_(domain), scale_(domain.scale()), shape_(shape), kind_(e)
{
  if (NDIM+1 > sizeof(size_t)*8)
    throw std::runtime_error("linear_interpolator<U,NDIM>::linear_interpolator: NDIM exceeds capacity of packed bit-width");

  if (transfer_ownership)
    data_.swap(y);
  else
    data_ = y;
}  


} // namespace linalg

#endif // __ND_matrix_template__h
