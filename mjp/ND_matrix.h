#if !defined(__ND_matrix__h)
#define __ND_matrix__h

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


/**
 * @file ND_matrix.h
 *   N-dimensional matrix and associated interpolation classes. 
 */
 // $Source: /usr/data0/leipzig_work/tmat_cvs/src/ND_matrix.h,v $
 
/**
 * @ingroup arbitrary_precision
 */

namespace linalg{
/**
 * @namespace linalg
 * @brief Support for fully-templated, arbitrary precision linear algebra.
 */


namespace BLAS{

/**
  * @brief Constant times a vector plus a vector: y += a*x: 
  *
  *   BLAS1: 
  *   @tparam  T1  source type
  *   @tparam  T2  destination type
  *   @tparam  S   scalar type
  *
  *   @param[in]  n  length of vectors 
  *     (for stride > 1: this is the number of elements to process, not the total length)
  *   @param[in]  a  scalar weight
  *   @param[in]  x  increment vector
  *   @param[in]  incx  stride for increment vector
  *     (negative stride processes vectors in reverse)
  *   @param[in/out]  y  accumulation vector
  *   @param[in]  incy  stride for accumulation vector
  */
template <class T1, class T2, class S>  
inline void axpy(size_t n, const S& a, const T1 *x, int incx, T2 *y, int incy)
{
  if (incx==1 && incy==1){
    const T1 *it_x(x);
    for(T2 *it_y = y, *it_yEnd = y + n;
        it_y != it_yEnd;
        ++it_y, ++it_x)
      *it_y += a * (*it_x);  
  }
  else{
    const T1 *it_x(incx > 0? x: x + -incx*n);
    for(T2 *it_y = (incy > 0? y: y + -incy*n), *it_yEnd = (incy > 0? y + n*incy: y-1);
        it_y != it_yEnd;
        it_y += incy, it_x += incx)
      *it_y += a * (*it_x);
  }
}

#if defined(GMM_USES_BLAS)
extern "C" void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern "C" void zaxpy_(...); // "..." follows usage in "gmm/gmm_blas_interface.h"

template < >
inline void axpy<double, double, double>(size_t n, const double& a, const double *x, int incx, double *y, int incy)
{
  int n_(static_cast<int>(n));
  daxpy_(&n_, &a, x, &incx, y, &incy); 
}

template < >
inline void axpy<std::complex<double>, std::complex<double>, std::complex<double> >
  (size_t n, const std::complex<double>& a, const std::complex<double> *x, int incx, std::complex<double> *y, int incy)
{
  int n_(static_cast<int>(n));
  zaxpy_(&n_, &a, x, &incx, y, &incy); 
}

template <class T, size_t NDIM>
inline void axpy(size_t n, const T& a, const ntuple<T,NDIM> *x, int incx, ntuple<T,NDIM> *y, int incy)
{
  assert(incx == incy == 1); // this specialization really doesn't make sense otherwise.
  size_t n_(n*NDIM);
  const T *x_(&(*x)[0]);
  T *y_(&(*y)[0]);
  axpy<T,T,T>(n_, a, x_, incx, y_, incy); 
}
#endif
 
} // namespace BLAS

// row-major index class (allowing arbitrary base and += sub-indices when type "N" is integer type):
// (IMPORTANT: this is a non-polymorphic class (no virtual functions) 
//     => static_cast<row_major_index<N,NDIM>&>(ntuple<N,NDIM>&) will work, without new allocation).
template <class N, size_t NDIM>
class row_major_index: public ntuple<N, NDIM>{

    // product of indices in [offset,NDIM)
    // (where product_(NDIM) \defas 1)
    // (templated to allow integer-based indices to work with size_t shape products)
    template <class N1>
    static N1 product_(const ntuple<N1,NDIM>& t, size_t offset=0);

  public:
    typedef ntuple<N, NDIM> base_class;
    typedef typename base_class::value_type value_type;
    typedef ntuple<size_t, NDIM> shape_type;

    row_major_index<N,NDIM> operator+(const row_major_index& other)const;
    row_major_index<N,NDIM> operator-(const row_major_index& other)const;

    row_major_index<N,NDIM>& operator+=(const row_major_index& other);
    row_major_index<N,NDIM>& operator-=(const row_major_index& other);

    // apply periodic B.C. as constrained by shape:
    row_major_index<N,NDIM>& mod_assign(const shape_type& shape);

    row_major_index<N,NDIM>& operator=(const row_major_index& other);

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
    template <class N1>
    row_major_index<N,NDIM>& operator=(const ntuple<N1,NDIM>& other);
#else
    // INTEL and PGI require body-definition at point-of-declaration:

    template <class N1>
    inline row_major_index<N,NDIM>& operator=(const ntuple<N1,NDIM>& other)
    { 
      return base_class::operator=(other); 
      return *this;
    }      
#endif

    size_t inverse(const shape_type& shape, value_type shift=0)const;

    /*
     * test wrap-around condition for index "x2" with respect to index "x1";
     *   initialize dimensions corresponding to wrapping with 0 or +-1, corresponding to wrap direction;
     *   return true if any dimensions wrap.
     */
    static bool wrap(const row_major_index& x1, const row_major_index& x2, const shape_type& shape,
                     row_major_index& wrap_dim);

    bool writeBinary(commUtil::abstractCommHandle* fp)const;
    bool readBinary(commUtil::abstractCommHandle* fp);
    size_t binarySize(void)const;

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
    /*
     * return the row_major_index corresponding to the position "p" within the hyper-cube "domain_cube"
     *   given the array "shape" and index-base "shift"
     */
    template <class R>
    static row_major_index<N,NDIM> ND_bin(const ntuple<R,NDIM>& p, const ntuple_interval<R,NDIM>& domain_cube,
                                  const shape_type& shape, value_type shift=0);
#else
    // INTEL and PGI require body-definition at point-of-declaration:

    /*
     * return the row_major_index corresponding to the position "p" within the hyper-cube "domain_cube"
     *   given the array "shape" and index-base "shift"
     */
    template <class R>
    row_major_index<N,NDIM> ND_bin(
      const ntuple<R,NDIM>& p, const ntuple_interval<R,NDIM>& domain_cube,
      const shape_type& shape, value_type shift)
    {
      row_major_index result(0,false);

      // "left_closure" \leftarrow  x \in [start, end)
      if (!(domain_cube.left_closure(p)))
         std::string("row_major_index<N,NDIM>::ND_bin: ntuple not in domain hypercube");

      for(size_t n = 0; n < NDIM; ++n){
        R L((domain_cube.end()[n] - domain_cube.start()[n])/integer<R>(shape[n])), 
          nx((p[n] - domain_cube.start()[n])/L);
        conv(result[n], floor(nx));
        result[n] += shift;  
      }

      return result;  
    }    
#endif

    bool in_domain(const shape_type& shape, value_type shift=0)const;


    row_major_index(value_type shift=0, bool clear=true);

    row_major_index(size_t n_linear, const shape_type& shape, value_type shift=0);

    row_major_index(const row_major_index& other);

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
    template <class N1>
    row_major_index(const ntuple<N1,NDIM>& other);  
#else
    // INTEL and PGI require body-definition at point-of-declaration:

    template <class N1>
    inline row_major_index(const ntuple<N1,NDIM>& other)
      :base_class(other)  
    { }    
#endif
}; // class row_major_index


/**
 * @brief Dereference an N-dimensional matrix.
 *
 *   @tparam U a linear container class with fixed-allocation value-type (e.g. scalar, ntuple, or tensor)
 *   @tparam NDIM  number of dimensions of this matrix
 *   @param[in]  u  container
 *   @param[in]  shape  dimension limits
 *   @param[in]  n  ntuple index of element to reference
 */    
template <class U, size_t NDIM>
typename gmm::linalg_traits<U>::reference
  ND_deref(U& u, const ntuple<size_t,NDIM>& shape, const ntuple<size_t,NDIM>& n);


/**
 * @brief Dereference an N-dimensional matrix.
 *
 *   @tparam U a linear container class with fixed-allocation value-type (e.g. scalar, ntuple, or tensor)
 *   @tparam NDIM  number of dimensions of this matrix
 *   @param[in]  u  container
 *   @param[in]  shape  dimension limits
 *   @param[in]  n  ntuple index of element to reference
 */    
template <class U, size_t NDIM>
const typename gmm::linalg_traits<U>::reference
  ND_deref(const U& u, const ntuple<size_t,NDIM>& shape, const ntuple<size_t,NDIM>& n);

    
/** 
 * @brief Weights for central-difference numerical derivative calculation.
 *
 *   param[in]  Np  number-of-points included in calculation of each derivative point
 *   param[in]  deriv_order  derivative order
 *   Np >= deriv_order+1, and must be odd
 */
template <class T>
std::vector<T> central_diff_weights(size_t Np=3, size_t deriv_order=1); 

/** 
 * @brief Central-difference numerical derivative calculation.
 *
 *   param[in]  src data points
 *   param[in]  dest output points
 *   param[in]  shape  dimension limits for src and dest  
 *   param[in]  dim dimension over which to take derivative
 *   param[in]  Np  number-of-points included in calculation of each derivative point
 *   param[in]  deriv_order  derivative order
 *   param[in]  boundary_condition  enum: 0: periodic_BC 
 *     - Np >= deriv_order+1, and must be odd
 *     - At present, only implemented B.C. is periodic
 *     .
 */
template <class U1, class U2, size_t NDIM>
void central_diff(U1& dest, const ntuple<size_t, NDIM>& shape, 
       const U2& src, size_t dim, size_t Np=3, size_t deriv_order=1, int BC_enum=0);

/**
 * @brief N-dimensional linear interpolation of values of arbitrary rank. 
 *
 *    - uses minimalist N-dimensional matrix semantics as (\<linear container\>, \<shape ntuple\>);
 *    - mesh corresponding to data instantiation is specified as ntuple_interval;
 *    - rank of value-type for N-dimensional matrix is arbitrary, but is a fixed-allocation type such as a scalar or ntuple<R,NDIM>.
 *    .
 */
template <class U, size_t NDIM>
class linear_interpolator{
  public:
  
  typedef typename gmm::linalg_traits<U>::value_type value_type;
  typedef typename tensor_traits<value_type>::scalar_type S;
  typedef typename numberTraits<S>::magnitudeType R;
  
  enum interp_kind {NEAREST=0, NEAREST_PERIODIC, LINEAR, LINEAR_PERIODIC};

  private:

    ntuple_interval<R,NDIM> domain_;
    ntuple<R,NDIM> scale_;
    
    U data_;
    ntuple<size_t,NDIM> shape_;
    
    interp_kind kind_;

    /**
     * @brief number of neighbor points participating in the interpolation
     */
    inline size_t N_neighbor_(void)const;
     
    /**
     * @brief index-tuple n0, and co-ordinate-tuple x0 of lowest-index corner of hypercube containing point x_i:
     */
    inline void base_index_(const ntuple<R,NDIM>& x_i, ntuple<size_t,NDIM>& n0, ntuple<R,NDIM>& x0)const;

    /**
     * @brief interpolation-kind includes periodic boundary-conditions.
     */
    inline bool periodic_BC_(void)const;
    
    /**
     * @brief index-tuple n of neighbor from base index-tuple n0 and packed-offset index nn.
     *
     * (neighbors are relative to lowest-index corner of hypercube containing x_i:
     *    offset bits of "nn", by dimension as:
     *      0: no offset 
     *      1: +1 offset)
     * @return  true if neighbor exists (note: if periodic B.C., neighbor always exists).
     */
    inline bool neighbor_(const ntuple<size_t,NDIM>& n0, size_t nn, ntuple<size_t,NDIM>& n)const;

    /**
     * @brief linear basis functions
     *
     *   @param[in] dx  N-dimensional offset from lowest-index corner of hypercube (i.e. containing x_i)
     *   @param[in] nn  index of basis-function: corresponds to packed neighbor offset-index  
     */
    inline S B_(const ntuple<R,NDIM>& dx, size_t nn)const;


  public:
  
  inline value_type operator()(const ntuple<R,NDIM>& x_i)const;
  
  void apply(const ntuple<R,NDIM>& x_i, value_type& y_i)const;  
  
#if !defined(__INTEL_COMPILER) && !defined(__PGI)
  /**
   * @brief interpolation over N-dimensional matrix of co-ordinates.
   *
   *   @tparam  U1  a linear container with value_type compatible with  ntuple<R,NDIM>
   *   @tparam  U2  a linear container with value_type compatible with linear_interpolator<U, NDIM>::value_type
   *   @param[in] x_i coordinate values
   *   @param[out] y_i  interpolated data values
   */
  template <class U1, class U2>
  void apply(const U1& x_i, U2& y_i)const;
#else
    // INTEL and PGI require body-definition at point-of-declaration:

    /**
     * @brief interpolation over N-dimensional matrix of co-ordinates.
     *
     *   @tparam  U1  a linear container with value_type compatible with  ntuple<R,NDIM>
     *   @tparam  U2  a linear container with value_type compatible with linear_interpolator<U, NDIM>::value_type
     *   @param[in] x_i coordinate values
     *   @param[out] y_i  interpolated data values
     */
    template <class U1, class U2>
    void apply(const U1& x_i, U2& y_i)const
    {
      // resize dest only if destination is _not_ a reference type:
      // see note at "simple_object_base::is" regarding type_info comparison by-name, to enable cross-module use with current GNU C++ ABI:
      if (0 == strcmp(typeid(typename gmm::linalg_traits<U2>::is_reference).name(), typeid(gmm::linalg_false).name()))
        gmm::linalg_traits<U2>::resize(dest, gmm::vect_size(x_i)); 
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
  
  linear_interpolator
    (const ntuple_interval<R,NDIM>& domain, const U& y, const ntuple<size_t,NDIM>& shape, interp_kind e);
    
  linear_interpolator
    (const ntuple_interval<R,NDIM>& domain, U& y, const ntuple<size_t,NDIM>& shape, interp_kind e, bool transfer_ownership);
}; 


} // namespace linalg

#include "ND_matrix_inline.h"

#if !defined(EXCLUDE_TEMPLATE_BODIES)
#include "ND_matrix_template.h"
#endif

#endif // __ND_matrix__h
