#if !defined(__linalgUtil__h)
#define __linalgUtil__h

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


#define __specialize_POD_binary__


// $Source: /usr/data0/leipzig_work/tmat_cvs/src/linalgUtil.h,v $

namespace linalg{

// most of the following are required in order to move "linalgUtil" methods to "namespace linalg":
using number::epsilon;
using number::zero;
using number::one;
using number::one_i;
using number::pi;
using number::integer;
using number::ratio;
using number::conv;
using number::real;
using number::imag;
using number::conj;
using number::fabs;
using number::max;
using number::min;
using number::random;
using number::numberTraits;

// return indices (as iterator offsets) of minimum and maximum elements of iterable container U:
template <class IT>
size_t min_index(IT itStart, IT itEnd);

template <class IT>
size_t max_index(IT itStart, IT itEnd);


/**
   @cond full_docs
   @brief STL-compatible functors implementing by-index comparison.
   
 (assume the appropriate definition of U::operator[](size_type)  )
 */
 
 
template <class U>
class less_dref{
    const U& u_;
  public:
  
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] < u_[n2]; }
    
    less_dref(const U& u)
      : u_(u)
    { }
};

template <class U>
class less_equal_dref{
    const U& u_;
  public: 
     
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] <= u_[n2]; }
    
    less_equal_dref(const U& u)
      : u_(u)
    { }
};

template <class U>
class greater_dref{
    const U& u_;
  public:   
   
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] > u_[n2]; }
    
    greater_dref(const U& u)
      : u_(u)
    { }
};

template <class U>
class greater_equal_dref{
    const U& u_;
  public:
  
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] >= u_[n2]; }
    
    greater_equal_dref(const U& u)
      : u_(u)
    { }
};

template <class U>
class equal_dref{
    const U& u_;
  public: 
     
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] == u_[n2]; }
    
    equal_dref(const U& u)
      : u_(u)
    { }
};

template <class U>
class not_equal_dref{
    const U& u_;
  public:
  
    bool operator()(size_t n1, size_t n2)
    { return u_[n1] != u_[n2]; }
    
    not_equal_dref(const U& u)
      : u_(u)
    { }
};

//! @endcond

//! index sort:
template <class U, class N>
std::vector<N> argsort(const U& u);

template <class U, class N>
void argsort(const U& u, std::vector<N>& dest);


//! product of vector elements (i.e. \f$\Pi v_k: k \in [beg, end) \f$):
template <class V>
typename gmm::linalg_traits<V>::value_type product(const V& v, size_t beg_=0, size_t end_=static_cast<size_t>(-1));

// linear index from N-dimensional row-major indices in shape (i.e. dimension limits):
//   here I assume indices are "size_t"
size_t row_major_index(const std::vector<size_t>& indices, const std::vector<size_t>& shape);

// N-dimensional row-major indices (constrained by shape) from linear index:  
void inverse_row_major_index(size_t n, const std::vector<size_t>& shape, std::vector<size_t>& dest);
std::vector<size_t> inverse_row_major_index(size_t n, const std::vector<size_t>& shape);

// increment a specified dimension of a set of indices subject to an interval-list constraint:
// (false return value => increment went out-of-range;
//    specifically: out-of-range return leaves highest dimension index at end-of-range, 
//    others at start-of-range positions)
bool increment_row_major_indices(const std::vector<size_t>& src, size_t ndim, 
                                 const std::vector<std::pair<size_t,size_t> >& intervals, 
                                 std::vector<size_t>& dest);


// gmm-interfaced test utilities (test value not by-reference => not efficient for mere types):
template <class U>
inline bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&));

template <class U>
inline bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&));

//  ... for C int-style bool test functions:
template <class U>
inline bool any(const U& u, int (*test_func) (const typename gmm::linalg_traits<U>::value_type&));

template <class U>
inline bool all(const U& u, int (*test_func) (const typename gmm::linalg_traits<U>::value_type&));


template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_vector);

template <class U>
inline bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_matrix);

template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::col_major);

template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::row_major);

template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_vector);

template <class U>
inline bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_matrix);

template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::col_major);

template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::row_major);


// in-place check for gmm-interfaced matrix or vector types U1, U2:
template <class U1, class U2>
bool sameOrigin(const U1& u1, const U2& u2);
	                           												

// convert for std::vector: same-type fall back:
template <class T>
void conv( std::vector<T>& vT1, const std::vector<T>& vT2);

// real for std::vector: same-type fall back:
template <class T>
void real( std::vector<T>& vT1, const std::vector<T>& vT2);

// imag for std::vector: same-type fall back:
template <class T>
void imag( std::vector<T>& vT1, const std::vector<T>& vT2);

template <class T, class U>
void conv( std::vector<T>& vT, const std::vector<U>& vU);

// real for std::vector:
template <class T, class U>
void real( std::vector<T>& vT, const std::vector<U>& vU);

// imag for std::vector:
template <class T, class U>
void imag( std::vector<T>& vT, const std::vector<U>& vU);

} // namespace linalg


namespace commUtil{

// binary read and write for std::vector:

template<class T>
bool writeBinary(abstractCommHandle *fp, const std::vector<T>& V);

template<class T>
bool readBinary(abstractCommHandle *fp, std::vector<T>& V );

template <class T>
size_t binarySize(const std::vector<T>& v);


// binary read and write for std::list:
template<class T>
bool writeBinary(abstractCommHandle *fp, const std::list<T>& U);

template<class T>
bool readBinary(abstractCommHandle *fp, std::list<T>& U);

template <class T>
size_t binarySize(const std::list<T>& u);


#if defined(__specialize_POD_binary__)
// specializations for contiguous POD types:
template<>
bool writeBinary<double>(abstractCommHandle *fp, const std::vector<double>& V);

template<>
bool readBinary<double>(abstractCommHandle *fp, std::vector<double>& V );


template<>
bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const std::vector<std::complex<double> >& V);

template<>
bool readBinary<std::complex<double> >(abstractCommHandle *fp, std::vector<std::complex<double> >& V );


template<>
bool writeBinary<long>(abstractCommHandle *fp, const std::vector<long>& V);

template<>
bool readBinary<long>(abstractCommHandle *fp, std::vector<long>& V );


template<>
bool writeBinary<size_t>(abstractCommHandle *fp, const std::vector<size_t>& V);

template<>
bool readBinary<size_t>(abstractCommHandle *fp, std::vector<size_t>& V );
#endif
} // namespace commUtil


namespace linalg{

// Do _not_ resize destination: allows to work correctly for vector references:
template <class V1, class V2, class V3>
void elementwiseDiv( const V1& v1, const V2& v2, V3& v3 )throw();

// Do _not_ resize destination: allows to work correctly for vector references:
template <class V1, class V2, class V3>
void elementwiseMult( const V1& v1, const V2& v2, V3& v3 )throw();

#if 0
template <class T1, class T2, class T3>
class multiplies{
public:
  inline T3 operator()(const T1& t1, const T2& t2)const
	  { return t1*t2; }
};

template <class T1, class T2, class T3>
class divides{
public:
  inline T3 operator()(const T1& t1, const T2& t2)const
	  { return t1/t2; }
};
#else
template <class T1, class T2, class T3>
T3 multiplies(const T1&, const T2&);
template <class T1, class T2, class T3>
T3 divides(const T1&, const T2&);
#endif

// Resize destination: will _not_ work correctly for vector references:
template <class V1, class V2>
void diff( const V1& v1, V2& v2, size_t nDiff=1 );

// for gmm::dense_matrix<T> only:
template <class T>
inline void transpose(const gmm::dense_matrix<T>& src, gmm::dense_matrix<T>& dest);

// Intended for vector or matrix "U" with implemented gmm::linalg_traits< U >::linalg_type definition.

// U <- U + delta
template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta );
template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta, gmm::abstract_vector );
template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta, gmm::abstract_matrix );

template <class U>
long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_=epsilon< typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type >() );
template <class U>
long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_, gmm::abstract_vector );
template <class U>
long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_, gmm::abstract_matrix );

template <class U>
void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val );
template <class U>
void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val, gmm::abstract_vector );
template <class U>
void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val, gmm::abstract_matrix );

// an extension for std::iota:
//   iota: < a vector of integers >
//   alternative iota: < a vector of small differences >
//
template <class T, class U>
void iota( typename U::iterator itStart, typename U::iterator itEnd, const T& initValue, const T& dt = one<T>() );

/*
 * These are convenient initialization routines, from simple arrays.
 * The number of elements copied corresponds to the existing _allocated_ destination vector or matrix.
 */
 
template <class U>
void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* val );
template <class U>
void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* val, gmm::abstract_vector );
template <class U>
void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* val, gmm::abstract_matrix );

// make unit vectors:
template <class V>
void e1( V& v );

template <class V>
void e2( V& v );

template <class V>
void e3( V& v );

// convert for gmm::dense_matrix:
template <class T, class U>
void conv( gmm::dense_matrix<T>& aT, const gmm::dense_matrix<U>& aU);

} // namespace linalg


namespace commUtil{

// binary read and write for gmm::dense_matrix:
// note: these immediately hand-off to the "namespace gmm" methods:

template <class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::dense_matrix<T>& M);

template <class T>
inline bool readBinary(abstractCommHandle *fp, gmm::dense_matrix<T>& M);

template <class T>
inline size_t binarySize(const gmm::dense_matrix<T>& M);


template <class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<T> >& M);

template <class T>
inline bool readBinary(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<T> >& M);

template <class T>
inline size_t binarySize(const gmm::row_matrix< gmm::slvector<T> >& M);



#if defined(__specialize_POD_binary__)
// specializations for contiguous POD types:
template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::dense_matrix<double>& M);

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::dense_matrix<double>& M);


template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::dense_matrix<std::complex<double> >& M);

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::dense_matrix<std::complex<double> >& M);


// The following are especially for factoredPropagator. It is very difficult to specialize 
//   these in a generic sense because of possibility of sub-vectors and sub-matrices with storage_type == abstract_skyline.
template <class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::slvector<T>& V);

template <class T>
inline bool readBinary(abstractCommHandle *fp, gmm::slvector<T>& V);

template <class T>
inline size_t binarySize(const gmm::slvector<T>& V);


template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::slvector<double>& V);

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::slvector<double>& V);


template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::slvector<std::complex<double> >& V);

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::slvector<std::complex<double> >& V);


template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<double> >& M);

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<double> >& M);

template <>
inline size_t binarySize<double>(const gmm::row_matrix< gmm::slvector<double> >& M);


template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<std::complex<double> > >& M);

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<std::complex<double> > >& M );

template <>
inline size_t binarySize<std::complex<double> >(const gmm::row_matrix< gmm::slvector<std::complex<double> > >& M);

#endif

// binary read and write for std::map:

template <class KEY, class DATA>
bool writeBinary(abstractCommHandle *fp, const std::map<KEY,DATA>& M);

template <class KEY, class DATA>
bool readBinary(abstractCommHandle *fp,  std::map<KEY,DATA>& M);

template <class KEY, class DATA>
size_t binarySize(const std::map<KEY,DATA>& M);


// binary read and write for std::pair:

template <class T1, class T2>
bool writeBinary(abstractCommHandle *fp, const std::pair<T1,T2>& p);

template <class T1, class T2>
bool readBinary(abstractCommHandle *fp, std::pair<T1,T2>& p);

template <class T1, class T2>
size_t binarySize(const std::pair<T1,T2>& p);

} // namespace commUtil


namespace linalg{

template <class T>
void saveVectorData(const std::vector<T>& V, const std::string& sFileName) throw();

template <class T>
void loadVectorData(std::vector<T>& V, const std::string& sFileName) throw();

#if 0 // omit for now
template <class domainType, class rangeType>
class unaryFunctionWrapper: public std::unary_function<const domainType&, rangeType>{
  public:
    unaryFunctionWrapper( rangeType (*foop)(const domainType&) )
      : foop_(foop)
      {}
    rangeType operator()(domainType& x) 
    { return (*foop_)(x); }  
  private:
    rangeType (*foop_)(const domainType&);
};
#endif

#if 1 // compiler work-around (gcc 4.1.1 not recognizing appropriate overload of "fabs", "real", "imag"); moved from numericalConstants.h
template <class T1, class T2>
void conv( std::complex<T1>& t1, const T2& t2); // T2 -> std::complex<T1>

template <class T1, class T2>
void conv( T1& t1, const std::complex<T2>& t2); // std::complex<T2> -> T1

template <class T1, class T2>
void conv( std::complex<T1>& t1, const std::complex<T2>& t2); // std::complex<T2> -> std::complex<T1>

template <class T>
const T max_abs(const T& a, const T& b);

template <class T>
const T min_abs(const T& a, const T& b);

template <class T>
const T max_abs(const T& a, const T& b, const T& c);

template <class T>
const T min_abs(const T& a, const T& b, const T& c);
#endif

// for 3D vectors:
// (sphere2cart and cart2sphere moved _here_ because of methods using std::vector).
#if 1 // moved from numericalConstants.h
// Note: here theta is zenith angle, phi is azimuth.
template <class T1, class T2>
void cart2sphere(const T1& x, const T1& y, const T1& z,
                 T2& r, T2& theta, T2& phi);
                 
template <class T1, class T2>
void sphere2cart(const T1& r, const T1& theta, const T1& phi,
                 T2& x, T2& y, T2& z);                 
#endif


template <class T1, class T2>
inline void cart2sphere(const std::vector<T1>& vSrc, 
                        std::vector<T2>& vDest);
                 
template <class T1, class T2>
inline void sphere2cart(const std::vector<T1>& vSrc,
                        std::vector<T2>& vDest);                

// (these methods not near sphere2cart and cart2sphere because they use ntuple).
template <class T1, class T2>
inline void cart2sphere(const ntuple<T1,3>& vSrc, 
                        ntuple<T2,3>& vDest);
                 
template <class T1, class T2>
inline void sphere2cart(const ntuple<T1,3>& vSrc,
                        ntuple<T2,3>& vDest);  
					
// Conversion of vector components:
// (these methods not near sphere2cart and cart2sphere because they use gmm).
template <class T>
void cartVect2sphere(const T& x, const T& y, const T& z,
                     const T& vx, const T& vy, const T& vz,
                     T& r, T& theta, T& phi,
										 T& vr, T& vtheta, T& vphi);
                 
template <class T>
void sphereVect2cart(const T& r, const T& theta, const T& phi,
                     const T& vr, const T& vtheta, const T& vphi,
                     T& x, T& y, T& z,
										 T& vx, T& vy, T& vz);                 

// for 3D vectors:
template <class T>
inline void cartVect2sphere(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV);
                 
template <class T>
inline void sphereVect2cart(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV);    
#if 0
// not implemented yet
//   requires explicit instantiation of a gmm::linalg_traits< ntuple<T,DIM> >

template <class T>
inline void cartVect2sphere(const ntuple<T,3>& vSrcP, const ntuple<T,3>& vSrcV, 
                            ntuple<T,3>& vDestP, ntuple<T,3>& vDestV);
                 
template <class T>
inline void sphereVect2cart(const ntuple<T,3>& vSrcP, const ntuple<T,3>& vSrcV, 
                            ntuple<T,3>& vDestP, ntuple<T,3>& vDestV);   
#endif																										


enum CSYS_KIND {NO_CSYS, CARTESIAN, CYLINDRICAL, SPHERICAL};



#if 1 // ========================== move from numericalFunctor.h: ============================================================
// ****************** Utility methods (in addition to those provided by class T) ************************:

  /**
   * @brief Solve the quadratic formula:
   * @tparam V1  container type of value_type real or complex
   * @tparam V2  container type of value_type complex (or value_type real \f$ \Leftrightarrow \f$ roots are somehow constrained to be real)
   * @param[in] vP  vector of polynomial coefficients
   * @param[out] vRoot  vector of roots
   * @param[out] D  discriminant
   */
  template <class V1, class V2>
  void quadraticFormula(const V1& vP, V2& vRoot, typename V1::value_type& D); 

  /**
   * @brief Calculate the discriminant for the quadratic formula:
   * @tparam V1  container type of value_type real or complex
   * @param[in] vP  vector of polynomial coefficients
   * @param[out] D  discriminant
   */
  template <class V1>
  void quadraticDiscriminant(const V1& vP, typename V1::value_type& D); 

  template <class T>    
  void cubicFormula( const std::vector<T>& vP, std::vector< std::complex<T> >& vRoot );
  template <class T>
  void cubicFormula( const std::vector< std::complex<T> >& vP, std::vector< std::complex<T> >& vRoot );

  template <class T>  
  T polyval( const std::vector<T>& vP, const T& arg ); 
  template <class T>
  std::complex<T> polyval( const std::vector< std::complex<T> >& vP, const std::complex<T>& arg );   

  /** 
    * @brief convert vector [v_r, v_theta, v_phi](r, theta, phi) to [v_minus, v_0, v_plus] in helical basis.
    * optionally use u = cos(theta) as polar argument.
    */
  template <class T>
  void sphericalToHelical( const T& v_r, const T& v_theta, const T& v_phi,
                           const T& r, const T& theta, const T& phi,
                           T& v_minus, T& v_0, T& v_plus, bool cosTheta=false )throw();    

  /** 
    * @brief convert vector [v_r, v_theta, v_phi](r, theta, phi) to [v_minus, v_0, v_plus] in helical basis.
    * optionally use u = cos(theta) as polar argument; special usage: additional parsing by m-value.
    */      
  template <class T>
  void sphericalToHelical( const T& v_r, const T& v_theta, const T& v_phi,
                           const T& r, const T& theta, const long& m,
                           T& v_minus, T& v_0, T& v_plus, bool cosTheta=false )throw();

// *******************************************************************************************************  
#endif // ========================== end, move from numericalFunctor.h: ============================================================


} // namespace linalg


namespace commUtil{
inline bool writeBinary(commUtil::abstractCommHandle* fp, const linalg::CSYS_KIND& e);

inline bool readBinary(commUtil::abstractCommHandle* fp, linalg::CSYS_KIND& e);

inline void write(std::ostream& os, const linalg::CSYS_KIND& e)throw();

inline void read(std::istream& is, linalg::CSYS_KIND& e)throw();

inline std::ostream& operator<<(std::ostream& os, const linalg::CSYS_KIND& e);

inline std::istream& operator>>(std::istream& is, linalg::CSYS_KIND& e);
} // namespace commUtil


namespace linalg{

// return an instance of a _uniform_ distribution for the specified co-ordinate (with the specified symmetry):
template <class R>
R randomCoord(CSYS_KIND eCSYS, size_t dimOffset)throw();

} // namespace linalg

// writeBinary, readBinary and binarySize for interfaced "gmm" classes.
// These are left in namespace "gmm" so that they won't conflict with the STL container methods;
//   inline methods then provide specific instantiations in "namespace commUtil" for the most widely-used gmm types.
namespace gmm{

template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u) throw();


template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector) throw();

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector, abstract_dense) throw();

// WARNING: the following almost certainly does _not_ work for sub-index types:
template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector, abstract_skyline) throw();



template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_matrix) throw();


template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, row_major) throw();

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, col_major) throw();


template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u) throw();

// allow "const reference" to temporary:
template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, const U& u) throw();


template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_false) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_modifiable) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_const) throw();


template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, abstract_dense) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, abstract_skyline) throw();



template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_false) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_modifiable) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_const) throw();


template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, row_major) throw();

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, col_major) throw();



// in consideration of variable-precision number types, the following "binarySize" method requires an instance of class U:
template <class U>
inline size_t binarySize(const U& u);

template <class U>
inline size_t binarySize(const U& u, abstract_vector);

template <class U>
inline size_t binarySize(const U& u, abstract_vector, abstract_dense);

template <class U>
inline size_t binarySize(const U& u, abstract_vector, abstract_skyline);

template <class U>
inline size_t binarySize(const U& u, abstract_matrix);

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, row_major, abstract_dense);

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, col_major, abstract_dense);

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, row_major, abstract_skyline);

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, col_major, abstract_skyline);

  
#if defined(__specialize_POD_binary__)
// specializations for contiguous POD types.
// (these correspond to existing specializations for writeBinary and readBinary).

#if 0 // no specialization required, size on media same as for non-POD case:
template<>
size_t binarySize<std::vector<double> >(const std::vector<double>& V);

template<>
size_t binarySize<std::vector<std::complex<double> > >(const std::vector<std::complex<double> >& V);

template<>
size_t binarySize<std::vector<long> >(const std::vector<long>& V);

template<>
size_t binarySize<std::vector<size_t> >(const std::vector<size_t>& V);
#endif

template <>
inline size_t binarySize<gmm::dense_matrix<double> >(const gmm::dense_matrix<double>& M);


template <>
inline size_t binarySize<gmm::dense_matrix<std::complex<double> > >(const gmm::dense_matrix<std::complex<double> >& M);

// The following are especially for factoredPropagator. It is very difficult to specialize 
//   these in a generic sense because of possibility of sub-vectors and sub-matrices with storage_type == abstract_skyline.
template <>
inline size_t binarySize<gmm::slvector<double> >(const gmm::slvector<double>& u);

template <>
inline size_t binarySize<gmm::slvector<std::complex<double> > >(const gmm::slvector<std::complex<double> >& u); 

#if 0 // no specialization required, size on media same as for non-POD case:
template <>
inline size_t binarySize<gmm::row_matrix< gmm::slvector<double> > >(const gmm::row_matrix< gmm::slvector<double> >& M);

template <>
inline size_t binarySize<gmm::row_matrix< gmm::slvector<std::complex<double> > > >(const gmm::row_matrix< gmm::slvector<std::complex<double> > >& M);
#endif

#endif

} // namespace gmm


#include "linalgUtil_inline.h"
#include "linalgUtil_template.h"

#endif // __linalgUtil__h
