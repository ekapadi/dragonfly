#if !defined(__linalgUtil_template__h)
#define __linalgUtil_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/linalgUtil_template.h,v $

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

// return indices (as iterator offsets) of minimum and maximum elements of iterable container U:
template <class IT>
size_t min_index(IT itStart, IT itEnd)
{
  IT itU = std::min_element(itStart, itEnd);
  if (itU == itEnd)
    throw std::runtime_error("min_index: empty range has no minimum");
  return static_cast<size_t>(itU - itStart); 
}

template <class IT>
size_t max_index(IT itStart, IT itEnd)
{
//  typedef typename std::iterator_traits<IT>::size_type size_type;
  
  IT itU = std::max_element(itStart, itEnd);
  if (itU == itEnd)
    throw std::runtime_error("min_index: empty range has no minimum");
  return static_cast<size_t>(itU - itStart); 
}


// index sort:
template <class U, class N>
inline std::vector<N> argsort(const U& u)
{
  std::vector<N> dest;
  argsort(u, dest);
  return dest;
}

template <class U, class N>
void argsort(const U& u, std::vector<N>& dest)
{
  dest.resize(u.size());
  _STL_EXT_NAMESPACE_::iota(dest.begin(), dest.end(), zero<N>());
  std::sort(dest.begin(), dest.end(), less_dref<U>(u));
}


// product of vector elements (i.e. \Pi):
template <class V>
typename gmm::linalg_traits<V>::value_type product(const V& v, size_t beg_, size_t end_)
{
  typedef typename gmm::linalg_traits<V>::value_type N;
  assert(static_cast<size_t>(-1) == beg_ || beg_ <= gmm::vect_size(v));
  assert(static_cast<size_t>(-1) == end_ || end_ <= gmm::vect_size(v));
  
  N val(one<N>());
  for(typename gmm::linalg_traits<V>::const_iterator 
    itV = (static_cast<size_t>(-1) == beg_? gmm::vect_const_end(v): gmm::vect_const_begin(v) + beg_),
    itVEnd = (static_cast<size_t>(-1) == end_? gmm::vect_const_end(v): gmm::vect_const_begin(v) + end_);
      itV != itVEnd;
      ++itV)
    val *= *itV;
  return val;
}

// gmm-interfaced test utilities:
template <class U>
inline bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&))
{ return any(u, test_func, typename gmm::linalg_traits<U>::linalg_type() ); }

template <class U>
inline bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&))
{ return all(u, test_func, typename gmm::linalg_traits<U>::linalg_type() ); }


//  ... for C int-style bool test functions:
template <class U>
inline bool any(const U& u, int (*test_func) (const typename gmm::linalg_traits<U>::value_type&))
{
 typedef bool (*test_func_type) (const typename gmm::linalg_traits<U>::value_type&);
 return any(u, reinterpret_cast<test_func_type>(test_func));
}

template <class U>
inline bool all(const U& u, int (*test_func) (const typename gmm::linalg_traits<U>::value_type&))
{
 typedef bool (*test_func_type) (const typename gmm::linalg_traits<U>::value_type&);
 return all(u, reinterpret_cast<test_func_type>(test_func));
}


template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_vector)
{
  bool val(false);
  for(typename gmm::linalg_traits<U>::const_iterator itU = gmm::vect_begin(u), itUEnd = gmm::vect_end(u);
      itU != itUEnd;
      ++itU)
    if (test_func(*itU)){
      val = true;
      break;
    }
  return val;   
}

template <class U>
inline bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_matrix)
{ 
  return any(u, test_func, gmm::abstract_matrix(), 
    typename gmm::principal_orientation_type<typename gmm::linalg_traits<U>::sub_orientation>::potype()); 
}

template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::col_major)
{
 bool val(false);
 for(typename gmm::linalg_traits<U>::const_col_iterator itCols = gmm::mat_col_begin(u), itColsEnd = gmm::mat_col_end(u);
     itCols != itColsEnd;
     ++itCols){

   typename gmm::linalg_traits<U>::const_sub_col_type vCol( gmm::linalg_traits<U>::col(itCols) );
   if (any(vCol, test_func)){
     val = true;
     break;   
   }
 }
 
 return val;
}    

template <class U>
bool any(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::row_major)
{
 bool val(false);
 for(typename gmm::linalg_traits<U>::const_row_iterator itRows = gmm::mat_row_begin(u), itRowsEnd = gmm::mat_row_end(u);
     itRows != itRowsEnd;
     ++itRows){

   typename gmm::linalg_traits<U>::const_sub_row_type vRow( gmm::linalg_traits<U>::row(itRows) );
   if (any(vRow, test_func)){
     val = true;
     break;   
   }
 }
 
 return val;
} 

template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_vector)
{
  bool val(true);
  for(typename gmm::linalg_traits<U>::const_iterator itU = gmm::vect_begin(u), itUEnd = gmm::vect_end(u);
      itU != itUEnd;
      ++itU)
    if (!test_func(*itU)){
      val = false;
      break;
    }
  return val;   
}

template <class U>
inline bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), gmm::abstract_matrix)
{ 
  return all(u, test_func, gmm::abstract_matrix(), 
    typename gmm::principal_orientation_type<typename gmm::linalg_traits<U>::sub_orientation>::potype()); 
}

template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type& t), 
           gmm::abstract_matrix, gmm::col_major)
{
 bool val(true);
 for(typename gmm::linalg_traits<U>::const_col_iterator itCols = gmm::mat_col_begin(u), itColsEnd = gmm::mat_col_end(u);
     itCols != itColsEnd;
     ++itCols){

   typename gmm::linalg_traits<U>::const_sub_col_type vCol( gmm::linalg_traits<U>::col(itCols) );
   if (!all(vCol, test_func)){
     val = false;
     break;   
   }
 }
 
 return val;
}
  
template <class U>
bool all(const U& u, bool (*test_func) (const typename gmm::linalg_traits<U>::value_type&), 
           gmm::abstract_matrix, gmm::row_major)
{
 bool val(true);
 for(typename gmm::linalg_traits<U>::const_row_iterator itRows = gmm::mat_row_begin(u), itRowsEnd = gmm::mat_row_end(u);
     itRows != itRowsEnd;
     ++itRows){

   typename gmm::linalg_traits<U>::const_sub_row_type vRow( gmm::linalg_traits<U>::row(itRows) );
   if (!all(vRow, test_func)){
     val = false;
     break;   
   }
 }
 
 return val;
} 



// in-place check for gmm-interfaced matrix or vector types U1, U2:
template <class U1, class U2>
inline bool sameOrigin(const U1& u1, const U2& u2)
{
 return
   reinterpret_cast<const void*>( gmm::linalg_traits<U1>::origin(u1) )
	 == reinterpret_cast<const void*>( gmm::linalg_traits<U2>::origin(u2) );
}


// Non-binary insertion and extraction for std::vector
//
//   These methods assume the existence of
//   insertion and extraction methods for type T
//   in the current scope.
//
#if 0
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& V)
{
 // Put it all on one line:
 os<<V.size()<<" ";
 for (typename std::vector<T>::const_iterator itV = V.begin(), itVEnd = V.end();
 
      itV != itVEnd;
      
      ++itV){
   os<<*itV<<" ";
   if (os.fail()) break;
 }
 return os;   
}

template<class T>
std::istream& operator>>(std::istream& is, std::vector<T>& V )
{
 typename std::vector<T>::size_type nSize;
 is>>nSize;
 V.clear();
 V.resize(nSize);
 for (typename std::vector<T>::iterator itV = V.begin(), itVEnd = V.end();
 
      itV != itVEnd;
      
      ++itV){
   is>>*itV;
   if (is.fail()) break;
 }
 return is;   
}
#endif


// convert for std::vector: same-type fall back:
template <class T>
void conv( std::vector<T>& vT, const std::vector<T>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 // ( kind of strange for "conv", but may occur in generic code )
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<T>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT)
   *itT = *itU;  
}

// real for std::vector: same-type fall back:
template <class T>
void real( std::vector<T>& vT, const std::vector<T>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<T>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT) 
   *itT = real(*itU);  
}

// imag for std::vector: same-type fall back:
template <class T>
void imag( std::vector<T>& vT, const std::vector<T>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<T>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT)
   *itT = imag(*itU);  
}

// convert for std::vector:
template <class T, class U>
void conv( std::vector<T>& vT, const std::vector<U>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 // ( kind of strange for "conv", but may occur in generic code )
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<U>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT)
   conv( *itT, *itU );  
}

// real for std::vector:
template <class T, class U>
void real( std::vector<T>& vT, const std::vector<U>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<U>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT) 
   conv( *itT, real(*itU) );  
}

// imag for std::vector:
template <class T, class U>
void imag( std::vector<T>& vT, const std::vector<U>& vU)
{
#if 0
 vT.clear();
#endif
 // _allow_ in-place:
 if (vT.size() != vU.size())
   vT.resize( vU.size(), zero<T>() );
 
 
 typename std::vector<T>::iterator itT = vT.begin();
 for(typename std::vector<U>::const_iterator itU = vU.begin(), itUEnd = vU.end();
     itU != itUEnd;
     ++itU,
     ++itT)
   conv( *itT, imag(*itU) );  
}

} // namespace linalg


namespace commUtil{

// Binary read and write for std::vector
//
//   These methods assume the existence of
//   bool writeBinary(abstractCommHandle *fp, const T& val);
//   bool readBinary(abstractCommHandle *fp, T& val);
//   in the current scope.
//

template<class T>
bool writeBinary(abstractCommHandle *fp, const std::vector<T>& V)
{
   bool rVal = true;
   
   // write size:
   size_t nSize = V.size();
   rVal = (rVal && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     for (typename std::vector<T>::const_iterator itV = V.begin(), itEnd = V.end();
     
          itV != itEnd;
          
          ++itV){
       rVal = (rVal && writeBinary(fp, *itV));
       if (!rVal)
         break;
     }
        
   return rVal;
}


template<class T>
bool readBinary(abstractCommHandle *fp, std::vector<T>& V)
{
   bool rVal = true;
   
   // read size
   size_t nSize;
   rVal = (rVal && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   V.resize( nSize ); 

   if (0 < nSize)
     for (typename std::vector<T>::iterator itV = V.begin(), itEnd = V.end();
     
          itV != itEnd;
          
          ++itV){
       rVal = (rVal && readBinary(fp, *itV));
       if (!rVal)
         break;
     }

   return rVal;
}


template<class T>
size_t binarySize(const std::vector<T>& V)
{
   size_t val(0);
   
   // write size:
   val += sizeof(size_t);
 
   // generic form: works correctly if vector elements do not have the same "binarySize":
   for (typename std::vector<T>::const_iterator itV = V.begin(), itEnd = V.end();     
        itV != itEnd;
        ++itV)
     val += binarySize(*itV);
        
   return val;
}

template <class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::slvector<T>& V)
{
 return gmm::writeBinary(fp, V);
}

template <class T>
inline bool readBinary(abstractCommHandle *fp, gmm::slvector<T>& V)
{
 return gmm::readBinary(fp, V);
}

template <class T>
inline size_t binarySize(const gmm::slvector<T>& V)
{
 return gmm::binarySize(V);
}

// binary read and write for std::list:
template<class T>
bool writeBinary(abstractCommHandle *fp, const std::list<T>& U)
{
   bool status = true;
   
   // write size:
   size_t nSize = U.size();
   status = (status && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     for (typename std::list<T>::const_iterator itU = U.begin(), itUEnd = U.end();
     
          status && (itU != itUEnd);
          
          ++itU)
       status = (status && writeBinary(fp, *itU));
        
   return status;
}


template<class T>
bool readBinary(abstractCommHandle *fp, std::list<T>& U)
{
   bool status = true;
   
   // read size
   size_t nSize(0);
   status = (status && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   // U.reserve(nSize); 

   T tmp /* (zero<T>()) */;  // do _not_ assume "T" is a numerical type
   for (size_t n = 0; status && (n < nSize); ++n){
     status = (status && readBinary(fp, tmp));
     if (status)
       U.push_back(tmp);
   }

   return status;
}

template<class T>
size_t binarySize(const std::list<T>& U)
{
   size_t val(0);
   
   // write size:
   val += sizeof(size_t);

   // generic form: works correctly if list elements don't have matching "binarySize":
   for (typename std::list<T>::const_iterator itU = U.begin(), itUEnd = U.end();
        itU != itUEnd;
        ++itU)
     val += binarySize(*itU);
        
   return val;
}

} // namespace commUtil


namespace linalg{

// Do _not_ resize destination: allows to work correctly for vector references:
template <class V1, class V2, class V3>
void elementwiseDiv( const V1& v1, const V2& v2, V3& v3 )
{
#if 0
 // Initialization may be required by some arbitrary precision implementations:
 if ( ( (void*)&v2 != (void*)&v3 ) // allow in-place
      && ( (void*)&v1 != (void*)&v3 ) )
   gmm::clear(v3); // initialize elements to zero.
#endif
 if( (gmm::vect_size(v1) != gmm::vect_size(v2) ) 
     || (gmm::vect_size(v2) != gmm::vect_size(v3) ))
   throw std::runtime_error("elementwiseMult: vector sizes don't match");
#if 0   
 std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                std::multiplies<typename gmm::linalg_traits<V1>::value_type>()); // value_type's must match
#else
 std::transform(gmm::vect_const_begin(v1), gmm::vect_const_end(v1), gmm::vect_const_begin(v2), gmm::vect_begin(v3),
                divides<typename gmm::linalg_traits<V1>::value_type, 
								        typename gmm::linalg_traits<V2>::value_type,
												typename gmm::linalg_traits<V3>::value_type>); 
#endif                
}

// Do _not_ resize destination: allows to work correctly for vector references:
template <class V1, class V2, class V3>
void elementwiseMult( const V1& v1, const V2& v2, V3& v3 )
{
#if 0
 // Initialization may be required by some arbitrary precision implementations:
 if ( ( (void*)&v2 != (void*)&v3 ) // allow in-place
      && ( (void*)&v1 != (void*)&v3 ) )
   gmm::clear(v3); // initialize elements to zero.
#endif
 if( (gmm::vect_size(v1) != gmm::vect_size(v2) ) 
     || (gmm::vect_size(v2) != gmm::vect_size(v3) ))
   throw std::runtime_error("elementwiseMult: vector sizes don't match");
#if 0   
 std::transform(v1.begin(), v1.end(), v2.begin(), v3.begin(),
                std::multiplies<typename gmm::linalg_traits<V1>::value_type>()); // value_type's must match
#else
 std::transform(gmm::vect_const_begin(v1), gmm::vect_const_end(v1), gmm::vect_const_begin(v2), gmm::vect_begin(v3),
                multiplies<typename gmm::linalg_traits<V1>::value_type, 
								           typename gmm::linalg_traits<V2>::value_type,
													 typename gmm::linalg_traits<V3>::value_type>); 
#endif                
}

#if 1
template <class T1, class T2, class T3>
T3 multiplies(const T1& t1, const T2& t2)
{ return t1*t2; }

template <class T1, class T2, class T3>
T3 divides(const T1& t1, const T2& t2)
{ return t1/t2; }
#endif

// for gmm::dense_matrix<T> only:
template <class T>
inline void transpose(const gmm::dense_matrix<T>& src, gmm::dense_matrix<T>& dest)
{
  if (&src != &dest) // allow in-place
	  dest.resize(gmm::mat_nrows(src), gmm::mat_ncols(src));
		
	gmm::copy(gmm::transposed(src), dest);	
}

// Works for vector or matrix "U" with implemented gmm::linalg_traits< U >::linalg_type definition.
// U <- U + delta
template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta )
{
  return offset( u, delta, typename gmm::linalg_traits< U >::linalg_type() );
}

template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta, gmm::abstract_vector )
{
  for(typename gmm::linalg_traits< U >::iterator 
        itU = gmm::vect_begin(u), itUEnd = gmm::vect_end(u);
      itU != itUEnd;
      ++itU)
    *itU += delta;      
}

template <class U>
void offset( U& u, const typename gmm::linalg_traits<U>::value_type& delta, gmm::abstract_matrix )
{
  for(typename gmm::linalg_traits< U >::row_iterator 
        itRow = gmm::mat_row_begin( u ), itRowEnd = gmm::mat_row_end( u );
      itRow != itRowEnd;
      ++itRow){
    typename gmm::linalg_traits< U >::sub_row_type 
      vRow( gmm::linalg_traits< U >::row(itRow) );
    for(typename gmm::linalg_traits< typename gmm::linalg_traits< U >::sub_row_type >::iterator 
          itV = gmm::vect_begin(vRow), itVEnd = gmm::vect_end(vRow);
        itV != itVEnd;
        ++itV)
      *itV += delta;      
  }  
}


// Works for vector or matrix "U" with implemented gmm::linalg_traits< U >::linalg_type definition.
template <class U>
inline long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_ )
{
  return numberNonzeroElements( u, eps_, typename gmm::linalg_traits< U >::linalg_type() );
}

template <class U>
long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_, gmm::abstract_vector )
{
  long nnzCount(0);

  for(typename gmm::linalg_traits< U >::const_iterator 
        itU = gmm::vect_const_begin(u), itUEnd = gmm::vect_const_end(u);
      itU != itUEnd;
      ++itU)
     if ( fabs(*itU) > fabs(eps_) )
       ++nnzCount;      
  
  return nnzCount;     
}

template <class U>
long numberNonzeroElements( const U& u, const typename gmm::number_traits<typename gmm::linalg_traits<U>::value_type>::magnitude_type& eps_, gmm::abstract_matrix )
{
  long nnzCount(0);
  for(typename gmm::linalg_traits< U >::const_row_iterator 
        itRow = gmm::mat_row_const_begin( u ), itRowEnd = gmm::mat_row_const_end( u );
      itRow != itRowEnd;
      ++itRow){
    typename gmm::linalg_traits< U >::const_sub_row_type 
      vRow( gmm::linalg_traits< U >::row(itRow) );
    for(typename gmm::linalg_traits< typename gmm::linalg_traits< U >::const_sub_row_type >::const_iterator 
          itV = gmm::vect_const_begin(vRow), itVEnd = gmm::vect_const_end(vRow);
        itV != itVEnd;
        ++itV)
       if ( fabs(*itV) > fabs(eps_) )
         ++nnzCount;      
  }
  
  return nnzCount;     
}

template <class U>
inline void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val )
{
  return fill( u, val, typename gmm::linalg_traits< U >::linalg_type() );
}

template <class U>
inline void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val, gmm::abstract_vector )
{
  std::fill( gmm::vect_begin(u), gmm::vect_end(u), val );
}

template <class U>
void fill( U& u, const typename gmm::linalg_traits<U>::value_type& val, gmm::abstract_matrix )
{
  for(typename gmm::linalg_traits< U >::row_iterator 
        itRow = gmm::mat_row_begin( u ), itRowEnd = gmm::mat_row_end( u );
      itRow != itRowEnd;
      ++itRow){
    typename gmm::linalg_traits< U >::sub_row_type 
      vRow( gmm::linalg_traits< U >::row(itRow) );
    std::fill( gmm::vect_begin(vRow), gmm::vect_end(vRow), val );
  }
}

template <class T, class U>
void iota( typename U::iterator itStart, typename U::iterator itEnd, const T& initValue, const T& dt )
{
 T t(initValue);
 for(typename U::iterator itU = itStart;
     itU != itEnd;
	   ++itU,
	   t += dt)
	 *itU = t;
}

template <class U>
inline void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* pVal )
{
  initFromArray( u, pVal, typename gmm::linalg_traits< U >::linalg_type() );
}

template <class U>
inline void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* pVal, gmm::abstract_vector )
{
  std::copy( pVal, pVal + gmm::vect_size(u), gmm::vect_begin(u) );
}

template <class U>
void initFromArray( U& u, const typename gmm::linalg_traits<U>::value_type* pVal, gmm::abstract_matrix )
{
  const typename gmm::linalg_traits<U>::value_type* itVal(pVal);
  const size_t ncols(gmm::mat_ncols(u));
  
  for(typename gmm::linalg_traits< U >::row_iterator 
        itRow = gmm::mat_row_begin( u ), itRowEnd = gmm::mat_row_end( u );
      itRow != itRowEnd;
      ++itRow,
      itVal += ncols){
    typename gmm::linalg_traits< U >::sub_row_type 
      vRow( gmm::linalg_traits< U >::row(itRow) );
    std::copy( itVal, itVal + ncols, gmm::vect_begin(vRow) );
  }
}

// Resize destination: will _not_ work correctly for vector references:
template <class V1, class V2>
void diff( const V1& v1, V2& v2, size_t nDiff )
{
 gmm::resize(v2, gmm::vect_size(v1)-nDiff );
 gmm::clear(v2); // initialize elements to zero.
 std::transform( (v1.begin()+nDiff), v1.end(), v1.begin(), v2.begin(),
                std::minus<typename gmm::linalg_traits<V1>::value_type>()); // value_type's must match
}

// unit vectors:
template <class V>
inline void e1( V& v )
{
 typedef typename gmm::linalg_traits<V>::value_type T;
 gmm::resize(v, 3);
 gmm::clear(v);
 v[0] = one<T>();
}

template <class V>
inline void e2( V& v )
{
 typedef typename gmm::linalg_traits<V>::value_type T;
 gmm::resize(v, 3);
 gmm::clear(v);
 v[1] = one<T>();
}

template <class V>
inline void e3( V& v )
{
 typedef typename gmm::linalg_traits<V>::value_type T;
 gmm::resize(v, 3);
 gmm::clear(v);
 v[2] = one<T>();
}


// convert for gmm::dense_matrix:
template <class T, class U>
void conv( gmm::dense_matrix<T>& aT, const gmm::dense_matrix<U>& aU)
{
 // resize the destination matrix:
 gmm::linalg_traits< gmm::dense_matrix<T> >::resize(aT, gmm::mat_nrows(aU), gmm::mat_ncols(aU));

 // iterate and convert elementwise:
 typename gmm::linalg_traits< gmm::dense_matrix<U> >::const_row_iterator itSrc = gmm::mat_row_const_begin(aU);
 for(typename gmm::linalg_traits< gmm::dense_matrix<T> >::row_iterator 
       itDest = gmm::mat_row_begin(aT), itDestEnd = gmm::mat_row_end(aT);
     itDest != itDestEnd;
     ++itSrc,
     ++itDest){
   typedef typename gmm::linalg_traits< gmm::dense_matrix<U> >::const_sub_row_type const_row_typeU;
   typedef typename gmm::linalg_traits< gmm::dense_matrix<T> >::sub_row_type row_typeT;
   const_row_typeU srcRow = gmm::linalg_traits< gmm::dense_matrix<U> >::row(itSrc);
   row_typeT destRow = gmm::linalg_traits< gmm::dense_matrix<T> >::row(itDest);

   typename gmm::linalg_traits< const_row_typeU >::const_iterator itSrcCol = gmm::vect_const_begin(srcRow);
   for(typename gmm::linalg_traits< row_typeT >::iterator 
         itDestCol = gmm::vect_begin(destRow),  itDestColEnd = gmm::vect_end(destRow);
       itDestCol != itDestColEnd;
       ++itSrcCol,
       ++itDestCol)
     conv(*itDestCol, *itSrcCol);           
 }   
}

} // namespace linalg


namespace commUtil{

// binary read and write for gmm::dense_matrix:

template<class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::dense_matrix<T>& M)
{
 return gmm::writeBinary(fp, M);
}


template<class T>
inline bool readBinary(abstractCommHandle *fp, gmm::dense_matrix<T>& M )
{
 return gmm::readBinary(fp, M);
}

template<class T>
inline size_t binarySize(const gmm::dense_matrix<T>& M )
{
 return gmm::binarySize(M);
}

// binary read and write for gmm::row_matrix< gmm::slvector<T> >:
template <class T>
inline bool writeBinary(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<T> >& M)
{ 
  return gmm::writeBinary(fp,M);
}

template <class T>
inline bool readBinary(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<T> >& M )
{
  return gmm::readBinary(fp,M);
}

template <class T>
inline size_t binarySize(const gmm::row_matrix< gmm::slvector<T> >& M)
{ 
  return gmm::binarySize(M);
}

// binary read and write for std::map:

template <class KEY, class DATA>
bool writeBinary(abstractCommHandle *fp, const std::map<KEY,DATA>& M)
{
 bool status(true);
 const size_t N(M.size());
 
 status = writeBinary(fp, N); 
 for (typename std::map<KEY,DATA>::const_iterator itM = M.begin(), itMEnd = M.end();
      status && (itM != itMEnd);
			++itM)
   status = (status && writeBinary( fp, (*itM).first ) && writeBinary( fp, (*itM).second ) );
 
 return status;	 
}

template <class KEY, class DATA>
bool readBinary(abstractCommHandle *fp,  std::map<KEY,DATA>& M)
{
 bool status(true);
 size_t N(0);
 
 status = readBinary(fp, N); 
 if (status){
   M.clear();
   M.reserve(N);
   for(size_t n=1; status && (n<=N); ++n){
     typename std::map<KEY,DATA>::value_type item;
	   status = (status && readBinary(fp, item));
		 if (status)		  
		   M.insert(item);
   }
 }
 
 return status;	 
}

template <class KEY, class DATA>
size_t binarySize(const std::map<KEY,DATA>& M)
{
 size_t val(0);
 
 // write size:
 val += sizeof(size_t);

 for (typename std::map<KEY,DATA>::const_iterator itM = M.begin(), itMEnd = M.end();
      itM != itMEnd;
			++itM)
   val += binarySize(*itM);
 
 return val;	 
}

// binary read and write for std::pair:

template <class T1, class T2>
bool writeBinary(abstractCommHandle *fp, const std::pair<T1,T2>& p)
{
 return
   writeBinary(fp, p.first)
	 && writeBinary(fp, p.second);
}

template <class T1, class T2>
bool readBinary(abstractCommHandle *fp, std::pair<T1,T2>& p)
{
 return
   readBinary(fp, p.first)
	 && readBinary(fp, p.second);
}

template <class T1, class T2>
inline size_t binarySize(const std::pair<T1,T2>& p)
{
 return binarySize(p.first) + binarySize(p.second);
}

} // namespace commUtil


namespace linalg{

template <class T>
void saveVectorData(const std::vector<T>& V, const std::string& sFileName)
{                    
    abstractCommHandle *fp = open(sFileName.c_str(),"wb");
    if (NULL==fp)
      throw std::runtime_error(std::string("Error: unable to open file ") + sFileName + " for write");

    if (!writeBinary( fp, V ))
      throw std::runtime_error(std::string("I/O error writing to file ") + sFileName);
		    
    close(fp);    
} 

template <class T>
void loadVectorData(std::vector<T>& V, const std::string& sFileName)
{                    
    abstractCommHandle *fp = open(sFileName.c_str(),"rb");
    if (NULL==fp)
      throw std::runtime_error(std::string("Error: unable to open file ") + sFileName + " for read");

    if (!readBinary( fp, V ))
      throw std::runtime_error(std::string("I/O error reading from file ") + sFileName);
		    
    close(fp);    
}
 
#if 1 // compiler work-around (gcc 4.1.1 not selecting overloaded "fabs", "real", "imag") moved from numericalConstants_template.h
template <class T1, class T2>
void conv( std::complex<T1>& t1, const T2& t2) // T2 -> std::complex<T1>
{
 T1 x(zero<T1>()), y(zero<T1>());
 conv(x, real(t2) ); 
 conv(y, imag(t2) );
 t1 = std::complex<T1>(x,y);
}


template <class T1, class T2>
void conv( T1& t1, const std::complex<T2>& t2) // std::complex<T2> -> T1
{
 conv(t1.real(), real(t2)); // this syntax _only_ makes sense for certain types T2...
 conv(t1.imag(), imag(t2));
}

template <class T1, class T2>
void conv( std::complex<T1>& t1, const std::complex<T2>& t2) // std::complex<T2> -> std::complex<T1>
{
 T1 x(zero<T1>()), y(zero<T1>());
 conv(x, real(t2) );
 conv(y, imag(t2) );
 t1 = std::complex<T1>(x,y);
}
#endif

#if 1 // compiler work-around (gcc 4.1.1 not recognizing appropriate overload of "fabs); moved from numericalConstants_template.h
template <class T>
inline const T max_abs(const T& a, const T& b)
{
 return T( fabs(a)>fabs(b)? fabs(a): fabs(b) );
}

template <class T>
inline const T min_abs(const T& a, const T& b)
{
 return T( fabs(a)<fabs(b)? fabs(a): fabs(b) );
}

template <class T>
inline const T max_abs(const T& a, const T& b, const T& c)
{
 const T d( max_abs(a,b) );
 return max_abs(d,c);
}

template <class T>
inline const T min_abs(const T& a, const T& b, const T& c)
{
 const T d( min_abs(a,b) );
 return min_abs(d,c);
}
#endif

#if 1 // moved from numericalConstants_template.h because of methods using std::vector:
/*
 * complex coordinate values:
 *   cartesian coordinates:  all of (x, y, z) may have a complex value
 *   spherical coordinates:  r may be complex, zenith and azimuth angles must be real
 */
 
//
// note: 
//   r = \tilde r \cdot \exp(i \phi_r), where \tilde r and \phi_r \in \Re
//   x = \tilde r \cdot \exp(i \phi_r) \sin(\theta) \cos(\phi) 
//   y = \tilde r \cdot \exp(i \phi_r) \sin(\theta) \sin(\phi)
// W.L.O.G: real zenith and azimuth angles:
//             \phi = atan2( \Re(y), \Re(x) )
//             \theta = acos( \Re(z/r) )
//
template <class T1, class T2>
void cart2sphere(const T1& x, const T1& y, const T1& z,
                   T2& r, T2& theta, T2& phi)
{
  typedef typename numberTraits<T1>::magnitudeType R1;
  typedef typename numberTraits<T2>::magnitudeType R2;

  // allow in-place:
  T2 r_(zero<R2>()); // allow complex r
  T2 tmp(zero<R2>());
  R2 theta_(zero<R2>()), 
     phi_(zero<R2>()),
     arg_(zero<R2>());

  // special case flags (zero coordinate cases):
  const unsigned long 
    zx( zero<T1>()==x? 0x1: 0x0 ), 
    zy( zero<T1>()==y? 0x2:0x0 ), 
    zz( zero<T1>()==z? 0x4: 0x0 ),
    test( zz | zy | zx );
    
  switch (test){
    case 0x0:
      // general case:

      // "branch" of sqrt defines "branch" of r:
      r_ = sqrt( x*x + y*y + z*z ); 

      tmp = z;
      tmp /= r_;
      theta_ = acos(real(tmp)); // see note above RE complex z,r

      phi_ = atan2(real(y), real(x)); // range: \phi \in [-\pi, \pi], see above note RE complex x,y
      if (zero<R2>() > phi_){
		    const R2 pi2(integer<R2>(2)*pi<R2>());
		    phi_ += pi2;
        if (phi_ >= pi2)    // this catches cases: <\Re(x) < 0>  or <\Re(y) < 0> &&  <|\Re(x)| < eps> or <|\Re(y)| < eps>
			    phi_ = zero<R2>();  
	    }

    break;

    case 0x1:
      // x == 0:

      r_ = sqrt(y*y + z*z);


      tmp = z;
      tmp /= r_;
      theta_ = acos(real(tmp)); // see note above RE complex z,r 

      arg_ = atan2(imag(y), real(y));
      if (fabs(arg_) < pi<R1>()/integer<R1>(2))
        phi_ = pi<R2>()/integer<R2>(2);
      else
        phi_ = ratio<R2,long>(3,2)*pi<R2>();

    break;
    case 0x2:
      // y == 0:

      r_ = sqrt(x*x + z*z);

      tmp = z;
      tmp /= r_;
      theta_ = acos(real(tmp)); // see note above RE complex z,r 

      arg_ = atan2(imag(x), real(x));
      if (fabs(arg_) < pi<R1>()/integer<R1>(2))
        phi_ = zero<R2>();
      else
        phi_ = pi<R1>();

    break;
    case 0x3:
      // x == 0  && y == 0:

      arg_ = atan2(imag(z), real(z));
      if (fabs(arg_) < pi<R1>()/integer<R1>(2)){
        r_ = z;
        theta_ = zero<R2>();
      }
      else{
        r_ = -z;
        theta_ = pi<R2>();
      }

      phi_ = zero<R2>();

    break;
    case 0x4:
      // z == 0:

      r_ = sqrt(x*x + y*y);

      theta_ = pi<R2>()/integer<R2>(2); 
      
      phi_ = atan2(real(y), real(x)); // range: \phi \in [-\pi, \pi], see above note RE complex x,y
      if (zero<R2>() > phi_){
		    const R2 pi2(integer<R2>(2)*pi<R2>());
		    phi_ += pi2;
        if (phi_ >= pi2)    // this catches cases: <\Re(x) < 0>  or <\Re(y) < 0> &&  <|\Re(x)| < eps> or <|\Re(y)| < eps>
			    phi_ = zero<R2>();  
	    }


    break;
    case 0x5:
      // x == 0 && z == 0:

      arg_ = atan2(imag(y), real(y));
      if (fabs(arg_) < pi<R1>()/integer<R1>(2)){
        r_ = y;
        theta_ = pi<R2>()/integer<R2>(2); 
        phi_ = pi<R2>()/integer<R2>(2);
      }
      else{
        r_ = -y;
        theta_ = pi<R2>()/integer<R2>(2);
        phi_ = ratio<R2,long>(3,2)*pi<R2>();
      }

    break;
    case 0x6:
      // y == 0 && z == 0:

      arg_ = atan2(imag(x), real(x));
      if (fabs(arg_) < pi<R1>()/integer<R1>(2)){
        r_ = x;
        theta_ = pi<R2>()/integer<R2>(2); 
        phi_ = zero<R2>();
      }
      else{
        r_ = -x;
        theta_ = pi<R2>()/integer<R2>(2);
        phi_ = pi<R2>();
      }


    break;
    case 0x7:
      // x == 0 && y == 0 && z == 0:

      r_ = zero<R2>();
      theta_ = zero<R2>();
      phi_ = zero<R2>();

    break;
    default:
      throw std::runtime_error("cart2sphere<T1,T2>: unrecognized special case");
  }

  // transfer to return values:
  r = r_;
  theta = theta_;
  phi = phi_; 
}

// 
// With respect to complex coordinates:
//   "sphere2cart" is simple, just use \Re(\theta), \Re(\phi):
//
template <class T1, class T2>
void sphere2cart(const T1& r_, const T1& theta_, const T1& phi_,
                   T2& x, T2& y, T2& z)
{
  typedef typename numberTraits<T1>::magnitudeType R1;

  // allow in-place
  T2 x_(zero<T2>()), y_(zero<T2>()), z_(zero<T2>()),
     rho(zero<T2>());

  const T1& r(r_);
  const R1 theta(real(theta_)), phi(real(phi_));
  assert(imag(theta_)==zero<R1>()  && imag(phi_)==zero<R1>());


  // special case flags (corresponding to cartesian-coordinate cases):
  const unsigned long 
    zx( (pi<T1>()/integer<T1>(2)==phi || ratio<T1,long>(3,2)*pi<T1>()==phi )? 0x1:0x0 ), 
    zy( (zero<T1>()==phi || pi<T1>()==phi )? 0x2: 0x0 ), 
    zz( theta==pi<R1>()/integer<R1>(2)? 0x4: 0x0 ),
    test( zz | zy | zx );
    
  switch (test){
    case 0x0:
      // general case:

      rho = r*sin(theta);
      x_ = rho*cos(phi);
      y_ = rho*sin(phi);
      z_ = r*cos(theta);
    
    break;

    case 0x1:
      // x == 0:

      if (pi<R1>()/integer<R1>(2) == phi){
        x_ = zero<T2>();
        y_ = r*sin(theta);
        z_ = r*cos(theta);
      }
      else{
        x_ = zero<T2>();
        y_ = -r*sin(theta);
        z_ = r*cos(theta);
      }

    break;
    case 0x2:
      // y == 0:

      if (zero<R1>() == phi){
        x_ = r*sin(theta);
        y_ = zero<T2>();
        z_ = r*cos(theta);
      }
      else{
        x_ = -r*sin(theta);
        y_ = zero<T2>();
        z_ = r*cos(theta);
      }

    break;
    case 0x3:
      // x == 0  && y == 0:

      if (zero<R1>() == theta){
        x_ = zero<T2>();
        y_ = zero<T2>();
        z_ = r;
      }
      else{
        x_ = zero<T2>();
        y_ = zero<T2>();
        z_ = -r;
      }      

    break;
    case 0x4:
      // z == 0:
      
      x_ = r*cos(phi);
      y_ = r*sin(phi);
      z_ = zero<T2>();

    break;
    case 0x5:
      // x == 0 && z == 0:

      if(pi<R1>()/integer<R1>(2) == phi){
        x_ = zero<T2>();
        y_ = r;
        z_ = zero<T2>();
      }
      else{
        x_ = zero<T2>();
        y_ = -r;
        z_ = zero<T2>();
      }

    break;
    case 0x6:
      // y == 0 && z == 0:

      if(zero<R1>() == phi){
        x_ = r;
        y_ = zero<T2>();
        z_ = zero<T2>();
      }
      else{
        x_ = -r;
        y_ = zero<T2>();
        z_ = zero<T2>();
      }

    break;
    case 0x7:
      // x == 0 && y == 0 && z == 0:

      x_ = zero<T2>();
      y_ = zero<T2>();
      z_ = zero<T2>();

    break;
    default:
      throw std::runtime_error("sphere2cart<T1,T2>: unrecognized special case");
  }

  // transfer to return values: 
  x = x_;
  y = y_;
  z = z_;
}
#endif

// for 3D vectors:
template <class T1, class T2>
inline void cart2sphere(const std::vector<T1>& vSrc, 
                 std::vector<T2>& vDest)
{
 assert(vSrc.size() == 3);
 vDest.resize(3);
 cart2sphere(vSrc[0], vSrc[1], vSrc[2], vDest[0], vDest[1], vDest[2]);
}								 
                 
template <class T1, class T2>
inline void sphere2cart(const std::vector<T1>& vSrc,
                 std::vector<T2>& vDest)  
{
 assert(vSrc.size() == 3);
 vDest.resize(3);
 sphere2cart(vSrc[0], vSrc[1], vSrc[2], vDest[0], vDest[1], vDest[2]);
}	

template <class T1, class T2>
inline void cart2sphere(const ntuple<T1,3>& vSrc, 
                 ntuple<T2,3>& vDest)
{
 cart2sphere(vSrc[0], vSrc[1], vSrc[2], vDest[0], vDest[1], vDest[2]);
}								 
                 
template <class T1, class T2>
inline void sphere2cart(const ntuple<T1,3>& vSrc,
                 ntuple<T2,3>& vDest)  
{
 sphere2cart(vSrc[0], vSrc[1], vSrc[2], vDest[0], vDest[1], vDest[2]);
}	
				         
// Conversion of vector components:
template <class T>
void cartVect2sphere(const T& x, const T& y, const T& z,
                     const T& vx, const T& vy, const T& vz,
                     T& r, T& theta, T& phi,
										 T& vr, T& vtheta, T& vphi)
{
 std::vector<T> 
   vSrcP(3,zero<T>()), vSrcV(3,zero<T>()),
   vDestP(3,zero<T>()), vDestV(3,zero<T>());
	  
 vSrcP[0] = x; vSrcP[1] = y; vSrcP[2] = z;
 vSrcV[0] = vx; vSrcV[1] = vy; vSrcV[2] = vz;
 
 cartVect2sphere(vSrcP, vSrcV, vDestP, vDestV);  

 r = vDestP[0]; theta = vDestP[1]; phi = vDestP[2];
 vr = vDestV[0]; vtheta = vDestV[1]; vphi = vDestV[2];  
}										 
                 
template <class T>
void sphereVect2cart(const T& r, const T& theta, const T& phi,
                     const T& vr, const T& vtheta, const T& vphi,
                     T& x, T& y, T& z,
										 T& vx, T& vy, T& vz)                 
{
 std::vector<T> 
   vSrcP(3,zero<T>()), vSrcV(3,zero<T>()),
   vDestP(3,zero<T>()), vDestV(3,zero<T>());
	  
 vSrcP[0] = r; vSrcP[1] = theta; vSrcP[2] = phi;
 vSrcV[0] = vr; vSrcV[1] = vtheta; vSrcV[2] = vphi;
 
 sphereVect2cart(vSrcP, vSrcV, vDestP, vDestV);
 
 x = vDestP[0]; y = vDestP[1]; z = vDestP[2];
 vx = vDestV[0]; vy = vDestV[1]; vz = vDestV[2];  
}	

// for 3D vectors:
template <class T>
inline void cartVect2sphere(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV)
{
  typedef typename numberTraits<T>::magnitudeType R;
  
  cart2sphere(vSrcP,vDestP);

  // see note at "cart2sphere" w.r.t. complex coordinates; at this point the angular coordinates will be real-valued:  
  const R st(sin(real(vDestP[1]))), ct(cos(real(vDestP[1]))),
	  sp(sin(real(vDestP[2]))), cp(cos(real(vDestP[2])));
	
	gmm::dense_matrix<T> aCartToSphere(3,3);  // could be of type "R", but leave as "T" to simplify "mult" below
	/*
  aCartToSphere = [ st*cp, st*sp,  ct ;...
                    ct*cp, ct*sp,  -st;...
                    -sp,   cp,     0           ];
	*/
	aCartToSphere(0,0) = st*cp; aCartToSphere(0,1) = st*sp; aCartToSphere(0,2) = ct;
	aCartToSphere(1,0) = ct*cp; aCartToSphere(1,1) = ct*sp; aCartToSphere(1,2) = -st;
	aCartToSphere(2,0) = -sp;   aCartToSphere(2,1) = cp;    aCartToSphere(2,2) = zero<T>();
	
	if (vDestV.size() != 3) // allows in-place (via gmm, in-place...)
	  vDestV.resize(3,zero<T>());
	gmm::mult(aCartToSphere, vSrcV, vDestV);									
}	
                 
template <class T>
inline void sphereVect2cart(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV)    
{
  typedef typename numberTraits<T>::magnitudeType R;
  
  sphere2cart(vSrcP,vDestP);

  // see note at "sphere2cart" w.r.t. complex coordinates; here the angular coordinates should be real-valued:  

  const R st(sin(real(vSrcP[1]))), ct(cos(real(vSrcP[1]))),
	  sp(sin(real(vSrcP[2]))), cp(cos(real(vSrcP[2])));
	
	gmm::dense_matrix<T> aSphereToCart(3,3); // could be of type "R", but leave as "T" to simplify "mult" below
	/*
    aSphereToCart = [ sin(t)*cos(p), cos(t)*cos(p),       -sin(p);...
                     sin(t)*sin(p), cos(t)*sin(p),        cos(p);...
                     cos(t),       -sin(t),             0];	
  */
	aSphereToCart(0,0) = st*cp; aSphereToCart(0,1) = ct*cp; aSphereToCart(0,2) = -sp;
	aSphereToCart(1,0) = st*sp; aSphereToCart(1,1) = ct*sp; aSphereToCart(1,2) = cp;
	aSphereToCart(2,0) = ct;    aSphereToCart(2,1) = -st;   aSphereToCart(2,2) = zero<T>();
	
	if (vDestV.size() != 3) // allows in-place (via gmm, in-place...)
	  vDestV.resize(3,zero<T>());

	gmm::mult(aSphereToCart, vSrcV, vDestV);	
}	



#if 1 // ========================================== moved from numericalFunctor_template.h: =====================================================
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
void quadraticFormula(const V1& vP, V2& vRoot, typename V1::value_type& D) 
{
  typedef typename V1::value_type T1;
  typedef typename V2::value_type T2;
  
  if ( vP.size() != 3 )
    throw std::runtime_error("quadraticFormula<V1,V2>: vP is not a 2nd degree polynomial");

  // Implementation note: could also assume that output parameter has already been allocated (in case of _reference_ output type);
  //   however, at the moment, I do not want to bring (or duplicate) the "gmm" linalg_traits mechanisms here, 
  //   and I'll assume that "V2" behaves as an STL container-type.
  vRoot.clear();   
  vRoot.resize(2,zero<T2>()); // two roots

  const T1 
    &a(vP[0]),
    &b(vP[1]),
    &c(vP[2]);

  // temporaries:
  T1 t1(zero<T1>());
  T2 t2(zero<T2>());

  D = integer<T1>(4);
  D *= -a*c;
  D += pow_n(b,2);   
  t2 = sqrt(T2(D)); // take the sqrt of "D" as type "T2"

  t1 = integer<T1>(2);
  t1 *= a;
  
  vRoot[0] = ( -b + t2 )/t1;
  vRoot[1] = ( -b - t2 )/t1;
}


/**
 * @brief Calculate the discriminant for the quadratic formula:
 * @tparam V1  container type of value_type real or complex
 * @param[in]  vP  vector of polynomial coefficients
 * @param[out]  D  discriminant
 */
template <class V1>
void quadraticDiscriminant(const V1& vP, typename V1::value_type& D)
{
  typedef typename V1::value_type T1;
  
  if ( vP.size() != 3 )
    throw std::runtime_error("quadraticDiscriminant<V1>: vP is not a 2nd degree polynomial");

  const T1 
    &a(vP[0]),
    &b(vP[1]),
    &c(vP[2]);

  T1 temp(zero<T1>());

  D = integer<T1>(4);
  D *= -a*c;
  D += pow_n(b,2);   
}
 

template <class T>
void cubicFormula( const std::vector<T>& vP, std::vector< std::complex<T> >& vRoot )
{
  if ( vP.size() != 4 )
    throw std::runtime_error("cubicFormula<T>: vP is not a 3rd degree polynomial");

  vRoot.clear();   
  vRoot.resize(3, std::complex<T>( zero<T>(), zero<T>() ) ); // three roots

  // x^3 + a2*x^2 + a1*x + a0 = 0
  T a2 = vP[1]/vP[0],
    a1 = vP[2]/vP[0],
    a0 = vP[3]/vP[0];

/* 
 ***********
  % From Mathworld: doesn't work for some reason:
  if 0
  Q = (3*a1 - a2^2)/9;
  R = (9*a2*a1 - 27*a0 - 2*a2^3)/54;
  D = (Q^3 + R^2);
  S = ( R + sqrt(D) )^(1/3);
  T = ( R - sqrt(D) )^(1/3);
  vR = [ ( (S+T) - a2/3 ),...
          ( -(1/2)*(S+T) + (1/2)*i*sqrt(3)*(S-T) - a2/3 ),...
          ( -(1/2)*(S+T) - (1/2)*i*sqrt(3)*(S-T) - a2/3 ) ...
    ];
  end
 ***********
 */

// From maple 'solve':
T a2_3 = pow_n(a2,3), a2_2 = pow_n(a2,2),
  a1_3 = pow_n(a1,3), a1_2 = pow_n(a1,2),
  a0_2 = pow_n(a0,2);

// ****** Perform calculations using explicit promotion ONLY (otherwise only double precision is achieved): *******     
 
T tIntFactor( zero<T>() );
T tInt( zero<T>() );
// T one( one<T>() ); // hold "exact" integer factors as a T.
T tRealPart( zero<T>() );


conv(tIntFactor,12);
tIntFactor *= a0*a2_3;
tRealPart += tIntFactor;

conv(tIntFactor,81);
tIntFactor *= a0_2;
tRealPart += tIntFactor;

conv(tIntFactor,-54);
tIntFactor *= a1*a2*a0;
tRealPart += tIntFactor;

conv(tIntFactor,-3);
tIntFactor *= a1_2*a2_2;
tRealPart += tIntFactor;

conv(tIntFactor,12);
tIntFactor *= a1_3;
tRealPart += tIntFactor;

std::complex<T> E10(tRealPart, zero<T>());  // (12*a1_3-3*a1_2*a2_2-54*a1*a2*a0+81*a0_2+12*a0*a2_3)

tRealPart = 0;

conv(tIntFactor,-8);
tIntFactor *= a2_3;
tRealPart += tIntFactor;

conv(tIntFactor,-108);
tIntFactor *= a0;
tRealPart += tIntFactor;

conv(tIntFactor,36);
tIntFactor *= a1*a2;
tRealPart += tIntFactor;

std::complex<T> E1(tRealPart, zero<T>()); // 36*a1*a2-108*a0-8*a2_3

conv(tInt,12);
E1 += tInt*sqrt(E10);
std::complex<T> z_tIntFactor( ratio<T,long>(1,3), zero<T>() );
E1 = pow(E1, z_tIntFactor ); // E1 = ( E1 + 12*E10^(1/2) )^(1/3)

std::complex<T> E2( ratio<T,long>(1,3)*a1 - ratio<T,long>(1,9)*a2_2 , zero<T>()); // (1/3*a1-1/9*a2_2);

tInt = 0;
std::complex<T> ISqrt3_4( tInt, sqrt( ratio<T,long>(3,4) ) ); // sqrt( -3/4 )

conv(tInt,6);
vRoot[0] = ratio<T,long>(1,6)*E1 - tInt*E2/E1 - ratio<T,long>(1,3)*a2; // ( 1/6*E1 -6.0*E2/E1 -1/3*a2 );

conv(tInt,3);
std::complex<T> E3( -ratio<T,long>(1,12)*E1 + tInt*E2/E1 -ratio<T,long>(1,3)*a2 );

conv(tInt,6);
std::complex<T> E4( ISqrt3_4*( ratio<T,long>(1,6)*E1 + tInt*E2/E1 ) );

vRoot[1] = E3 + E4; // ( -1/12*E1 +3*E2/E1 +ISqrt3_4*(  1/6*E1 +6*E2/E1) -1/3*a2 );
vRoot[2] = E3 - E4; // ( -1/12*E1  +3*E2/E1 -ISqrt3_4*(  1/6*E1 +6*E2/E1) -1/3*a2 );   

}


template <class T>
void cubicFormula( const std::vector< std::complex<T> >& vP, std::vector< std::complex<T> >& vRoot )
{
  if ( vP.size() != 4 )
    throw std::runtime_error("cubicFormula<T>: vP is not a 3rd degree polynomial");

  vRoot.clear();   
  vRoot.resize(3, std::complex<T>( zero<T>(), zero<T>() ) ); // three roots

  // x^3 + a2*x^2 + a1*x + a0 = 0
  std::complex<T>
    a2 = vP[1]/vP[0],
    a1 = vP[2]/vP[0],
    a0 = vP[3]/vP[0];

/* 
 ***********
  % From Mathworld: doesn't work for some reason:
  if 0
  Q = (3*a1 - a2^2)/9;
  R = (9*a2*a1 - 27*a0 - 2*a2^3)/54;
  D = (Q^3 + R^2);
  S = ( R + sqrt(D) )^(1/3);
  T = ( R - sqrt(D) )^(1/3);
  vR = [ ( (S+T) - a2/3 ),...
          ( -(1/2)*(S+T) + (1/2)*i*sqrt(3)*(S-T) - a2/3 ),...
          ( -(1/2)*(S+T) - (1/2)*i*sqrt(3)*(S-T) - a2/3 ) ...
    ];
  end
 ***********
 */

// From maple 'solve':
std::complex<T>
  a2_3 = pow_n(a2,3), a2_2 = pow_n(a2,2),
  a1_3 = pow_n(a1,3), a1_2 = pow_n(a1,2),
  a0_2 = pow_n(a0,2);

// ****** Perform calculations using explicit promotion ONLY (otherwise only double precision is achieved): *******     
  
T tInt( zero<T>() );
// T one( one<T>() ); // hold "exact" integer factors as a T.

std::complex<T> ztTermPart( zero<T>(), zero<T>() ), ztFactor( zero<T>(), zero<T>() );

conv(tInt,12);
ztFactor = tInt;
ztFactor *= a0*a2_3;
ztTermPart += ztFactor;

conv(tInt,81);
ztFactor = tInt;
ztFactor *= a0_2;
ztTermPart += ztFactor;

conv(tInt,-54);
ztFactor = tInt;
ztFactor *= a1*a2*a0;
ztTermPart += ztFactor;

conv(tInt,-3);
ztFactor = tInt;
ztFactor *= a1_2*a2_2;
ztTermPart += ztFactor;

conv(tInt,12);
ztFactor = tInt;
ztFactor *= a1_3;
ztTermPart += ztFactor;

std::complex<T> E10(ztTermPart);  // (12*a1_3-3*a1_2*a2_2-54*a1*a2*a0+81*a0_2+12*a0*a2_3)

ztTermPart = zero<T>();

conv(tInt,-8);
ztFactor = tInt;
ztFactor *= a2_3;
ztTermPart += ztFactor;

conv(tInt,-108);
ztFactor = tInt;
ztFactor *= a0;
ztTermPart += ztFactor;

conv(tInt,36);
ztFactor = tInt;
ztFactor *= a1*a2;
ztTermPart += ztFactor;

std::complex<T> E1(ztTermPart); // 36*a1*a2-108*a0-8*a2_3

conv(tInt,12);
E1 += tInt*sqrt(E10);
std::complex<T> z_tIntFactor( ratio<T,long>(1,3), zero<T>() );
E1 = pow(E1, z_tIntFactor ); // E1 = ( E1 + 12*E10^(1/2) )^(1/3)

std::complex<T> E2( ratio<T,long>(1,3)*a1 - ratio<T,long>(1,9)*a2_2 ); // (1/3*a1-1/9*a2_2);

tInt = 0;
std::complex<T> ISqrt3_4( tInt, sqrt( ratio<T,long>(3,4) ) ); // sqrt( -3/4 )

conv(tInt,6);
vRoot[0] = ratio<T,long>(1,6)*E1 - tInt*E2/E1 - ratio<T,long>(1,3)*a2; // ( 1/6*E1 -6.0*E2/E1 -1/3*a2 );

conv(tInt,3);
std::complex<T> E3( -ratio<T,long>(1,12)*E1 + tInt*E2/E1 -ratio<T,long>(1,3)*a2 );

conv(tInt,6);
std::complex<T> E4( ISqrt3_4*( ratio<T,long>(1,6)*E1 + tInt*E2/E1 ) );

vRoot[1] = E3 + E4; // ( -1/12*E1 +3*E2/E1 +ISqrt3_4*(  1/6*E1 +6*E2/E1) -1/3*a2 );
vRoot[2] = E3 - E4; // ( -1/12*E1  +3*E2/E1 -ISqrt3_4*(  1/6*E1 +6*E2/E1) -1/3*a2 );   

}


template <class T>
T polyval( const std::vector<T>& vP, const T& arg ) 
{
  T tPow = one<T>(); 
  T val = zero<T>(); 

  // Evaluate the polynomial (start with the constant term vP(end) and reverse iterate through vP):  
  for (typename std::vector<T>::const_reverse_iterator itP = vP.rbegin(), itPEnd = vP.rend();
       itP != itPEnd;
       ++itP,
       tPow *= arg){
     val += (*itP)*tPow;
  }
  
  return val;
}  
  
template <class T>
std::complex<T> polyval( const std::vector< std::complex<T> >& vP, const std::complex<T>& arg ) 
{
  T tInt = one<T>();
  std::complex<T> ztPow( zero<T>(), zero<T>() ); 

  ztPow = tInt; // arg^0
  
  std::complex<T> val( zero<T>(), zero<T>() );

  // Evaluate the polynomial (start with the constant term vP(end) and reverse iterate through vP):  
  for (typename std::vector< std::complex<T> >::const_reverse_iterator itP = vP.rbegin(), itPEnd = vP.rend();
       itP != itPEnd;
       ++itP,
       ztPow *= arg){
     val += (*itP)*ztPow;
  }
  
  return val;
}  

// convert vector [v_r, v_theta, v_phi](r, theta, phi) to [v_minus, v_0, v_plus] in helical basis.  
//   implementation note: "inline" so that this function can be efficient in iterator loops, and so that
//   requirement for vector and matrix versions doesn't need to be designed.  I.E. Assume that this function
//   is only used inside of other functors (such as the general tmatrix<T>).  If this changes, this routine
//   should be re-implemented as a functor.
//
#if defined(__ICC)
#pragma warning(disable:869) // remark #869: parameter "r" was never referenced
#endif
template <class T>
inline void sphericalToHelical( const T& v_r, const T& v_theta, const T& v_phi,
                                const T& r, const T& theta, const T& phi,
                                T& v_minus, T& v_0, T& v_plus, bool cosTheta )
{
 const T 
   ct( cosTheta? theta: cos(theta) ), st( cosTheta? sqrt(one<T>()-sqr(ct)): sin(theta) ),
   exp_p( exp(one_i<T>()*phi) ), exp_m( exp(-one_i<T>()*phi ) ),
   sqrt2_2( sqrt(integer<T>(2))/integer<T>(2) );

 v_minus = sqrt2_2 * exp_p * (st*v_r  + ct*v_theta  + one_i<T>()*v_phi); 
 v_0 =     ct*v_r  - st*v_theta;
 v_plus =  -sqrt2_2 * exp_m * (st*v_r  + ct*v_theta  - one_i<T>()*v_phi);

}

// Same thing, but parsed by m-value:                                   
template <class T>
inline void sphericalToHelical( const T& v_r, const T& v_theta, const T& v_phi,
                                const T& r, const T& theta, const long& m,
                                T& v_minus, T& v_0, T& v_plus, bool cosTheta )
{
 const T 
   ct( cosTheta? theta: cos(theta) ), st( cosTheta? sqrt(one<T>() - sqr(ct)): sin(theta) ),
   sqrt2_2( sqrt(integer<T>(2)) / integer<T>(2) );

 switch (m){
   case -1:
     v_minus = zero<T>();
     v_0 = zero<T>();
     v_plus =  -sqrt2_2 * (st * v_r  + ct * v_theta  - one_i<T>() * v_phi);
   break;
   
   case 0:
     v_minus = zero<T>();
     v_0 =     ct * v_r  - st * v_theta;
     v_plus = zero<T>();
   break;
   
   case 1:
     v_minus = sqrt2_2 * (st * v_r  + ct * v_theta  + one_i<T>() * v_phi); 
     v_0 = zero<T>();
     v_plus = zero<T>();
   break;
   
   default:
     throw std::runtime_error("sphericalToHelical: m-value out of domain {-1, 0, 1} ");
   // break;
 }
} 
#if defined(__ICC)
#pragma warning(default:869)                              
#endif
  // *******************************************************************************************************  
#endif // ========================================== end, moved from numericalFunctor_template.h: =====================================================


	
	
// return a _uniform_ distribution for the specified co-ordinate variable (with the specified symmetry):
template <class R>
R randomCoord(CSYS_KIND eCSYS, size_t dimOffset)
{
 R val(zero<R>());

 R r(random<R>()); // r \in [0,1], uniform distribution; 

 switch (eCSYS){
 
	 case NO_CSYS:
		 val = r; // rand \in [0,1]
	 break;

	 case CARTESIAN: // (x, y, z, ... )
		 val = r - ratio<R>(1,2); // rand \in [-0.5, 0.5]
	 break;

	 case CYLINDRICAL: // (\rho, \phi, z, ...)
		 switch (dimOffset){
			 case 0:
				 val = sqrt(r);  // rand \in [0, 1]					 
			 break;

			 case 1:
				 val = r*( integer<R>(2)*pi<R>() ); // rand \in [0, 2*pi]
			 break;

			 case 2:
				 val = r - ratio<R>(1,2);  // rand \in [-0.5, 0.5]
			 break; 

			 default:
				 throw std::runtime_error("linalgUtil::randomCoord: cylindrical symmetry with co-ordinate offset out of range [0,2]");
			 break; 
		 }			 
	 break;

	 case SPHERICAL: // (r,\theta, \phi)

		 switch (dimOffset){
			 case 0:
				 val = cbrt(r);  // rand \in [0, 1]					 
			 break;

			 case 1:
				 val = asin(r) * pi<R>(); // rand \in [0,pi]
			 break;

			 case 2:
				 val = r*( integer<R>(2)*pi<R>() ); // rand \in [0, 2*pi]
			 break;

			 default:
				 throw std::runtime_error(std::string("linalgUtil::randomCoord: spherical symmetry with co-ordinate offset out of range [0,2]\n")
						               + " (hyperspheres _not_ implemented because DIM !=> <max number of co-ords>)");
			 break; 
		 }			 
	 break;

	 default:
		 throw std::runtime_error("linalgUtil::randomCoord: unknown co-ordinate symmetry");
	 break;
 }
 
 return val;
}

#if 0
// OK: original reference, and its idea for 3D Box-Mueller transform appears to *not* be correct.
//   This reference was: Aravind Kalaiah and Amitabh Varshney, "Statistical Point Geometry", Eurographics Symposium on Geometry Processing (2003), 
//     L. Kobbelt and P. Schroeder and H. Hoppe (Editors). 

// return an instance of an isotropic 3D _normal_ distribution with the given mean and variance (cartesian form):
template <class R>
void normalDist(const ntuple<R,3>& mean, const R& variance, ntuple<R,3>& p)
{
 // This is _exact_, but not optimized for speed (or stability).  OK for arbitrary precision, for now...
 // 3D Box-Mueller transform:
 // ( see G.E.P. Box and M.E. M\"{u}ller, "A note on the generation of random normal deviates." Ann. Math. Stat., 28:610-611, 1958
 
 // first calculate the distribution with zero mean and unit variance:
 R r0(random<R>()), r1(random<R>()), r2(integer<R>(2)*random<R>() - one<R>()), tau(zero<R>()), rTmp1(zero<R>()), rTmp2(zero<R>());
 if (r0 < epsilon<R>()) 
   r0 = epsilon<R>(); 
 tau = sqrt(-integer<R>(2)*log(r0));	// radius 
 rTmp1 = one<R>();
 rTmp1 -= sqr(r2);
 rTmp2 = sqrt(rTmp1); // sin(<zenith angle>)
 rTmp2 *= tau;        // transverse radius
 rTmp1 = cos(integer<R>(2)*pi<R>()*r1);
 p[0] = rTmp2;
 p[1] = rTmp2;
 p[0] *= rTmp1; // x
 rTmp2 = one<R>();
 rTmp2 -= sqr(rTmp1);
 rTmp1 = sqrt(rTmp2); // sin(<azimuth angle>)
 p[1] *= rTmp1; // y;
 p[2] = tau;
 p[2] *= r2;    // z;
 
 // correct for specified mean and variance:
 p *= sqrt(integer<R>(2)*variance);
 p += mean;
}
#endif
								 
} // namespace linalg


// writeBinary, readBinary for interfaced "gmm" classes.
// This is located in namespace "gmm" so it won't conflict with the namespace TMatrix methods.
namespace gmm{

template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u)
{
 return gmm::writeBinary(fp, u, typename linalg_traits<U>::linalg_type());
}

template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector)
{
 using commUtil::writeBinary;
 bool status(true);
 status = ( status && writeBinary(fp, vect_size(u)) );
 
 return 
   status
   && writeBinary(fp, u, abstract_vector(), typename linalg_traits<U>::storage_type());
}

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector, abstract_dense)
{
 using commUtil::writeBinary;
 bool status(true);

 // iterate and write elementwise:
 for(typename linalg_traits<U>::const_iterator 
       itU = vect_const_begin(u), itUEnd = vect_const_end(u);
     status && (itU != itUEnd);
     ++itU)
	 status = (status && writeBinary(fp, *itU) );	 

 return status;    
}

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_vector, abstract_skyline)
{
 using commUtil::writeBinary;
 bool status(true);
 const size_t 
   size_(linalg_traits<U>::size(u)),
   data_size_( /* linalg_traits<U>:: */ nnz(u)),
   shift_(linalg_traits<U>::origin(u)->first()); 
   
 status = ( status && writeBinary(fp, size_) ); // redundant write (to outer above) but simplest to _leave_ in (for purposes of read-back verify).
 status = ( status && writeBinary(fp, data_size_) );
 status = ( status && writeBinary(fp, shift_) ); 
  
 // iterate and write elementwise:
 for(typename linalg_traits<U>::const_iterator 
       itU = vect_const_begin(u), itUEnd = vect_const_end(u);
     status && (itU != itUEnd);
     ++itU)
	 status = (status && writeBinary(fp, *itU) );	 

 return status;    
}

template <class U>
inline bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, abstract_matrix)
{
 using commUtil::writeBinary;
 bool status(true);
 
 // write the size information:
 status = (status && writeBinary(fp, mat_nrows(u) ));
 status = (status && writeBinary(fp, mat_ncols(u) ));

 return 
   status
   && gmm::writeBinary(fp, u, 
                  typename principal_orientation_type
                    <typename linalg_traits<U>::sub_orientation>::potype());
}

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, row_major)
{
 bool status(true);
 
 // iterate and write by row vector:
 for(typename linalg_traits<U>::const_row_iterator 
       itRow = mat_row_const_begin(u), itRowEnd = mat_row_const_end(u);
     status && (itRow != itRowEnd);
     ++itRow)
   status = (status && writeBinary(fp, linalg_traits<U>::row(itRow)));   
 
 return status;   
}

template <class U>
bool writeBinary(commUtil::abstractCommHandle *fp, const U& u, col_major)
{
 bool status(true);
 
 // iterate and write by column vector:
 for(typename linalg_traits<U>::const_col_iterator 
       itCol = mat_col_const_begin(u), itColEnd = mat_col_const_end(u);
     status && (itCol != itColEnd);
     ++itCol)
   status = (status && writeBinary(fp, linalg_traits<U>::col(itCol)));     

 return status;   
}

// _allow_ references => do _not_ necessarily resize destinations:
template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u)
{
 return readBinary(fp, u, typename linalg_traits<U>::linalg_type(), typename linalg_traits<U>::is_reference());
}

// allow "const reference" to temporary (this is the gmm reference type, which may actually be modifiable!).
template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, const U& u)
{
 return readBinary(fp, linalg_const_cast(u));
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_false)
{
 using commUtil::readBinary;
 bool status(true);
 
 size_t N(0);
 status = ( status && readBinary(fp, N) );
 
 resize( u, N );

 return
   status
   && readBinary(fp, u, abstract_vector(), typename linalg_traits<U>::storage_type());
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_modifiable)
{
 using commUtil::readBinary;
 bool status(true);
 
 size_t N(0);
 status = ( status && readBinary(fp, N) );
 
 if (N != vect_size(u))
   throw std::runtime_error("gmm::readBinary(abstractCommHandle *, U&, abstract_vector, linalg_modifiable): size of destination doesn't match read size");

 return
   status
   && gmm::readBinary(fp, u, abstract_vector(), typename linalg_traits<U>::storage_type());	 
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, linalg_const)
{
 throw std::runtime_error("gmm::readBinary: can't modify a const reference");
 return false;
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, abstract_dense)
{
 using commUtil::readBinary;
 bool status(true);

 // iterate and read elementwise:
 for(typename linalg_traits<U>::iterator 
       itU = vect_begin(u), itUEnd = vect_end(u);
     status && (itU != itUEnd);
     ++itU)
	 status = (status && readBinary(fp, *itU) );	 

 return status;     
}

template <class U>
bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_vector, abstract_skyline)
{
 using commUtil::readBinary;
 bool status(true);
 
 typedef typename linalg_traits<U>::value_type T;
 
 size_t 
   size_(0),
   data_size_(0),
   shift_(0); 
   
 status = ( status && readBinary(fp, size_) );
 status = ( status && readBinary(fp, data_size_) );
 status = ( status && readBinary(fp, shift_) ); 
 
 if (status){ 
   // This is the only efficient (hopefully) and _legal_ way to set the "shift" (apart from re-writing gmm code).
   // (here it is assumed that "swap" swaps the data pointers).
   slvector<T> u_(size_, data_size_, shift_);

   // iterate and read elementwise:
   for(typename linalg_traits<slvector<T> >::iterator 
         itU = vect_begin(u_), itUEnd = vect_end(u_);
       status && (itU != itUEnd);
       ++itU)
	   status = (status && readBinary(fp, *itU) );
     
   if (status)
     std::swap(*linalg_traits<U>::origin(u), u_);   	 
 }

 return status;    
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_false)
{
 using commUtil::readBinary;
 bool status(true);
 
 // read size information:
 size_t nrows(0), ncols(0);
 status = (status && readBinary(fp, nrows ));
 status = (status && readBinary(fp, ncols ));
 
 resize(u, nrows, ncols); 
 	
 return 
   status
	 && readBinary(fp, u, 
        typename principal_orientation_type <typename linalg_traits<U>::sub_orientation>::potype() );	    
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_modifiable)
{
 using commUtil::readBinary;
 bool status(true);
 
 // read size information:
 size_t nrows(0), ncols(0);
 status = (status && readBinary(fp, nrows ));
 status = (status && readBinary(fp, ncols ));
  
 // don't resize _reference_ types:
 if ( (mat_nrows(u) != nrows) || (mat_ncols(u) != ncols) )
   throw std::runtime_error("gmm::readBinary(abstractCommHandle*, U&, abstract_matrix, linalg_modifiable): size of destination doesn't match read size");
	 	
 return 
   status
	 && readBinary(fp, u, 
        typename principal_orientation_type <typename linalg_traits<U>::sub_orientation>::potype());	 
}

template <class U>
inline bool readBinary(commUtil::abstractCommHandle *fp, U& u, abstract_matrix, linalg_const)
{
 throw std::runtime_error("gmm::readBinary: can't modify a const reference");
 return false;
}

template <class U>
bool readBinary(commUtil::abstractCommHandle *fp, U& u, row_major)
{
 bool status(true);

 // iterate and read by row vector:
 for(typename linalg_traits< U >::row_iterator 
       itRow = mat_row_begin(u), itRowEnd = mat_row_end(u);
     status && (itRow != itRowEnd);
     ++itRow)
   status = (status && readBinary(fp, linalg_traits<U>::row(itRow)));  
 
 return status;   
}

template <class U>
 bool readBinary(commUtil::abstractCommHandle *fp, U& u, col_major)
{
 bool status(true);

 // iterate and read by col vector:
 for(typename linalg_traits< U >::col_iterator 
       itCol = mat_col_begin(u), itColEnd = mat_col_end(u);
     status && (itCol != itColEnd);
     ++itCol)
   status = (status && readBinary(fp, linalg_traits<U>::col(itCol)));  
 
 return status;   
}


// in consideration of variable-precision number types, the following "binarySize" method requires an instance of class U:
template <class U>
inline size_t binarySize(const U& u)
{
 return gmm::binarySize(u, typename linalg_traits<U>::linalg_type());
}

template <class U>
inline size_t binarySize(const U& u, abstract_vector)
{
  return binarySize(u, abstract_vector(), typename linalg_traits<U>::storage_type());
}

template <class U>
inline size_t binarySize(const U& u, abstract_vector, abstract_dense)
{    
  using commUtil::binarySize;
  
  // assume: all elements of vector have same size: 
  size_t N_elt(vect_size(u)), N(0);

  if (N_elt > 0)
    N = N_elt*binarySize(u[0]); 

  return N + sizeof(size_t); // account for size_t vect_size in written info.
}

template <class U>
inline size_t binarySize(const U& u, abstract_vector, abstract_skyline)
{    
  using commUtil::binarySize;
  
  // assume: all elements of vector have same size: 
  size_t N_elt(vect_const_end(u)-vect_const_begin(u)), N(0);

  if (N_elt > 0)
    N = N_elt*binarySize(*vect_const_begin(u)); 

  return N + 4*sizeof(size_t); // account for size_t: vect_size, vect_size(redundant), data_size, shift in written info.
}   

template <class U>
inline size_t binarySize(const U& u, abstract_matrix)
{
  return binarySize(u, abstract_matrix(),
                    typename principal_orientation_type
                      <typename linalg_traits<U>::sub_orientation>::potype(),
                       typename linalg_traits<U>::storage_type()); 
}


template <class U>
inline size_t binarySize(const U& u, abstract_matrix, row_major, abstract_dense)
{    
  using commUtil::binarySize;
  
  // assume: all elements of matrix have same size: 
  size_t N_row(mat_nrows(u)), N(0);

  if (N_row > 0)
    N = N_row*binarySize(linalg_traits<U>::row(mat_row_const_begin(u)));

  return N + 2*sizeof(size_t);  // account for size_t: mat_nrows, and mat_ncols in written info.
}

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, col_major, abstract_dense)
{    
  using commUtil::binarySize;
  
  // assume: all elements of matrix have same size: 
  size_t N_col(mat_ncols(u)), N(0);

  if (N_col > 0)
    N = N_col*binarySize(linalg_traits<U>::col(mat_col_const_begin(u)));

  return N + 2*sizeof(size_t);  // account for size_t: mat_nrows, and mat_ncols in written info.
}

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, row_major, abstract_skyline)
{    
  using commUtil::binarySize;
  
  // with skyline vectors, matrix rows may have differing sizes: 
  size_t N(0);

  for(typename linalg_traits<U>::const_row_iterator itRow = mat_row_const_begin(u), itRowEnd = mat_row_const_end(u);
      itRow != itRowEnd;
      ++itRow)
    N += binarySize(linalg_traits<U>::row(itRow));

  return N + 2*sizeof(size_t);  // account for size_t: mat_nrows, and mat_ncols in written info.
}

template <class U>
inline size_t binarySize(const U& u, abstract_matrix, col_major, abstract_skyline)
{    
  using commUtil::binarySize;
  
  // with skyline vectors, matrix columns may have differing sizes: 
  size_t N(0);

  for(typename linalg_traits<U>::const_col_iterator itCol = mat_col_const_begin(u), itColEnd = mat_col_const_end(u);
      itCol != itColEnd;
      ++itCol)
    N += binarySize(linalg_traits<U>::col(itCol));

  return N + 2*sizeof(size_t);  // account for size_t: mat_nrows, and mat_ncols in written info.
}

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
inline size_t binarySize<gmm::dense_matrix<double> >(const gmm::dense_matrix<double>& M)
{
 using commUtil::binarySize;
  
 size_t N_elt(mat_nrows(M)*mat_ncols(M)), N(0);
 if (N_elt > 0)
   N = N_elt * binarySize(M(0,0));
 return N + 2*sizeof(size_t); // account for size_t: mat_nrows, and mat_ncols in written info.
}

template <>
inline size_t binarySize<gmm::dense_matrix<std::complex<double> > >(const gmm::dense_matrix<std::complex<double> >& M)
{
 using commUtil::binarySize;
  
 size_t N_elt(mat_nrows(M)*mat_ncols(M)), N(0);
 if (N_elt > 0)
   N = N_elt * binarySize(M(0,0));
 return N + 2*sizeof(size_t); // account for size_t: mat_nrows, and mat_ncols in written info.
}

// The following are especially for factoredPropagator. It is very difficult to specialize 
//   these in a generic sense because of possibility of sub-vectors and sub-matrices with storage_type == abstract_skyline.
template <>
inline size_t binarySize<gmm::slvector<double> >(const gmm::slvector<double>& u)
{    
  using commUtil::binarySize;
  
  // assume: all elements of vector have same size: 
  size_t N_elt(vect_const_end(u)-vect_const_begin(u)), N(0);

  if (N_elt > 0)
    N = N_elt*binarySize(*vect_const_begin(u)); 

  return N + 3*sizeof(size_t); // account for size_t: vect_size, data_size, shift in written info (note: this case: _no_ redundant vect_size).
} 

template <>
inline size_t binarySize<gmm::slvector<std::complex<double> > >(const gmm::slvector<std::complex<double> >& u)
{    
  using commUtil::binarySize;
  
  // assume: all elements of vector have same size: 
  size_t N_elt(vect_const_end(u)-vect_const_begin(u)), N(0);

  if (N_elt > 0)
    N = N_elt*binarySize(*vect_const_begin(u)); 

  return N + 3*sizeof(size_t); // account for size_t: vect_size, data_size, shift in written info (note: this case: _no_ redundant vect_size).
} 

#if 0 // no specialization required, size on media same as for non-POD case:
template <>
inline size_t binarySize<gmm::row_matrix< gmm::slvector<double> > >(const gmm::row_matrix< gmm::slvector<double> >& M);

template <>
inline size_t binarySize<gmm::row_matrix< gmm::slvector<std::complex<double> > > >(const gmm::row_matrix< gmm::slvector<std::complex<double> > >& M);
#endif

#endif
            

} // namespace gmm


#endif // __linalgUtil_template__h
