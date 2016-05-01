#if !defined(__ntuple_template__h)
#define __ntuple_template__h

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



template <class T, size_t DIM>
inline T& ntuple<T,DIM>::operator[](size_t n)
{ 
 assert( n < DIM );
 return v_[n]; 
}

template <class T, size_t DIM>
inline const T& ntuple<T,DIM>::operator[](size_t n)const
{ 
 assert( n < DIM );
 return v_[n]; 
}



// allow use with STL iterator methods:
template <class T, size_t DIM>
inline T* ntuple<T,DIM>::begin()
{ return &(v_[0]); }

template <class T, size_t DIM>
inline const T* ntuple<T,DIM>::begin()const
{ return &(v_[0]); }

template <class T, size_t DIM>
inline T* ntuple<T,DIM>::end()
{ return &(v_[DIM]); }

template <class T, size_t DIM>
inline const T* ntuple<T,DIM>::end()const
{ return &(v_[DIM]); }
		

template <class T, size_t DIM>
inline size_t ntuple<T,DIM>::size(void)const
{ return DIM; }

template <class T, size_t DIM>
inline void ntuple<T,DIM>::resize(size_t n)
{ 
  if (n != DIM)
    throw std::runtime_error("ntuple<T,DIM>::resize: size not equal to fixed ntuple dimension");
} 

template <class T, size_t DIM>
typename numberTraits<T>::magnitudeType ntuple<T,DIM>::dist(const ntuple& other)const
{
 using number::sqrNorm;
 R val(zero<R>());
 for (size_t n=0; n<DIM; ++n)
	 val += sqrNorm(v_[n] - other[n]);
 return sqrt(val);	 
}
		

template <class T, size_t DIM>
inline bool ntuple<T,DIM>::operator==(const ntuple<T,DIM>& other)const
{
  bool test(true);
  for (size_t n=0; n<DIM; ++n)
    if (!(v_[n] == other[n])){
      test = false;
      break;
    }
  return test;  
}

template <class T, size_t DIM>
inline bool ntuple<T,DIM>::operator!=(const ntuple<T,DIM>& other)const
{
  return !operator==(other);
}

// inner product
template <class T, size_t DIM>
T ntuple<T,DIM>::operator*(const ntuple& other)const
{
 T val(zero<T>());
 for (size_t n=0; n<DIM; ++n)
	 val += v_[n] * other[n];
 return val;	 
}

    // magnitude
template <class T, size_t DIM>		
inline typename numberTraits<T>::magnitudeType ntuple<T,DIM>::abs(void)const
{
 return sqrt(absSqr());	 
}
		
template <class T, size_t DIM>		
typename numberTraits<T>::magnitudeType ntuple<T,DIM>::absSqr(void)const
{
 using number::sqrNorm;
 R val(zero<R>());
 for (size_t n=0; n<DIM; ++n)
	 val += sqrNorm(v_[n]);
 return val;	 
}

template <class T, size_t DIM>		
inline size_t ntuple<T,DIM>::numberNonzeroElements(const R& eps)const
{
 size_t rval(0);
 for(size_t n=0; n<DIM; ++n)
   if (fabs(v_[n]) > eps) 
     ++rval;
 return rval;
}

template <class T, size_t DIM>		
inline size_t ntuple<T,DIM>::nnz(const R& eps)const
{ return numberNonzeroElements(eps); }

// bit-mask corresponding to non-zero elements:
template <class T, size_t DIM>	
inline size_t ntuple<T,DIM>::mask_nz(void)const
{
 size_t bits(0);
 for(size_t n=0; n < DIM; ++n)
   if (v_[n] != zero<T>()) bits |= (1U<<n);
 return bits; 
}
		
template <class T, size_t DIM>		
void ntuple<T,DIM>::normalize(void)
{
	R norm_(abs());
	assert(norm_ > epsilon<R>());
	for (size_t n=0; n<DIM; ++n)
		v_[n]/= norm_;
}

template <class T, size_t DIM>		
void ntuple<T,DIM>::invert(bool spherical)
{
  using number::mod;
  if (!spherical)
	  for (size_t n=0; n<DIM; ++n)
		  v_[n] = -v_[n];   
	else{
	  assert(DIM == 3);
    v_[1] = pi<T>() - v_[1];
    v_[2] = mod(pi<T>() + v_[2], integer<T>(2)*pi<T>());    
  }  
}

template <class T, size_t DIM>		
void ntuple<T,DIM>::clear(void)
{
	for (size_t n=0; n<DIM; ++n)
		v_[n] = zero<T>();
}
		
template <class T, size_t DIM>		
inline ntuple<T,DIM>& ntuple<T,DIM>::operator+(void) { return *this; }

template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::operator-(void)const 
{
 ntuple<T,DIM> val(*this);
 val *= -one<T>();
 return val;
} 
							
template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::operator+(const ntuple<T,DIM>& other)const
{
 ntuple<T,DIM> val(*this);
 val += other;
 return val;
}

template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::operator-(const ntuple<T,DIM>& other)const
{
 ntuple<T,DIM> val(*this);
 val -= other;
 return val;
}
		
template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::operator*(const T& r)const
{
 ntuple<T,DIM> val(*this);
 val *= r;
 return val;
}

template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::operator/(const T& r)const
{
 ntuple<T,DIM> val(*this);
 val /= r;
 return val;
}

template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::mod(const T& r)const
{
 ntuple<T,DIM> val(*this);
 val.mod_assign(r);
 return val;
}

template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::operator+=(const ntuple<T,DIM>& other)
{		 
 for(size_t n=0; n<DIM; ++n)
	 v_[n] += other[n];

 return *this;	 
}
		
template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::operator-=(const ntuple<T,DIM>& other)
{
 for(size_t n=0; n<DIM; ++n)
	 v_[n] -= other[n];

 return *this;	 
}

template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::operator*=(const T& r)
{
 for(size_t n=0; n<DIM; ++n)
	 v_[n] *= r;

 return *this;	 
}

template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::operator/=(const T& r)
{
 for(size_t n=0; n<DIM; ++n)
	 v_[n] /= r;

 return *this;	 
}	

// set all elements to scalar value:
// (implementation note: default assigment op (i.e. ntuple& operator=(const ntuple& other)) is automatically generated by compiler)
template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::operator=(const T& r)
{
  for(size_t n=0; n<DIM; ++n)
    v_[n] = r;
  return *this;
}

template <class T, size_t DIM>		
ntuple<T,DIM>& ntuple<T,DIM>::mod_assign(const T& r)
{
 using number::mod;
 for(size_t n=0; n<DIM; ++n)
	 v_[n] = mod(v_[n], r);

 return *this;	 
}
	
template <class T, size_t DIM>
T ntuple<T,DIM>::dot(const ntuple<T,DIM>& other)const
{
 T val(zero<T>());
 for(size_t n=0; n<DIM; ++n)
   val += v_[n]*other.v_[n];
 return val;
}

template <class T, size_t DIM>
template <class S>
T ntuple<T,DIM>::dot(const ntuple<S,DIM>& other)const
{
 T val(zero<T>());
 for(size_t n=0; n<DIM; ++n)
   val += v_[n]*other[n];
 return val;
}

template <class T, size_t DIM>		
ntuple<T,DIM> ntuple<T,DIM>::cross(const ntuple<T,DIM>& other)const
{
 if (DIM != 3)
   throw std::runtime_error("ntuple<T,DIM>::cross: not implemented for DIM != 3");
 
 return
 ntuple<T,DIM>(
   v_[1]*other.v_[2] - v_[2]*other.v_[1],
	 -v_[0]*other.v_[2] + v_[2]*other.v_[0],
	 v_[0]*other.v_[1] - v_[1]*other.v_[0] 
 );
}
		

template <class T, size_t DIM>		
template <class S>
ntuple<T,DIM> ntuple<T,DIM>::cross(const ntuple<S,DIM>& other)const
{
 if (DIM != 3)
   throw std::runtime_error("ntuple<T,DIM>::cross: not implemented for DIM != 3");
 
 return
 ntuple<T,DIM>(
   v_[1]*other[2] - v_[2]*other[1],
	 -v_[0]*other[2] + v_[2]*other[0],
	 v_[0]*other[1] - v_[1]*other[0] 
 );
}

/*
 * ntuple<T,DIM> unit()
 *   ntuple of unit vector (i.e. normalized ntuple, _not_ ntuple with unit elements)
 */	
template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::unit(void)
{ return filled(one<T>() / sqrt(integer<T>(DIM))); }

/*
 * ntuple<T,DIM> filled(const T& t)
 *   ntuple with constant elements)
 */
template <class T, size_t DIM>		
inline ntuple<T,DIM> ntuple<T,DIM>::filled(const T& t)
{
 ntuple<T,DIM> val;  
 for(size_t n=0; n<DIM; ++n)
	 val[n] = t;
 return val;	 
}

template <class T, size_t DIM>
inline void ntuple<T,DIM>::sort(ntuple<T,DIM>& p)
{
 std::sort( p.begin(), p.end() );
}
    
#if 0 // --------------------------- obsolete: use linalgUtil.h: "argsort" ------------------------------------
template <class T, size_t DIM>
void ntuple<T,DIM>::index_sort(const ntuple<T,DIM>& p, ntuple<size_t,DIM>& vn)
{
 ntuple<T,DIM> p_(p); // copy it
 
 _STL_EXT_NAMESPACE_::iota(vn.begin(), vn.end(), 0);
	 
 typedef linalg::ghostIterator< T*, size_t* > combined_data_iterator;
 // Sort both p_ and vn based on ordering of p_:
 std::sort( combined_data_iterator(p_.begin(), vn.begin()), 
            combined_data_iterator(p_.end(), vn.end()) );
}
#endif // --------------------------- end: obsolete ----------------------------------------------------------

template <class T, size_t DIM>
void ntuple<T,DIM>::apply( T (*func)(const T&), const ntuple<T,DIM>& src, ntuple<T,DIM>& dest)
{
 for(size_t n=0; n<DIM; ++n)
   dest.v_[n] = func(src.v_[n]); 
}
		
template <class T, size_t DIM>		
bool ntuple<T,DIM>::writeBinary(std::ostream &out)const
{
 for(size_t n=0; out && n<DIM; ++n)
   out.write(reinterpret_cast<const char*>(&v_[n]), sizeof(T));

 return out.good();	 
}

template <class T, size_t DIM>		
bool ntuple<T,DIM>::readBinary(std::istream &in)
{
 for(size_t n=0; in && n<DIM; ++n)
   in.read(reinterpret_cast<char *>(&v_[n]), sizeof(T));
   
 return in.good();	 
}

template <class T, size_t DIM>		
size_t ntuple<T,DIM>::binarySize(void)
{ return DIM * sizeof(T); }

template <class T, size_t DIM>		
void ntuple<T,DIM>::write(std::ostream& os)const
{
 // width for "(" DIM*< <digit> " "|")" >
 size_t itemWidth( (os.width() - 1)>static_cast<int>(DIM)? ((os.width() - 1)/DIM - 1): 0 );

 os << "(";
 for(size_t n=0; os && (n<DIM); ++n){
   os.width(itemWidth);
   os << v_[n];
	 os << (n<(DIM-1)?", ":")");  // 06.2009: new-format is (<number>, <number>, <number>) -- see if this breaks anything...
 }
}

template <class T, size_t DIM>		
void ntuple<T,DIM>::read(std::istream& is)
{
 char ch;
 while(is){
   is.get(ch);
	 if (std::isspace(ch)) continue;
	 if (ch == '(')
	   break;
	 else{
	   is.setstate(std::ios::failbit);
		 break;
	 }
 }

 for(size_t n=0; is && (n<DIM); ++n){
   is >> v_[n];
   if (n < DIM-1)
	 while(is){
		is.get(ch);
		if (std::isspace(ch)) continue;
    // allow alternative format: "(" <number>, <number>, <number> ")":
		if (ch == ',')
			break;
    // allow alternative format: "(" <number> <number> <number> ")":
    if (std::isdigit(ch) || (ch=='+') || (ch=='-')){
      is.putback(ch);
      break;
    }  
		else{
			is.setstate(std::ios::failbit);
			break;
		}
	 }
 }

 while(is){
	is.get(ch);
	if (std::isspace(ch)) continue;
	if (ch == ')')
		break;
	else{
		is.setstate(std::ios::failbit);
		break;
	}
 }  	 
}

template <class T, size_t DIM>		
void ntuple<T,DIM>::debug_print(void)const
{
 write(std::cout);
 std::cout<<std::endl;
}

template <class T, size_t DIM>		
ntuple<T,DIM>::ntuple(bool clear)
{
 if (clear)
   for(size_t n=0; n<DIM; ++n)
	   v_[n] = zero<T>();
}

// allow "ntuple" to function as a gmm abstract_vector:
template <class T, size_t DIM>		
ntuple<T,DIM>::ntuple(size_t DIM_)
{
 assert(DIM_ == DIM);
 for(size_t n=0; n<DIM; ++n)
	 v_[n] = zero<T>();
}

template <class T, size_t DIM>		
ntuple<T,DIM>::ntuple(const T& e1)
{
 assert(DIM >= 1);
 v_[0] = e1;
 for(size_t n=1; n<DIM; ++n)
	 v_[n] = zero<T>();
}
		
template <class T, size_t DIM>		
ntuple<T,DIM>::ntuple(const T& e1, const T& e2)
{
 assert(DIM >= 2);
 v_[0] = e1;
 v_[1] = e2;
 for(size_t n=2; n<DIM; ++n)
	 v_[n] = zero<T>();
}

template <class T, size_t DIM>		
ntuple<T,DIM>::ntuple(const T& e1, const T& e2, const T& e3)
{
 assert(DIM >= 3);
 v_[0] = e1;
 v_[1] = e2;
 v_[2] = e3;
 for(size_t n=3; n<DIM; ++n)
	 v_[n] = zero<T>();
}

#if !defined(__ICC) && !defined(__PGI) // GNU: place body _here_, otherwise at primary declaration:							
template <class T, size_t DIM>		
template <class T1>
ntuple<T,DIM>::ntuple(const T1& e1, const T1& e2, const T1& e3)
{
 assert(DIM >= 3);
 conv(v_[0],e1);
 conv(v_[1],e2);
 conv(v_[2],e3);
 for(size_t n=3; n<DIM; ++n)
	 v_[n] = zero<T>();
}
#endif

template <class T, size_t DIM>		
inline ntuple<T,DIM> operator*(const T& r, const ntuple<T,DIM>& p)
{
 ntuple<T,DIM> val(p);
 val *= r;
 return val;
}

template <class T, size_t DIM>		
inline typename numberTraits<T>::magnitudeType dist(const ntuple<T,DIM>& p1, const ntuple<T,DIM>& p2)
{ return p1.dist(p2); }

template <class T, size_t DIM>		
inline typename numberTraits<T>::magnitudeType absSqr(const ntuple<T,DIM>& p)
{ return p.absSqr(); }


template <class T, size_t DIM>		
inline size_t numberNonzeroElements(const ntuple<T,DIM>& p, const T& eps)
{
  return p.numberNonzeroElements(eps);
}

// true if _all_ values are _exactly_ zero:
template <class T, size_t DIM>
inline bool isZero(const ntuple<T,DIM>& p, const T& eps)
{
 bool val(true);

 for(size_t n=0; val && (n<DIM); ++n)
   if (fabs(p[n]) > eps){
     val = false;
	   break;
   }

 return val;	 
}

// true if _any_ values are NAN:
template <class T, size_t DIM>
inline bool isNAN(const ntuple<T,DIM>& p)
{
 bool val(false);

 for(size_t n=0; !val && (n<DIM); ++n)
 if (isnan(p[n])){
   val = true;
	 break;
 }	 

 return val;	 
}

// in analogy to std::min_element returning an iterator, these methods return indices:
template <class T, size_t DIM>
inline size_t min_element(const ntuple<T,DIM>& p)
{
 size_t val(0);
 T rMin(p[0]);
 for(size_t n=0; n<DIM; ++n)
 if (p[n] < rMin){
   rMin = p[n];
	 val = n;
 }
 
 return val;
}

template <class T, size_t DIM>
inline size_t max_element(const ntuple<T,DIM>& p)
{
 size_t val(0);
 T rMax(p[0]);
 for(size_t n=0; n<DIM; ++n)
 if (p[n] > rMax){
   rMax = p[n];
	 val = n;
 }
 
 return val;
}

template <class T, size_t DIM>
inline void sort(ntuple<T,DIM>& p)
{
 ntuple<T,DIM>::sort(p);
}

template <class T, size_t DIM>
inline void index_sort(const ntuple<T,DIM>& p, ntuple<size_t,DIM>& vn)
{
 ntuple<T,DIM>::index_sort(p, vn);
}

template <class T1, class T2, size_t DIM>
inline void conv(ntuple<T1,DIM>& dest, const ntuple<T2,DIM>& src)
{
 for(size_t n=0; n<DIM; ++n)
   conv(dest[n], src[n]);
}

template <class T, size_t DIM>
inline bool writeBinary(std::ostream &out, const ntuple<T,DIM>& p) { return p.writeBinary(out); }

template <class T, size_t DIM>
inline bool readBinary(std::istream &in, ntuple<T,DIM>& p) { return p.readBinary(in); }

template <class T, size_t DIM>
inline std::ostream& operator<<(std::ostream& os, const ntuple<T,DIM>& p) { p.write(os); return os; }

template <class T, size_t DIM>
inline std::istream& operator>>(std::istream& is, ntuple<T,DIM>& p) { p.read(is); return is; }


// -------------------------- ntuple_interval: -----------------------------------------------------------------
// interval class for ntuple:

template <class T, size_t DIM>
inline const ntuple<T,DIM>& ntuple_interval<T,DIM>::start(void)const
{ return start_; }

template <class T, size_t DIM>
inline ntuple<T,DIM>& ntuple_interval<T,DIM>::start(void)
{ return start_; }

template <class T, size_t DIM>
inline const ntuple<T,DIM>& ntuple_interval<T,DIM>::end(void)const
{ return end_; }

template <class T, size_t DIM>
inline ntuple<T,DIM>& ntuple_interval<T,DIM>::end(void)
{ return end_; }

template <class T, size_t DIM>
inline const T& ntuple_interval<T,DIM>::left_epsilon(void)const
{ return left_epsilon_; }

template <class T, size_t DIM>
inline const T& ntuple_interval<T,DIM>::right_epsilon(void)const
{ return right_epsilon_; }

/* test if point p is contained in associated sets:
 *
 *  rules:  
 *    - R^n == <interior> U <exterior>
 *    - <closure> == <interior> U <boundary>     
 *    - <closure> == <left-closure> U <right-closure>
 *    .
 */
template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::interior(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!((p[n] >= start_[n] + left_epsilon_) && (p[n] <= end_[n] - right_epsilon_))){
      test = false;
      break;
    }
  return test; 
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::boundary(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!(((p[n] >= start_[n] - left_epsilon_) && (p[n] < start_[n] + left_epsilon_))
       || ((p[n] > end_[n] - right_epsilon_) && (p[n] <= end_[n] + right_epsilon_)))){
      test = false;
      break;
    }
  return test; 
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::exterior(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!((p[n] < start[n] - left_epsilon_) || (p[n] > end_[n] + right_epsilon_))){
      test = false;
      break;
    }
  return test; 
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::closure(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!((p[n] >= start_[n] - left_epsilon_) && (p[n] <= end_[n] + right_epsilon_))){
      test = false;
      break;
    }
  return test; 
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::left_closure(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!((p[n] >= start_[n] - left_epsilon_) && (p[n] <= end_[n] - right_epsilon_))){
      test = false;
      break;
    }
  return test; 
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::right_closure(const ntuple<T,DIM>& p)const
{
  bool test(true);
  for(size_t n = 0; n < DIM; ++n)
    if (!((p[n] >= start_[n] + left_epsilon_) && (p[n] <= end_[n] + right_epsilon_))){
      test = false;
      break;
    }
  return test; 
}


template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::operator==(const ntuple_interval<T,DIM>& other)const
{
  return (start_ == other.start_) && (end_ == other.end_) 
    && (left_epsilon_ == other.left_epsilon_) && (right_epsilon_ == other.right_epsilon_);
}

template <class T, size_t DIM>
inline bool ntuple_interval<T,DIM>::operator!=(const ntuple_interval<T,DIM>& other)const
{
  return !operator==(other);
}

template <class T, size_t DIM>
bool ntuple_interval<T,DIM>::writeBinary(std::ostream &out)const
{
  bool status(false);
  
  status = (status && start_.writeBinary(out));
  status = (status && end_.writeBinary(out));
  status = (status && out.write(reinterpret_cast<const char *>(&left_epsilon_), sizeof(T)));
  status = (status && out.write(reinterpret_cast<const char *>(&right_epsilon_), sizeof(T)));
  return status;
}

template <class T, size_t DIM>
bool ntuple_interval<T,DIM>::readBinary(std::istream &in)
{
  bool status(false);
  
  status = (status && start_.readBinary(in));
  status = (status && end_.readBinary(in));
  status = (status && in.read(reinterpret_cast<char *>(&left_epsilon_), sizeof(T)));
  status = (status && in.read(reinterpret_cast<char *>(&right_epsilon_), sizeof(T)));
  return status;
}

template <class T, size_t DIM>
size_t ntuple_interval<T,DIM>::binarySize(void)const
{
  size_t val(0);
  
  val += start_.binarySize();
  val += end_.binarySize();
  val += sizeof(T);
  val += sizeof(T);
  return val;
}

template <class T, size_t DIM>
void ntuple_interval<T,DIM>::write(std::ostream& os, const ntuple<T,DIM>* p)const
{
  using number::epsilon;
  using std::setw;
  using std::left;
  using std::right;
  
  start_.write(os);
  os<<" -> ";
  end_.write(os);
  os<<"(left_eps, right_eps: "<<left_epsilon_<<", "<<right_epsilon_<<")";
  if (p != NULL){
    os<<"\n"
      <<"  normalized boundary differences for ntuple "<<*p<<"(units of epsilon<T>()):\n"
      <<"  DIM start           end             \n";  // widths:2, DIM:4, start:16, end:16
    for(size_t n = 0; n < DIM; ++n){
      os<<"  "
        <<setw(4)<<left<<n
        <<setw(16)<<left<</* fabs */((*p)[n] - start_[n])/epsilon<T>()
        <<setw(16)<<left<</* fabs */ ((*p)[n] - end_[n])/epsilon<T>()
        <<right<<"\n";              
    }  
  }
}

/**
 * @brief  Calculate scale tuple corresponding to this ntuple-interval (i.e. as end() - start()).
 */
template <class T, size_t DIM>
inline ntuple<T,DIM> ntuple_interval<T,DIM>::scale(void)const
{ return end() - start(); }

#if !defined(__PGI) && !defined(__INTEL_COMPILER)
/**
  * @brief Calculate extent of a container of ntuple<T,DIM> as ntuple_interval<T,DIM>
  * Assumes typename U::value_type is ntuple<T,DIM>
  */
template <class T, size_t DIM>
template <class U>
ntuple_interval<T,DIM> ntuple_interval<T,DIM>::extent(const U& u, const T& left_eps, const T& right_eps)
{
  // empty container has zero extent:
  ntuple_interval val(ntuple<T,DIM>(true), ntuple<T,DIM>(true), left_eps, right_eps);

  if (!u.empty()){
    val.start() = *u.begin();
    val.end() = *u.begin();
    for(typename U::const_iterator itU = u.begin(), itUEnd = u.end();
        itU != itUEnd;
        ++itU)
    for(size_t n = 0; n < DIM; ++n){
      const T& test((*itU)[n]);
      if (test < val.start()[n])
        val.start()[n] = test;
      if (test > val.end()[n])
        val.end()[n] = test;
    }
  }

  return val;
}

/**
  * @brief Calculate extent of a container of ntuple<T,DIM> as ntuple_interval<T,DIM>
  * Assumes typename U::value_type is ntuple<T,DIM>
  *   (this signature assumes left_epsilon == right_epsilon)
  */
template <class T, size_t DIM>
template <class U>
inline ntuple_interval<T,DIM> ntuple_interval<T,DIM>::extent(const U& u, const T& eps)
{ return ntuple_interval<T,DIM>::template extent<U>(u, eps, eps); }
    
#endif

template <class T, size_t DIM>
ntuple_interval<T,DIM>::ntuple_interval(const ntuple<T,DIM>& start, const ntuple<T,DIM>& end, const T& left_eps, const T& right_eps)
  : start_(start), end_(end), left_epsilon_(left_eps), right_epsilon_(right_eps)
{ }

//! this signature assumes left_eps == right_eps:
template <class T, size_t DIM>
ntuple_interval<T,DIM>::ntuple_interval(const ntuple<T,DIM>& start, const ntuple<T,DIM>& end, const T& eps)         
  : start_(start), end_(end), left_epsilon_(eps), right_epsilon_(eps)
{ }

template <class T, size_t DIM>
ntuple_interval<T,DIM>::ntuple_interval(void)
  : start_(true), end_(true), left_epsilon_(zero<T>()), right_epsilon_(zero<T>())
{ }

template <class T, size_t DIM>
inline bool writeBinary(std::ostream &out, const ntuple_interval<T,DIM>& p)
{ return p.writeBinary(out); }

template <class T, size_t DIM>
inline bool readBinary(std::istream &in, ntuple_interval<T,DIM>& p)
{ return p.readBinary(in); }

template <class T, size_t DIM>
inline std::ostream& operator<<(std::ostream& os, const ntuple_interval<T,DIM>& p)
{ 
  p.write(os);
  return os;
}
      
// --------------------------------------------------------------------------------------------------------------


} // namespace linalg

#endif // __ntuple_template__h
