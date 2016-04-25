#if !defined(__numberTraits__h)
#define __numberTraits__h

namespace number{


// -------------- allow this module to be independent of the TMatrix namespace: ------------------------
#if !defined(__numericalConstants__h)

template <class T> 
class numberTraits{
public:

#if 1
  typedef T complexType;   // complex type if T is real.
  typedef T magnitudeType; // magnitude type if T is complex.
  typedef T integralType;  // integer type of corresponding range.
#endif

};

// explicit specialization for double std::complex<double>, and std::complex<T>
template <> 
class numberTraits<double>{
public:
  typedef std::complex<double> complexType;   // complex type if T is real.  
  typedef double magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.

};

#if 1
template <> 
template <class T>
class numberTraits< std::complex<T> >{
public:
  typedef std::complex<T> complexType;   
  typedef T magnitudeType; 
  typedef typename numberTraits<T>::integralType integralType;  

};
#endif

#if defined(__USE_MERE)
// explicit specializations for mere::C, mere::R, mere::Z
template <> 
class numberTraits<mere::C>{
public:
  typedef mere::C complexType;   // complex type if T is real. 
  typedef mere::R magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
template <> 
class numberTraits<mere::R>{
public:
  typedef mere::C complexType;   // complex type if T is real.
  typedef mere::R magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
template <> 
class numberTraits<mere::Z>{
public:

  typedef mere::Z magnitudeType; // magnitude type if T is complex.
  typedef mere::Z integralType;  // integer type of corresponding range.

};
#endif

#if 0
template <> 
class numberTraits< std::complex<double> >{
public:
  typedef std::complex<double> complexType; // complex type if T is real.
  typedef double magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.
};
#endif
// *******************************

// explicit instantiation for long and unsigned long
template <> 
class numberTraits<long>{
public:

  typedef long magnitudeType; // magnitude type if T is complex.
  typedef long integralType;  // integer type of corresponding range.

};

template <> 
class numberTraits<unsigned long>{
public:

  typedef unsigned long magnitudeType; // magnitude type if T is complex.
  typedef unsigned long integralType;  // integer type of corresponding range.

};


template <class T, class N>
inline T integer(const N& n)
{ return static_cast<T>(n); }

template <class T, class N>
inline T ratio(const N& numerator, const N& denominator)
{ return static_cast<T>(numerator) / static_cast<T>(denominator); }

template <class T>
inline T zero(void)
{ return T(0); }

template <class T>
inline T one(void)
{ return T(1); }

template <class T>
inline T epsilon(void)
{ return std::numeric_limits<T>::epsilon(); }

template <class T0, class T1>
inline T0 conv(const T1& t1)
{ return static_cast<T0>(t1); }


#if defined(__USE_MERE)
template <>
mere::C integer<C,long>(const long& n);

template <>
const mere::C& zero<mere::C>(void);

template <>
const mere::C& epsilon<mere::C>(void);

template <>
mere::R integer<R,long>(const long& n);

template <>
const mere::R& zero<mere::R>(void);

template <>
const mere::R& epsilon<mere::R>(void);

#else
template <>
inline std::complex<double> integer<std::complex<double>,long>(const long& n)
{ return std::complex<double>(static_cast<double>(n)); }

template <>
inline  std::complex<double> zero<std::complex<double> >(void)
{ return std::complex<double>(0.0); }

template <>
inline std::complex<double> epsilon<std::complex<double> >(void)
{ return std::complex<double>(epsilon<double>()); }
#endif


} // namespace number
#endif // !defined(__numericalConstants__h)


#endif // __numberTraits__h
