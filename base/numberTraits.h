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

template <class T>
inline T pi(void)
{ throw std::runtime_error("number::pi<T>: not implemented"); }

// Default versions of the following template functions
//   assume *real* number class "T".

template <class T>
inline T one_i(void)
{ throw std::runtime_error("number::one_i<T>: not implemented"); }

template <class T>
inline typename numberTraits<T>::magnitudeType real(const T& t)
{ return t; }

template <class T>
inline typename numberTraits<T>::magnitudeType imag(const T& t)
{ return zero<T>(); }

template <class T>
inline T conj(const T& t)
{ return real<T>(t); }

template <class T0, class T1>
inline T0 conv(const T1& t1)
{ return static_cast<T0>(t1); }

template <class T>
inline T sqr(const T& t)
{ return t * t; }

template <class T>
inline typename numberTraits<T>::magnitudeType sqrNorm(const T& t)
{ return sqr(t); }

template <class T>
inline typename numberTraits<T>::magnitudeType fabs(const T& t)
{ throw std::runtime_error("number::fabs: *generic* absolute value not implemented"); }

template <class T>
inline T mod(const T& x, const T& y)
{ throw std::runtime_error("number::mod: *generic* modulus not implemented"); }

template <class T>
inline const T& min(const T& t0, const T& t1)
{ return std::min(t0, t1); }

template <class T>
inline const T& max(const T& t0, const T& t1)
{ return std::max(t0, t1); }

// Instance from *uniform* distribution: t \in [0, 1].
template <class T>
inline T random(void)
{ throw std::runtime_error("number::random: *generic* \"random\" not implemented"); }

template <class T>
inline T pow_n(const T& x, long n)
{ throw std::runtime_error("number::pow_n: not implemented"); }

template <>
inline long mod(const long& x, const long& y)
{ return x % y; }

#if defined(__USE_MERE)
template <>
mere::C integer<C,long>(const long& n);

template <>
const mere::C& zero<mere::C>(void);

template <>
const mere::C& epsilon<mere::C>(void);

template <>
inline mere::R sqrNorm(const mere::C& c)
{ return c * mere::conj(c); }

template <>
inline mere::C mod(const mere::C& x, const mere::C& y);

template <>
inline mere::C pow_n(const mere::C& c, long n)
{ return mere::pow_n(c, n); }

template <>
inline mere::C& one_i<mere::C>(void);

template <>
inline mere::R real<mere::C>(const mere::C& c);

template <>
inline mere::R imag<mere::C>(const mere::C& c);

template <>
inline mere::C conj<mere::C>(const mere::C& c);

template <>
inline mere::R fabs(const mere::C& c);
{ return mere::fabs(c); }

template <>
mere::R integer<R,long>(const long& n);

template <>
const mere::R& zero<mere::R>(void);

template <>
const mere::R& epsilon<mere::R>(void);

template <>
inline const mere::R& pi(void);

template <>
inline mere::R mod(const mere::R& x, const mere::R& y);

template <>
inline mere::R pow_n(const mere::R& r, long n)
{ return mere::pow_n(r, n); }

template <>
inline mere::R fabs(const mere::R& r);
{ return mere::fabs(r); }

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

template <>
inline double sqrNorm(const std::complex<double>& c)
{ return sqr(std::real(c)) + sqr(std::imag(c)); }

template <>
inline double fabs(const std::complex<double>& c)
{ return std::real(std::abs(c)); }

template <>
inline std::complex<double> one_i<std::complex<double> >(void)
{ return std::complex<double>(0.0, 1.0); }

template <>
inline double real(const std::complex<double>& c)
{ return std::real(c); }

template <>
inline double imag(const std::complex<double>& c)
{ return std::imag(c); }

template <>
inline std::complex<double> conj<std::complex<double> >(const std::complex<double>& c)
{ return std::conj(c); }

template <>
inline double mod(const double& x, const double& y)
{ return std::fmod(x, y); }

template <>
inline double pow_n(const double& r, long n)
{ return std::pow(r, static_cast<double>(n)); }

template <>
inline double pi(void)
{ return 3.14159265358979323846; }

template <>
inline double fabs(const double& r)
{ return std::fabs(r); }

#endif


} // namespace number

#if defined(__commUtil__h)
#if !defined(__USE_MERE)
// ----------- Specialize binary I/O for std::complex<double>, which is *not* a POD type: --------------

namespace commUtil{


template<>
inline bool writeBinary(abstractCommHandle *fp, const std::complex<double>& c)
{ return writeBinary_POD_dispatch_(fp, c, std::true_type()); }

template<>
inline bool readBinary(abstractCommHandle *fp, std::complex<double>& c)
{ return readBinary_POD_dispatch_(fp, c, std::true_type()); }


template <>
inline size_t binarySize(const std::complex<double>& c)
{ return binarySize_POD_dispatch_(c, std::true_type()); }


} // namespace commUtil

// -----------------------------------------------------------------------------------------------------
#endif
#endif

#endif // !defined(__numericalConstants__h)


#endif // __numberTraits__h
