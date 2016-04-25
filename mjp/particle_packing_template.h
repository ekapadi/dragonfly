#if !defined(__particle_packing_template__h)
#define __particle_packing_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/particle_packing_template.h,v $

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
template <class T, class N>
T ratio(const N& numerator, const N& denominator)
{
  T val(integer<T,N>(numerator));
  val /= integer<T,N>(denominator);
  return val;
}
#endif

// allow this module to be independent of namespace TMatrix:
//   privately define static "readBinary_", "writeBinary_", and "binarySize_" methods for use within this class.

#include "binaryIO_stub_template.h"
#undef _FULLY_QUALIFIED_CLASS_NAME_
#undef _OUTER_TEMPLATE_PREFIX_

#if 1
// definitions in "binaryIO_stub.h" do not compile, for some reason:
template <class R, size_t NDIM>   
inline size_t system<R,NDIM>::binarySize_(const ntuple<size_t,NDIM>& t)
{ return NDIM * binarySize_(t[0]); }

template <class R, size_t NDIM>   
inline size_t system<R,NDIM>::binarySize_(const ntuple<long,NDIM>& t)
{ return NDIM * binarySize_(t[0]); }

template <class R, size_t NDIM>   
inline size_t system<R,NDIM>::binarySize_(const ntuple<R,NDIM>& t)
{ return NDIM * binarySize_(t[0]); }

template <class R, size_t NDIM>   
inline size_t system<R,NDIM>::binarySize_(const ntuple_interval<R,NDIM>& t2)    
{ 
  size_t val(0);
  val += binarySize_(t2.start());
  val += binarySize_(t2.end());
  val += binarySize_(t2.epsilon());
  return val;
}
#endif

#if !defined(__numericalFunctor__h)
    
/**
 * @brief Solve the quadratic formula:
 * @tparam V1: container type of value_type real or complex
 * @tparam V2: container type of value_type complex (or value_type real \f$ \Leftrightarrow \f$ roots are somehow constrained to be real)
 * @param[in]: vP: vector of polynomial coefficients
 * @param[out]: vRoot: vector of roots
 * @param[out]: D: discriminant
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

#endif

// ---------------------------------- end: section for TMatrix namespace independence -------------------------

#if defined(__USE_PTHREAD)

// pair lock and unlock methods:

/*
 * blocks until pair can be acquired simultaneously (keyed to addresses of p1 and p2);
 *   returns pair of keys to pass to matching "unlock_pair" method
 */
template <class R, size_t NDIM> 
std::pair<size_t, size_t> system<R,NDIM>::lock_pair_(void *p0, void *p1)
{
  bool pair_locked(true);
  std::pair<size_t, size_t> keys;
  do {
    pair_locked = true;
    if (pthread_mutex_lock(&pair_mutex_))
      throw std::runtime_error("system<R,NDIM>::lock_pair: pthread_mutex_lock error return");
      
    pair_locked = (pair_locked && particle0_mutex_.trylock(p0, keys.first));
    pair_locked = (pair_locked && particle1_mutex_.trylock(p1, keys.second));
    
    if (!pair_locked){
      // wait for another thread to release a pair:
      if (pthread_cond_wait(&pair_cond_, &pair_mutex_))
        throw std::runtime_error("system<R,NDIM>::lock_pair: pthread_cond_signal error return");        
    }
      
    if (pthread_mutex_unlock(&pair_mutex_))
      throw std::runtime_error("system<R,NDIM>::lock_pair: pthread_mutex_unlock error return");      
  } while(!pair_locked);
  
  return keys;
}

/*
 * uses pair of keys from matching "lock_pair" method to release aquired pair.
 */
template <class R, size_t NDIM> 
void system<R,NDIM>::unlock_pair_(const std::pair<size_t, size_t>& keys)
{
  if (pthread_mutex_lock(&pair_mutex_))
    throw std::runtime_error("system<R,NDIM>::unlock_pair: pthread_mutex_lock error return");
    
  particle0_mutex_.unlock(keys.first);
  particle1_mutex_.unlock(keys.second);
  
  if (pthread_mutex_unlock(&pair_mutex_))
    throw std::runtime_error("system<R,NDIM>::unlock_pair: pthread_mutex_unlock error return"); 
         
  if (pthread_cond_broadcast(&pair_cond_))
    throw std::runtime_error("system<R,NDIM>::unlock_pair: pthread_cond_broadcast error return");        
}  

#endif
// ------------------------------------------------------------

template <class R, size_t NDIM>   
void system<R,NDIM>::parameters::update_cache(bool derived, bool to_cache)
{
  // usage note: _local_ copies of parameters are the "cache"

  base_class::update_cache(true, to_cache); 

  if (to_cache){
    L = base_class::template get_named_parm<R>("L");
    
    // "N_particle" is a mandatory parameter:
    N_particle = static_cast<size_t>(base_class::template get_named_parm<Z>("N_particle"));

    #if 0
    if (base_class::has_named_parm("t_max"))
      t_max = base_class::template get_named_parm<R>("t_max");
    else
      t_max = zero<R>();
    #endif

    if (base_class::has_named_parm("NSTEP"))
      NSTEP = base_class::template get_named_parm<Z>("NSTEP");
    else
      NSTEP = static_cast<size_t>(-1);

    if (base_class::has_named_parm("sticking_probability"))
      sticking_probability = base_class::template get_named_parm<R>("sticking_probability");
    else
      sticking_probability = zero<R>();
      
    #if 0
    particle_per_cell = static_cast<size_t>(base_class::template get_named_parm<Z>("particle_per_cell"));
    N_seed = static_cast<size_t>(base_class::template get_named_parm<Z>("N_seed"));
    #endif
  }
  else{
    base_class::template set_named_parm<R>("L", L, true);
    base_class::template set_named_parm<Z>("N_particle", static_cast<Z>(N_particle), true);
    #if 0
    base_class::template set_named_parm<R>("t_max", static_cast<R>(t_max), true);
    #endif
    base_class::template set_named_parm<Z>("NSTEP", static_cast<Z>(NSTEP), true);
    base_class::template set_named_parm<R>("sticking_probability", sticking_probability, true);
    #if 0
    base_class::template set_named_parm<Z>("particle_per_cell", static_cast<Z>(particle_per_cell), true);
    base_class::template set_named_parm<Z>("N_seed", static_cast<Z>(N_seed), true);  
    #endif
  }
  
  // most-derived class sets currency flag:
  if (!derived)
    base_class::set_cache_current(true);
}
        

template <class R, size_t NDIM>   
void system<R,NDIM>::parameters::valid_check(void)const
{
  // cache update should have happened prior to this method
  if (!base_class::cache_current())
    throw std::runtime_error("system<R,NDIM>::parameters::valid_check: cached parameters not current");
  if ((N_particle < 1)
     #if 0 
      || (particle_per_cell < 1)
     #endif
      || (L < zero<R>())
      || (base_class::has_named_parm("t_max") && (base_class::template get_named_parm<R>("t_max") < zero<R>()))
      || (sticking_probability < zero<R>())
      || (sticking_probability > integer<R>(1)))
    throw std::runtime_error("system<R,NDIM>::parameters::valid_check: invalid parameters");
}

template <class R, size_t NDIM>   
typename python_util::options_map<typename system<R,NDIM>::C>* system<R,NDIM>::parameters::clone(void)const
{ return new parameters(*this); }

// ------- just wrap base_class::copy -----------------------------
template <class R, size_t NDIM>   
inline void system<R,NDIM>::parameters::copy(const parameters& other)
{
  base_class::copy(other);
}
// ----------------------------------------------------------------

template <class R, size_t NDIM>   
typename system<R,NDIM>::parameters& system<R,NDIM>::parameters::operator=(const parameters& other)
{
  copy(other);
  return *this;
}

template <class R, size_t NDIM>   
typename system<R,NDIM>::parameters& system<R,NDIM>::parameters::operator=(const python_util::options_map<R>& other)
{
  base_class::copy(other);
  return *this;
}

template <class R, size_t NDIM>   
void system<R,NDIM>::parameters::write(std::ostream& os)const
{
  os<<"system<R,NDIM::parameters: \n";
  base_class::template write<C,R,Z>(os);
  os<<"\n";
}

        
/** 
 * @brief Initialize from a simple_object_base*.
 */
template <class R, size_t NDIM>   
inline void system<R,NDIM>::parameters::extract(const python_util::simple_object_base* src)
{ 
  python_util::extract(*reinterpret_cast<base_class*>(this), src);
  update_cache();
  valid_check();
}

template <class R, size_t NDIM>   
system<R,NDIM>::parameters::~parameters(void)
{ }

template <class R, size_t NDIM>   
system<R,NDIM>::parameters::parameters(bool derived)
  : base_class(true),  // true => allocate simple_object_base*
    L(integer<R>(1)),
    N_particle(1),
    #if 0
    t_max(zero<R>()),
    #endif
    NSTEP(static_cast<size_t>(-1)),
    sticking_probability(zero<R>())
    #if 0
    , particle_per_cell(1),
    N_seed(0)
    #endif
{
  update_cache(derived, false); // transfer from cache to options_map; most-derived sets cache currency flag 
}

template <class R, size_t NDIM>   
system<R,NDIM>::parameters::parameters(const parameters& other)
{
  copy(other);
}

template <class R, size_t NDIM>   
system<R,NDIM>::parameters::parameters(const python_util::options_map<C>& other)
{
  base_class::copy(other);
}



#if !defined(__ND_matrix__h)         
// ---------------------- moved to namespace linalg: (ND_matrix.h): -----------------------------

template <class R, size_t NDIM>   
template <class N1>     
inline N1 system<R,NDIM>::row_major_index::product_(const ntuple<N1,NDIM>& t, size_t offset)
{
 N1 P(static_cast<N1>(1));
 for(size_t n = offset; n < NDIM; ++n)
   P *= t[n];
 return P;
}
  
template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index system<R,NDIM>::row_major_index::operator+(const row_major_index& other)const
{ 
  row_major_index val(*this);
  val += other;
  return val;
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index system<R,NDIM>::row_major_index::operator-(const row_major_index& other)const
{ 
  row_major_index val(*this);
  val -= other;
  return val;
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index& system<R,NDIM>::row_major_index::operator+=(const row_major_index& other)
{ 
  base_class::operator+=(other); 
  return *this;
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index& system<R,NDIM>::row_major_index::operator-=(const row_major_index& other)
{ 
  base_class::operator-=(other); 
  return *this;
}

// apply periodic B.C. as constrained by shape:
template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index& system<R,NDIM>::row_major_index::mod_assign(const shape_type& shape)
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

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::row_major_index& system<R,NDIM>::row_major_index::operator=(const row_major_index& other)
{ 
  base_class::operator=(other); 
  return *this;
}

template <class R, size_t NDIM>        
template <class N1>
inline typename system<R,NDIM>::row_major_index& system<R,NDIM>::row_major_index::operator=(const ntuple<N1,NDIM>& other)
{ 
  return base_class::operator=(other); 
  return *this;
}

/*
 * row_major_index::inverse(shape, shift=0)
 *   linear index from index-tuple constrained by shape
 */
 
template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::row_major_index::inverse(const shape_type& shape, value_type shift)const
{ 
  size_t n_linear(0);
  for(size_t n = 0; n < NDIM; ++n)
    n_linear += static_cast<size_t>(base_class::operator[](n) - shift) * product_(shape, n+1);
    
  return n_linear;  
}

/*
 * test wrap-around condition for index "x2" with respect to index "x1";
 *   initialize dimensions corresponding to wrapping with 0 or +-1, corresponding to wrap direction;
 *   return true if any dimensions wrap.
 */
template <class R, size_t NDIM>        
bool system<R,NDIM>::row_major_index::wrap(const row_major_index& x1, const row_major_index& x2, const shape_type& shape,
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

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::row_major_index::writeBinary(commUtil::abstractCommHandle* fp)const
{ return base_class::writeBinary(fp); }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::row_major_index::readBinary(commUtil::abstractCommHandle* fp)
{ return base_class::readBinary(fp); }

template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::row_major_index::binarySize(void)const
{ return base_class::binarySize(); }

/*
 * return the row_major_index corresponding to the position "p" within the hyper-cube "domain_cube"
 *   given the array "shape" and index-base "shift"
 */
template <class R, size_t NDIM>        
typename system<R,NDIM>::row_major_index system<R,NDIM>::row_major_index::ND_bin(
  const ntuple<R,NDIM>& p, const ntuple_interval<R,NDIM>& domain_cube,
  const shape_type& shape, value_type shift)
{
  row_major_index result(0,false);
  
  // "left_closure" \leftarrow  x \in [start, end)
  if (!(domain_cube.left_closure(p)))
    throw std::runtime_error("row_major_index::ND_bin: ntuple not in domain hypercube");
    
  for(size_t n = 0; n < NDIM; ++n){
    R L((domain_cube.end()[n] - domain_cube.start()[n])/integer<R>(shape[n])), 
      nx((p[n] - domain_cube.start()[n])/L);
    conv(result[n], floor(nx));
    result[n] += shift;  
  }
  
  return result;  
}  

#if 0  // completely screws-up the compiler for some reason (in unrelated code sections):
/*
 * return the sub-hypercube corresponding to the row_major_index within the hyper-cube "domain_cube"
 *   given the array "shape" and index-base "shift"
 */
template <class R, size_t NDIM>        
ntuple_interval<R,NDIM> system<R,NDIM>::row_major_index::ND_edges(
  const row_major_index& indices, const ntuple_interval<R,NDIM>& domain_cube,
  const shape_type& shape, value_type shift, const R& epsilon_)
{
  ntuple<R,NDIM> start_, end_;

  if (!indices.in_domain(shape, shift))
    throw std::runtime_error("row_major_index::ND_edges: indices out-of-range for specified shape and base shift");

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
template <class R, size_t NDIM>        
inline bool system<R,NDIM>::row_major_index::in_domain(const shape_type& shape, value_type shift)const
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

template <class R, size_t NDIM>        
inline system<R,NDIM>::row_major_index::row_major_index(value_type shift, bool clear)
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
template <class R, size_t NDIM>        
inline system<R,NDIM>::row_major_index::row_major_index(size_t n_linear, const shape_type& shape, value_type shift)
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


template <class R, size_t NDIM>        
inline system<R,NDIM>::row_major_index::row_major_index(const row_major_index& other)
  :base_class(other)
{ }

template <class R, size_t NDIM>        
template <class N1>
inline system<R,NDIM>::row_major_index::row_major_index(const ntuple<N1,NDIM>& other)
  :base_class(other)  
{ }

// ---------------------- end: moved to namespace linalg: (ND_matrix.h): -----------------------------
#endif



template <class R, size_t NDIM>        
inline const ntuple<R,NDIM>& system<R,NDIM>::state::position(void)const
{ return position_; }

template <class R, size_t NDIM>        
inline ntuple<R,NDIM>& system<R,NDIM>::state::position(void)
{ return position_; }

template <class R, size_t NDIM>        
inline const ntuple<R,NDIM>& system<R,NDIM>::state::velocity(void)const
{ return velocity_; }

template <class R, size_t NDIM>        
inline ntuple<R,NDIM>& system<R,NDIM>::state::velocity(void)
{ return velocity_; }

template <class R, size_t NDIM>        
inline const R& system<R,NDIM>::state::time(void)const
{ return time_; }

template <class R, size_t NDIM>        
inline R& system<R,NDIM>::state::time(void)
{ return time_; }

#if defined(state_partner_lists)
template <class R, size_t NDIM>        
const typename system<R,NDIM>::particle_list_type& system<R,NDIM>::state::partners(void)const
{ return partners_; }

template <class R, size_t NDIM>        
typename system<R,NDIM>::particle_list_type& system<R,NDIM>::state::partners(void)
{ return partners_; }

template <class R, size_t NDIM>        
void system<R,NDIM>::state::add_partner(typename system<R,NDIM>::particle *p) 
{ 
  // add particle only if it isn't already in the partners-list:
  typename particle_list_type::const_iterator itP(std::find(partners_.begin(), partners_.end(), p));
  if (itP == partners_.end()) 
    partners_.push_back(p);
}
#endif

template <class R, size_t NDIM>        
inline void system<R,NDIM>::state::init(const ntuple<R,NDIM>& x, const ntuple<R,NDIM>& v, const R& t)
{
  position_ = x;
  velocity_ = v;
  time_ = t;
  #if defined(state_partner_lists)
  partners_.clear();
  #endif
  event_type_ = 0;
}

        
template <class R, size_t NDIM>        
inline int  system<R,NDIM>::state::event_type(void)const
{ return event_type_; }

template <class R, size_t NDIM>        
inline void  system<R,NDIM>::state::set_event_type(int type)
{ event_type_ = type; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::state& system<R,NDIM>::state::operator=(const state& other)
{   
  // assert(this != &other); 
  position_ = other.position_;
  velocity_ = other.velocity_;
  time_ = other.time_;

  #if defined(state_partner_lists)
  partners_.clear();
  partners_ = other.partners_;
  #endif
  
  event_type_ = other.event_type_;
  
  return *this;
}

template <class R, size_t NDIM>        
inline void system<R,NDIM>::state::write(std::ostream& os)const
{
  size_t prec_save(os.precision(10));
  os<<"state: \n"
    <<"  position: ";
    position_.write(os);
    os<<"\n"
    <<"  velocity: ";
    velocity_.write(os);
    os<<"\n"
    <<"  time: "<<time_;
  #if defined(state_partner_lists)
  os<<"\n";
  write_partners(os);
  #endif
  os<<"\n";
  os<<"  event type: ";
  system<R,NDIM>::event::write(os, static_cast<typename system<R,NDIM>::event::event_kind>(event_type_));
  os<<"\n";
  os.precision(prec_save);
}

#if defined(state_partner_lists)
template <class R, size_t NDIM>        
inline void system<R,NDIM>::state::write_partners(std::ostream& os)const
{
  os<<"partners: ";
  if (!partners_.empty()){
    bool first(true);
    for(typename particle_list_type::const_iterator itP = partners_.begin(), itPEnd = partners_.end();
        itP != itPEnd;
        ++itP){
      if (first){
        os<<" ";
        first = false;
      }
      else
        os<<", ";
      os<<(*itP)->linear_index();
    }    
  }
  else
    os<<" none";
}
#endif

template <class R, size_t NDIM>        
inline  bool system<R,NDIM>::state::writeBinary(commUtil::abstractCommHandle* fp)const
{
  bool status(true);
  status = (status && position_.writeBinary(fp));
  status = (status && velocity_.writeBinary(fp));
  status = (status && system<R,NDIM>::writeBinary_(fp, time_)); 
  // partners_ is not read or written (pointer-based data); 
  //   it will be reconstructed next "step" cycles. 
  status = (status && system<R,NDIM>::writeBinary_(fp, event_type_));
  return status;
}

template <class R, size_t NDIM>        
inline  bool system<R,NDIM>::state::readBinary(commUtil::abstractCommHandle* fp)
{
  bool status(true);
  status = (status && position_.readBinary(fp));
  status = (status && velocity_.readBinary(fp));
  status = (status && system<R,NDIM>::readBinary_(fp, time_));
  #if defined(state_partner_lists)
  // partners_ is not read (it is pointer-based data):
  partners_.clear();  
  #endif
  status = (status && system<R,NDIM>::readBinary_(fp, event_type_));
  return status;
}

template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::state::binarySize(void)const
{
  size_t val(0);
  val += system<R,NDIM>::binarySize_(position_);
  val += system<R,NDIM>::binarySize_(velocity_);
  val += system<R,NDIM>::binarySize_(time_);  
  val += system<R,NDIM>::binarySize_(event_type_);  
  return val;
}

template <class R, size_t NDIM>        
inline system<R,NDIM>::state::state(void)
  : position_(), velocity_(), time_(0.0), event_type_(0)
{ }

template <class R, size_t NDIM>        
inline system<R,NDIM>::state::state(const state& other)
  : position_(other.position_), velocity_(other.velocity_), time_(other.time_), event_type_(0)
{ }
  
template <class R, size_t NDIM>        
inline system<R,NDIM>::state::state(const ntuple<R,NDIM>& x, const ntuple<R,NDIM>& v, const R& t)         
  : position_(x), velocity_(v), time_(t), event_type_(0)
{ }    
      

template <class R, size_t NDIM>        
typename system<R,NDIM>::particle::particle_kind system<R,NDIM>::particle::kind(void)const
{ return ABSTRACT_PARTICLE; }


/*
 * stick to another particle (modifies state of _both_ particles):
 *   sets particle velocities to zero (and sets "frozen_" flags)
 *   (where applicable (i.e. _not_ in jam case), particles have been moved to kinematic coordinates of event,
 *      prior to this method)
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::particle::stick(event* E)
{
  // note: stick-propagation goes in both directions:
  //  => _either_ of particle0_ or particle1_ may be "this".
  assert((E->particle1() != NULL)
         && ((E->particle0() == this) || (E->particle1() == this))); 

  particle &other(*((E->particle0() == this)? E->particle1(): E->particle0()));
  state
    &s1(current_state()),
    &s2(other.current_state());

  #if defined(_debug_print_)
  // *** DEBUG ***
  cout<<"STICK EVENT: "<<id_str()<<" <--> "<<other.id_str()<<endl;;
  #endif

  // transfer information to new states:
  rotate_state_buffer();
  other.rotate_state_buffer();

  // extrinsic state:
  current_state() = s1;
  other.current_state() = s2;
  current_state().velocity().clear();
  other.current_state().velocity().clear();
  #if defined(state_partner_lists)
  current_state().partners().clear();
  other.current_state().partners().clear();
  #endif
  current_state().set_event_type(event::STICK_EVENT);
  other.current_state().set_event_type(event::STICK_EVENT);  

  // intrinsic state:
  freeze();
  other.freeze();
}



template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::state& system<R,NDIM>::particle::current_state(void)const
{ 
  #if 0
  return state_[current_state_]; 
  #else
  // *** DEBUG ***
  return state_.back();
  #endif
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::state& system<R,NDIM>::particle::current_state(void)
{ 
  #if 0
  return state_[current_state_]; 
  #else
  // *** DEBUG ***
  return state_.back();
  #endif
}

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::state& system<R,NDIM>::particle::previous_state(void)const
{
  #if 0
  // note: using "%" operator for the following is not correct: 
  return state_[(current_state_==0? (state_buf_size_-1): current_state_-1)]; 
  #else
  // *** DEBUG ***
  assert(state_.size() > 1);
  return *(state_.end() - 2);
  #endif
}

template <class R, size_t NDIM>        
inline void system<R,NDIM>::particle::rotate_state_buffer(void)
{ 
  #if 0
  current_state_ = (current_state_ + 1) % state_buf_size_; 
  #else
  // *** DEBUG ***
  state_.push_back(state());
  if (state_.size() > state_buf_size_)
    state_.pop_front();
  #endif
}

// support for examination of entire state buffer:
template <class R, size_t NDIM>        
inline typename system<R,NDIM>::particle::state_buf_iterator 
  system<R,NDIM>::particle::state_buf_begin(void)const
{
  #if 0
  return state_;
  #else
  // *** DEBUG ***
  return state_.begin();
  #endif
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::particle::state_buf_iterator 
  system<R,NDIM>::particle::state_buf_end(void)const
{
  #if 0
  return state_ + state_buf_size_;
  #else
  // *** DEBUG ***
  return state_.end();
  #endif
}


template <class R, size_t NDIM>        
void system<R,NDIM>::particle::write_state_buf(std::ostream& os)const
{
  os<<"state history: \n";
  int nS(0);
  for(state_buf_iterator itS = state_buf_end() - 1, itSEnd = state_buf_begin() - 1;
      itS != itSEnd;
      --itS, --nS){
    if (0 != nS)
      os<<"\n";
    os<<"------ "<<nS<<": ------";
    (*itS).write(os);      
  }
  os<<"\n-----------------------------";    
}



template <class R, size_t NDIM>        
inline long system<R,NDIM>::particle::linear_index(void)const
{ return linear_index_; }

template <class R, size_t NDIM>        
inline long& system<R,NDIM>::particle::linear_index(void)
{ return linear_index_; }


template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell* system<R,NDIM>::particle::owner(void)const
{ return owner_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell*& system<R,NDIM>::particle::owner(void)
{ return owner_; }


template <class R, size_t NDIM>        
inline bool system<R,NDIM>::particle::frozen(void)const
{ return frozen_; }

template <class R, size_t NDIM>        
inline void system<R,NDIM>::particle::freeze(bool flag)
{ frozen_ = flag; }

// id-string as "p<n>,c<n>":
template <class R, size_t NDIM>        
std::string system<R,NDIM>::particle::id_str(const particle* p, const cell* c)
{  
  std::ostringstream oss;
  oss<<"p<";
  if (p != NULL)
    oss<<p->linear_index();
  else
    oss<<"NULL";
  oss<<">,c<";
  if (c != NULL)
    oss<<c->linear_index();
  else
    oss<<"NULL";
  oss<<">";
  if (p->frozen())
    oss<<"_frozen";
  return oss.str();
}

template <class R, size_t NDIM>        
inline std::string system<R,NDIM>::particle::id_str(void)const
{ return id_str(this, owner()); }


// update any _intrinsic_ attributes associated with system time change.
// (does _not_ modify particle "state" (here considered extrinsic))
template <class R, size_t NDIM>        
void system<R,NDIM>::particle::move(const R& t)
{ }

// overlap condition test:
template <class R, size_t NDIM>        
bool system<R,NDIM>::particle::overlap(const particle* pother)const
{ return false; }

/**
 * @brief Velocity magnitude corresponding to a specified kinetic temperature.
 *   This virtual method allows any method using temperature to initialize state velocity to be completely virtual.
 *   - for spheres, the zero-radius case uses density and "dr" (instead of radius) to calculate an effective velocity;
 *   - for abstract particles, mass is taken as unity.
 *   .
 */
template <class R, size_t NDIM>        
R system<R,NDIM>::particle::kinetic_velocity(const R& T)const
{ return sqrt(integer<R>(2) * T); }  

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::particle& system<R,NDIM>::particle::operator=(const particle& other)
{
  #if 0
  std::copy(&(other.state_[0]), &(other.state_[0])+state_buf_size_, &(state_[0]));
  current_state_ = other.current_state_;
  #else
  // *** DEBUG ***
  state_ = other.state_;
  #endif
  linear_index_ = other.linear_index_;
  owner_ = other.owner_;
  frozen_ = other.frozen_;
  return *this;
}

template <class R, size_t NDIM>        
typename system<R,NDIM>::particle* system<R,NDIM>::particle::clone(void)const
{ return new particle(*this); }

// methods to allow binary read and write from pointer to base-class:            
template <class R, size_t NDIM>        
bool system<R,NDIM>::particle::writeBinaryVirtual(commUtil::abstractCommHandle* fp, const particle* p)
{
  bool status(true);
  status = (status && system<R,NDIM>::writeBinary_(fp, static_cast<long>(p->kind())));
  status = (status && p->writeBinary(fp));
  return status;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::particle::readBinaryVirtual(commUtil::abstractCommHandle* fp, particle*& p)
{
  bool status(true);
  long nKind(0);
  status = (status && system<R,NDIM>::readBinary_(fp, nKind));
  if (status)
    switch (static_cast<particle_kind>(nKind)){
      case ABSTRACT_PARTICLE:
      p = new typename system<R,NDIM>::particle();
      break;
      case SPHERE_PARTICLE:
      p = new typename mjp_system<R,NDIM>::sphere();
      break;
      default:
      status = false;
      break;
    }
  status = (status && p->readBinary(fp));
  return status;  
}

template <class R, size_t NDIM>        
size_t system<R,NDIM>::particle::binarySizeVirtual(const particle* p)
{
  size_t val(0);
  val += system<R,NDIM>::binarySize_(static_cast<long>(p->kind()));
  val += p->binarySize();
  return val;  
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::particle::writeBinary(commUtil::abstractCommHandle* fp)const
{
  bool status(true);
  #if 0
  for(size_t s = 0; status && (s < state_buf_size_); ++s)
    status = (status && state_[s].writeBinary(fp));
  status = (status && writeBinary_(fp,current_state_));
  #else
  // *** DEBUG ***
  // define this explicitly, to continue independent implementation of read/writeBinary_ methods:
  size_t size_(state_.size());
  status = (status && writeBinary_(fp,size_));
  for(typename std::deque<state>::const_iterator itS = state_.begin(), itSEnd = state_.end();
      status && (itS != itSEnd);
      ++itS)
    status = (status && (*itS).writeBinary(fp));  
  #endif
  status = (status && system<R,NDIM>::writeBinary_(fp, linear_index_));  
  // note: "owner_" is not written
  status = (status && system<R,NDIM>::writeBinary_(fp, frozen_));    
  return status;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::particle::readBinary(commUtil::abstractCommHandle* fp)
{
  bool status(true);
  #if 0
  for(size_t s = 0; status && (s < state_buf_size_); ++s)
    status = (status && state_[s].readBinary(fp));
  status = (status && readBinary_(fp,current_state_));
  #else
  // *** DEBUG ***
  // define this explicitly, to continue independent implementation of read/writeBinary_ methods:
  size_t size_(0);
  status = (status && readBinary_(fp,size_));
  if (status){
    state_.clear();
    state_.resize(size_);
    for(typename std::deque<state>::iterator itS = state_.begin(), itSEnd = state_.end();
        status && (itS != itSEnd);
        ++itS)
      status = (status && (*itS).readBinary(fp));  
  }     
  #endif
  status = (status && system<R,NDIM>::readBinary_(fp, linear_index_));  
  owner_ = NULL;
  status = (status && system<R,NDIM>::readBinary_(fp, frozen_));      
  return status;
}

template <class R, size_t NDIM>        
size_t system<R,NDIM>::particle::binarySize(void)const
{
  bool val(0);
  #if 0
  for(size_t s = 0; s < state_buf_size_; ++s)
    val += state_[s].binarySize();
  #else
  val += sizeof(size_t);
  for(typename std::deque<state>::const_iterator itS = state_.begin(), itSEnd = state_.end();
      itS != itSEnd;
      ++itS)
    val += (*itS).binarySize();  
  #endif
  val += sizeof(long); // linear_index_
  // note: owner_ is not written
  val += system<R,NDIM>::binarySize_(frozen_);
  return val;  
}


template <class R, size_t NDIM>        
system<R,NDIM>::particle::~particle(void)
{ }

template <class R, size_t NDIM>        
inline system<R,NDIM>::particle::particle(void)
  : 
  #if 0
  current_state_(0), 
  #else
  state_(state_buf_size_), 
  #endif
  linear_index_(LONG_MAX), owner_(NULL), frozen_(false)
{ }

template <class R, size_t NDIM>        
inline system<R,NDIM>::particle::particle(const particle& other)
{ operator=(other); }

template <class R, size_t NDIM>        
inline system<R,NDIM>::particle::particle(const state& s, long linear_index, cell* owner)
  : 
  #if 0
  current_state_(0), 
  #else
  state_(state_buf_size_),
  #endif
  linear_index_(linear_index), owner_(owner), frozen_(false)
{
  current_state() = s;
}

// usage note: a cell _owns_ its particles, 
//   transfer of particle pointer transfers ownership:
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::particle_list_type& system<R,NDIM>::cell::particles(void)const
{ return particles_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::particle_list_type& system<R,NDIM>::cell::particles(void)
{ return particles_; }

// a cell does _not_ own its neighbors:
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell::neighbor_list_type& system<R,NDIM>::cell::neighbors(void)const
{ return neighbors_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell::neighbor_list_type& system<R,NDIM>::cell::neighbors(void)
{ return neighbors_; }

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell::neighbor_list_type& system<R,NDIM>::cell::positive_neighbors(void)const
{ return positive_neighbors_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell::neighbor_list_type& system<R,NDIM>::cell::positive_neighbors(void)
{ return positive_neighbors_; }

template <class R, size_t NDIM>        
inline const std::vector<int>& system<R,NDIM>::cell::wrap(void)const
{ return wrap_; }

template <class R, size_t NDIM>        
inline std::vector<int>& system<R,NDIM>::cell::wrap(void)
{ return wrap_; }

template <class R, size_t NDIM>        
inline const std::vector<typename system<R,NDIM>::row_major_index>& system<R,NDIM>::cell::wrap_dim(void)const
{ return wrap_dim_; }

template <class R, size_t NDIM>        
inline std::vector<typename system<R,NDIM>::row_major_index>& system<R,NDIM>::cell::wrap_dim(void)
{ return wrap_dim_; }


template <class R, size_t NDIM>        
inline const size_t system<R,NDIM>::cell::linear_index(void)const
{ return linear_index_; }

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::row_major_index& system<R,NDIM>::cell::indices(void)const
{ return indices_;}

        
// neighbor from relative index:
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell* system<R,NDIM>::cell::neighbor(const row_major_index& nx)const
{ return neighbors_[nx.inverse(system<R,NDIM>::cell::offset_shape(), -1)]; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell* system<R,NDIM>::cell::neighbor(const row_major_index& nx)
{ return neighbors_[nx.inverse(system<R,NDIM>::cell::offset_shape(), -1)]; }

// face neighbor from face index pair:
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell* system<R,NDIM>::cell::face_neighbor(size_t dimension, size_t face)const
{
  row_major_index nx_offset(0,false);
  nx_offset.clear(); // init to zero 
  nx_offset[dimension] = (face==0? -1: 1);
  return neighbor(nx_offset);
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell* system<R,NDIM>::cell::face_neighbor(size_t dimension, size_t face)
{
  row_major_index nx_offset(0,false);
  nx_offset.clear(); // init to zero 
  nx_offset[dimension] = (face==0? -1: 1);
  return neighbor(nx_offset);
}



template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::cell::N_face(void)
{ return 2*NDIM;}

template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::cell::N_neighbor(void)
{ return _STL_EXT_NAMESPACE_::power(3, NDIM) - 1;}

template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::cell::N_positive_neighbor(void)
{ return _STL_EXT_NAMESPACE_::power(2, NDIM) - 1;}

// static constant (and associated method):
template <class R, size_t NDIM>        
typename system<R,NDIM>::shape_type system<R,NDIM>::cell::offset_shape_(shape_type::filled(3));

template <class R, size_t NDIM>        
const typename system<R,NDIM>::shape_type& system<R,NDIM>::cell::offset_shape(void)
{ return system<R,NDIM>::cell::offset_shape_; }


// transfer a particle to another cell:
template <class R, size_t NDIM>        
void system<R,NDIM>::cell::transfer_particle(particle* p, cell* pother)
{
  typename particle_list_type::iterator itP = std::find(particles_.begin(), particles_.end(), p);
  if (itP == particles_.end())
    throw std::runtime_error("cell::transfer_particle: particle not owned by cell");
  pother->particles_.push_back(*itP);
  (*itP)->owner() = pother;
  particles_.erase(itP);  
}

// assignment or clone does _not_ transfer neighbor-lists (and associated boundary wrapping info).
// (these can only be constructed by parent "cell_array")template <class R, size_t NDIM>        
template <class R, size_t NDIM>        
typename system<R,NDIM>::cell* system<R,NDIM>::cell::clone(void)const
{ return new cell(*this); }

template <class R, size_t NDIM>        
typename system<R,NDIM>::cell& system<R,NDIM>::cell::operator=(const cell& other)
{
  free_pointers_();
  particles_.resize(other.particles_.size(), NULL);
  typename particle_list_type::const_iterator itOtherP = other.particles_.begin();
  for(typename particle_list_type::iterator itP = particles_.begin(), itPEnd = particles_.end();
      itP != itPEnd;
      ++itP, ++itOtherP)
    *itP = (*itOtherP)->clone();
  linear_index_ = other.linear_index_;
  indices_ = other.indices_;
  return *this;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::cell::writeBinary(commUtil::abstractCommHandle* fp)const
{
 // Implementation note: neighbors_, positive_neighbors_, wrap_, wrap_dim_ are _not_ written.
 //   These will be reconstructed by the parent "cell_array" (if applicable).
 bool status(true);
 status = (status && system<R,NDIM>::writeBinary_(fp, particles_.size()));
 for(typename particle_list_type::const_iterator itP = particles_.begin(), itPEnd = particles_.end();
     status && (itP != itPEnd);
     ++itP)
   status = (status && system<R,NDIM>::particle::writeBinaryVirtual(fp, *itP)); 
 status = (status && system<R,NDIM>::writeBinary_(fp, linear_index_));
 status = (status && indices_.writeBinary(fp));
 return status;    
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::cell::readBinary(commUtil::abstractCommHandle* fp)
{
 // (see comment at "cell::writeBinary".)
 bool status(true);
 size_t N(0);
 status = (status && system<R,NDIM>::readBinary_(fp, N));
 if (status){
   free_pointers_();
   particles_.resize(N,NULL);
 }
 for(typename particle_list_type::iterator itP = particles_.begin(), itPEnd = particles_.end();
     status && (itP != itPEnd);
     ++itP){
   status = (status && system<R,NDIM>::particle::readBinaryVirtual(fp, *itP));
   if (status)
     (*itP)->owner() = this;
 }
 status = (status && system<R,NDIM>::readBinary_(fp, linear_index_));
 status = (status && indices_.readBinary(fp));
 return status;    
}

template <class R, size_t NDIM>        
size_t system<R,NDIM>::cell::binarySize(void)const
{
  size_t val(0);
  val += sizeof(size_t);
  for(typename particle_list_type::const_iterator itP = particles_.begin(), itPEnd = particles_.end();
      itP != itPEnd;
      ++itP)
    val += system<R,NDIM>::particle::binarySizeVirtual(*itP);
  val += system<R,NDIM>::binarySize_(linear_index_);
  val += indices_.binarySize();
  return val;  
}

template <class R, size_t NDIM>        
void system<R,NDIM>::cell::free_pointers_(void)
{
  // delete owned particles
  for(typename particle_list_type::iterator itP = particles_.begin(), itPEnd = particles_.end();
      itP != itPEnd;
      ++itP)
    if (*itP != NULL)
      delete *itP;    
  particles_.clear();
}

template <class R, size_t NDIM>        
system<R,NDIM>::cell::~cell(void)
{
  // delete owned particles
  free_pointers_();     
}

template <class R, size_t NDIM>        
inline system<R,NDIM>::cell::cell(void)
{ }  

// copy does _not_ initialize neighbor-lists (and associated boundary wrapping info).
// (these can only be initialized by parent "cell_array")
template <class R, size_t NDIM>        
system<R,NDIM>::cell::cell(const cell& other)
{ 
  operator=(other);
}

template <class R, size_t NDIM>        
inline system<R,NDIM>::cell::cell(const size_t linear_index, const row_major_index& indices)
  : linear_index_(linear_index), indices_(indices)
{
  // other attributes are initialized by parent "cell_array" 
}                 

    
template <class R, size_t NDIM>        
inline const size_t& system<R,NDIM>::cell_array::N1(void)const
{ return N1_; }

template <class R, size_t NDIM>        
inline const size_t& system<R,NDIM>::cell_array::N_cell(void)const
{ return data_.size(); }

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::shape_type& system<R,NDIM>::cell_array::shape(void)const
{ return shape_; }

// usage note: a cell_array _owns_ its cells,
//   transfer of cell pointer transfers ownership
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell_array::cell_list_type& system<R,NDIM>::cell_array::data(void)const
{ return data_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell_array::cell_list_type& system<R,NDIM>::cell_array::data(void)
{ return data_; }


// allow usage as container to be (somewhat) transparent:
template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell_array::iterator system<R,NDIM>::cell_array::begin(void)
{ return data_.begin(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell_array::const_iterator system<R,NDIM>::cell_array::begin(void)const
{ return data_.begin(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell_array::iterator system<R,NDIM>::cell_array::end(void)
{ return data_.end(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell_array::const_iterator system<R,NDIM>::cell_array::end(void)const
{ return data_.end(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell* system<R,NDIM>::cell_array::operator[](size_t linear_index)
{ 
  if (linear_index >= data_.size())
    throw std::runtime_error("cell_array::operator[]: linear index out-of-range");
  return data_[linear_index]; 
}

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell* system<R,NDIM>::cell_array::operator[](size_t linear_index)const
{ 
  if (linear_index >= data_.size())
    throw std::runtime_error("cell_array::operator[]: linear index out-of-range");
  return data_[linear_index]; 
}

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::cell* system<R,NDIM>::cell_array::operator[](const row_major_index& indices)
{ return operator[](indices.inverse(shape())); }

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::cell* system<R,NDIM>::cell_array::operator[](const row_major_index& indices)const
{ return operator[](indices.inverse(shape())); }




template <class R, size_t NDIM>        
bool system<R,NDIM>::cell_array::writeBinary(commUtil::abstractCommHandle* fp)const
{
  bool status(true);
  status = (status && system<R,NDIM>::writeBinary_(fp, N1_));
  status = (status && shape_.writeBinary(fp));
  size_t size_(data_.size());
  status = (status && system<R,NDIM>::writeBinary_(fp,size_));
  for(const_iterator itD = data_.begin(), itDEnd = data_.end();
      status && (itD != itDEnd);
      ++itD)
    status = (status && (*itD)->writeBinary(fp));  
  return status;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::cell_array::readBinary(commUtil::abstractCommHandle* fp)
{
  bool status(true);
  status = (status && system<R,NDIM>::readBinary_(fp, N1_));
  status = (status && shape_.readBinary(fp));
  size_t size_(0);
  status = (status && system<R,NDIM>::readBinary_(fp, size_));
  if (status){
    free_pointers_();
    data_.resize(size_, NULL);
  }
  for(iterator itD = data_.begin(), itDEnd = data_.end();
      status && (itD != itDEnd);
      ++itD){
    *itD = new cell();
    status = (status && (*itD)->readBinary(fp));
  }
  if (status)
    init_neighbors_();    
  return status;
}

template <class R, size_t NDIM>        
size_t system<R,NDIM>::cell_array::binarySize(void)const
{
  size_t val(0);
  val += system<R,NDIM>::binarySize_(N1_);
  val += shape_.binarySize();
  val += sizeof(size_t); // data_ size field
  for(const_iterator itD = data_.begin(), itDEnd = data_.end();
      itD != itDEnd;
      ++itD)
    val += (*itD)->binarySize();  
  return val;  
}

template <class R, size_t NDIM>        
void system<R,NDIM>::cell_array::free_pointers_(void)
{
  for(iterator itD = data_.begin(), itDEnd = data_.end();
      itD != itDEnd;
      ++itD)
    if (*itD != NULL)
      delete *itD;
  data_.clear();
} 

template <class R, size_t NDIM>        
typename system<R,NDIM>::cell_array& system<R,NDIM>::cell_array::operator=(const cell_array& other)
{
  free_pointers_();
  N1_ = other.N1_;
  shape_ = other.shape_;
  
  data_.resize(other.data_.size());
  const_iterator itOtherD = other.data_.begin();
  for(iterator itD = data_.begin(), itDEnd = data_.end();
      itD != itDEnd;
      ++itD, ++itOtherD)
    *itD = (*itOtherD)->clone();
  
  // rebuild cell-neighbor lists (and associated boundary wrapping information):
  init_neighbors_();
    
  return *this;
} 

template <class R, size_t NDIM>        
system<R,NDIM>::cell_array::~cell_array(void)
{
  free_pointers_();
} 

/*
 * initialize cell neighbor-lists and boundary wrapping information:
 *   definitions:
 *     neighbor_: _all_ adjacent cells (for convenience, cell itself is at central index location)
 *     positive_neighbor_: group of adjacent cells "owned" by cell 
 *       (i.e. which don't overlap with analogous groups from adjacent cells))
 */ 
template <class R, size_t NDIM>        
void system<R,NDIM>::cell_array::init_neighbors_(void)
{
  for(typename std::vector<typename system<R,NDIM>::cell*>::iterator itD = data_.begin(), itDEnd = data_.end();
      itD != itDEnd;
      ++itD){
      
    cell *cell_(*itD);
    cell_->neighbors().clear();
    cell_->positive_neighbors().clear();
    cell_->wrap().clear();
    cell_->wrap_dim().clear();
    
    const index &nx(cell_->indices());
    
    // offset from cell at central position:
    //   shift will be -1 to produce offset sub-indices: -1, 0, 1
    shape_type offset_shape(shape_type::filled(3));     
    for(size_t nn = 0; nn < cell::N_neighbor()+1; ++nn){
      // note: one more than "formal" N_neighbor:
      //  initialize neighbor-index (0,0,...) to cell itself
      index nx_offset(nn, offset_shape, -1),
        nx_neighbor(nx_offset);
      nx_neighbor += nx;
      
      // application of modulus => periodic boundary conditions for the global hypercube:
      nx_neighbor.mod_assign(shape_);

      // add pointer to neighboring cell to the neighbor list:
      cell *neighbor(data_[nx_neighbor.inverse(shape_)]);
      cell_->neighbors().push_back(neighbor); 
      
      // add neighbors in the positive wrapped half-space corresponding to the last (fastest changing) dimension 
      //   to the positive_neighbors_ list:
      if (nx_offset[NDIM-1] >= 0){
        cell_->positive_neighbors().push_back(neighbor);
        index wrap_dim;
        // record details about wrapped neighbors:
        bool wrap = index::wrap(nx, nx_neighbor, shape_, wrap_dim);
        cell_->wrap().push_back(static_cast<int>(wrap));
        cell_->wrap_dim().push_back(wrap_dim);        
      }
    } 
  } 
}

/**
 * @brief Initialize the cell array using a "suggested" value for the number-of-cells.
 *   the actual number of cells will be the value corresponding to the closest (>=) integer
 *   that is <some other integer>^NDIM.
 * For correct functioning of periodic-boundary-condition wrap-around, the minimum number of
 *   cells in the linear dimension must be >= 3; this insures that there is only one way to
 *   reach a cell from a neighboring cell.
 */ 
template <class R, size_t NDIM>        
void system<R,NDIM>::cell_array::init(size_t N_cell)
{
  using std::max;
  // allow re-initialization:
  free_pointers_();

  // adjust number of cells so that
  //   hypercube side length is an integer >= 3:
  N1_ = max(static_cast<size_t>(3), static_cast<size_t>(ceil(pow(double(N_cell),(1.0/double(NDIM))))));
  size_t size_ = _STL_EXT_NAMESPACE_::power(N1_, NDIM);
  std::fill(shape_.begin(), shape_.end(), N1_);
  
  // initialize cell array:
  data_.resize(size_, NULL);
  size_t cell_number(0);
  for(typename std::vector<cell*>::iterator itD = data_.begin(), itDEnd = data_.end();
      itD != itDEnd;
      ++itD, ++cell_number)
    *itD = new cell(cell_number, row_major_index(cell_number, shape_));  
  
  // initialize cell neighbor-lists and boundary wrapping information:
  init_neighbors_();
}

template <class R, size_t NDIM>        
system<R,NDIM>::cell_array::cell_array(void)
 : N1_(0)
{ }

template <class R, size_t NDIM>        
system<R,NDIM>::cell_array::cell_array(const cell_array& other)
  : N1_(0)
{ 
  operator=(other);
}


/*
 * initialize the cell array using a "suggested" value for the number-of-cells.
 *   the actual number of cells will be the value corresponding to the closest (>=) integer
 *   that is <some other integer>^NDIM
 */      
template <class R, size_t NDIM>        
system<R,NDIM>::cell_array::cell_array(size_t N_cell)    
  : N1_(0) 
{
  init(N_cell);
}


template <class R, size_t NDIM>        
void system<R,NDIM>::event::write(std::ostream& os, event_kind kind_)
{
  switch(kind_){
    case NULL_EVENT:
    os<<"NULL";
    break;
    
    case TIMESTEP_EVENT:
    os<<"TIMESTEP";
    break;
    
    case MOVE_EVENT:
    os<<"MOVE";
    break;
    
    case CELL_EXIT_EVENT:
    os<<"CELL EXIT";
    break;
    
    case COLLISION_EVENT:
    os<<"COLLISION";
    break;
    
    case JAM_EVENT:
    os<<"JAM";
    break;
    
    case STICK_EVENT:
    os<<"STICK";
    break;
    
    default:
    throw std::runtime_error("event::write: unknown event kind");
    break;
  }
}


template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event::event_kind system<R,NDIM>::event::kind(void)const
{ return kind_; }

template <class R, size_t NDIM>        
inline size_t system<R,NDIM>::event::id(void)const
{ return id_; }

template <class R, size_t NDIM>        
const R& system<R,NDIM>::event::time(void)const
{ return time_; }

// usage note: event does _not_ own its cell or particle[s]
//   for consistency with other classes these methods return pointers.
template <class R, size_t NDIM>        
const typename system<R,NDIM>::cell* system<R,NDIM>::event::event_cell(void)const
{ return cell_; }

template <class R, size_t NDIM>        
typename system<R,NDIM>::cell* system<R,NDIM>::event::event_cell(void)
{ return cell_; }

template <class R, size_t NDIM>        
const typename system<R,NDIM>::particle* system<R,NDIM>::event::particle0(void)const
{ return particle0_; }

template <class R, size_t NDIM>        
typename system<R,NDIM>::particle* system<R,NDIM>::event::particle0(void)
{ return particle0_; }

template <class R, size_t NDIM>        
const typename system<R,NDIM>::particle* system<R,NDIM>::event::particle1(void)const
{ return particle1_; }

template <class R, size_t NDIM>        
typename system<R,NDIM>::particle* system<R,NDIM>::event::particle1(void)
{ return particle1_; }

template <class R, size_t NDIM>        
const size_t& system<R,NDIM>::event::dimension(void)const
{ return dimension_; }

template <class R, size_t NDIM>        
const size_t& system<R,NDIM>::event::face(void)const
{ return face_; }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event::is_cell_exit(void)const
{ return kind_ == CELL_EXIT_EVENT; }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event::is_collision(void)const
{ return kind_ == COLLISION_EVENT; }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event::is_jam(void)const
{ return kind_ == JAM_EVENT; }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event::is_stick(void)const
{ return kind_ == STICK_EVENT; }

/*
 * duplicate events have same event_kind and same particles (or NULL particle*):
 *  (particle ordering is ignored)
 */ 
template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event::duplicate(const event* other)const
{
  bool duplicate(false);
  if (kind_ == other->kind_){
    const particle *p0(particle0_), *p1(particle1_);
    // ignore order when comparing particles:
    if ((p0 != NULL) && (p1 != NULL) && (p0->linear_index() > p1->linear_index()))
      std::swap(p0, p1);

    const particle *p0_(other->particle0_), *p1_(other->particle1_);
    if ((p0_ != NULL) && (p1_ != NULL) && (p0_->linear_index() > p1_->linear_index()))
      std::swap(p0_, p1_);

    if ((((p0_ == NULL) && (p0 == NULL))
         || ((p0_ != NULL) && (p0 != NULL) && (p0_->linear_index() == p0->linear_index())))
        &&
        (((p1_ == NULL) && (p1 == NULL))
         || ((p1_ != NULL) && (p1 != NULL) && (p1_->linear_index() == p1->linear_index()))))
      duplicate = true;
  }

  return duplicate;
}



template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event* system<R,NDIM>::event::clone(void)const
{ return new event(*this); }

template <class R, size_t NDIM>        
typename system<R,NDIM>::event& system<R,NDIM>::event::operator=(const event& other)
{
  id_ = other.id_;  // note "id_" counter is associated with initialization, not copy
  
  kind_ = other.kind_;
  time_ = other.time_;
  cell_ = other.cell_;
  particle0_ = other.particle0_;
  particle1_ = other.particle1_;
  dimension_ = other.dimension_;
  face_ = other.face_;
  jam_ = other.jam_;
  return *this;
}


template <class R, size_t NDIM>        
void system<R,NDIM>::event::write(std::ostream& os)
{
  os<<"event("<<id_<<"): ";
  write(os,kind_);
  os<<"\n"
    <<"  time: "<<time_<<"\n"
    <<"  cell number: "<<cell_->linear_index()<<"\n";
    
  os<<"  particle0: ";    
  if (particle0_ != NULL){  
    os<<particle0_->linear_index();
    #if defined(state_partner_lists)
    if (!particle0_->current_state().partners().empty()){
      os<<" (";
      particle0_->current_state().write_partners(os);
      os<<")";
    }
    #endif
    os<<"\n";
  }  
  else
    os<<"NULL\n";
    
  os<<"  particle1: ";
  if (particle1_ != NULL){  
    os<<particle1_->linear_index();
    #if defined(state_partner_lists)
    if (!particle1_->current_state().partners().empty()){
      os<<"(";
      particle1_->current_state().write_partners(os);
      os<<")";
    }
    #endif
    os<<"\n";
  }  
  else
    os<<"NULL\n";
    
  os<<"  dimension: "<<dimension_<<"\n"
    <<"  face: "<<face_<<"\n";
}

/*
 * static counter for class to generate unique ids:
 * (note: only "init" increments this counter ("operator=" and "clone" do _not_):
 */
template <class R, size_t NDIM>        
size_t system<R,NDIM>::event::event_counter_ = 0;
           
template <class R, size_t NDIM>        
void system<R,NDIM>::event::init(event_kind kind,
                                 const R& time, cell* c, particle *particle0, 
                                 particle *particle1,
                                 size_t dimension, size_t face)
{
  if (kind != NULL_EVENT){
    ++event_counter_; // static  
    id_ = event_counter_;
  }
  else
    id_ = static_cast<size_t>(-1); // NULL_EVENT marker
    
  kind_ = kind;
  time_ = time;
  cell_ = c;
  particle0_ = particle0;
  particle1_ = particle1;
  dimension_ = dimension;
  face_ = face;
}
 
template <class R, size_t NDIM>        
system<R,NDIM>::event::~event(void)
{ }

template <class R, size_t NDIM>        
system<R,NDIM>::event::event(void)
 : id_(static_cast<size_t>(-1)),
   kind_(NULL_EVENT), 
   time_(0.0), cell_(NULL), particle0_(NULL), particle1_(NULL), 
   dimension_(ULONG_MAX), face_(ULONG_MAX)
{ }
 
template <class R, size_t NDIM>        
system<R,NDIM>::event::event(const event& other)
{ operator=(other); }

template <class R, size_t NDIM>        
system<R,NDIM>::event::event(
      event_kind kind,
      const R& time, cell* c, particle *particle0, 
      particle *particle1,
      size_t dimension, size_t face)
{
  init(kind, time, c, particle0, particle1, dimension, face);
}
                  


// free owned data:
template <class R, size_t NDIM>        
void system<R,NDIM>::event_list::free_pointers_(void)
{
  for(typename list_type::iterator itE = data_.begin(), itEEnd = data_.end();
      itE != itEEnd;
      ++itE)
    if (*itE != NULL)
      delete *itE;
  data_.clear();
}      

// allow (somewhat) transparent usage as container-type:
template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event_list::const_iterator system<R,NDIM>::event_list::begin(void)const
{ return data_.begin(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event_list::iterator system<R,NDIM>::event_list::begin(void)
{ return data_.begin(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event_list::const_iterator system<R,NDIM>::event_list::end(void)const
{ return data_.end(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event_list::iterator system<R,NDIM>::event_list::end(void)
{ return data_.end(); }  

template <class R, size_t NDIM>        
inline void system<R,NDIM>::event_list::push_back(event *E)
{ data_.push_back(E); }

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::event*& system<R,NDIM>::event_list::back(void)const
{ return data_.back(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event*& system<R,NDIM>::event_list::back(void)
{ return data_.back(); }

// push an event pointer onto front of list (transfers ownership of pointer-object): 
template <class R, size_t NDIM>        
inline void system<R,NDIM>::event_list::push_front(event *E)
{ data_.push_front(E); }

// push a non-duplicate event object onto front of list (pointer-object is cloned if transfer occurs): 
template <class R, size_t NDIM>        
inline void system<R,NDIM>::event_list::push_front_unique(const event *E)
{    
  bool duplicate(false);
  for(const_iterator itE = data_.begin(), itEEnd = data_.end();
      itE != itEEnd;
      ++itE)
    if (E->duplicate(*itE)){
        duplicate = true;
        break;
    }  
  
  if (!duplicate)   
    data_.push_front(E->clone()); 
}

template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::event*& system<R,NDIM>::event_list::front(void)const
{ return data_.front(); }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event*& system<R,NDIM>::event_list::front(void)
{ return data_.front(); }

template <class R, size_t NDIM>        
inline bool system<R,NDIM>::event_list::empty(void)const
{ return data_.empty(); }

// usage note: event_list _owns_ its events,
//   transfer of event pointer transfers ownership:
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::event_list::list_type& system<R,NDIM>::event_list::data(void)const
{ return data_; }

template <class R, size_t NDIM>        
inline typename system<R,NDIM>::event_list::list_type system<R,NDIM>::event_list::data(void)
{ return data_; }

template <class R, size_t NDIM>        
typename system<R,NDIM>::event_list& system<R,NDIM>::event_list::operator=(const event_list& other)
{
  free_pointers_();
  data_.resize(other.data_.size(), NULL);
  
  typename list_type::const_iterator itOtherE = other.data_.begin();
  for(typename list_type::iterator itE = data_.begin(), itEEnd = data_.end();
      itE != itEEnd;
      ++itE, ++itOtherE)
    *itE = (*itOtherE)->clone();

  return *this;
}


#if 0 // ------------ event is fundamentally a pointer-based structure               ----------------
// ------------   => implementing these methods doesn't really make any sense: ---------------- 

template <class R, size_t NDIM>        
bool system<R,NDIM>::event_list::writeBinary(commUtil::abstractCommHandle* fp)const

template <class R, size_t NDIM>        
bool system<R,NDIM>::event_list::readBinary(commUtil::abstractCommHandle* fp)

template <class R, size_t NDIM>        
size_t system<R,NDIM>::event_list::binarySize(void)const
#endif // -------------------------------------------------------------------------------------------

template <class R, size_t NDIM>        
void system<R,NDIM>::event_list::write(std::ostream& os)const
{
  os<<"----- event_list: -----\n";
  for(const_iterator itE = data_.begin(), itEEnd = data_.end();
      itE != itEEnd;
      ++itE)
    (*itE)->write(os);
  os<<"\n";  
}

// remove list-position corresponding to next finite-time event (and free its object):
template <class R, size_t NDIM>        
inline void system<R,NDIM>::event_list::remove_next_event(void)
{
  assert(!data_.empty() && !data_.back()->is_jam());
  delete data_.back(); 
  data_.pop_back();  
}

template <class R, size_t NDIM>        
inline void system<R,NDIM>::event_list::clear(void)
{
  free_pointers_();
}

template <class R, size_t NDIM>        
system<R,NDIM>::event_list::~event_list(void)
{
  free_pointers_();
}

template <class R, size_t NDIM>        
system<R,NDIM>::event_list::event_list(void)
{ }

template <class R, size_t NDIM>        
system<R,NDIM>::event_list::event_list(const event_list& other)    
{
  operator=(other);
}

        
template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::parameters& system<R,NDIM>::get_parameters(void)const
{
  assert(pparam_ != NULL);
  return *pparam_;
}

template <class R, size_t NDIM>        
void system<R,NDIM>::set_parameters(const parameters& param)
{
  if (pparam_ != NULL)
    delete pparam_;
  pparam_ = reinterpret_cast<parameters*>(param.clone());
  pparam_->valid_check();
}        


#if defined(event_history)
// *** DEBUG ***

// -------------- system<R,NDIM>: static variable: --------------------
template <class R, size_t NDIM>
typename system<R,NDIM>::event_list system<R,NDIM>::event_history_;
  
template <class R, size_t NDIM>        
inline void system<R,NDIM>::write_event_history(std::ostream& os)
{ 
  size_t prec_save(os.precision(18));
  event_history_.write(os); 
  os.precision(prec_save);
}

#endif



template <class R, size_t NDIM>        
inline const typename system<R,NDIM>::parameters* system<R,NDIM>::pparam(void)const
{ return pparam_; }

// protected initialization method 
//  (initialization relevent to local class only, does not call virtual methods or super-class methods):           
template <class R, size_t NDIM>        
void system<R,NDIM>::init_(void)
{
  const parameters &param = get_parameters();
  cells_.init(static_cast<size_t>(ceil(double(param.template get_named_parm<Z>("N_particle"))
               /double(param.template get_named_parm<Z>("particle_per_cell")))));

  // note: "apply_directed" will override the following initialization of "t_max_", if it is called:
  if (param.has_named_parm("t_max"))
    set_t_max(param.template get_named_parm<R>("t_max"));
  else
    set_t_max(zero<R>());  
}
           
#if 0
/*
* create wrapped-clone of particle, using wrapping information
*   from cell_array neighbor lists:
*/
template <class R, size_t NDIM>        
typename system<R,NDIM>::particle* system<R,NDIM>::wrap_particle(const particle* src, bool wrap, const row_major_index& wrap_dim)const
{
  const parameters &param(get_parameters());
  particle *dest(src->clone());
  if (wrap){
    // use negative linear-index to indicate wrap;
    //   don't allow multiple-wrapping:
    assert(src->linear_index() >= 0); 
    dest->linear_index() *= -1;

    // image position across periodic B.C. (i.e. extend virtually):
    // (note: can't be done as one operation with current ntuple methods as elements of wrap-dim are "long")
    for(size_t n = 0; n < NDIM; ++n)
      dest->current_state().position()[n] += integer<R>(wrap_dim[n]) * param.L;  
  }
  
  return dest;
}
#endif

/* 
 * hooks for initialization and update of per-particle system statistics:
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::init_statistics(void)const
{ }

template <class R, size_t NDIM>        
void system<R,NDIM>::update_statistics(const particle* p)const
{ } 

/*
 * hook for application of system temperature control:
 * (note: this is a per-step hook, rather than a per-particle hook)
 */
template <class R, size_t NDIM>        
R system<R,NDIM>::thermostat_v_factor(void)const
{ return one<R>(); }    

/**
* @brief Offset that was applied at "init_state" to obtain positive-definite positions from incoming state.
*/ 
template <class R, size_t NDIM>        
inline const ntuple<R,NDIM>& system<R,NDIM>::position_offset(void)const
{ return position_offset_; }

/**
* @brief Set current position offset (used by "init_state" to obtain positive definite positions):
*/ 
template <class R, size_t NDIM>        
inline void system<R,NDIM>::set_position_offset(const ntuple<R,NDIM>& offset)
{ position_offset_ = offset; } 


/**
 * @brief Set the maximum-time limit value used by the completion test method (see "complete").
 *   This is a protected method: end-user applications set the parameters "t_max" attribute, and cannot use this method.
 */
template <class R, size_t NDIM>        
inline void system<R,NDIM>::set_t_max(const R& t)const
{
  // note: t_max_ == 0.0 indicates that it is not being applied as a time-limit:
  assert(t >= zero<R>());
  t_max_ = t;
}

/**
 * @brief Time limit value used by the completion test method (see "complete").
 */
template <class R, size_t NDIM>        
inline const R& system<R,NDIM>::t_max(void)const { return t_max_; }
 
                 
// either  all of "jam_particles", "move_particles", and "process_event"
//   and/or "step_" *must* be defined (in the latter case, the former may remain stubs)
template <class R, size_t NDIM>        
bool system<R,NDIM>::jam_particles(event* e)
{ return true; }

/*
 * implementation note:
 *  At present, "move_particles" should NOT require a mutex lock on the particle:
 *    as this method is threaded by _cell_, and particles are not shared between cells,
 *    there is no _valid_ situation where more than one thread of this method would access a given
 *    particle (read _or_ write) at the same time. 
 */
template <class R, size_t NDIM>        
bool system<R,NDIM>::move_particles(cell* c, event* E)
{
  bool status(true);
  const parameters &param __attribute_unused__ (get_parameters());  
  
  // move all particles in cell to kinematic coordinates of event:
  for(typename cell::particle_list_type::iterator itP = c->particles().begin(), itPEnd = c->particles().end();
      itP != itPEnd;
      ++itP){
    // implement extrinsic _and_ intrinsic move (i.e. call system::move, rather than particle::move):  
    move(*itP, E->time());    
  }
  
  return status;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::process_event(event* e)
{ return true; }

template <class R, size_t NDIM>        
bool system<R,NDIM>::step_(event_list& events)
{
  bool status(true);

  // processing JAM events and FINITE-TIME events now done in separate cycles.  
  assert(events.empty() 
         || !events.front()->is_jam() 
         || events.back()->is_jam()); 
  
  #if defined(__USE_PTHREAD)
  typename event_list::iterator it_jam_end(events.end());
  // last event is next finite-time event, all preceeding events will be particle-jam events:
  if (!events.empty() && !events.back()->is_jam())
    --it_jam_end;
  parallelUtil::parallel_for< system<R,NDIM>, bool, bool (system<R,NDIM>::*)(event*), typename event_list::iterator> 
    (const_cast<system<R,NDIM>*>(this), &system<R,NDIM>::jam_particles, events.begin(), it_jam_end);
  
  if (!events.empty() && !events.back()->is_jam()){
  
    parallelUtil::parallel_for< system<R,NDIM>, bool, bool (system<R,NDIM>::*)(cell*, event*), typename cell_array::iterator, event*>  
      (const_cast<system<R,NDIM>*>(this), &system<R,NDIM>::move_particles, cells_.begin(), cells_.end(), events.back()); 
    
    process_event(events.back());
    
  }  
  #else
  typename event_list::iterator it_jam_end(events.end());
  // last event is next finite-time event, all preceeding events will be particle-jam events:
  if (!events.empty() && !events.back()->is_jam())
    --it_jam_end;  
    
  for(typename event_list::iterator itE = events.begin(), itEEnd = it_jam_end;
      itE != itEEnd;
      ++itE)
    jam_particles(*itE);
    
  if (!events.empty() && !events.back()->is_jam()){  
    for(typename cell_array::iterator itC = cells_.begin(), itCEnd = cells_.end();
        itC != itCEnd;
        ++itC)
      move_particles(*itC, events.back());

    process_event(events.back());  
  }
  #endif
  return status;
}
     
/**
 * @brief protected apply method to be used by the "apply" method when a target-distribution is specified.
 * @param[in] arg: 
 *    - @b position[opt]: list of particle positions (cartesian coordinates)
 *    - @b velocity[opt]: when corresponding positions are supplied, a list of particle vector velocities
 *    .
 *    (General system-state initialization will have been extracted from "arg" prior to this method).
 * @param[out] val:
 *    - @b position: evolved particle positions
 *    - @b velocity: evolved particle vector velocities
 *    .    
 *   (at successful return from this method, "val" should be ready for "apply" to use as its own return-value)
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::apply_directed_(void)const
{
  using python_util::simple_object_base;
  using python_util::extract;
  using TMatrix::inDomain; 

  // implementation note: both data and gradient must use "dense_vector_ref" type of interpolator,
  //   otherwise initialization of the gradient _vector_ will copy all of the gradient component arrays.
  #if 0 // --- moved to outer scope: ---
  typedef linalg::linear_interpolator<gmm::dense_vector_ref<const R*>,NDIM> interpolator;
  #endif 
   
  const parameters& param(get_parameters());
  const object_map& dist_param(param.get_named_object("target_distribution")->template as<object_map>());
  
  /* ******************************************
  default_target_distribution_parameter = {
    'coord_interval':((0.0, 0.0, 0.0), tuple([default_parameters['L'] for n in range(default_parameters['NDIM'])])), 
    'data':NP.array([],dtype=float), 'data_shape':tuple([20 for n in range(default_parameters['NDIM'])]), 'nearest':False,
    'dt_range':(0.001, 0.5), 'T_B':0.1}
  */
  
  // extract target distribution parms:
  const simple_object_base *pobj(simple_object_base::get_named_object(dist_param, "data"));
  const R* data_ptr(pobj->template ptr<R>());
  const size_t data_len(pobj->size());
  gmm::dense_vector_ref<const R*> dist_data(data_ptr, data_ptr + data_len); 
    
  pobj = simple_object_base::get_named_object(dist_param, "data_shape");
  ntuple<size_t,NDIM> data_shape;
  extract<ntuple<size_t,NDIM> >(data_shape, pobj);
  
  pobj = simple_object_base::get_named_object(dist_param, "coord_interval");
  ntuple_interval<R,NDIM> coord_interval;
  extract<ntuple<R,NDIM> >(coord_interval.start(), pobj->template as<object_list>()[0]);
  extract<ntuple<R,NDIM> >(coord_interval.end(), pobj->template as<object_list>()[1]);
  
  const bool nearest(simple_object_base::template get_named_parm<bool>(dist_param, "nearest")); 
  
  pobj = simple_object_base::get_named_object(dist_param, "dt_range");
  const R 
    &dt_min((pobj->template as<object_list>()[0])->template as<R>()),
    &dt_max((pobj->template as<object_list>()[1])->template as<R>()),
    &T_B(simple_object_base::template get_named_parm<R>(dist_param, "T_B")),
    &T_0(param.template get_named_parm<R>("T_0")),
    &t_max(param.template get_named_parm<R>("t_max"));

  pobj = simple_object_base::get_named_object(dist_param, "multistep_factors");
  const R 
    &step_contract_f((pobj->template as<object_list>()[0])->template as<R>()),
    &step_expand_f((pobj->template as<object_list>()[1])->template as<R>());

  interpolator target_interp(coord_interval, dist_data, data_shape, (nearest?interpolator::NEAREST_PERIODIC: interpolator::LINEAR_PERIODIC));
  
  // distribution gradient calculation and interpolation:
  std::vector<std::vector<R> > dist_gradient(NDIM);
  std::vector<interpolator> target_gradient_interp;
  target_gradient_interp.reserve(NDIM);
  for(size_t nd = 0; nd < NDIM; ++nd){
    linalg::central_diff(dist_gradient[nd], data_shape, dist_data, nd);
    target_gradient_interp.push_back(
      interpolator(coord_interval, 
        gmm::dense_vector_ref<const R*>(&(*dist_gradient[nd].begin()), &(*dist_gradient[nd].begin()) + dist_gradient[nd].size()),data_shape,
                    (nearest?interpolator::NEAREST_PERIODIC: interpolator::LINEAR_PERIODIC)));
  }
  
  // evolution time-step vars:
  R dt_step(sqrt(dt_min * dt_max)), // start at the geometric mean of the time-step range
    fitness(-inv(epsilon<R>())),    // start at negative max: distributions are not necessarily positive definite, or normalized
    max_f(fitness),
    last_t(zero<R>());
  
  // gradient-directed evolution to t_max:
  #if 0
  for(R t = min(dt_step,t_max); t <= t_max; t += max(min(dt_step, t_max - t), epsilon<R>())){
  #endif
  R t(zero<R>());
  int fsa_state(0); // 0:start, 1: expand, 2:contract
  do{
    fitness = evaluate_fitness_(target_interp);
    clog(logger::LOG_VERBOSE)<<"system<R,NDIM>::apply_directed_: evolution to time: "<<t<<", fitness: "<<fitness<<endl;
  
    if (fitness > max_f){
      if (1 == fsa_state){  
        // time-step success, expand dt_step:
        dt_step *= min(step_expand_f, dt_max/dt_step);
      }        
      else
      if (0 == fsa_state){
        fsa_state = 1;
      }  
      else
        fsa_state = 0;  
    }
    else{
      if (2 == fsa_state){
        // time-step failure, contract dt_step:
        dt_step *= max(step_contract_f, dt_min/dt_step);
      }        
      else
      if (0 == fsa_state){
        fsa_state = 2;
      }  
      else
        fsa_state = 0;  
    }
  
    max_f = max(fitness, max_f);
    last_t = t;
    
    // evolve system state:
    t += min(dt_step, t_max - t);

    clog(logger::LOG_VERBOSE)<<" --- dt: "<<(t - last_t)<<" --- "<<endl;
    
    // (regarding const_cast, it is a work-around to use "system" as a functor with "const" apply method)
    const_cast<system<R,NDIM>*>(this)->direct_velocities_(target_gradient_interp, T_0, T_B);

    set_t_max(t); // time-step completion condition.
    for(size_t nstep = 0; nstep < param.NSTEP; ++nstep)
      if (!const_cast<system<R,NDIM>*>(this)->step()) break;
  }
  while((t < t_max) && inDomain(dt_step, dt_min, dt_max, true));    

  if (inDomain(dt_step, dt_min, dt_max, true)){
    clog(logger::LOG_VERBOSE)<<"system<R,NDIM>::apply_directed_: evolution to time: "<<t<<", fitness: "<<fitness<<"\n"
       <<"  directed-evolution sequence completed."<<endl;  
  }
  else{
    clog(logger::LOG_QUIET)<<"system<R,NDIM>::apply_directed_: evolution to time: "<<t<<", fitness: "<<fitness<<"\n"
       <<"  WARNING: directed-evolution sequence terminated prematurely: dt_step clamped at boundary of range."<<endl;     
    std::cerr<<"system<R,NDIM>::apply_directed_: evolution to time: "<<t<<", fitness: "<<fitness<<"\n"
       <<"  WARNING: directed-evolution sequence terminated prematurely: dt_step clamped at boundary of range."<<endl;     
  }   
}  


/**
 * @brief Fitness evaluator based on particle spatial distribution with respect to a target distribution, 
 *    specified as an N-dimensional interpolator.
 */
template <class R, size_t NDIM>        
R system<R,NDIM>::evaluate_fitness_(const interpolator& target_interp)const
{
  R rval(zero<R>());
  size_t N_particle(0);
  
  // evaluate the particle-number normalized sum of the target distribution evaluated at the particle positions:
  for(typename std::vector<cell*>::const_iterator itC = cells_.data().begin(), itCEnd = cells_.data().end();
      itC != itCEnd;
      ++itC)
    for(typename cell::particle_list_type::const_iterator itP0 = (*itC)->particles().begin(), itP0End = (*itC)->particles().end();
        itP0 != itP0End;
        ++itP0, ++N_particle){    
      const ntuple<R,NDIM>& x((*itP0)->current_state().position()); 
      rval += target_interp(x);
    }        

  rval /= integer<R>(N_particle);
  return rval;     
} 

/**
 * @brief Re-initialize particle velocities based on the gradient of a target distribution.
 *   @param[in] target_gradient_interp  vector of interpolators for each component of the target-distribution gradient;
 *   @param[in] T  kinetic temperature of directed velocity component (randomized in the half-space);
 *   @param[in] T_B  Brownian temperature: kinetic temperature of non-directed component (completely randomized).
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::direct_velocities_(const std::vector<interpolator>& target_gradient_interp, const R& T, const R& T_B) 
{
  const parameters &param __attribute_unused__ (get_parameters());
  const R 
    v_scale(sqrt(T)), // (inv(sqrt(integer<R>(param.N_particle)))),
    v_B_scale(sqrt(T_B));
  
  for(typename std::vector<cell*>::iterator itC = cells_.data().begin(), itCEnd = cells_.data().end();
      itC != itCEnd;
      ++itC){
    for(typename cell::particle_list_type::iterator itP0 = (*itC)->particles().begin(), itP0End = (*itC)->particles().end();
        itP0 != itP0End;
        ++itP0){
      particle *p0(*itP0);  
      const state &current_state(p0->current_state());
      p0->rotate_state_buffer();
      state &new_state(p0->current_state());
      new_state = current_state; // leave everything the same, except for velocity (and state event marker).
    
      // note: this should apply a "drift" velocity, 
      //   and at present, it is not quite correct 
      //   (with respect to randomization and later normalization of the directed component).
      // I will leave this alone to test that the initial output matches that from the python version, and
      //   after that, it should be fixed.
      
      const ntuple<R,NDIM>& x(current_state.position());
      ntuple<R,NDIM> v_directed, v_random; 
      for(size_t nd = 0; nd < NDIM; ++nd){
        R gradient = target_gradient_interp[nd](x);
        v_directed[nd] = gradient * random<R>();    // randomize in the half-sphere: *** NOT correct! ***
        v_random[nd] = random<R>() - ratio<R>(1,2); // completely random component
      }
      
      // v_directed.normalize(); // *** NOT correct! ***
      #if 0
      v_directed *= p0->kinetic_velocity(T);
      #else
      v_directed *= v_scale;
      #endif
      
      // v_random.normalize(); // *** NOT correct! ***
      #if 0
      v_random *= p0->kinetic_velocity(T_B);
      #else
      v_random *= v_B_scale;
      #endif
      
      new_state.velocity() = v_directed + v_random;
      new_state.set_event_type(event::TIMESTEP_EVENT);
    }        
  }   
} 

template <class R, size_t NDIM>        
inline const size_t system<R,NDIM>::DIM(void)const
{ return NDIM; }

template <class R, size_t NDIM>        
inline const size_t system<R,NDIM>::N_particle(void)const
{ return get_parameters().N_particle; }

template <class R, size_t NDIM>        
inline const R& system<R,NDIM>::L(void)const
{ return get_parameters().L; }


// length corresponding to edge of cell-cube:
template <class R, size_t NDIM>        
R system<R,NDIM>::cell_edge_length(void)const
{ return L()/integer<R>(cells_.N1()); }

// center position of cell-cube:
template <class R, size_t NDIM>        
ntuple<R,NDIM> system<R,NDIM>::cell_center(const cell* cell_)const
{
  const row_major_index &nx(cell_->indices());
  const R& cell_L1(cell_edge_length());
  
  ntuple<R,NDIM> center;
  for(size_t n = 0; n < NDIM; ++n)
    center[n] = (integer<R>(nx[n]) + ratio<R,long>(1,2)) * cell_L1;
  return center;
}

// ntuple-interval representing cell-cube:
template <class R, size_t NDIM>        
ntuple_interval<R,NDIM> system<R,NDIM>::cell_cube(const cell* cell_)const
{
  const parameters &param(get_parameters());
  const row_major_index &nx(cell_->indices());
  const R& cell_L1(cell_edge_length());

  ntuple<R,NDIM> start_, end_;
  for(size_t n = 0; n < NDIM; ++n){
    start_[n] = integer<R>(nx[n]) * cell_L1;
    end_[n] = start_[n] + cell_L1;
  }
  
  // left closure:  x \in [start_, end_)  ; with resolution determined by the _system_ hypercube dimension
  return ntuple_interval<R,NDIM>(start_, end_, zero<R>(), epsilon<R>()*param.L); 
}

// ntuple-interval representing entire system domain:
template <class R, size_t NDIM>        
ntuple_interval<R,NDIM> system<R,NDIM>::system_cube(void)const
{
  const parameters &param(get_parameters());
  ntuple<R,NDIM> start_, end_(ntuple<R,NDIM>::filled(param.L));
  start_.clear();
  
  return ntuple_interval<R,NDIM>(start_, end_, zero<R>(), epsilon<R>()*param.L);  
}

/*
 * test for particle cell-exit event:
 *   returns true if particle will exit cell
 *   dimension, face: face-pair corresponding to cube face of exit
 *   time: entry-time into the other cell (i.e. exit-time + epsilon)
 */
template <class R, size_t NDIM>        
bool system<R,NDIM>::cell_exit(const particle* particle_, const cell* cell_, R& time, size_t& dimension, size_t& face)const    
{
  bool test(false);

  if (!particle_->frozen()){
    const state &s1(particle_->current_state());
    const ntuple<R,NDIM>
      &p1(s1.position()),
      &v1(s1.velocity());
    const R& t1(s1.time());

    #if 0 // --------------- *** DEBUG *** :off ------------------------------
    // only deal with cube-exit case:
    const ntuple_interval<R,NDIM> cube(cell_cube(cell_));
    
    #if 1 || defined(_debug_print_)
    // *** DEBUG ***
    if (!cube.closure(p1)){
      cerr<<"ERROR: "<<particle_->id_str()<<"@ "<<p1<<", NOT IN CELL-CUBE: \n";
      cube.write(cout, &p1);
      cerr<<"cell edge dimension (units of epsilon<T>()): "<<cell_edge_length()/epsilon<R>();
      cerr<<"\n";
      particle_->write_state_buf(cerr);
      #if defined(state_partner_lists)
      if (!particle_->current_state().partners().empty()){
        cerr<<"\n   PARTNER "<<particle_->current_state().partners().back()->id_str()<<":\n";
        particle_->current_state().partners().back()->write_state_buf(cerr);
      }
      #endif
      cerr<<endl;
      throw std::runtime_error("-- ABORTING --");
    }
    #endif
    #if 0 // *** DEBUG *** : off
    assert(cube.closure(p1));  
    #endif
    
    const R delta(integer<R>(10)*epsilon<R>()); // test implied end location s.t. it is in interior of adjacent cell cube.  
    
    #else // ------------------ *** DEBUG *** :on --------------------------
    // only deal with cube-exit case:
    const ntuple_interval<R,NDIM> cube(cell_cube(cell_));
    
    if (!cube.closure(p1)){
      cerr<<"ERROR: "<<particle_->id_str()<<"@ "<<p1<<" NOT IN CELL-CUBE: \n";
      cube.write(cout, &p1);
      cerr<<"cell edge dimension (units of epsilon<T>()): "<<cell_edge_length()/epsilon<R>();
      cerr<<"\n";
      particle_->write_state_buf(cerr);
      #if defined(state_partner_lists) && 0
      if (!particle_->current_state().partners().empty()){
        cerr<<"\n   PARTNER "<<particle_->current_state().partners().back()->id_str()<<":\n";
        particle_->current_state().partners().back()->write_state_buf(cerr);
      }
      #endif
      cerr<<endl;
      throw std::runtime_error("-- ABORTING --");
    }    
    const R delta(zero<R>());       
    
    #endif // ------------------ end: *** DEBUG *** 
    
    R dt(inv(epsilon<R>())),  // init to max time.
      dt_(zero<R>());

    // with respect to cell-transit, exit events are determined by particle center position,
    //   particle radius is ignored:
    for(size_t n = 0; n < NDIM; ++n){
      const R 
        &x(p1[n]),
        &v(v1[n]);

      bool test_(false);
      size_t face_(0);

      #if 1 // ------------- new version -----------------------
      if (v > epsilon<R>() * ((cube.end()[n] + delta) - x)){
        test_ = true;
        dt_ = ((cube.end()[n] + delta) - x)/v;
        face_ = 1;          
      }
      else
      if (fabs(v) > epsilon<R>() * (x - (cube.start()[n] - delta))){
        test_ = true;
        dt_ = (x - (cube.start()[n] - delta))/fabs(v);
        face_ = 0;        
      }
      #endif

      #if 0 // ------- previous version: ----------------
      if (fabs(v) > epsilon<R>()){
        test_ = true;
        if (v > zero<R>()){
          dt_ = ((cube.end()[n] + delta) - x)/v;
          face_ = 1;
        }  
        else{
          dt_ = (x - (cube.start()[n] - delta))/fabs(v);
          face_ = 0;
        }
      }
      #endif
      
      if (test_ && (dt_ < dt)){
        test = true;
        dt = dt_;
        dimension = n;
        face = face_;
      }  
    }

    if (test){
      time = t1 + dt;
      // In the following, _allow_ dt == 0 if t1 <= eps.  This indicates this is the first event-step after
      //   a system re-initialization; in this case, it is not an error if the previous system ended 
      //   directly after a cell-exit event, and the velocity vectors have been changed to return to the just-exited cell.
      if (!(dt > zero<R>()) && (t1 > epsilon<R>()))
        throw std::runtime_error("system<R,NDIM>::cell_exit: zero-time exit event detected");  
    }
  } 
      
  return test;
}


/*
 * move particle (deals with both extrinsic (i.e. kinematic-state related) and intrinsic changes):
 *   (calls particle.move(t) to update intrinsic change)
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::move(particle* p, const R& t)
{
  const parameters &param(get_parameters());
  
  // transfer information to new state:
  state &s0(p->current_state());
  p->rotate_state_buffer();
  state &s1(p->current_state());
  s1 = s0;
  
  R dt(t - s0.time());
  
  if (dt >= epsilon<R>()){ 
  
    // intrinsic state update:
    // (must happen _first_, as extrinsic update modifies state->time_):
    p->move(t);  

    // extrinsic state update:

    s1.position() += s0.velocity() * dt;
        
    // apply periodic B.C.:
    s1.position().mod_assign(param.L);
    
    // apply thermostat:
    s1.velocity() *= thermostat_v_factor();
      
    s1.time() = t;
    s1.set_event_type(event::MOVE_EVENT);
  }
  else
  if (t > epsilon<R>())
    // note: dt == 0 and t > eps => first event-step after system re-initialization, which may be zero-time event (see comments at cell_exit).
    std::cerr<<"WARNING: system<R,NDIM>::move: dt < epsilon"<<std::endl;       
}

/*
 * implement a cell-exit event, including any required position wrap-around:
 * (note: presently transfer only occurs through one face at a time (corner transfers are zero cross-section events...))
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::exit_cell(particle* p, size_t dimension, size_t face)
{
  const parameters &param __attribute_unused__(get_parameters());
  assert((face==0 || face==1) && (0 <= dimension <= NDIM));

  state &s0(p->current_state());
  p->rotate_state_buffer();
  state &s1(p->current_state());
  s1 = s0;
  
  cell *src(p->owner()), *dest(p->owner()->face_neighbor(dimension, face));

  // transfer from src particle-list to dest particle-list:
  src->transfer_particle(p, dest);

  // adjust position s.t. transfer is exact:
  //   (see row_major_index::ND_bin: "left_closure" \leftarrow  x \in [start, end))
  ntuple_interval<R,NDIM> cube(cell_cube(dest));
  if (face == 0)
    s1.position()[dimension] = cube.end()[dimension] - cube.right_epsilon();
  else
    s1.position()[dimension] = cube.start()[dimension];
  
  // mark the event in the state:
  s1.set_event_type(event::CELL_EXIT_EVENT);
}
       
/**
 * @brief Evolve system through next event.
 * Return false if completion criteria are satisfied.
 */
template <class R, size_t NDIM>        
bool system<R,NDIM>::step(void)
{
  #if 0 
  // *** DEBUG ***
  static bool PROCESS_JAM_EVENTS(true);
  #else
  // *** DEBUG ***
  static size_t JAM_ONLY_COUNT = 0;
  #endif
  
  bool complete_(false);
  event_list events;  
  
  #if defined(__USE_PTHREAD)  
  // initialize system per-particle statistics measurement:
  // (note: next_event accumulates per-particle statistics (used by "complete"))
  init_statistics();
    
  // event_list[s] are organized as [ <jam event 1>, <jam event 1>, ..., <next finite-time event>]:
  std::vector<event_list> vl_events 
    = parallelUtil::parallel_for< system<R,NDIM>, event_list, event_list (system<R,NDIM>::*)(const cell*)const, typename cell_array::iterator> 
        (const_cast<system<R,NDIM>*>(this), &system<R,NDIM>::next_event, cells_.begin(), cells_.end());
  
  if (!(complete_ = complete())){
    // combine thread event_lists:
    bool valid_next_event(false);
    events.push_back(new event(event::NULL_EVENT, inv(epsilon<R>()),NULL,NULL)); // position for next finite-time event (init next->time_ to max)
    event *next(events.back()); 

    for(typename std::vector<event_list>::const_iterator itVL =  vl_events.begin(), itVLEnd = vl_events.end();
        itVL != itVLEnd;
        ++itVL){
        
      for(typename event_list::const_iterator itL = (*itVL).begin(), itLEnd = (*itVL).end();
          itL != itLEnd;
          ++itL){
        event *E(*itL);
        if (E->is_jam())
          // put the event in the jam-list, if it doesn't duplicate an event already there:
          events.push_front_unique(E);
        else
        if (E->time() < next->time()){
          *next = *E; // copy event object
          valid_next_event = true;
        }       
      }        
    }

    if (!valid_next_event){
      // note: jam events may still exist:
      events.remove_next_event();
    }

    // clean-up
    vl_events.clear();
  }
  #else
  bool valid_next_event(false); // true when any finite-time events actually exist
  events.push_back(new event(event::NULL_EVENT, inv(epsilon<R>()),NULL,NULL)); // position for next finite-time event (init next->time_ to max)
  event *next(events.back());                               

  // initialize system per-particle statistics measurement:
  init_statistics();

  for(typename std::vector<cell*>::iterator itC = cells_.data().begin(), itCEnd = cells_.data().end();
      itC != itCEnd;
      ++itC){
    // get event list for cell:
    // (jam events always added to start of list; last event is next finite-time event)
    event_list events_ = next_event(*itC);
    for(typename event_list::const_iterator itE_ = events_.begin(), itE_End = events_.end();
        itE_ != itE_End;
        ++itE_){
      const event *E_(*itE_);
      if (E_->is_jam())
        // put the event in the jam-list, if it doesn't duplicate an event already there:
        events.push_front_unique(E_);
      else
      if (E_->time() < next->time()){
        // copy event object: each event_list owns its event data:
        *next = *E_;  
        valid_next_event = true;
      }         
    }  
  }

  if (!valid_next_event){
    // note: jam events may still exist:
    events.remove_next_event();
  }

  complete_ = complete();
  #endif

  if (!complete_){
    #if 0
    // *** DEBUG ***
    // check "next_event" after _any_ modification of kinematic state:
    if (PROCESS_JAM_EVENTS){
      if (!events.empty() && !events.back()->is_jam())
        events.remove_next_event();
      // allow "step_" to process _only_ JAM events.
    }
    else{
      if (!events.empty() && events.front()->is_jam()){

        typename event_list::iterator 
          it_jam_begin(events.begin()),
          it_jam_end(events.end());
        if (!events.front()->is_jam())
          ++it_jam_begin;
        if (!events.back()->is_jam())
          --it_jam_end;
        events.data().erase(it_jam_begin, it_jam_end);  
      }
    }
    PROCESS_JAM_EVENTS = !PROCESS_JAM_EVENTS;
    #else
    // *** DEBUG ***
    // if there are jam events, process only the jam events:
    assert(!events.empty());
    if (events.front()->is_jam()){
      if(!events.back()->is_jam())
         events.remove_next_event();
      ++JAM_ONLY_COUNT;      
      if (JAM_ONLY_COUNT > 10){
        #if defined(event_history)
        cout<<"EVENT LIST:\n";
        write_event_history(cout);
        cout<<endl;
        typename event_list::list_type::reverse_iterator itE = event_history_.data().rbegin(), itEEnd = event_history_.data().rend();
        // for some reason, checking itE != itEEnd doesn't work in the next lines...
        for(size_t ne = 0; 
            /* (itE != itEEnd) && */ (ne < JAM_ONLY_COUNT); 
            ++ne, ++itE){
          event *E(*itE);
          if (E->is_jam()){
            // mark for special display:
            E->particle0()->current_state().velocity().clear();
            E->particle1()->current_state().velocity().clear();            
          }          
        }
        #endif
        throw std::runtime_error("ABORT: JAM_ONLY_COUNT exceeds 10");
      }  
    }
    else
      JAM_ONLY_COUNT = 0; // reset the count   
    #endif

    step_(events);
  }
  
  return !complete_;
}

// completion test will be "OR" chain down through derived classes:
template <class R, size_t NDIM>        
bool system<R,NDIM>::complete(void)const
{ return false; }

// ================== apply related methods: ======================

/**
 * @brief Check for consistency between a given argument and the present value of the functor parameters.
 * This method does not check for state consistency itself (e.g. lack of particle overlap); this is assumed.
 */
template <class R, size_t NDIM>        
void system<R,NDIM>::arg_valid_check(const object_map& arg)const
{
  using python_util::simple_object_base;
  const parameters& param(get_parameters());
  
  if (simple_object_base::has_named_parm(arg, "position")){
    const simple_object_base* pobj = simple_object_base::get_named_object(arg, "position");
    
    // allow zero-size as parameter "place-holder" (i.e. to name the possible parameter for end-user usage info):
    // (for position/velocity vectors: corresponding simple_objects are RANK=1 number-type ("size()" returns number of scalar entries)).
    if (!pobj->is_empty() && (pobj->size() != param.N_particle * NDIM))
      throw std::runtime_error("system<R,NDIM>::arg_valid_check: number of entries in position list doesn't match number of particles");
      
    if (simple_object_base::has_named_parm(arg, "velocity")){
      pobj = simple_object_base::get_named_object(arg, "velocity");
      if (!pobj->is_empty() && (pobj->size() != param.N_particle * NDIM))
        throw std::runtime_error("system<R,NDIM>::arg_valid_check: number of entries in velocity list doesn't match number of particles");
    }
  }
  else{
    if (simple_object_base::has_named_parm(arg, "velocity"))
      throw std::runtime_error("system<R,NDIM>::arg_valid_check: velocity list provided without corresponding position list");
  }
}

#if 0
// base-class: pure virtual method

/**
 * @brief generic apply method
 * @param[in] arg: 
 *    - @b position[opt]: list of particle positions (cartesian coordinates)
 *    - @b velocity[opt]: when corresponding positions are supplied, a list of particle vector velocities
 *    .
 * When particle positions are not supplied, values for number of particles "N_particle", initial temperature "T_0", radial growth rate "dr"
 *   will be taken from parameters to build the initial system state prior to evolution.
 * @param[out] val:
 *    - @b position: evolved particle positions
 *    - @b velocity: evolved particle vector velocities
 *    .    
 */
template <class R, size_t NDIM>        
bool system<R,NDIM>::apply(const python_util::generic_object arg, python_util::generic_object val)const=0;
#endif

/// bool return version compatible with "classic" functor implementation:
template <class R, size_t NDIM>        
bool system<R,NDIM>::setParameters(const parameters& param)
{
  bool status(true);
  try{
    set_parameters(param);
  }
  catch (std::string& msg){
    setStatusString(msg);
    status = false;
  }
  return status;
}

template <class R, size_t NDIM>        
inline void system<R,NDIM>::setStatusString(const std::string& msg)const
{ status_ = msg; }

template <class R, size_t NDIM>        
inline const std::string& system<R,NDIM>::getStatusString(void)const
{ return status_; }

// =======================================================

template <class R, size_t NDIM>
system<R,NDIM>& system<R,NDIM>::operator=(const system& other)
{
  set_parameters(other.get_parameters());
  cells_ = other.cells_;
  t_max_ = other.t_max_;
  return *this;
}

template <class R, size_t NDIM>        
bool system<R,NDIM>::writeBinary(commUtil::abstractCommHandle *fp)const
{
  bool status(true);
  // allow other usage of output file, which requires knowledge of NDIM:
  status = (status && writeBinary_(fp, NDIM));
  status = (status && get_parameters().writeBinary(fp));
  status = (status && cells_.writeBinary(fp));
  status = (status && writeBinary(fp, t_max_));
  return status;
}
    
template <class R, size_t NDIM>        
bool system<R,NDIM>::readBinary(commUtil::abstractCommHandle *fp)
{
  bool status(true);
  size_t NDIM_(0); // ignored input value (present for other output file usage)
  status = (status && readBinary_(fp, NDIM_));
  
  // read correct derived-class of parameters:
  parameters* pparam(reinterpret_cast<parameters*>(get_parameters().clone()));
  status = (status && pparam->readBinary(fp));
  if (status)
    set_parameters(*pparam);
  delete pparam;
  
  status = (status && cells_.readBinary(fp));
  status = (status && readBinary(fp, t_max_));
  return status;
}

template <class R, size_t NDIM>        
size_t system<R,NDIM>::binarySize(void)const
{
  size_t val(0);
  val += sizeof(size_t); // include NDIM in size calc.
  val += get_parameters().binarySize();
  val += cells_.binarySize();
  val += binarySize(t_max_);
  return val;
}


template <class R, size_t NDIM>        
system<R,NDIM>::~system(void)
{
  if (pparam_ != NULL){
    delete pparam_;
    pparam_ = NULL;
  }     
}

template <class R, size_t NDIM>        
system<R,NDIM>::system(bool initialize)
  : 
#if defined(__USE_PTHREAD)  
    mutex_(parallelUtil::CONST_PTHREAD_MUTEX_INITIALIZER),
    pair_mutex_(parallelUtil::CONST_PTHREAD_MUTEX_INITIALIZER),
    pair_cond_(parallelUtil::CONST_PTHREAD_COND_INITIALIZER),
    particle0_mutex_(parallelUtil::OMP_NUM_THREADS()),
    particle1_mutex_(parallelUtil::OMP_NUM_THREADS()),
#endif
    pparam_(initialize? new parameters(): NULL),
    t_max_(zero<R>())
{ }

template <class R, size_t NDIM>        
system<R,NDIM>::system(const system& other)
  : 
#if defined(__USE_PTHREAD)    
    mutex_(parallelUtil::CONST_PTHREAD_MUTEX_INITIALIZER),
    pair_mutex_(parallelUtil::CONST_PTHREAD_MUTEX_INITIALIZER),
    pair_cond_(parallelUtil::CONST_PTHREAD_COND_INITIALIZER),
    particle0_mutex_(parallelUtil::OMP_NUM_THREADS()),
    particle1_mutex_(parallelUtil::OMP_NUM_THREADS()),
#endif
    pparam_(NULL),
    t_max_(zero<R>())
{
  operator=(other); 
}


template <class R, size_t NDIM>   
void mjp_system<R,NDIM>::parameters::update_cache(bool derived, bool to_cache)
{
  // usage note: _local_ copies of parameters are the "cache"

  base_class::update_cache(true, to_cache); 

  if (to_cache){
    r_max = base_class::template get_named_parm<R>("r_max");
    T_max = base_class::template get_named_parm<R>("T_max");
  }
  else{
    base_class::template set_named_parm<R>("r_max", r_max, true);
    base_class::template set_named_parm<R>("T_max", T_max, true);
  }
  
  // most-derived class sets currency flag:
  if (!derived)
    root_class::set_cache_current(true);
}

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::parameters::valid_check(void)const
{
  using python_util::simple_object_base;
  
  base_class::valid_check();
  if(
      #if 0 // disable checks on non-cached parameters:
      || ((base_class::sticking_probability > zero<R>()) && (dr > zero<R>()))
      || (v < zero<R>())
      || (dr < zero<R>())
      #endif
      (r_max < zero<R>())
      || (T_max < zero<R>()))
    throw std::runtime_error("mjp_system<R,NDIM>::parameters::valid_check: invalid parameters");
  if (root_class::has_named_parm("density")){
    const simple_object_base* parm(root_class::get_named_object("density"));
    if (parm->template is<object_list>() && parm->template as<object_list>().size() != base_class::N_particle)
      throw std::runtime_error("mjp_system<R,NDIM>::parameters::valid_check: number of entries in density list doesn't match number of particles"); 
  }
  if (root_class::has_named_parm("dr")){
    const simple_object_base* parm(root_class::get_named_object("dr"));
    if (parm->template is<object_list>() && parm->template as<object_list>().size() != base_class::N_particle)
      throw std::runtime_error("mjp_system<R,NDIM>::parameters::valid_check: number of entries in dr list doesn't match number of particles"); 
  }    
}

template <class R, size_t NDIM>
typename python_util::options_map<typename mjp_system<R,NDIM>::C>* mjp_system<R,NDIM>::parameters::clone(void)const
{ return new parameters(*this); }

template <class R, size_t NDIM>
inline void mjp_system<R,NDIM>::parameters::copy(const parameters& other)
{
  base_class::copy(other);
}

template <class R, size_t NDIM>
typename mjp_system<R,NDIM>::parameters& mjp_system<R,NDIM>::parameters::operator=(const parameters& other)
{
  copy(other);
  return *this;
}

template <class R, size_t NDIM>   
typename mjp_system<R,NDIM>::parameters& mjp_system<R,NDIM>::parameters::operator=(const python_util::options_map<R>& other)
{
  root_class::copy(other);
  return *this;
}

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::parameters::write(std::ostream& os)const
{
  os<<"mjp_system<R,NDIM>::parameters: \n";
  root_class::template write<C,R,Z>(os);
  os<<"\n";
}

/** 
 * @brief Initialize from a simple_object_base*.
 */
template <class R, size_t NDIM>   
inline void mjp_system<R,NDIM>::parameters::extract(const python_util::simple_object_base* src)
{ 
  python_util::extract(*reinterpret_cast<root_class*>(this), src);
  update_cache();
  valid_check();
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::parameters::parameters(bool derived)
  : base_class(true), // derived => most-derived sets cache currency flag
    r_max(zero<R>()),
    T_max(zero<R>())
    #if 0
    , v(zero<R>()), dr(zero<R>())
    #endif 
{
  update_cache(derived, false);  // transfer from cache to options_map; most-derived sets cache currency flag 
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::parameters::parameters(const parameters& other)
{ copy(other); }


template <class R, size_t NDIM>   
mjp_system<R,NDIM>::parameters::parameters(const python_util::options_map<C>& other)
{
  python_util::options_map<C>::copy(other);
}

template <class R, size_t NDIM>
inline const typename mjp_system<R,NDIM>::parameters& mjp_system<R,NDIM>::get_parameters(void)const
{
  #if !defined(NDEBUG)
  const parameters& param(*dynamic_cast<const parameters*>(base_class::pparam()));
  assert(&param != NULL);
  #else
  const parameters& param(*reinterpret_cast<const parameters*>(base_class::pparam()));  
  #endif
  return param;
}


template <class R, size_t NDIM>
inline void mjp_system<R,NDIM>::set_parameters(const parameters& param)
{  base_class::set_parameters(static_cast<const typename base_class::parameters&>(param)); }



/*
 * Reset all sphere intrinsic "dr" attributes to zero:
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::freeze_growth(void)
{
  // set radial growth-velocity to zero:
  const parameters &param(get_parameters());
  assert(param.dr == zero<R>());
  
  // set each sphere's copy of same parameter to zero:
  for(typename cell_array::iterator itC = base_class::cells_.begin(), itCEnd = base_class::cells_.end();
      itC != itCEnd;
      ++itC){
    cell *cell_(*itC);  
    for(typename cell::particle_list_type::iterator itP = cell_->particles().begin(), itPEnd = cell_->particles().end();
        itP != itPEnd;
        ++itP){
      #if defined(NDEBUG)  
      sphere *particle_(dynamic_cast<sphere*>(*itP));
      assert(particle_ != NULL);
      #else
      sphere *particle_(reinterpret_cast<sphere*>(*itP));
      #endif
      particle_->dr() = zero<R>();  
    }  
  }
}

/*
 * Prepare system for DLA sequence: 
 *  -- reset radial growth-rate to zero
 *  -- fix seed particles
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::init_DLA(void)
{
  // set radial growth-velocity to zero:
  const parameters &param(get_parameters());
  
  // ---------------- note: parameters initialization block duplicates that in "init_state", although many will not be used here: ------------------

  // required non-cached parameters:
  const size_t 
    N_seed(static_cast<size_t>(param.template get_named_parm<Z>("N_seed")));
  
  // required vector (or constant) parameters:
  // (if vector, shall have N_particle entries)
  std::vector<R> dr, density;

  dr.reserve(param.N_particle);
  if (typeid(typename parameters::object_list) != param.named_parm_type("dr")){
    const R dr_(param.template get_named_parm<R>("dr"));
    dr.resize(param.N_particle);
    std::fill(dr.begin(), dr.end(), dr_);
  }
  else
    param.template get_named_parm<std::vector<R> >("dr", dr);

  typename std::vector<R>::const_iterator it_dr_max(std::max_element(dr.begin(), dr.end()));  
  if (it_dr_max == dr.end())
    throw std::runtime_error("mjp_system<R,NDIM>::init_DLA: max_element error return");  
  assert((*it_dr_max) == zero<R>());
  
  assert(N_seed > 0);
  
  double P_seed(double(N_seed)/double(param.N_particle));
  const size_t max_tries(2);
  
  #if defined(event_history)
  // *** DEBUG ***
  base_class::event_history_.clear();
  #endif
  
  // set each sphere's copy of same parameter to zero:
  size_t N_planted(0), tries(0);
  while ((N_planted < N_seed) && (tries++ < max_tries)){
    // the init_state loop is _also_ used to acquire system starting statistics.
    // clear attributes of system-statistics (including non-single-step attributes):
    system_stats_.clear(false);
    
    for(typename cell_array::iterator itC = base_class::cells_.begin(), itCEnd = base_class::cells_.end();
        itC != itCEnd;
        ++itC){
      cell *cell_(*itC);  
      for(typename cell::particle_list_type::iterator itP = cell_->particles().begin(), itPEnd = cell_->particles().end();
          itP != itPEnd;
          ++itP){
        #if defined(NDEBUG)  
        sphere *particle_(dynamic_cast<sphere*>(*itP));
        assert(particle_ != NULL);
        #else
        sphere *particle_(reinterpret_cast<sphere*>(*itP));
        #endif
        // set intrinsic radial growth-velocity attribute to zero:
        particle_->dr() = zero<R>();  

        // randomly freeze "N_seed" particles:
        if ((N_planted < N_seed) && (random<double>() < P_seed)){
          particle_->current_state().velocity().clear();
          particle_->freeze();
          ++N_planted;
        }

        update_statistics(particle_);
      }  
    }
  }
  
  if (N_planted < N_seed)
    throw std::runtime_error("mjp_system<R,NDIM>::init_DLA: seed particle initialization fails");
}

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::sphere::radius(void)const
{ return radius_; }

template <class R, size_t NDIM>
inline R& mjp_system<R,NDIM>::sphere::radius(void)
{ return radius_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::sphere::dr(void)const
{ return dr_; }

template <class R, size_t NDIM>
inline R& mjp_system<R,NDIM>::sphere::dr(void)
{ return dr_; }


template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::sphere::density(void)const
{ return density_; }

template <class R, size_t NDIM>
inline R& mjp_system<R,NDIM>::sphere::density(void)
{ return density_; }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::volume(void)const
{ return volume(radius_); }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::hypervolume(void)const
{ 
  return hypervolume(radius_); 
}

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::volume(const R& radius)
{ return ratio<R>(4,3) * pi<R>() * pow_n(radius, 3); }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::hypervolume(const R& radius)
{ 
  return C_n_ * pow_n(radius, static_cast<long>(NDIM)); 
}

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::mass(void)const
{ return mass(density_, radius_); }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::hypermass(void)const
{ return  hypermass(density_, radius_); }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::mass(const R& density, const R& radius)
{ return  density * volume(radius); }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::sphere::hypermass(const R& density, const R& radius)
{ return  density * hypervolume(radius); }

/*
 * kinetic temperature:
 */
template <class R, size_t NDIM>
inline const R mjp_system<R,NDIM>::sphere::T(void)const
{ return ratio<R>(1,2) * mass() * base_class::current_state().velocity().absSqr(); }

// static constant
template <class R, size_t NDIM>
const R mjp_system<R,NDIM>::sphere::C_n_ = pow(pi<R>(), ratio<R>(static_cast<long>(NDIM),2L)) / tgamma(ratio<R>(static_cast<long>(NDIM),2L) + one<R>()); 

/**
 * @brief Calculate time for maximum radius limit condition.
 *   @param[in] p0  sphere to test
 *   @param[out] time time when maximum radius will be reached
 *   @return  true if time is finite
 */
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::radius_limit(const sphere* p0, R& time)const
{
  bool rval(false);
  const parameters &param(get_parameters());
  
  if ((param.r_max > zero<R>()) && (p0->dr() > zero<R>())){
    const state &s0(p0->current_state());
    R dt((param.r_max - p0->radius())/p0->dr());
    time = s0.time() + dt; 
    rval = true;
  }
  
  return rval;
}

/*
 * test for pair collision between two spheres:
 *   returns true if collision will actually occur
 *   time: absolute time of collision
 *   jam: spheres are touching with no possible relative motion
 *        (or have motion that will not be resolved using pair-only dynamics)
 */
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::collision(const sphere* particle0_, const sphere* particle1_, R& time, bool& jam)const
{
  const parameters &param(get_parameters());

  const state 
    &s0(particle0_->current_state()), 
    &s1(particle1_->current_state());
  
  // default return values: no collision occurs, relative-time = 0:     
  bool test(false);
  time = s0.time();
  jam = false;

  #if 0
  // *** DEBUG ***
  cout<<"testing: "<<particle0_->id_str()<<" <--> "<<particle1_->id_str()<<", t = "<<time<<"\n"
      <<"   cell indices: "<<particle0_->owner()->indices()<<", "<<particle1_->owner()->indices()<<endl;
  #endif    


  // calculate wrapped-position (particle 1):
  ntuple<R,NDIM> wrap_position1(s1.position());
  row_major_index wrap_dim(0,false);
  if (row_major_index::wrap(particle0_->owner()->indices(), particle1_->owner()->indices(), base_class::cells_.shape(), wrap_dim)){
    for(size_t n = 0; n < NDIM; ++n)
      wrap_position1[n] += integer<R>(wrap_dim[n]) * param.L;
  } 

  if (!(particle0_->frozen() && particle1_->frozen())){

#if 0 // ------------------------ obsolete: incorrect: tests only coincidence along original p0 -- p1 normal direction: --------------------------

    // test implied collision location s.t. it happens before actual location (particle OVERLAP protection):
    const R delta(integer<R>(10)*epsilon<R>());  
  
    R dt(zero<R>()), v01(zero<R>()), v01_eff(zero<R>());

    // calculate normal vector: particle0_ --> particle1_:
    ntuple<R,NDIM> n01(wrap_position1);
    n01 -= s0.position();
    const R d01(n01.abs());
    assert(d01 > epsilon<R>());
    n01 /= d01;

    // relative-velocity and effective-relative-velocity: center of mass system:
    v01 = (s1.velocity() - s0.velocity()).dot(n01);
    v01_eff = v01 - (particle0_->dr() + particle1_->dr());

    // surface-surface distance:
    R r01_0(d01 - (particle0_->radius() + particle1_->radius()));

    /*
     * four cases:
     *   (1) particles have no relative motion;
     *   (2) particles are moving away from each other
     *   (3) particles are moving towards each other
     *   (4) particles are jammed together (special "not quite physical" case)
     *
     * solve the equation:
     *   t: 
     *     r01(t) = |(s0.position_ + s0.velocity_*t) - (s1.position_ + s1.velocity_*t)| -  ((r0 + dr0*t) + (r1 + dr1*t))
     *            = 0
     *     \rightarrow
     *  r01(0) = (v_01:i - dr0 - dr1)*t   # i: intrinsic (i.e. without radial expansion)
     *          = v_{0 1:eff}*t
     */


    if (r01_0 >= zero<R>()){

      if (v01_eff < zero<R>()){

        #if 1 
        // finite collision time:
        if (fabs(v01_eff) > r01_0 * epsilon<R>()){

          // kinematic collision time:
          dt = -r01_0/v01_eff;

          // adjust time by transit-time corresponding to overlap-padding space:
          R t_delta(zero<R>());
          if (fabs(v01_eff) > delta * epsilon<R>())
            t_delta = -delta/v01_eff;

          #if 1 // works OK, but possibly non-physical:
          if (dt > t_delta){
            dt -= t_delta/integer<R>(2);
            time += dt;
            test = true;         
            #if 0
            // *** DEBUG ***
            cout<<"case 1: "
              <<"r01_0, v01_eff: "<<r01_0<<", "<<v01_eff<<endl;
            #endif
          }  
          else
          if ((particle0_->dr() + particle1_->dr()) - fabs(v01) > integer<R>(10) * epsilon<R>()){
            // particle-jam: radial growth > relative velocity:
            #if 0
            // *** DEBUG ***
            cout<<"JAM EVENT"<<endl;
            #endif
            jam = true;
            test = true;
            #if 0
            // *** DEBUG ***
            cout<<"case 2: "
              <<"r01_0, v01_eff: "<<r01_0<<", "<<v01_eff<<endl;
            #endif
           }
          else{
            time += dt;
            test = true;         
            #if 0
            // *** DEBUG ***
            cout<<"case 3: "
              <<"r01_0, v01_eff: "<<r01_0<<", "<<v01_eff<<"\n"
              <<(particle0_->frozen()?"this FROZEN ":"")<<(particle1_->frozen()?", other FROZEN":"")<<"\n";
            s0.write(cout);
            cout<<"\n";
            s1.write(cout);
            cout<<endl;  
            #endif
           } 
          #endif

        } 
        #endif

      }

    }  
    else{
      const R R01(particle0_->radius() + particle1_->radius());
      
      cout<<"WARNING: OVERLAP( "
          <<"( "<<particle0_->id_str()<<" with "<<particle1_->id_str()<<" ), by "
          <<(d01 - R01)/epsilon<R>()<<" eps"
          <<": \n"
          <<"    | "<<s0.position()<<" - "<<wrap_position1<<" | = "<<d01
          <<" IS LESS THAN "<<particle0_->radius()<<" + "<<particle1_->radius()<<" = "<<R01<<endl;
    
      #if 1
      // *** DEBUG ***
      
      #if defined(event_history)
      cout<<"EVENT HISTORY:\n";
      system<R,NDIM>::write_event_history(cout);
      cout<<endl;
      #endif
      
      // mark particles so that they can be identified in system display (zero v special):
      const_cast<state&>(s0).velocity() *= integer<R>(0);
      const_cast<state&>(s1).velocity() *= integer<R>(0);
      throw std::runtime_error("ABORT");
      #endif
    }


    #if 0
    assert dist(s0.position(), wrap_position1) >= (particle0_->radius() + particle1_->radius())    
    #endif

    if (test){
      #if 1 || defined(_debug_print_)
      // *** DEBUG ***
      cout<<"collision time "
          <<"( "<<particle0_->id_str()<<" -> "<<particle1_->id_str()<<" )";
      if (!jam)
        cout<<", dt: "<<time<<", "<<dt<<"( = "<<dt/epsilon<R>()<<" eps)\n";
      else
        cout<<": *** PARTICLE-JAM ***\n";
      cout<<"  r0, r1, r01_0, v01, v01_eff: "
        <<particle0_->radius()<<", "<<particle1_->radius()<<", " 
        <<r01_0<<", "<<(s1.velocity() - s0.velocity()).dot(n01)<<", "<<v01_eff<<endl;
      #endif

      if (!jam && (dt < epsilon<R>())){
        cout<<"WARNING: sphere::collision: dt < epsilon"<<endl;
      }    
    }

#else // ------------------------ end obsolete version  ------------------------------------------------------------------------------------------


    // test implied collision location s.t. it happens before actual location (particle OVERLAP protection):
    const R delta(integer<R>(10)*epsilon<R>());

    #if 1
    // *** DEBUG ***
    // OVERLAP test:
    const R 
      d01((wrap_position1 - s0.position()).abs()),
      R01(particle0_->radius() + particle1_->radius());
    if (d01 - R01 < zero<R>()){

      cerr<<"WARNING: OVERLAP( "
          <<"( "<<particle0_->id_str()<<" with "<<particle1_->id_str()<<" ), by "
          <<(d01 - R01)/epsilon<R>()<<" eps"
          <<": \n"
          <<"    | "<<s0.position()<<" - "<<wrap_position1<<" | = "<<d01
          <<" IS LESS THAN "<<particle0_->radius()<<" + "<<particle1_->radius()<<" = "<<R01<<"\n"
          <<" -------- particle states: ----------------\n"
          <<" --- "<<particle0_->id_str()<<": --- \n";
      particle0_->write_state_buf(cerr);
      cerr<<"\n --- "<<particle1_->id_str()<<": --- \n";
      particle1_->write_state_buf(cerr);
      cerr<<"\n ----------------------------------- "<<endl;
      
      #if 1
      // *** DEBUG ***
      
      #if defined(event_history)
      cout<<"EVENT HISTORY:\n";
      system<R,NDIM>::write_event_history(cout);
      cout<<endl;
      #endif
      
      // mark particles so that they can be identified in system display (zero v special):
      const_cast<state&>(s0).velocity() *= integer<R>(0);
      const_cast<state&>(s1).velocity() *= integer<R>(0);
      throw std::runtime_error("ABORT");
      #endif    
    }
    #else
    assert dist(s0.position(), wrap_position1) >= (particle0_->radius() + particle1_->radius())    
    #endif
      
    /*
     * four cases:
     *   (1) particles have no relative motion;
     *   (2) particles are moving away from each other
     *   (3) particles are moving towards each other
     *   (4) particles are jammed together (special "not quite physical" case)
     *
     */

    R dt(zero<R>());
    // kinematic collision test result will be accepted provided it satisfies certain constraints 
    if (collision(
              s0.position(), particle0_->radius(), particle0_->dr(), s0.velocity(),
              wrap_position1, particle1_->radius(), particle1_->dr(), s1.velocity(),
              dt)){
      
      // Use magnitude of relative velocity to test particle-jam case,  
      //   and to generate an overlap-prevention padding space:
      R v_r((s1.velocity() - s0.velocity()).abs());

      // finite collision time:
      if (dt > epsilon<R>()){
       
        // adjust time by transit-time corresponding to overlap-padding space:
        R t_delta(zero<R>());
        if (v_r > delta * epsilon<R>())
          t_delta = delta/v_r;

        #if 1 // works OK, but possibly non-physical:
        if (dt > t_delta){
          dt -= t_delta/integer<R>(2);
          
          // transfer to return values:
          time += dt;
          test = true;         
        }  
        else
        if ((particle0_->dr() + particle1_->dr()) - v_r > integer<R>(10) * epsilon<R>()){
          // particle-jam: radial growth > relative velocity:
          #if 0
          // *** DEBUG ***
          cout<<"JAM EVENT"<<endl;
          #endif
          
          // transfer to return values:
          // (do not modify time)
          jam = true;
          test = true;
        }
        else{
          // very small relative collision time, but not obviously a particle-jam:
        
          // transfer to return values:
          time += dt;
          test = true;         
        } 
        #endif

      } 

    }
 
    if (test){
      #if defined(_debug_print_)
      // *** DEBUG ***
      cout<<"collision time "
          <<"( "<<particle0_->id_str()<<" -> "<<particle1_->id_str()<<" )";
      if (!jam)
        cout<<", dt: "<<time<<", "<<dt<<"( = "<<dt/epsilon<R>()<<" eps)\n";
      else
        cout<<": *** PARTICLE-JAM ***\n";
      cout<<"  r0, r1: "
        <<particle0_->radius()<<", "<<particle1_->radius()<<"\n"
        <<"---  s0: ---\n";
      particle0_->current_state().write(cout);
      cout<<"--- s1: ---\n";
      particle1_->current_state().write(cout);
      cout<<"-----------"<<endl;
      #endif
      if (!jam && (dt < epsilon<R>())){
        cout<<"WARNING: sphere::collision: dt < epsilon"<<endl;
      }    
    }


#endif // ----------------------- end new version  -----------------------------------------------------------------------------------------------

  }
        
  return test;    
}


/**
 * @brief Test sphere-sphere pair collision.
 *   This method only tests collision kinematics.
 *   returns true if collision will occur
 * @param[in]  p0: position of first particle
 * @param[in]  r0: radius of first particle
 * @param[in]  dr0: radial growth velocity of first particle
 * @param[in]  v0: velocity of first particle
 * @param[in]  p1: position of second particle
 * @param[in]  r1: radius of second particle
 * @param[in]  dr1: radial growth velocity of second particle
 * @param[in]  v1: velocity of second particle     
 * @param[out] dt: relative time of collision
 */
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::collision(
                           const ntuple<R,NDIM>& p0, const R& r0, const R& dr0, const ntuple<R,NDIM>& v0,
                           const ntuple<R,NDIM>& p1, const R& r1, const R& dr1,const ntuple<R,NDIM>& v1,
                           R& dt)
{
  bool test(false); // dt modified only in case of collision
  R t1(zero<R>()), t2(zero<R>()); // temporaries
  
  // relative position and velocity:
  ntuple<R,NDIM> p_r(p1), v_r(v1);
  p_r -= p0;
  v_r -= v0; 
  
  // combined radii, and combined radial growth velocity:
  R R_(r0), dR_(dr0);
  R_ += r1;
  dR_ += dr1;

  // coefficient and root vectors for quadratic formula:
  std::vector<R> vP(3, zero<R>());
  R D(zero<R>()); // discriminant
  std::vector<C> vRoot(2, zero<C>()); // must be complex
  
  t1 = v_r.absSqr();
  t1 -= sqr(dR_);
  vP[0] = t1;
  
  t1 = p_r.dot(v_r);
  t1 *= integer<R>(2);  
  t2 = R_;
  t2 *= dR_;
  t2 *= integer<R>(2);
  t1 -= t2;
  vP[1] = t1;
  
  t1 = p_r.absSqr();
  t1 -= sqr(R_);
  vP[2] = t1;
  
  quadraticFormula<std::vector<R>, std::vector<C> >(vP, vRoot, D);

  #if 0
  // *** DEBUG ***
  cout<<"  dt (roots): "<<vRoot[0]<<", "<<vRoot[1]<<", D: "<<D<<"\n"
      <<"  p0, r0, dr0, v0:\n"
      <<"    "<<p0<<", "<<r0<<", "<<dr0<<", "<<v0<<"\n"
      <<"  p1, r1, dr1, v1:\n"
      <<"    "<<p1<<", "<<r1<<", "<<dr1<<", "<<v1<<endl;
  #endif
      
  if (zero<R>() <= D){
    // D==0: duplicate real roots OK
    const R max_(max(real(vRoot[0]), real(vRoot[1])));
    
    if (zero<R>() <= max_){
      // positive root exists => future collision

      // two positive roots possible (e.g. if one particle is stationary):
      const R min_(min(real(vRoot[0]), real(vRoot[1])));
      if (zero<R>() <= min_)
        dt = min_;
      else
        dt = max_;

      test = true;
    }
  }

  return test;
}                           
                  
/*
 * sphere pair-collision (modifies _both_ spheres as appropriate):
 *   finite-time collisions and
 *   zero-time particle-jam events are treated uniformly
 *   (where applicable (i.e. _not_ in jam case), particles have been moved to kinematic coordinates of event,
 *      prior to this method)
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::collide(event* E)
{
  const parameters &param(get_parameters());

  #if 1
  // *** DEBUG ***
  static size_t ncollide(0);
  ++ncollide;
  #endif
  
  assert(E->particle1() != NULL);
  #if !defined(NDEBUG)
  sphere *particle0_(dynamic_cast<sphere*>(E->particle0()));
  assert(particle0_ != NULL);
  sphere *particle1_(dynamic_cast<sphere*>(E->particle1()));
  assert(particle1_ != NULL);
  #else
  sphere *particle0_(reinterpret_cast<sphere*>(E->particle0()));
  sphere *particle1_(reinterpret_cast<sphere*>(E->particle1()));
  #endif
    
  if (!(particle0_->frozen() && particle1_->frozen())){

    const R 
      m0 __attribute_unused__ (particle0_->mass()), 
      m1 __attribute_unused__ (particle1_->mass());

    // apply all changes to new states:
    state *s_(&particle0_->current_state());
    particle0_->rotate_state_buffer();
    particle0_->current_state() = *s_;

    s_ = &particle1_->current_state();
    particle1_->rotate_state_buffer();
    particle1_->current_state() = *s_;
        
    state
      &s0(particle0_->current_state()),
      &s1(particle1_->current_state());

    // NEED _CURRENT_ wrap-position and normal, _NOT_ position and normal at time of collision calculation:
    ntuple<R,NDIM> wrap_position(particle1_->current_state().position());
    row_major_index wrap_dim(0,false);
    if (row_major_index::wrap(E->event_cell()->indices(), E->particle1()->owner()->indices(), base_class::cells_.shape(), wrap_dim)){
      for(size_t n = 0; n < NDIM; ++n)
        wrap_position[n] += integer<R>(wrap_dim[n]) * param.L;
    }    

    // note: parallel unit vector calculation is still valid in jamming case: particle surfaces are touching, not centers:
    ntuple<R,NDIM> n01(wrap_position);
    n01 -= s0.position();
    n01.normalize();  


    if (!E->is_jam()){
      #if 1
      // *** DEBUG *** ======== KEEP THIS TEST ON FOR A WHILE... =============
      const R& r0(particle0_->radius()), &r1(particle1_->radius());
      if (fabs((wrap_position - s0.position()).abs() - (r0 + r1)) > integer<R>(10000)*epsilon<R>()){
        cerr<<"ERROR (collision #"<<ncollide<<"): COLLIDING PARTICLES DON'T TOUCH ("<<particle0_->id_str()<<" -> "<<particle1_->id_str()<<")\n"
          <<"  r0, r1, (r0 + r1), dx, |dx - (r0 + r1)|, |dx - (r0 + r1)|/eps: "
          <<r0<<", "<<r1<<", "
          <<(r0 + r1)<<", "<<(wrap_position - s0.position()).abs()<<", "
          <<fabs((wrap_position - s0.position()).abs() - (r0 + r1))<<", "
          <<fabs((wrap_position - s0.position()).abs() - (r0 + r1))/epsilon<R>()<<"\n"
          <<"  |v1 - v0|: "<<(s1.velocity() - s0.velocity()).abs()<<endl;
      }
      #endif
      
      // exact relative particle placement (this should eliminate any placement error and associated overlap glitches):
      // (note: 2*eps seems to work for most usage (at least in 2D), but relaxed to 10*eps may be better)
      s1.position() = s0.position();
      s1.position() += n01 * (particle0_->radius() + particle1_->radius() + integer<R>(10) * epsilon<R>());    
      s1.position().mod_assign(param.L);
    }


    if (particle0_->frozen() || particle1_->frozen())
      particle0_->stick(E);
    else{    

      // separate out components parallel(p), and perpendicular(s) to collision direction 
      //   (perpendicular components will be unmodified):
      ntuple<R,NDIM> 
        v0_p(n01),
        v0_s(s0.velocity()),
        v1_p(n01),
        v1_s(s1.velocity());
      v0_p *= s0.velocity().dot(n01);
      v0_s -= v0_p;
      v1_p *= s1.velocity().dot(n01);
      v1_s -= v1_p;

      // retain perpendicular component:
      ntuple<R,NDIM>
       v0(v0_s),
       v1(v1_s);

      // add explicit kinetic energy increment for radial growth:
      // (also added momentum ratio here to be consistent with momentum exchange below)
      #if 0
      // this version works with non-heterogenous dr (or heterogenous + incorrect momenta):
      v0 -= n01 * (integer<R>(2) * particle1_->dr());
      v1 += n01 * (integer<R>(2) * particle0_->dr());
      #else

        #if 1
          v0 += n01 * (particle0_->dr() * (m0-m1)/(m0+m1) - particle1_->dr() * (integer<R>(2) *m1)/(m0+m1)) * integer<R>(2);
          v1 += n01 * (-particle1_->dr() * (m1-m0)/(m0+m1) + particle0_->dr() * (integer<R>(2) *m0)/(m0+m1)) * integer<R>(2);
        #endif
        
      #endif

  #if 0 // ------- *** DEBUG *** : off
      if (!E->is_jam()){
        // exchange parallel momenta components
        v1_p *= m1;
        v1_0 /= m0;
        v0 += v1_p;
        
        v0_p *= m0;
        v0_p /= m1;
        v1 += v0_p;
      }
      else{  // THE FOLLOWING IS A CROCK:
        // zero relative parallel component
        // (this case deals with continuous growth-induced collisions, in glancing case;
        //    here the problem seems to be the back-and-forth cycle between particles, thus
        //    the reason for zeroing the relative velocity (apart from the "dr" part))
        ntuple<R,NDIM> 
          v_p(v0_p); 
        v_p += v1_p;
        v_p /= integer<R>(2);
        v0 += v_p;
        v1 += v_p;

        #if defined(_debug_print_)
        // *** DEBUG ***
        cout<<"GLANCE, v1_p, v2_p, v12_eff: "<<(v_p - n01 * particle1_->dr())<<", "<<v_p + n01 * particle0_->dr()<<", "<<(v1-v0).dot(n01) - (particle0_->dr() + particle1_->dr())<<endl;
        #endif
      }
  #else // *** DEBUG *** :on --------------------
        // exchange parallel momenta components
        #if 0
        // non-heterogenous r:
        v0 += v1_p;
        v1 += v0_p;
        #else
         // if (!E->is_jam()){
          v0 += v0_p * (m0-m1)/(m0+m1) + v1_p * (integer<R>(2) *m1)/(m0+m1);
          v1 += v1_p * (m1-m0)/(m0+m1) + v0_p * (integer<R>(2) *m0)/(m0+m1);
         //}
        #if 0
          else{
            ntuple<R,NDIM> 
              v_p(v0_p * (m0/(m0 + m1)));   
            v_p += v1_p * (m1/(m0 + m1));
            v0 += v_p;
            v1 += v_p;        
          }
        #endif
        #endif

        #if 0
        if (E->is_jam()){
          const R dR(particle0_->dr() + particle1_->dr());
          cout<<"\nJAM: initial relative parallel velocity (<= 0): "<<(-dR + (s1.velocity() - s0.velocity()).dot(n01))<<"\n"
              <<"   non-kinematic relative v: "<<-dR<<"\n"
              <<"   non-kinematic exchange increment: "
                <<(particle0_->dr() * (m0-m1)/(m0+m1) - particle1_->dr() * (integer<R>(2) *m1)/(m0+m1))
                  - (-particle1_->dr() * (m1-m0)/(m0+m1) + particle0_->dr() * (integer<R>(2) *m0)/(m0+m1))<<"\n"
              <<"  new relative parallel velocity (>= 0): "<<(-dR + (v1 - v0).dot(n01))<<endl;

          cout<<"initial kinematic relative v: "<<(v1_p.dot(n01) - v0_p.dot(n01))<<"\n"
              <<"  kinematic exchange increment: "
                <<-(v0_p.dot(n01) * (m0-m1)/(m0+m1) + v1_p.dot(n01) * (integer<R>(2) *m1)/(m0+m1))
                  + (v1_p.dot(n01) * (m1-m0)/(m0+m1) + v0_p.dot(n01) * (integer<R>(2) *m0)/(m0+m1))
                <<endl;
          cout<<"  final effective relative v: "<<(v1 - v0).dot(n01) - dR<<endl;
          cout<<"  final kinematic relative v: "<<(v1 - v0).dot(n01)<<endl;
          #if 0
          cout<<"p_r0, p_r1, pk_0, pk1:\n"
            <<particle0_->dr()*m0<<", "<<-particle1_->dr()*m1<<", "<<s0.velocity().dot(n01)*m0<<", "<<s1.velocity().dot(n01)*m1<<endl;
          #endif
        }
        #endif
  #endif

      // transfer velocities to new states:
      particle0_->current_state().velocity() = v0; // here _actual_ position is used (s0.position()), _not_ wrap_position!

      particle1_->current_state().velocity() = v1; // here _actual_ position is used (s0.position()), _not_ wrap_position!

      #if defined(state_partner_lists)
      if (E->is_jam()){
        // particle-jam adds to partner lists:
        particle0_->current_state().add_partner(particle1_);
        particle1_->current_state().add_partner(particle0_);      
        particle0_->current_state().set_event_type(event::JAM_EVENT);
        particle1_->current_state().set_event_type(event::JAM_EVENT);        
      }
      else{
        // collision resets partner lists:
        particle0_->current_state().partners().clear();
        particle1_->current_state().partners().clear();    
        #if 1
        // *** DEBUG ***
        // temporarily use partner list to retain "last collision" information:
        particle0_->current_state().add_partner(particle1_);
        particle1_->current_state().add_partner(particle0_);              
        #endif
        
        particle0_->current_state().set_event_type(event::COLLISION_EVENT);
        particle1_->current_state().set_event_type(event::COLLISION_EVENT);        
       }
      #else
      if (E->is_jam()){
        particle0_->current_state().set_event_type(event::JAM_EVENT);
        particle1_->current_state().set_event_type(event::JAM_EVENT);        
      }
      else{
        particle0_->current_state().set_event_type(event::COLLISION_EVENT);
        particle1_->current_state().set_event_type(event::COLLISION_EVENT);        
      }      
      #endif
    }
  }
  else
    throw std::runtime_error("mjp_system<R,NDIM>::collide: both particles frozen");  
}
        
        
template <class R, size_t NDIM>
typename system<R,NDIM>::particle::particle_kind mjp_system<R,NDIM>::sphere::kind(void)const
{ return base_class::SPHERE_PARTICLE; }


// update any _intrinsic_ attributes associated with system time change.
// (does _not_ modify particle "state" (here considered extrinsic))
template <class R, size_t NDIM>        
void mjp_system<R,NDIM>::sphere::move(const R& t)
{
  // grow:
  radius_ +=  dr_ * (t - base_class::current_state().time());
}


// overlap condition test:
template <class R, size_t NDIM>        
bool mjp_system<R,NDIM>::sphere::overlap(const particle* pother)const
{ 
  #if !defined(NDEBUG)
  const sphere &other(*dynamic_cast<const sphere*>(pother));
  assert(&other != NULL);
  #else
  const sphere &other(*reinterpret_cast<const sphere*>(pother));
  #endif
  const state 
    &s1(base_class::current_state()), 
    &s2(other.base_class::current_state());
  // surface-surface distance:
  R r12_0(dist(s1.position(), s2.position()) - (radius() + other.radius()));
  return (r12_0 < zero<R>());
}

/**
 * @brief Velocity magnitude corresponding to a specified kinetic temperature.
 *   This virtual method allows any method using temperature to initialize state velocity to be completely virtual.
 *   - for spheres, the zero-radius case uses density and "dr" (instead of radius) to calculate an effective velocity;
 *   - for abstract particles, mass is taken as unity.
 *   .
 */
template <class R, size_t NDIM>        
R mjp_system<R,NDIM>::sphere::kinetic_velocity(const R& T)const  
{
  R mass_(mass()); // use "mass" for kinetic calculations, not "hypermass"
  if (mass_ < epsilon<R>()){
    mass_ = sphere::mass(density(), dr());
    if (mass_ < epsilon<R>())
      throw std::runtime_error("mjp_system<R,NDIM>::sphere::kinetic_velocity: particle density and radius parameters imply infinite velocity");
  }
  
  return sqrt(integer<R>(2)*T/mass_);
}

template <class R, size_t NDIM>
typename mjp_system<R,NDIM>::sphere& mjp_system<R,NDIM>::sphere::operator=(const sphere& other)
{
  base_class::operator=(static_cast<const base_class&>(other));
  radius_ = other.radius_;
  dr_ = other.dr_;
  density_ = other.density_;
  return *this;
}

template <class R, size_t NDIM>
typename system<R,NDIM>::particle* mjp_system<R,NDIM>::sphere::clone(void)const
{ return new sphere(*this); }

template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::sphere::writeBinary(commUtil::abstractCommHandle* fp)const
{
  bool status(true);
  status = (status && (base_class::writeBinary(fp)));
  status = (status && (system<R,NDIM>::writeBinary_(fp, radius_)));
  status = (status && (system<R,NDIM>::writeBinary_(fp, dr_)));
  status = (status && (system<R,NDIM>::writeBinary_(fp, density_)));
  return status;
}

template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::sphere::readBinary(commUtil::abstractCommHandle* fp)
{
  bool status(true);
  status = (status && (base_class::readBinary(fp)));
  status = (status && (system<R,NDIM>::readBinary_(fp, radius_)));
  status = (status && (system<R,NDIM>::readBinary_(fp, dr_)));
  status = (status && (system<R,NDIM>::readBinary_(fp, density_)));
  return status;
}

template <class R, size_t NDIM>
size_t mjp_system<R,NDIM>::sphere::binarySize(void)const
{ 
  size_t val(0);
  val += base_class::binarySize();
  val += system<R,NDIM>::binarySize_(radius_);
  val += system<R,NDIM>::binarySize_(dr_);
  val += system<R,NDIM>::binarySize_(density_);
  return val;
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::sphere::~sphere(void)
{ }

template <class R, size_t NDIM>
mjp_system<R,NDIM>::sphere::sphere(void)
 : base_class(), radius_(zero<R>()), dr_(zero<R>()), density_(one<R>())
{ }

template <class R, size_t NDIM>
mjp_system<R,NDIM>::sphere::sphere(const sphere& other) 
 : base_class(), radius_(zero<R>()), dr_(zero<R>()), density_(one<R>())
{ operator=(other); }

template <class R, size_t NDIM>
mjp_system<R,NDIM>::sphere::sphere(const state& s, long linear_index, cell* owner, const R& r, const R& dr, const R& density)       
 : base_class(s, linear_index, owner), radius_(r), dr_(dr), density_(density)
{ }


template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::r_min(void)const
{ return r_min_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::r_max(void)const
{ return r_max_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::dr_min(void)const
{ return dr_min_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::dr_max(void)const
{ return dr_max_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::v_min(void)const
{ return v_min_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::v_max(void)const
{ return v_max_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::t_min(void)const
{ return t_min; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::t_max(void)const
{ return t_max_; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::T_min(void)const
{ return T_min; }

template <class R, size_t NDIM>
inline const R& mjp_system<R,NDIM>::statistics::T_max(void)const
{ return T_max_; }

template <class R, size_t NDIM>
inline R mjp_system<R,NDIM>::statistics::T_mean(void)const
{
  R T(T_mean_);
  T /= (N_step_update_ > 0? integer<R>(N_step_update_): one<R>()); 
  return T; 
}

template <class R, size_t NDIM>
const R& mjp_system<R,NDIM>::statistics::hypervolume(void)const
{ return hypervolume_; }

template <class R, size_t NDIM>
inline const size_t& mjp_system<R,NDIM>::statistics::N_frozen(void)const
{ return N_frozen_; }

template <class R, size_t NDIM>
inline const size_t& mjp_system<R,NDIM>::statistics::N_step_update(void)const
{ return N_step_update_; }

template <class R, size_t NDIM>
inline const size_t& mjp_system<R,NDIM>::statistics::N_total_update(void)const
{ return N_total_update_; }
        
// called to update mutable per-particle system statistics:   
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::statistics::update(const particle* p)
{
  // place NULL test here to simplify external usage:
  if (p != NULL){
    #if !defined(NDEBUG)
    const sphere* particle0_(dynamic_cast<const sphere*>(p));
    assert(particle0_ != NULL);
    #else
    const sphere* particle0_(reinterpret_cast<const sphere*>(p));
    #endif
    const state &s1(particle0_->current_state());
    // particle radius and growth velocity:
    const R &r(particle0_->radius()), &dr(particle0_->dr());
    if (r < r_min_)
      r_min_ = r;
    else
    if (r > r_max_)
      r_max_ = r;
    if (dr < dr_min_)
      dr_min_ = dr;
    else
    if (dr > dr_max_)
      dr_max_ = dr;    
    // velocity magnitude:
    const R v(s1.velocity().abs());
    if (v < v_min_)
      v_min_ = v;
    else
    if (v > v_max_)
      v_max_ = v;
    // time:
    const R &t(s1.time());
    if (t < t_min_)
      t_min_ = t;
    else
    if (t > t_max_)
      t_max_ = t;
     
    const R T(particle0_->T());
    if (T < T_min_)
      T_min_ = T;
    else
    if (T > T_max_)
      T_max_ = T;
    // T_mean_ holds accumulation "N_step_update_" times:  
    T_mean_ += T;
      
    // hypervolume_ holds total volume of all sphere:
    hypervolume_ += particle0_->hypervolume();
      
    if (particle0_->frozen())
      ++N_frozen_; 
      
    ++N_step_update_;
    ++N_total_update_;   
  }  
}

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::statistics::write(std::ostream& os)const
{
  os<<"system statistics: \n"
    <<"  r:           "<<r_min_<<" -> "<<r_max_<<"\n"
    <<"  dr:          "<<dr_min_<<" -> "<<dr_max_<<"\n"    
    <<"  |v|:         "<<v_min_<<" -> "<<v_max_<<"\n"
    <<"  t:           "<<t_min_<<" -> "<<t_max_<<"\n"
    <<"  T:           "<<T_min_<<" -> "<<T_max_<<", mean: "<<T_mean()<<"\n"
    <<"  hypervolume: "<<hypervolume_<<"\n"
    <<"  N_frozen: "<<N_frozen_<<"\n"
    <<"  N_update (step):    "<<N_step_update_<<endl
    <<"  N_update (total):   "<<N_total_update_<<endl;
}

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::statistics::clear(bool single_step)
{
  // note: all values are positive definite w.r.t. time,
  //   except "v", which may be positive semidefinite.
  r_min_ = inv(epsilon<R>());
  r_max_ = zero<R>();
  dr_min_ = inv(epsilon<R>());
  dr_max_ = zero<R>();
  v_min_ = inv(epsilon<R>());
  v_max_ = zero<R>();
  t_min_ = inv(epsilon<R>());
  t_max_ = zero<R>();
  T_min_ = inv(epsilon<R>());
  T_max_ = zero<R>();
  T_mean_ = zero<R>();
  hypervolume_ = zero<R>();
  N_frozen_ = 0;
  
  N_step_update_ = 0;
  if (!single_step)
    N_total_update_ = 0;
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::statistics::statistics(void)
{ clear(false); }


/**
 * @brief Initialize system state from parameters or optional incoming state.
 *   @param[in] incoming_state:
 *      - @b position: list of particle positions
 *      - @b velocity[opt]: list of particle velocities
 *      - @b radius[opt]:   list of particle radii.
 *      .
 */   
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::init_state(const python_util::simple_object_base* incoming_state)
{
  using TMatrix::seedRandom;
  using python_util::simple_object_base;
  const parameters &param(get_parameters());

  // ------------ initialize the random number generator: -----------------------------
  // parameters "random_seed" is allowed to have zero length, to be used as a placeholder: 
  std::string seed;
  if (param.has_named_parm("random_seed"))
    param.template get_named_parm<std::string>("random_seed", seed);
    
  if (seed.size() > 0)
    seedRandom<R>(seed);
  else
    seedRandom<R>(); // seed using /dev/urandom 
  // ----------------------------------------------------------------------------------
  
  // assemble required state values, including scalar, and vector parameters, and values from incoming state (if any):
    
  /*
   * Implementation note: many incoming values, from both parameters and incoming state are optional, or may be single-value, or may be vector.
   * Due to this complexity, all values will be transferred (or generated) to std::vector prior to transfer to the system.
   */
     
  // values from parameters:
  // (if vector, parameters::valid_check method has checked number of entries at "set_parameters")
  std::vector<R> dr, density;    
  
  dr.reserve(param.N_particle);
  if (param.get_named_object("dr")->is_scalar()){
    const R dr_(param.template get_named_parm<R>("dr"));
    dr.resize(param.N_particle);
    std::fill(dr.begin(), dr.end(), dr_);
  }
  else
    param.template get_named_parm<std::vector<R> >("dr", dr);
  
  #if 0 // --- omit, for the moment... ---
  typename std::vector<R>::const_iterator it_dr_max(std::max_element(dr.begin(), dr.end()));  
  if (it_dr_max == dr.end())
    throw std::runtime_error("mjp_system<R,NDIM>::init_state: max_element error return");
  const R dr_max(*it_dr_max);
  if ((*it_dr_max) >= v)
    clog(logger::LOG_NORMAL)<<"WARNING radial growth velocity >= kinematic velocity"<<endl;
  #endif

  density.reserve(param.N_particle);
  if (param.get_named_object("density")->is_scalar()){
    // see comment at "simple_object_base::is": type_info must be compared by-name for cross-module compatibility.
    const R density_(param.template get_named_parm<R>("density"));
    density.resize(param.N_particle);
    std::fill(density.begin(), density.end(), density_);
  }
  else
    param.template get_named_parm<std::vector<R> >("density", density);

  // values from incoming state:
  std::vector<R> radius;
  std::vector<ntuple<R,NDIM> > position, velocity;
  const ntuple_interval<R,NDIM> system_cube_(base_class::system_cube());    
   
  // check consistency of incoming state:
  if (NULL != incoming_state)
    arg_valid_check(incoming_state->as<object_map>());
  
  bool zero_radius(false); // note whether or not incoming radii were specified.  
  if ((NULL == incoming_state)  // empty "radius" list allowed as "place-holder" parm: 
      || !simple_object_base::has_named_parm(incoming_state->template as<object_map>(), "radius")
      || simple_object_base::get_named_object(incoming_state->template as<object_map>(), "radius")->is_empty()){

    radius.resize(param.N_particle);
    std::fill(radius.begin(), radius.end(), zero<R>());
    zero_radius = true;
  }
  else
    simple_object_base::get_named_parm(incoming_state->as<object_map>(), "radius", radius);

  if ((NULL == incoming_state) // empty "position" list allowed as "place-holder" parm:
      || !simple_object_base::has_named_parm(incoming_state->template as<object_map>(), "position")
      || simple_object_base::get_named_object(incoming_state->template as<object_map>(), "position")->is_empty()){
    // randomly place N_particle spheres:

    ntuple<R,NDIM> p_;
    for(size_t np = 0; np < param.N_particle; ++np){
      // random position vector (within hypercube):
      for(size_t n = 0; n < NDIM; ++n){
        R &x(p_[n]);
        x = random<R>(); // in [0,1]
        R L(system_cube_.end()[n] - system_cube_.start()[n]);
        x *= L;
        x += system_cube_.start()[n];
      }
      position.push_back(p_);
    }
    // retain offset to allow re-centering of output positions:
    ntuple<R,NDIM> offset_;
    for(size_t n = 0; n < NDIM; ++n)
      offset_[n] = -(system_cube_.end()[n] - system_cube_.start()[n])/integer<R>(2);
    base_class::set_position_offset(offset_);   
  }
  else{  
    simple_object_base::get_named_parm(incoming_state->as<object_map>(), "position", position);
    // offset the positions so that they are positive-definite:
    ntuple_interval<R,NDIM> extent_(ntuple_interval<R,NDIM>::extent(position));
    #if 0
    for(typename std::vector<ntuple<R,NDIM> >::iterator itP = position.begin(), itPEnd = position.end();
        itP != itPEnd;
        ++itP)
      (*itP) -= extent_.start();
    #endif
    base_class::set_position_offset(extent_.start());
  }  
       
  if ((NULL == incoming_state)  // empty "velocity" list allowed as "place-holder" parm:
      || !simple_object_base::has_named_parm(incoming_state->template as<object_map>(), "velocity")
      || simple_object_base::get_named_object(incoming_state->template as<object_map>(), "velocity")->is_empty()){
    
    const R& T_0(param.template get_named_parm<R>("T_0"));    
    ntuple<R,NDIM> v_;
    
    // in case where no incoming radii were specified, calculate v_magnitude as if
    //   the maximum mass particles (using "dr" instead of "r") have a mass  of 1.0:
    // (note: at the moment kinetics is calculated using the particle _mass_, not the _hypermass_).
    R m_default_scale(zero<R>());
    for(size_t np = 0; np < param.N_particle; ++np)
      m_default_scale = max(sphere::mass(density[np], dr[np]), m_default_scale);
    if (m_default_scale < epsilon<R>())
      throw std::runtime_error("mjp_system<R,NDIM>::init_state: given specified density and dr, effective velocity will be infinite");
    
    for(size_t np = 0; np < param.N_particle; ++np){
      
      R v_magnitude(zero_radius? 
          sqrt(integer<R>(2) * (sphere::mass(density[np], dr[np])/m_default_scale) * T_0)
          : sqrt(integer<R>(2) * sphere::mass(density[np], radius[np]) * T_0));
    
      // random velocity vector (of magnitude param.v):
      for(size_t n = 0; n < NDIM; ++n)
        v_[n] = random<R>() - ratio<R,long>(1,2);
      v_.normalize();
      v_ *= v_magnitude;

      velocity.push_back(v_); 
    }  
  }
  else 
    simple_object_base::get_named_parm(incoming_state->as<object_map>(), "velocity", velocity);  
        
  // initialize the cell-array, and other base-class attributes:
  base_class::init_(); 
  
  #if defined(event_history)
  // *** DEBUG ***
  base_class::event_history_.clear();
  #endif
       
  // the init_state loop is _also_ used to acquire system starting statistics.
  // clear attributes of system-statistics (including non-single-step attributes):
  system_stats_.clear(false);

  
  // place the particles in the cells:
  for(size_t sphere_number = 0; sphere_number < param.N_particle; ++sphere_number){

    // calculate correct cell:
    const row_major_index cell_indices_(row_major_index::ND_bin(position[sphere_number], system_cube_, base_class::cells_.shape()));
    cell *cell_(base_class::cells_[cell_indices_]);

    particle *particle_(
      new sphere(
        state(position[sphere_number], velocity[sphere_number], zero<R>()), 
        sphere_number, cell_, radius[sphere_number], dr[sphere_number], density[sphere_number]));

    cell_->particles().push_back(particle_);      
    update_statistics(particle_);    
  }
}


/**
 * @brief Extract system state to external representation.
 *   @param[out] dest:
 *      - @b NDIM: system dimensions
 *      - @b parameters: system parameters (as object map)
 *      - @b position: list of particle positions
 *      - @b velocity: list of particle velocities
 *      - @b radius:   list of particle radii.
 *      - @b density:   list of particle densities.
 *      .
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::extract_state(python_util::simple_object_base *dest)const
{
  using python_util::simple_object_base;
  using python_util::simple_object;
  using python_util::new_simple_object;
  const parameters &param(get_parameters());
  
  // transfer kinematic information to dense structures to take advantage of simple_object contiguous-memory functionality:

  // (implementation note: at the moment this requires 1 unnecessary copy at simple-object conversion, which can be bypassed by using the special
  //    simple_object array constructors and gmm::dense_vector_ref, but this would be premature optimization and I want to test the more normal
  //    code pathways).

  std::vector<ntuple<R,NDIM> > position(param.N_particle), velocity(param.N_particle);
  std::vector<R> radius(param.N_particle), density(param.N_particle);

  // offset subtracted at "init_state" to produce positive-definite  positions from incoming state:
  const ntuple<R,NDIM>& offset_ __attribute_unused__ (base_class::position_offset());
  
  typename std::vector<ntuple<R,NDIM> >::iterator it_position = position.begin(), it_velocity = velocity.begin();
  typename std::vector<R>::iterator it_radius = radius.begin(), it_density = density.begin();
  for(typename base_class::cell_array::const_iterator itCells = base_class::cells_.begin(), itCellsEnd = base_class::cells_.end(); 
      itCells != itCellsEnd;
      ++itCells)
  for(typename base_class::particle_list_type::const_iterator itP = (*itCells)->particles().begin(), itPEnd = (*itCells)->particles().end();
      itP != itPEnd;
      ++itP,
      ++it_position, ++it_velocity, ++it_radius, ++it_density){
    const sphere& sphere_(*dynamic_cast<const sphere*>(*itP));
    const state& state_(sphere_.current_state());
    
    *it_position = state_.position() /* + offset_ */;
    *it_velocity = state_.velocity();
    *it_radius = sphere_.radius();      
    *it_density = sphere_.density();      
  }        
  
  // transfer to return-value:
  dynamic_cast<simple_object<object_map>*>(dest)->clear(); 
  // transfer ownership of all simple_object_base pointers:
  #if 1
  // *** DEBUG ***
  simple_object_base::set_named_object(dest->as<object_map>(), "RANK", new_simple_object<Z>(static_cast<Z>(MPI::COMM_WORLD.Get_rank())), true); 
  simple_object_base::set_named_object(dest->as<object_map>(), "REFCOUNT_TEST", new_simple_object<std::string>(std::string("REFCOUNT_TEST")), true); 
  #endif
  simple_object_base::set_named_object(dest->as<object_map>(), "NDIM", new_simple_object<Z>(static_cast<Z>(NDIM)), true); 
  simple_object_base::set_named_object(dest->as<object_map>(), "parameters", get_parameters().object_pointer()->clone(), true);     
  simple_object_base::set_named_object(dest->as<object_map>(), "position", new_simple_object<std::vector<ntuple<R,NDIM> > >(position), true); 
  simple_object_base::set_named_object(dest->as<object_map>(), "velocity", new_simple_object<std::vector<ntuple<R,NDIM> > >(velocity), true);
  simple_object_base::set_named_object(dest->as<object_map>(), "radius", new_simple_object<std::vector<R> >(radius), true);
  #if 0
  // *** DEBUG *** this should already be in "parameters":
  simple_object_base::set_named_object(dest->as<object_map>(), "density", new_simple_object<std::vector<R> >(density), true);
  #endif
}

 
/*
 * next_event *always* is defined:
 * (note: this is per-cell < <zero-time jam event 1>, <zero-time jam event 2>,...,<finite-time soonest event> >,
 *    where there may be many jam events, but only one soonest finite-time event)
 */
template <class R, size_t NDIM>
typename system<R,NDIM>::event_list mjp_system<R,NDIM>::next_event(const cell* c)const
{
  const parameters &param(get_parameters());

  // possibility exists that there will be no events:
  //   - no particles in cell
  //   - cell contains stationary particles
  // additional possibility of no finite-time events, 
  //   however with some zero-time particle-jam events (i.e. no "next event"):
  
  bool valid_next_event(false); // true when any finite-time events actually exist
  
  const cell* cell0(c);
  event_list result;
  
  if (!cell0->particles().empty()){
    // result.back() will be next finite-time event(if such exists):
    result.push_back(new event(event::NULL_EVENT, inv(epsilon<R>()),NULL,NULL)); // position for next finite-time event (init next->time_ to max)
    event *next(result.back());    

    #if defined(_debug_print_)
    // *** DEBUG ***
    cout<<"scanning cell: "<<cell0->linear_index()<<endl;
    #endif    
    for(typename cell::particle_list_type::const_iterator itP0 = cell0->particles().begin(), itP0End = cell0->particles().end();
        itP0 != itP0End;
        ++itP0){
      #if !defined(NDEBUG)  
      const sphere* particle0_(dynamic_cast<const sphere*>(*itP0));
      assert(particle0_ != NULL);
      #else
      const sphere* particle0_(reinterpret_cast<const sphere*>(*itP0));
      #endif  

      // event parameters (modified by "cell_exit" and "collision" methods):
      R time_(zero<R>());
      size_t dimension_(0), face_(0);
      bool jam_(false);

      // update system per-particle statistics:

      #if defined(__USE_PTHREAD)
      if (pthread_mutex_lock(&(this->base_class::mutex_)))
        throw std::runtime_error("mjp_system<R,NDIM>::next_event: pthread_mutex_lock error return");
      #endif

      // update system per-particle statistics:
      update_statistics(particle0_);

      #if defined(__USE_PTHREAD)
      if (pthread_mutex_unlock(&(this->base_class::mutex_)))
        throw std::runtime_error("mjp_system<R,NDIM>::next_event: pthread_mutex_unlock error return");
      #endif
      
      // don't evolve beyond maximum allowed particle radius:
      if (radius_limit(particle0_, time_)){
        valid_next_event = true;

        if (time_ < next->time()){
          // cast cell* and particle* from const: event will eventually be used to modify the system state:
          // re-init to timestep only event:
          next->init(
            event::TIMESTEP_EVENT,
            time_, 
            const_cast<cell*>(cell0), 
            const_cast<sphere*>(particle0_));
        }          
      }
             
      // check cell-exit events:
      if (cell_exit(particle0_, cell0, time_, dimension_, face_)){
        valid_next_event = true;
        if (time_ < next->time()){
          // cast cell* and particle* from const: event will eventually be used to modify the system state:
          next->init(
            event::CELL_EXIT_EVENT,
            time_, 
            const_cast<cell*>(cell0), 
            const_cast<sphere*>(particle0_), 
            NULL, dimension_, face_);
        }
      }

      // check pair collision events (these will possibly include jammed particles):
      // ----------- OPTIMIZATION: look only at the half of the neighbors of each cell that haven't already been checked: ------------        
      size_t nn(0);
      for(typename cell::neighbor_list_type::const_iterator itN = cell0->positive_neighbors().begin(),
            itNEnd = cell0->positive_neighbors().end();
          itN != itNEnd;
          ++itN, ++nn){
        const cell *cell1(*itN); 
        #if defined(_debug_print_)
        // *** DEBUG ***
        cout<<"    scanning cell "<<cell0->linear_index()<<" neighbor "<<nn<<"(= cell "<<cell1->linear_index()<<")"<<endl;
        #endif 

        // note: cell1 now may be cell0 itself
        typename cell::particle_list_type::const_iterator itP1(cell1->particles().begin()),
          itP1End = cell1->particles().end();
        if (cell1 == cell0){
          // within the same cell, only check particles later in the particle-list:
          itP1 = std::find(cell1->particles().begin(), cell1->particles().end(), particle0_);
          assert(itP1 != cell1->particles().end());
          ++itP1;
        }
        for(  ;itP1 != itP1End;
               ++itP1){
          #if !defined(NDEBUG)  
          const sphere* particle1_(dynamic_cast<const sphere*>(*itP1));
          assert(particle1_ != NULL);
          #else
          const sphere* particle1_(reinterpret_cast<const sphere*>(*itP1));
          #endif  
          if (particle1_ != particle0_){
            if (collision(particle0_, particle1_, time_, jam_)){
              if (jam_){
                // save _all_ particle-jam events to the front of the list:
                // (at this phase, duplicate checking is _not_ performed (that's done later, during step))
                // (const_cast: see previous comment at cell-exit event check)
                result.push_front(
                  new event(event::JAM_EVENT, 
                            time_, 
                            const_cast<cell*>(cell0), 
                            const_cast<sphere*>(particle0_), const_cast<sphere*>(particle1_)));
              }              
              else
              if (time_ < next->time()){
                valid_next_event = true;
                typename event::event_kind event_kind_(event::COLLISION_EVENT);
                #if 1
                if (param.sticking_probability > zero<R>()){
                  // convert random fraction of COLLISION_EVENT into STICK_EVENT:
                  R test(random<R>());
                  if (test < param.sticking_probability)
                    event_kind_ = event::STICK_EVENT;
                }
                #endif
                next->init(
                  event_kind_,
                  time_, 
                  const_cast<cell*>(cell0), 
                  const_cast<sphere*>(particle0_), const_cast<sphere*>(particle1_));

                  #if 0
                  // *** DEBUG ***
                  // OVERLAP glitch:
                    sphere 
                      *p0(reinterpret_cast<sphere*>(particle0_->clone())), 
                      *p1(reinterpret_cast<sphere*>(wrap_particle1->clone()));
                    const_cast<mjp_system<R,NDIM>*>(this)->base_class::move(p0, time_);
                    const_cast<mjp_system<R,NDIM>*>(this)->base_class::move(p1,time_);
                    if (p0->overlap(p1)){
                      const state &s0_(p0->current_state()), &s1_(p1->current_state());
                      cout<<"next collision: imminent OVERLAP detected:\n"                          
                          <<" (dt = "<<(time_ - particle0_->current_state().time())/epsilon<R>()<<" eps): distance: "
                          <<dist(s0_.position(), s1_.position()) - (p0->radius() + p1->radius())<<endl; 
                    }
                    delete p0;
                    delete p1;
                  #endif

              }    
            }
          }
        }
      }
    }

    if (valid_next_event){
      // don't evolve beyond maximum allowed timestep:
      if ((base_class::t_max() > zero<R>()) && (next->time() > base_class::t_max())){
        // re-init to timestep only event:
        next->init(
          event::TIMESTEP_EVENT,
          base_class::t_max(), 
          next->event_cell(), 
          next->particle0());             
      }     
      #if defined(_debug_print_)
      // *** DEBUG ***
      cout<<"next: \n";
      next->write(cout);
      cout<<endl;
      #endif
    }
    else{
      // remove the position for "next-event"
      // (there may still be particle-jam events)
      result.remove_next_event(); 
    }
    
    #if defined(_debug_print_)
    // *** DEBUG ***
    size_t N_jam(result.data().size());
    if (!result.data().empty() && !result.back()->is_jam())
      --N_jam;
    if (N_jam > 0)
      cout<<"****** "<<N_jam<<" JAM EVENTS ******"<<endl;
    #endif  
  }
  #if defined(_debug_print_)
  else
    // *** DEBUG ***
    cout<<"cell contains NO particles"<<endl;;
  #endif  
  
  return result;
}


// EITHER  ALL of "jam_particles", "move_particles", and "process_event"
//   OR "step_" *must* be defined (in the latter case, the former may remain stubs, or not, as required)
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::jam_particles(event* E)
{
  bool status(true);
  const parameters &param __attribute_unused__ (get_parameters());
  
  assert(E->is_jam() && (E->particle0() != NULL) && (E->particle1() != NULL));
  
  #if defined(__USE_PTHREAD)
  std::pair<size_t, size_t> keys = base_class::lock_pair_(E->particle0(), E->particle1());
  #endif
  
  // unified treatment of jam and collision events:
  //   (except that jam-events processed at zero-time)
  collide(E);

  #if defined(__USE_PTHREAD)
  base_class::unlock_pair_(keys);  
  #endif

  #if defined(event_history)
  // *** DEBUG *** ----------------------------------------------------------------------------
  #if defined(__USE_PTHREAD)
  if (pthread_mutex_lock(&(this->base_class::mutex_)))
    throw std::runtime_error("mjp_system<R,NDIM>::jam_particles: pthread_mutex_lock error return");
  #endif
  
  base_class::event_history_.push_back(E->clone());
  
  #if defined(__USE_PTHREAD)
  if (pthread_mutex_unlock(&(this->base_class::mutex_)))
    throw std::runtime_error("mjp_system<R,NDIM>::jam_particles: pthread_mutex_unlock error return");
  #endif
  #endif // ------------------------------------------------------------------------------------
  
  return status;   
}

#if 0 // -------- not special: moved up to base_class: ---------------------------------------------------------
/*
 * implementation note:
 *  At present, "move_particles" should NOT require a mutex lock on the particle:
 *    as this method is threaded by _cell_, and particles are not shared between cells,
 *    there is no _valid_ situation where more than one thread of this method would access a given
 *    particle (read _or_ write) at the same time. 
 */
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::move_particles(cell* c, event* E)
{
  bool status(true);
  const parameters &param __attribute_unused__ (get_parameters());  
  
  // move all particles in cell to kinematic coordinates of event:
  for(typename cell::particle_list_type::iterator itP = c->particles().begin(), itPEnd = c->particles().end();
      itP != itPEnd;
      ++itP){
    // implement extrinsic _and_ intrinsic move (i.e. call system::move, rather than particle::move):  
    base_class::move(*itP, E->time());    
  }
  
  return status;
}
#endif

template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::process_event(event* E)
{
  bool status(true);
  const parameters &param __attribute_unused__ (get_parameters());
  
  switch (E->kind()){
    case event::TIMESTEP_EVENT:
    // mark the state: TIMESTEP_EVENT is a special type of MOVE_EVENT:
    E->particle0()->current_state().set_event_type(event::TIMESTEP_EVENT);
    break;
    
    case event::CELL_EXIT_EVENT:
    #if 0
    E->event_cell()->transfer_particle(E->particle0(), E->event_cell()->face_neighbor(E->dimension(), E->face()));
    #else
    base_class::exit_cell(E->particle0(), E->dimension(), E->face());
    #endif
    break;
    
    case event::COLLISION_EVENT:
    // unified treatment of jam and collision events:
    //   (except that collision events processed after time-step 
    //    and associated move to kinematic coordinates of event)
    collide(E);
    break;
     
    case event::STICK_EVENT:
    (E->particle0())->stick(E);
    break;
    
    case event::JAM_EVENT:  
    default:
    throw std::runtime_error("mjp_system<R,NDIM>::process_event: unrecognized or non-finite-time event");
    break;
  }
  
  #if defined(event_history)
  // *** DEBUG *** ----------------------------------------------------------------------------
  #if defined(__USE_PTHREAD)
  if (pthread_mutex_lock(&(this->base_class::mutex_)))
    throw std::runtime_error("mjp_system<R,NDIM>::process_event: pthread_mutex_lock error return");
  #endif
  
  base_class::event_history_.push_back(E->clone());
  
  #if defined(__USE_PTHREAD)
  if (pthread_mutex_unlock(&(this->base_class::mutex_)))
    throw std::runtime_error("mjp_system<R,NDIM>::process_event: pthread_mutex_unlock error return");
  #endif
  #endif // ------------------------------------------------------------------------------------
  
  
  return status;   
}

/* 
 * hooks for initialization and update of per-particle system statistics:
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::init_statistics(void)const
{ system_stats_.clear(); }

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::update_statistics(const particle* p)const
{ system_stats_.update(p); }


/*
 * hook for application of system temperature control:
 * (note: this is a per-step hook, rather than a per-particle hook)
 */
template <class R, size_t NDIM>
R mjp_system<R,NDIM>::thermostat_v_factor(void)const
{
  R v_factor(one<R>());
  const parameters &param(get_parameters());  
  const R system_T(system_stats().T_mean());
  if (param.T_max > zero<R>()){
    // cooling only (allow heat-up):
    if (system_T > param.T_max)
      v_factor = sqrt(param.T_max/system_T);
  }
  return v_factor;
}

// accumulate statistics during "move_particles" for use by "complete":
template <class R, size_t NDIM>
inline const typename mjp_system<R,NDIM>::statistics& mjp_system<R,NDIM>::system_stats(void)const
{ return system_stats_; }

template <class R, size_t NDIM>
R mjp_system<R,NDIM>::packing_fraction(void)const
{ return system_stats_.hypervolume()/pow_n(get_parameters().L, static_cast<long>(NDIM)); }

template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::complete(void)const
{
  bool test(base_class::complete());
  
  if (!test){
    const parameters &param(get_parameters());    
    const statistics& stats(system_stats());

    // check _error_ conditions (or rather warnings) prior to by-request completion tests:
    if (stats.r_max() * integer<R>(2) > base_class::cell_edge_length()){
      // this could actually be an exception throw...
      clog(logger::LOG_QUIET)
        <<"WARNING: maximum particle diameter "<<stats.r_max() * integer<R>(2)
        <<" exceeds cell edge dimension "<<base_class::cell_edge_length()<<endl;
      test = true;
    }
    else  
    if ((base_class::t_max() > zero<R>()) && (stats.t_max() >= base_class::t_max())){
      clog(logger::LOG_VERBOSE)<<"mjp_system<R,NDIM>::complete: maximum particle time "<<stats.t_max()
          <<" exceeds maximum time allowed for current evolution phase "<<base_class::t_max()<<endl;
      test = true;
    }
    else
    if ((param.r_max > zero<R>()) && (stats.r_max() >= param.r_max)){
      clog(logger::LOG_VERBOSE)<<"mjp_system<R,NDIM>::complete: maximum particle radius "<<stats.r_max()
          <<" exceeds maximum radius allowed by parameters "<<param.r_max<<endl;
      test = true;
    } 
    else
    if (stats.N_frozen() == param.N_particle){
      clog(logger::LOG_VERBOSE)<<"mjp_system<R,NDIM>::complete: all particles frozen"<<endl;    
      test = true;
    }
  }  
    
  #if defined(_debug_print_)
  cout<<"--- state summary (time: "<<stats.t_max()<<"): ---\n";
  stats.write(cout);
  cout<<endl;
  #endif
    
  return test;
}


// ================== apply related methods: ======================

/**
 * @brief Check for consistency between a given argument and the present value of the functor parameters.
 * This method does not check for state consistency itself (e.g. lack of particle overlap); this is assumed.
 */
template <class R, size_t NDIM>
void mjp_system<R,NDIM>::arg_valid_check(const object_map& arg)const
{
  using python_util::simple_object_base;
  const parameters& param(get_parameters());
  
  base_class::arg_valid_check(arg);
  if (simple_object_base::has_named_parm(arg, "radius")){
    // allow empty "radius" list as placeholder parameter:
    const simple_object_base *pobj = simple_object_base::get_named_object(arg, "radius");
    if (!pobj->is_empty() &&
        (pobj->size() != param.N_particle))
      throw std::runtime_error("mjp_system<R,NDIM>::arg_valid_check: number of entries in radius list doesn't match number of particles");
  }
}

/// wrapper "apply" method:
template <class R, size_t NDIM>
inline bool mjp_system<R,NDIM>::apply(const python_util::generic_object<C,R,Z>& arg, python_util::generic_object<C,R,Z>& val)const
{ 
  /*
   * Implementation note:
   *   this whole method is a Kluge:  when parallelFunctor<...>::workspace::result is created using a RANGE type of "generic_object"
   *   it does not have the information to initialize the object-pointer to *simple_object<object_map>.
   * This seems to suggest that "generic_object" needs to be re-worked so that it initializes appropriately.
   */
  if (NULL == val.ptr())
    val.clone_from(new python_util::simple_object<object_map>(), true); // transfer object ownership

  return apply(arg.ptr(), val.ptr()); 
}

/**
 * @brief generic apply method
 * @param[in] arg: 
 *    - @b position[opt]: list of particle positions (cartesian coordinates)
 *    - @b velocity[opt]: when corresponding positions are supplied, a list of particle vector velocities
 *    - @b radius[opt]: when corresponding positions are supplied, a list of particle radii
 *    .
 * When particle positions are not supplied, values for number of particles "N_particle", initial temperature "T_0", radial growth rate "dr"
 *   will be taken from parameters to build the initial system state prior to evolution.
 * @param[out] val:
 *    - @b position: evolved particle positions
 *    - @b velocity: evolved particle vector velocities
 *    - @b radius:   evolved particle radii
 *    .    
 */
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::apply(const python_util::simple_object_base *arg, python_util::simple_object_base *val)const
{
  bool status(true);
  int log_level_save(statusUtil::clog.output_level());
    
  try{
    /*
     * Implementation note: temporary work-around: for use as functor, methods which modify the system-state
     *   will be called from "const_cast<mjp_system<R,NDIM>*>(this)->" rather than redeclared as "const" with "mutable" system attributes.
     * In the longer term, this should be analyzed and remedied as appropriate (e.g. by using a workspace structure).
     */

    using python_util::simple_object_base;
    const parameters& param(get_parameters());
    
    // allow any previously-specified increase in logger output_level to override the parameters value:
    int log_level(0);
    conv(log_level, param.template get_named_parm<Z>("log_level"));  
    if (log_level > log_level_save)
      statusUtil::clog.output_level(log_level);
      
    // "arg_valid_check" called by "init_state":
    const_cast<mjp_system<R,NDIM>*>(this)->init_state(arg);

    if (param.has_named_parm("target_distribution")){
      base_class::apply_directed_();
    }
    else{
      size_t nstep;
      for(nstep = 0; nstep < param.NSTEP; ++nstep)
        if (!const_cast<mjp_system<R,NDIM>*>(this)->base_class::step()) break;

      if (nstep == param.NSTEP)
        clog(logger::LOG_VERBOSE)<<"mjp_system<R,NDIM>::apply: "<<param.NSTEP<<" event steps completed"<<endl;
    }

    // transfer system state to return value:
    extract_state(val);
  }
  catch(std::string& msg){
    base_class::setStatusString(msg);
    status = false;
  }
    
  // restore incoming log output-level:
  statusUtil::clog.output_level(log_level_save);
  
  return status;
}

// =======================================================


template <class R, size_t NDIM>
void mjp_system<R,NDIM>::write(std::ostream& os)const
{
  for(typename cell_array::const_iterator itC = base_class::cells_.begin(), itCEnd = base_class::cells_.end();
      itC != itCEnd;
      ++itC){
    const cell *cell_(*itC);  
    os<<"--- cell: "<<cell_->linear_index()<<"(";
    cell_->indices().write(os);
    os<<"): ---\n";
    
    for (typename cell::particle_list_type::const_iterator itP = cell_->particles().begin(), itPEnd = cell_->particles().end();
         itP != itPEnd;
         ++itP){
      const particle *particle_(*itP);
      os<<"particle: "<<particle_->linear_index()<<(particle_->frozen()?": *** FROZEN ***:":":")<<"\n";
      particle_->current_state().write(os);
      os<<"\n";        
    }        
  }    
}

template <class R, size_t NDIM>
void mjp_system<R,NDIM>::debug_print(void)const
{ 
  write(std::cout);
  std::cout<<std::endl;
}

template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::writeBinary(commUtil::abstractCommHandle *fp)const
{
  bool status(true);
  // note: parameters written during base_class::writeBinary:
  status = (status && base_class::writeBinary(fp));
  return status;
}
    
template <class R, size_t NDIM>
bool mjp_system<R,NDIM>::readBinary(commUtil::abstractCommHandle *fp)
{
  bool status(true);
  status = (status && base_class::readBinary(fp));
  return status;
}
    

template <class R, size_t NDIM>
size_t mjp_system<R,NDIM>::binarySize(void)const
{
  size_t val(0);
  val += base_class::binarySize();
  return val;
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>& mjp_system<R,NDIM>::operator=(const mjp_system& other)
{
  base_class::operator=(static_cast<base_class&>(other));
  return *this;
}


template <class R, size_t NDIM>
mjp_system<R,NDIM>::~mjp_system(void)
{ }

template <class R, size_t NDIM>
mjp_system<R,NDIM>::mjp_system(bool initialize)
  : system<R,NDIM>(false)
{
  if (initialize){
    // establish correct virtual pointer for base_class::pparam_:
    set_parameters(parameters());
  }    
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::mjp_system(const mjp_system& other)
  : system<R,NDIM>(false)
{
  operator=(other);
}

template <class R, size_t NDIM>
mjp_system<R,NDIM>::mjp_system(size_t N_particle, const R& L, const R& v, const R& dr)
  : system<R,NDIM>(false)
{
  // full re-initialization sequence:
  parameters param;
  param.N_particle = N_particle;
  param.L = L;
  param.v = v;
  param.dr = dr;
  set_parameters(param);

  init_state();
}  


    
} // namespace particle_packing

#endif // __particle_packing_template__h
