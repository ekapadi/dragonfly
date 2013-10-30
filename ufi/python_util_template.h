#if !defined(__python_util_template__h)
#define __python_util_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util_template.h,v $

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


namespace python_util{



/**
 * @brief template-function wrapper for external function to be called from a python extension module.
 *
 * @param[in] self: functor-state, interpreted as its parameters
 * @param[in] args: arguments
 * @return  functor::apply(args) with specified parameters
 *
 * @tparam FUNCTOR the functor class
 * @tparam RVAL the return-value class
 * @tparam ARG  the argument class
 * @tparam PARAM the parameters class
 *
 * "self", "args", and return-value are extracted-from and inserted-into the python objects by passing through
 *    an intermediate conversion to "simple_object_base*"; this means that any argument and return-value types
 *    that have defined "new_simple_object" and "extract" methods (used by "simple_object_base*") are supported.
 *
 * FUNCTOR shall implement the following methods:
 *   - "bool setParameters(const PARAM&)", 
 *   - "bool apply(const ARG&, RVAL&)", 
 *   - "std::string getStatusString(void)const"
 *   .
 * Error return: python-style (NULL) with python error-state set. 
 */
template <class FUNCTOR, class RVAL, class ARG, class PARAM>
PyObject*  wrapped_external_functor(PyObject *self, PyObject *args)
{
  bool status(true);
  PyObject* python_rval(NULL);

  #if 0
  // *** DEBUG ***
  bool debug_wait(true);
  while (debug_wait);
  #endif
  
  try{  

  /*
   * Implementation note (in order to avoid confusion):
   *   the following methods have very different signatures
   *   (eventually these should be unified).
   * For clarity the explicit namespace qualification will be left in,
   *   for the PyObject* methods.
   *
   * WARNING: the default extract/insert methods for class simple_object_base must be
   *   initialized before this function is called (e.g. using): 
   *     simple_object_base::set_default_extractor(extract<C,R,Z>);
   *     simple_object_base::set_default_inserter(insert<C,R,Z>); 
   *  
   */
  #if 0
  // insert/extract to and from simple_object_base*:
  using python_util::new_simple_object;
  using python_util::extract;
  
  // insert/extract to and from PyObject*:
  using implementation_module::insert;
  using implementation_module::extract;
  #endif
  
  FUNCTOR functor;
  PARAM param;
  ARG arg;
  RVAL rval;
  
  const simple_object_base* pconst_obj(NULL);

  // extract the parameters:
  pconst_obj = implementation_module::template extract<simple_object_base*>(self);
  extract(param, pconst_obj);
  delete pconst_obj;
  
  // extract the argument:
  pconst_obj = implementation_module::template extract<simple_object_base*>(args);
  extract(arg, pconst_obj);
  delete pconst_obj;
  
  // functor must allocate "rval" as required, there is no "correct" (i.e. fully generic) way to do it here:
  if (!functor.setParameters(param) || !functor.apply(arg, rval)){
    // do not have the information here to set exception-type to anything besides runtime-error:
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(functor.getStatusString().c_str()));
    status = false;
  }
  
  if (status){
    simple_object_base* pobj = new_simple_object(rval);
    
    // use the read/write version of the "insert" function (the ownership of the simple_object data will be passed to python):
    python_rval = implementation_module::template insert<simple_object_base* const>(pobj);
    delete pobj;    
  }

  }
  catch(std::exception& X){
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(X.what()));
    return NULL;
  }
  catch(std::string& X){
    PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(X.c_str()));
    return NULL;
  }
  catch(...){
    PyErr_SetString(PyExc_RuntimeError, "unknown EXCEPTION type");
    return NULL;
  }
  
  return python_rval;
}


// ------------------------------ simple_object_base: -----------------------------------------

template <class T>
inline const T& simple_object_base::as_(void)const
{ 
  #if 0
  // strict aliasing:
  const generic_ptr<T>* gptr(dynamic_cast<const generic_ptr<T>*>(pg_));
  assert(is<T>()); // *** DEBUG ***  this is OK, when type_info comparison occurs "by-name"
  assert(NULL != gptr);
  return *(gptr->ptr);
  #else
  // strict aliasing:  *** DEBUG *** 
  // ------- duplicate virtual-tables (cross-module RTTI problems...): just use "reinterpret_cast"
  const generic_ptr<T>* gptr(reinterpret_cast<const generic_ptr<T>*>(pg_));
  assert((NULL != gptr->ptr) && (size() == 1)); /* at present NULL ptr represents empty array object: this may not be a good idea. */
  return *(gptr->ptr);  
  #endif
}

template <class T>
inline T& simple_object_base::as_(void)
{ 
  #if 0
  // strict aliasing:
  generic_ptr<T>* gptr(dynamic_cast<generic_ptr<T>*>(pg_));
  assert(NULL != gptr);
  return *(gptr->ptr);
  #else
  // strict aliasing:  *** DEBUG *** 
  // ------- duplicate virtual-tables (cross-module RTTI problems...): just use "reinterpret_cast"
  generic_ptr<T>* gptr(reinterpret_cast<generic_ptr<T>*>(pg_));
  return *(gptr->ptr);  
  #endif
}

template <class T>
inline T*& simple_object_base::ptr_(void)
{ 
  #if 0
  // strict aliasing:
  generic_ptr<T>* gptr(dynamic_cast<generic_ptr<T>*>(pg_));
  assert(NULL != gptr);
  return gptr->ptr;
  #else
  // strict aliasing:  *** DEBUG *** 
  // ------- duplicate virtual-tables (cross-module RTTI problems...): just use "reinterpret_cast"
  generic_ptr<T>* gptr(reinterpret_cast<generic_ptr<T>*>(pg_));
  return gptr->ptr;  
  #endif
}

template <class T>
inline const T* simple_object_base::ptr_(void)const
{ 
  #if 0
  // strict aliasing:
  const generic_ptr<T>* gptr(dynamic_cast<const generic_ptr<T>*>(pg_));
  assert(NULL != gptr);
  return gptr->ptr;
  #else
  // strict aliasing:  *** DEBUG *** 
  // ------- duplicate virtual-tables (cross-module RTTI problems...): just use "reinterpret_cast"  
  const generic_ptr<T>* gptr(reinterpret_cast<const generic_ptr<T>*>(pg_));
  return gptr->ptr;  
  #endif
}



template <class T>
/* inline */ bool simple_object_base::is(void)const
{
  #if 1
  // *** DEBUG *** RTTI cross-module problems:
  return (0 == strcmp(typeid(T).name(), type().name()));
  #else
  return typeid(T) == type(); // doesn't work  with "-O2" without extensive header-file / template-body separation
                              // (and I'm not sure this is worthwhile, at this time...)
  #endif
}

// --------- prevent multiple instantiations: ------------------
#if defined(__USE_MERE)
extern template bool simple_object_base::is<mere::C>(void)const;
extern template bool simple_object_base::is<mere::R>(void)const;
extern template bool simple_object_base::is<mere::Z>(void)const;
#endif
extern template bool simple_object_base::is<std::complex<double> >(void)const;
extern template bool simple_object_base::is<double>(void)const;
extern template bool simple_object_base::is<long>(void)const;

extern template bool simple_object_base::is<bool>(void)const;
extern template bool simple_object_base::is<std::string>(void)const;

extern template bool simple_object_base::is<simple_object_base::object_map>(void)const;
extern template bool simple_object_base::is<simple_object_base::object_list>(void)const;
// -------------------------------------------------------------



template <class T>
inline const T& simple_object_base::as(void)const throw(std::string)
{
  if (!is<T>())
    throw std::string("simple_object_base::as: object does not have requested type: ") + typeid(T).name() + std::string(", defined as: ") + type().name();
  return as_<T>();  
}

template <class T>
inline T& simple_object_base::as(void) throw(std::string)
{
  if (!is<T>())
    throw std::string("simple_object_base::as: object does not have requested type: ") + typeid(T).name() + std::string(", defined as: ") + type().name();
  return as_<T>(); 
}

template <class T>
inline const T* simple_object_base::ptr(void)const throw(std::string)
{
  if (!is<T>())
    throw std::string("simple_object_base::ptr: object does not have requested type: ") + typeid(T).name() + std::string(", defined as: ") + type().name();
  return ptr_<T>();  
}

template <class T>
inline T* simple_object_base::ptr(void) throw(std::string)
{
  if (!is<T>())
    throw std::string("simple_object_base::ptr: object does not have requested type: ") + typeid(T).name() + std::string(", defined as: ") + type().name();
  if (size()>1 && !own_data())
    throw std::string("simple_object_base::ptr: non-const pointer access to RANK>0 number objects requires data ownership.");
  return ptr_<T>(); 
}    
// ----------------------- object_map utility methods: ----------------------------------------


// single_entity types U (see simple_object_traits<T>):
// (note: use of this form with list_entity or object_entity types throws exception)
template <class U>
const U& simple_object_base::get_named_parm(const object_map& map, const std::string& key) throw(std::string)
{
  const U* p(NULL);
  // note on form: dispatcher-type partial specializations only allowed at namespace scope.
  if (simple_object_traits<U>::is_single()){
    object_map::const_iterator itP = map.find(key);
    if (itP != map.end())
      p = &((*itP).second->as<U>());
    else  
      throw std::string("get_named_parm<U>: no entry found for specified key: ") + key;
  }
  else
    throw std::string("get_named_parm<U>: no reference available on target type \"U\"; \n"
                      "  list_entity or map_entity requires conversion");
  return *p;
}

template <class U>
void simple_object_base::set_named_parm(object_map& map, const std::string& key, const U& val) throw(std::string)
{
  // Notes on form: 
  //   -- dispatcher-type partial specializations only allowed at namespace scope;
  //   -- one cannot assume that a map-entry pointer is initialized to NULL;
  //   -- no attempt is made to re-use allocation.

  std::pair<object_map::iterator, bool> test = map.insert(std::pair<std::string, simple_object_base*>(key,NULL));
  object_map::iterator &itP(test.first);
  if (!test.second && ((*itP).second != NULL))
    delete (*itP).second;
     
  (*itP).second = new_simple_object(val);   
}


// container types U (see simple_object_traits<T>):
template <class U>
void simple_object_base::get_named_parm(const object_map& map, const std::string& key, U& val) throw(std::string)
{
  using python_util::extract;

  // note on form: dispatcher-type partial specializations only allowed at namespace scope.
  object_map::const_iterator itP = map.find(key);
  const simple_object_base *psrc(NULL);
  if (itP != map.end())
    psrc = (*itP).second;
  else  
    throw std::string("get_named_parm<U>: no entry found for specified key: ") + key;
  
  extract(val, psrc); 
}


// ----------------------- I/O methods (primary justification for this class) -----------------

// initialize a simple_object_base pointer from a stream representation as python syntax 
// (Allows arbitrary python code, including preamble and postscript. Note: name of python object is "object".):
template <class C, class R, class Z>
simple_object_base* simple_object_base::read_repn_virtual(const std::string& src, const std::string* ppreamble, const std::string* ppostscript)
{
  simple_object_base* dest(NULL);
  implementation_module::extract_simple_object<C,R,Z>(dest, src, ppreamble, ppostscript);
  return dest;
}

// This implementation is a bit odd, but the immediate alternative was passing filename strings:
template <class C, class R, class Z>
simple_object_base* simple_object_base::read_repn_virtual(std::istream& src, 
                                                          std::istream* ppreamble, 
                                                          std::istream* ppostscript)
{
  using scanner_util::load_from_stream;
  simple_object_base *dest(NULL);
  
  // transfer src std::istream to std::string:
  std::string src_, preamble_, postscript_;
  load_from_stream(src_, src);
  if (ppreamble != NULL)
    load_from_stream(preamble_, *ppreamble);
  if (ppostscript != NULL)
    load_from_stream(postscript_, *ppostscript);

  dest = read_repn_virtual<C,R,Z>(src_, &preamble_, &postscript_);
  return dest;
}

// output the object to a stream representation (as python syntax):
template <class C, class R, class Z>
void simple_object_base::write_repn(std::ostream &dest)const    
{
  std::string dest_;
  implementation_module::insert_simple_object<C,R,Z>(dest_,this);
  dest<<dest_;
}
    
// output the object to "pickle" binary representation:
template <class C, class R, class Z>
void simple_object_base::write_pickle(const std::string& filename)const
{
  implementation_module::write_pickle<C,R,Z>(filename,this);
}    

// construct object from "pickle" binary representation:
template <class C, class R, class Z>
simple_object_base* simple_object_base::read_pickle_virtual(const std::string& filename)
{
  simple_object_base *dest(NULL);
  implementation_module::read_pickle<C,R,Z>(dest, filename);
  return dest;
}   

// output the object to a stream representation (in syntax appropriate for display to end-user):
// (note: for possible debugging use this method should not require the embedded python interpreter.)
template <class C, class R, class Z>
void simple_object_base::write(std::ostream &dest)const
{
  // closely follows implementation from: python_util::implementation_module::insert_simple_object:
  const simple_object_base* src(this);
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  

  //  PyObject *src_(src.ptr());
  if (src->is<bool>())
    dest<<(src->as<bool>()?"true":"false");
  else
  if (src->is<C>())
    dest<<src->as<C>();
  else
  if (src->is<R>())
    dest<<src->as<R>();
  else
  if (src->is<Z>())
    dest<<src->as<Z>();
  else
  if (src->is<std::string>())
    dest<<"\""<<src->as<std::string>()<<"\"";
  else
  if (src->is<LIST>()){
    const LIST &l(src->as<LIST>());
    dest<<"(";
    bool first_element(true);
    for(LIST::const_iterator itl = l.begin(), itlEnd = l.end();
        itl != itlEnd;
        ++itl){
      if (!first_element)
        dest<<", ";
      else
        first_element = false;    
      (*itl)->write<C,R,Z>(dest);
    }
    dest<<")\n"; // end-of-line after LIST.    
  }
  else
  if (src->is<MAP>()){
    const MAP &m(src->as<MAP>());
    dest<<"{";
    bool first_element(true);
    for(MAP::const_iterator itm = m.begin(), itmEnd = m.end();
        itm != itmEnd;
        ++itm){
      if (!first_element)
        dest<<", ";
      else
        first_element = false;   
      dest<<"\""<<(*itm).first<<"\":";
      (*itm).second->write<C,R,Z>(dest);
    }
    dest<<"}\n"; // end-of-line after MAP.   
  }
  else  
    throw std::string("simple_object_base::write: object type not in { <complex>, <real>, <integer>, std::string, object_list, object_map }");
}        

// --------------------------------------------------------------------------------------------

// ---------------------------------------- dispatcher classes: ----------------------------------

// type-specific extraction from simple_object_base* :
template <class T>
void extract(T& t, const simple_object_base* pobject) throw(std::string)
{
  dispatch_extract(t, pobject, typename simple_object_traits<T>::entity_type());
}

template <class T, class ET>
inline void dispatch_extract(T& t, const simple_object_base *pobject, ET) throw(std::string)
{ 
  t.extract(pobject);
}

template <class T>
inline void dispatch_extract(T& t, const simple_object_base *pobject, single_entity) throw(std::string)
{ extract_single(t, pobject); }

template <class U>
inline void dispatch_extract(U& u, const simple_object_base *pobject, list_entity) throw(std::string)
{ 
  // using "typename gmm::linalg_traits<U>::value_type" is also a possibility here, but at this point it isn't known
  //   whether or not U is a container of numbers:
  extract_list(u, pobject, 
    typename simple_object_traits<
      typename linalg::tensor_traits<typename U::value_type>::scalar_type>::number_type()); 
}

template <class M>
inline void dispatch_extract(M& m, const simple_object_base *pobject, map_entity) throw(std::string)
{ extract_map(m, pobject); }

// extract single_entity types with number promotion as T is C <-- C, R, Z;  T is R <-- R, Z;  T is Z <-- Z
template <class T>
inline void extract_single(T& dest, const simple_object_base *pobject) throw(std::string)
{ dispatch_extract_single(dest, pobject, typename simple_object_traits<T>::number_type()); }

template <class T, class NT>
inline void dispatch_extract_single(T& dest, const simple_object_base *pobject, NT) throw(std::string)
{
  if (pobject->size() > 1)
    throw std::runtime_error("extract: SCALAR number cannot be extracted from RANK > 0 number object"); 
  extract_number(dest, pobject); 
}

template <class T>
inline void dispatch_extract_single(T& dest, const simple_object_base *pobject, abstract_null_entity) throw(std::string)
{ dest = pobject->as<T>(); }

// default: extract container from a RANK=1 number object: "NT" is simple_object_traits<typename U::value_type>::number_type
template <class U, class NT>
inline void extract_list(U& dest, const simple_object_base *pobject, NT) throw(std::string)
{
  // Implementation note: _allow_ lists of numbers to be stored in dense RANK=1 format, 
  //   or alternatively as object lists:
  
  if (!pobject->is<simple_object_base::object_list>()){
    // here gmm methods could also be used: dest is known to be a container of number-type objects (possibly with RANK>0).

    // the number-type to extract:
    typedef typename linalg::tensor_traits<typename U::value_type>::scalar_type S;
    
    size_t N_scalar(pobject->size());
    if (N_scalar % linalg::tensor_traits<typename U::value_type>::N_components != 0)
      throw std::runtime_error("python_util::extract_list: RANK > 0 destination:\n"
        "  number of components in destination object incompatible with number of elements in simple_object");
    size_t len(N_scalar / linalg::tensor_traits<typename U::value_type>::N_components);
    dest.resize(len);
    
    std::copy(reinterpret_cast<const typename U::value_type*>(pobject->ptr<S>()), 
              reinterpret_cast<const typename U::value_type*>(pobject->ptr<S>()) + len, dest.begin());
  }
  else{
    // _allow_ lists of numbers to be expressed alternatively as object lists:
    extract_list(dest, pobject, abstract_null_entity());
  }  
}

// special case: not a number-type: extract container from an object-list:
template <class U>
inline void extract_list(U& dest, const simple_object_base *pobject, abstract_null_entity) throw(std::string)
{
  typedef simple_object_base::object_list object_list;
  const object_list &src(pobject->template as<object_list>());
  
  // assume STL access methods (source is a generic list, gmm methods don't really make sense here):
  dest.clear();
  dest.resize(src.size());

  object_list::const_iterator itSrc = src.begin();
  for(typename U::iterator itDest = dest.begin(), itDestEnd = dest.end();
      itDest != itDestEnd;
      ++itDest, ++itSrc)
    extract(*itDest, *itSrc); 
}

template <class M>
inline void extract_map(M& dest, const simple_object_base *pobject) throw(std::string)
{
  typedef simple_object_base::object_map object_map;
  const object_map &src(pobject->template as<object_map>());
  // assume STL-type container access methods:
  dest.clear();
  for(typename object_map::const_iterator itSrc = src.begin(), itSrcEnd = src.end();
      itSrc != itSrcEnd;
      ++itSrc)
    extract(dest[(*itSrc).first], (*itSrc).second); 
}

// -------------------- specializations of extract _bypassing_ dispatcher: -------------------------------

// python object, returns new reference:
template < >
inline void extract<PyObject*>(PyObject*& dest, const simple_object_base* pobject) throw(std::string)
{
  // nomenclature: _insert_ into the python-object, _extract_ from the simple_object_base* (what this method is doing):
  // use the default method for the class simple_object_base:
  dest = simple_object_base::insert(pobject);
}

// simple_object_base* itself (with clone()):
template < >
inline void extract<simple_object_base*>(simple_object_base*& dest, const simple_object_base* pobject) throw(std::string)
{
  assert(NULL == dest);
  dest = pobject->clone();
}

// generic_object:
template <class C, class R, class Z>
inline void extract(generic_object<C,R,Z>& dest, const simple_object_base* pobject) throw(std::string)
{
  dest = pobject; // generic_object::operator=(const simple_object_base* other) performs clone
}

// extensible_parameters_base:
template <class T>
inline void extract(extensible_parameters_base<T>& dest, const simple_object_base* pobject) throw(std::string)
{
  dest.clone_from(pobject);
}

// options_map<T>:
template <class T>
inline void extract(options_map<T>& dest, const simple_object_base* pobject) throw(std::string)
{
  dest.clone_from(pobject);
}

// -----------------------------------------------------------------------------------------------------

// generic allocation of simple_object_base*
template <class T>
simple_object_base* new_simple_object(const T& t) throw(std::string)
{ return new_object_dispatch(t, typename simple_object_traits<T>::entity_type()); }

// ------------- specializations of "new_simple_object" _bypassing_ dispatcher : --------------------------------

// simple_object_base* itself (using clone()):
inline simple_object_base* new_simple_object(const simple_object_base* ptr) throw(std::string)
{ return ptr->clone(); }

// generic_object (using clone()):
template <class C, class R, class Z>
inline simple_object_base* new_simple_object(const generic_object<C,R,Z>& t) throw(std::string)
{ return t.ptr()->clone(); }

// extensible_parameters_base (using clone()):
template <class T>
inline simple_object_base* new_simple_object(const extensible_parameters_base<T>& t) throw(std::string)
{ return t.object_pointer()->clone(); }

// options_map<T> (using clone()):
template <class T>
inline simple_object_base* new_simple_object(const options_map<T>& t) throw(std::string)
{ return t.object_pointer()->clone(); }

// --------------------------------------------------------------------------------------------------------------

template <class T, class ET>
inline simple_object_base* new_object_dispatch(const T& t, ET) throw(std::string)
{ throw std::string("new_object_dispatch: unrecognized simple_object_traits<T>::entity_type"); }

template <class T>
simple_object_base* new_object_dispatch(const T& t, single_entity) throw(std::string)
{ return new_single(t); }

template <class U>
simple_object_base* new_object_dispatch(const U& u, list_entity) throw(std::string)
{ 
  // using "typename gmm::linalg_traits<U>::value_type" is also a possibility here, but at this point it isn't
  //   known whether or not U is a container of numbers: 
  return new_list(u, typename simple_object_traits<
                       typename linalg::tensor_traits<typename U::value_type>::scalar_type>::number_type()); 
}

template <class M>
simple_object_base* new_object_dispatch(const M& m, map_entity) throw(std::string)
{ return new_map(m); }

template <class T>
inline simple_object_base* new_single(const T& t) throw(std::string)
{ return new simple_object<T>(t); }

template <class U, class NT>
inline simple_object_base* new_list(const U& src, NT) throw(std::string)
{ 
  // at this point gmm::linalg_traits could also be used: U is known to be a container of number-type (possibly with RANK>0):
  
  // the number type for the simple_object:
  typedef typename linalg::tensor_traits<typename U::value_type>::scalar_type S;
  
  // create a RANK=1 simple_object, of number-type:
  size_t 
    len(src.size()), 
    N_scalar(len * linalg::tensor_traits<typename U::value_type>::N_components);    
  simple_object_base *dest = new simple_object<S>(NULL, N_scalar, true, false);
  
  std::copy(src.begin(), src.end(), reinterpret_cast<typename U::value_type*>(dest->ptr<S>()));
  
  return dest; 
}

template <class U>
inline simple_object_base* new_list(const U& src, abstract_null_entity) throw(std::string)
{ 
  // not a number-type: create an object-list:
  return new simple_object<simple_object_base::object_list>(src); 
}

template <class M>
inline simple_object_base* new_map(const M& src) throw(std::string)
{ return new simple_object<simple_object_base::object_map>(src); }



// promotion to more-general number classes:
template <class T>
inline void extract_number(T& t, const simple_object_base* pobject) throw(std::string)
{ dispatch_extract_number(t, pobject, typename simple_object_traits<T>::number_type()); }

template <class T, class NT>
inline void dispatch_extract_number(T& t, const simple_object_base* pobject, NT) throw(std::string)
{ throw std::string("python_util::extract_number: unrecognized extraction type"); }

template <class T>
inline void dispatch_extract_number(T& t, const simple_object_base* pobject, complex_number_entity) throw(std::string)
{ extract_complex_number(t, pobject); }

template <class T>
inline void dispatch_extract_number(T& t, const simple_object_base* pobject, real_number_entity) throw(std::string)
{ extract_real_number(t, pobject); }

template <class T>
inline void dispatch_extract_number(T& t, const simple_object_base* pobject, integral_number_entity) throw(std::string)
{ extract_integral_number(t, pobject); }

template <class T>
inline void extract_complex_number(T& t, const simple_object_base* pobject) throw(std::string)
{
  if (pobject->template is<typename simple_object_traits<T>::C>())
    t = pobject->template as<typename simple_object_traits<T>::C>();
  else
  if (pobject->template is<typename simple_object_traits<T>::R>())
    t = pobject->template as<typename simple_object_traits<T>::R>();
  else
  if (pobject->template is<typename simple_object_traits<T>::Z>())
    t = pobject->template as<typename simple_object_traits<T>::Z>();
  else
    throw std::string("python_util::extract_complex_number: object is not <complex number> or <real number> or <integral number>");  
}

template <class T>
inline void extract_real_number(T& t, const simple_object_base* pobject) throw(std::string)
{
  if (pobject->template is<typename simple_object_traits<T>::R>())
    t = pobject->template as<typename simple_object_traits<T>::R>();
  else
  if (pobject->template is<typename simple_object_traits<T>::Z>())
    t = pobject->template as<typename simple_object_traits<T>::Z>();
  else
    throw std::string("python_util::extract_real_number: object is not <real number> or <integral number>");  
}

template <class T>
inline void extract_integral_number(T& t, const simple_object_base* pobject) throw(std::string)
{
  if (pobject->template is<typename simple_object_traits<T>::Z>())
    t = pobject->template as<typename simple_object_traits<T>::Z>();
  else
    throw std::string("python_util::extract_integral_number: object is not <integral number>");  
}

template < >
inline void extract_integral_number<size_t>(size_t& t, const simple_object_base* pobject) throw(std::string)
{
  if (pobject->is<simple_object_traits<size_t>::Z>())
    t = static_cast<size_t>(pobject->as<simple_object_traits<size_t>::Z>());
  else
    throw std::string("python_util::extract_integral_number: object is not <integral number>");  
}

// ---------------------------------------- end: dispatcher classes: ----------------------------------


// ------------------------- simple_object<T>: -----------------------------------------

// allocation control methods:
template <class T>
void simple_object<T>::free_pointers_(void)
{  
  // note: if size is zero (pre-allocation, e.g. during readBinaryVirtual): ptr_<T>() will be NULL:
  //  (at present NULL ptr _also_ represents an empty RANK>0 number-object: this may not be a good idea.)
  assert((size_ > 0)  || (NULL == ptr_<T>()));
  
  if (own_data_ && (NULL != ptr_<T>())){
      // Implementation note: allocate with "array" new regardless of size.
      //   Otherwise there will be a bug for _generic_ deallocation when the
      //   data-ownership is transferred away from this object for arrays of size == 1.
      delete[] &as_<T>();
  }
       
  ptr_<T>() = NULL;
  size_ = 0;
}

template <class T>
void simple_object<T>::resize_(size_t new_size, bool new_own_data)
{
  // dissallow resize to zero:
  assert(new_size > 0);

  if ((own_data_ != new_own_data) || (size_ != new_size)){
    free_pointers_();
    if (new_own_data){
      // Implementation note: allocate with "array" new regardless of size.
      //   Otherwise there will be a bug for _generic_ deallocation when the
      //   data-ownership is transferred away from this object for arrays of size == 1.
      ptr_<T>() = new T[new_size];
    }
    size_ = new_size;
    own_data_ = new_own_data;
  }  
}

// Note: explicit instantiantions of simple_object<object_list> and simple_object<object_map> mean that
//   more general pointer-object ownership does not need to be dealt with in these methods.

template <class T>
inline bool simple_object<T>::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{
  using commUtil::readBinary;
  bool status(true);
  
  free_pointers_();
  
  size_t len(0);
  status = (status && readBinary(fp, len));
  
  if (status && (len > 0)){
    // binary-read => data-ownership:
    resize_(len, true);
    
    for(T *it = ptr_<T>(), *itEnd = it + size_;
        status && (it != itEnd);
        ++it)
      status = (status && readBinary(fp, *it));  
  }
    
  return status; 
}

template <class T>
inline bool simple_object<T>::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string)      
{
  using commUtil::writeBinary;
  bool status(true);
  
  status = (status && writeBinary(fp, size_));
  
  for(const T *it = ptr_<T>(), *itEnd = it + size_;
      status && (it != itEnd);
      ++it)
    status = (status && writeBinary(fp, *it));  
  
  return status; 
}

template <class T>
inline size_t simple_object<T>::binarySize(void)const throw(std::string)
{ 
  using commUtil::binarySize;
  size_t val(0);
  
  val += sizeof(size_t);
  val += size_ * binarySize(*ptr_<T>());
  return val; 
}

template <class T>
const std::type_info& simple_object<T>::type(void)const
{
  return typeid(T);
}

/**
 * @brief Allow contiguous storage arrays (type method will return type of value-type).
 *   IMPLEMENTATION NOTES:
 *     - allowing size > 1 will complicate extraction for target list types, but these _target_ types in general
 *         are single value-type container types, so it shouldn't be too difficult to implement;
 *     - this method shall return 1 for simple_object<object_map>, simple_object<object_list> specializations,
 *       which will also _hide_ the associated (*T, len, own_data=true) constructor signature;
 *     - all allocation management is at simple_object.
 */
template <class T>
size_t simple_object<T>::size(void)const
{ return size_; }

    
/**
 * @brief Data ownership state for the simple_object.
 */
template <class T>
bool simple_object<T>::own_data(void)const
{ return own_data_; }
    
/**
 * @brief Reset data-ownership flag associated with the simple_object.
 *   - At present, this is allowed for all single-entity simple_object 
 *    (i.e. not simple_object<object_list> or simple_object<object_map>).
 *   - Calling release_data for non single-entity objects, or when the 
 *     ownership flag is not set is an error.
 *   .
 */
template <class T>
void simple_object<T>::release_data(void)const throw(std::runtime_error)
{
  if (own_data_)
    own_data_ = false;
  else
    throw std::runtime_error("simple_object<T>::release_data: ownership flag not set");
}
 
template <class T>
simple_object_base* simple_object<T>::clone(void)const
{
  return new simple_object(*this);
}

template <class T>
simple_object<T>& simple_object<T>::operator=(const simple_object<T>& other) throw(std::runtime_error)
{
  resize_(other.size_, other.own_data_); 
  if (own_data_){
    if (1 == size_)  
      as_<T>() = other.as_<T>();
    else
      std::copy(other.ptr_<T>(), other.ptr_<T>() + size_, ptr_<T>());  
  }
  else
    ptr_<T>() = other.ptr_<T>();
    
  return *this;
}

template <class T>
simple_object<T>& simple_object<T>::operator=(const T& other) throw(std::runtime_error)
{
  resize_(1, true); // scalar assignment => data-ownership
  as_<T>() = other;
  return *this;
}


template <class T>
simple_object<T>::~simple_object(void)
{
  // Implementation note: at present: size_==0 represents an empty RANK>0 number object: this may not be a good idea.
  //   (I would prefer always having the pointer non-NULL once the object exists).
  assert((NULL != ptr_<T>()) || (size_ == 0)); 
  
  free_pointers_();
}

template <class T>
simple_object<T>::simple_object(void)
  : simple_object_base(new simple_object_base::template generic_ptr<T>()), size_(0), own_data_(true)
{
  // defer any allocation until the real size is known (at first assignment or binary-read):
  #if 0
  ptr_<T>() = new T[1];
  #endif 
}

template <class T>
simple_object<T>::simple_object(const simple_object<T>& other)
  : simple_object_base(new simple_object_base::template generic_ptr<T>()), size_(other.size_), own_data_(other.own_data_)
{
  // allow (len == 0) as an alternative to the void constructor:
  // note: ptr<T>() is initialized to NULL by "generic_ptr<T>" constructor.
  
  if (own_data_){
    if (size_ > 0){
      // all allocation with "array" new:
      ptr_<T>() = new T[size_];
      std::copy(other.ptr_<T>(), other.ptr_<T>() + size_, ptr_<T>());
    }
  }
  else
    ptr_<T>() = const_cast<simple_object<T>&>(other).ptr_<T>();  
}

template <class T>
simple_object<T>::simple_object(const T& t)
  : simple_object_base(new simple_object_base::template generic_ptr<T>()), size_(1), own_data_(true)
{
  // all allocation with "array" new:
  ptr_<T>() = new T[1];
  *ptr_<T>() = t;
}

/**
 * @brief Construct simple object from contiguous in-memory data.
 *   @param[in] ptr  data pointer
 *   @param[in] len  number of data elements
 *   @param[in] own_data  true => allocate, otherwise reference the data
 *   @param[in] copy_data  copy from the in-memory data
 */
template <class T>
simple_object<T>::simple_object(const T* ptr, size_t len, bool own_data, bool copy_data)    
  : simple_object_base(new simple_object_base::template generic_ptr<T>()), size_(len), own_data_(own_data)
{
  // allow (len == 0) as an alternative to the void constructor:
  // note: ptr<T>() is initialized to NULL by "generic_ptr<T>" constructor.
  if (size_ > 0){
    if (own_data_){
      ptr_<T>() = new T[size_];
      if (copy_data)
        std::copy(ptr, ptr + size_, ptr_<T>());
    }
    else
      ptr_<T>() = const_cast<T*>(ptr);
  }    
}

// special constructors for simple_object<object_list> and simple_object<object_map> explicit instantiations:
// (note: "size_" and "own_data_" attributes should not exist for these specializations).
template <class U>
simple_object<simple_object_base::object_list>::simple_object(const U& other)
  : simple_object_base(new simple_object_base::template generic_ptr<simple_object_base::object_list>())
{
  ptr_<object_list>() = new object_list(other.size(), NULL);
  object_list &self(as_<object_list>());
  typename U::const_iterator itOther(other.begin());
  for(object_list::iterator itSelf = self.begin(), itSelfEnd = self.end();
      itSelf != itSelfEnd;
      ++itSelf, ++itOther)
    *itSelf = new_simple_object(*itOther); // possibly recursive
}

template <class M>
simple_object<simple_object_base::object_map>::simple_object(const M& other)
  : simple_object_base(new simple_object_base::template generic_ptr<simple_object_base::object_map>())
{
  ptr_<object_map>() = new object_map();
  object_map &self(as_<object_map>());
  for(typename M::const_iterator itOther = other.begin(), itOtherEnd = other.end();
      itOther != itOtherEnd;
      ++itOther)
    self[(*itOther).first] = new_simple_object((*itOther).second); // possibly recursive
}

// --------------------------------- extensible_parameters_base<T>: ----------------------------------------

template <class T>
void extensible_parameters_base<T>::init(const std::string* pname, const std::vector<T>* pcoeff)
{
  assert(pobject_ == NULL);
  pobject_ = new simple_object<object_list>();
  object_list &l(pobject_->as<object_list>());
  l.resize(3,NULL); // [<name>, <coefficients>, <named-parm map>]: where coefficients is list of T*
  
  // set name:
  l[0] = new simple_object<std::string>(pname == NULL? "": pname->c_str());

  // set coefficients:
  if (pcoeff == NULL)
    l[1] = new simple_object<object_list>();
  else
    set_coeff(*pcoeff);
    
  // set named-parm map:
  l[2] = new simple_object<object_map>();  
}

template <class T>
bool extensible_parameters_base<T>::enforce_usage(void)const
{
  bool test(true);
  if (pobject_ != NULL && pobject_->is<object_list>()){
    const object_list& l(pobject_->as<object_list>());
    if (!(l.size() >= 3)
        || !(l[0]->is<std::string>())
        || !(l[1]->is<object_list>())
        || !(l[2]->is<object_map>()))
      test = false;
  }
  else
    test = false;
  return test;   
}

// named-parm map:
template <class T>
inline simple_object_base::object_map& extensible_parameters_base<T>::parm_map(void)
{ return (pobject_->as<object_list>())[2]->as<object_map>(); }

template <class T>
inline const simple_object_base::object_map& extensible_parameters_base<T>::parm_map(void)const
{ return (pobject_->as<object_list>())[2]->as<object_map>(); }
 

template <class T>
inline const simple_object_base* extensible_parameters_base<T>::object_pointer(void)const
{ return pobject_; }


template <class T>
void extensible_parameters_base<T>::clone_from(const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
{ 
  if (pobject_ != NULL)
    delete pobject_;
    
  if (!transfer_ownership)
    pobject_ = pobject->clone();
  else
    pobject_ = const_cast<simple_object_base*>(pobject); 
     
  if (!enforce_usage())
    throw std::string("extensible_parameters_base<T>::clone_from: object as cloned does not have correct structure: \n"
                      "   (<name>, <coefficients list>, <named parm map> ...)"); 
}

template <class T>
extensible_parameters_base<T>* extensible_parameters_base<T>::clone(void)const throw(std::string)
{ return new extensible_parameters_base(*this); }

template <class T>
inline std::string& extensible_parameters_base<T>::name(void) throw(std::string)
{ return (pobject_->as<object_list>())[0]->as<std::string>(); }

template <class T>
inline const std::string& extensible_parameters_base<T>::name(void)const throw(std::string)
{ return (pobject_->as<object_list>())[0]->as<std::string>(); }


// U is any iterable container type of value_type convertable to T via direct assignment as T = U::value_type:
template <class T>
template <class U>
void extensible_parameters_base<T>::set_coeff(const U& src) throw(std::string)
{
  // see comment on "usage": object will be a list, and must contain at least two elements, the second of which is a list of T* (possibly empty).
  assert(enforce_usage());
  
  // delete the pointer to the coefficients-list and re-create the list (easiest way to correctly delete simple_object_base*):       
  object_list &l(pobject_->as<object_list>());
  if (l[1] != NULL) // allow use in "init" method where l[1] == NULL
    delete l[1];
  l[1] = new simple_object<object_list>();
  object_list &coeff(l[1]->as<object_list>());
  coeff.resize(src.size(),NULL);
  
  typename U::const_iterator it_src(src.begin());
  for(object_list::iterator it_coeff = coeff.begin(), it_coeffEnd = coeff.end();
      it_coeff != it_coeffEnd;
      ++it_coeff, ++it_src)
    *it_coeff = new simple_object<T>(*it_src);           
}

// U is any iterable container type of value_type convertable to T via direct assignment as U::value_type = T:
template <class T>
template <class U>
void extensible_parameters_base<T>::get_coeff(U& dest)const throw(std::string)
{
  // see comment on "usage": object will be a list, and must contain at least two elements, the second of which is a list of T* (possibly empty).
  assert(enforce_usage());
        
  const object_list &coeff((pobject_->as<object_list>())[1]->as<object_list>());
  dest.clear();
  dest.resize(coeff.size()); // assume U::value_type() initializes the value_type correctly.
  object_list::const_iterator it_coeff(coeff.begin());
  for(typename U::iterator it_dest = dest.begin(), it_destEnd = dest.end();
      it_dest != it_destEnd;
      ++it_dest, ++it_coeff)
    extract_number(*it_dest, *it_coeff); // extract number types, with promotions C <-- {C, R, Z}, R <-- {R, Z}, Z <-- Z  
}

    
// arbitrary named parameters:

// test if a parm exists in the map:
template <class T>
inline bool extensible_parameters_base<T>::has_named_parm(const std::string& key)const
{ return simple_object_base::has_named_parm(parm_map(), key); }

template <class T>
inline const std::type_info& extensible_parameters_base<T>::named_parm_type(const std::string& key)const throw(std::string) 
{ return simple_object_base::named_parm_type(parm_map(), key); }

// single_entity types U (see simple_object_traits<T>):
// (note: use of this form with list_entity or object_entity types throws exception)
template <class T>
template <class U>
inline const U& extensible_parameters_base<T>::get_named_parm(const std::string& key)const throw(std::string)
{ return simple_object_base::template get_named_parm<U>(parm_map(), key); }

template <class T>
template <class U>
inline void extensible_parameters_base<T>::set_named_parm(const std::string& key, const U& val) throw(std::string)
{ simple_object_base::template set_named_parm<U>(parm_map(), key, val); }

// container types U (see simple_object_traits<T>):
template <class T>
template <class U>
inline void extensible_parameters_base<T>::get_named_parm(const std::string& key, U& val)const throw(std::string)
{ simple_object_base::template get_named_parm<U>(parm_map(), key, val); }

// object pointers themselves:
template <class T>
inline const simple_object_base* extensible_parameters_base<T>::get_named_object(const std::string& key)const throw(std::string)
{ return simple_object_base::get_named_object(parm_map(), key); }

template <class T>
inline void extensible_parameters_base<T>::set_named_object(const std::string& key, const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
{ simple_object_base::set_named_object(parm_map(), key, pobject, transfer_ownership); }


template <class T>
void extensible_parameters_base<T>::dump_parm_map(void)const
{
  using std::cout;
  using std::endl;
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  cout<<"dump of parm map: \n";
  if (pobject_ != NULL){
    for(object_map::const_iterator itM = parm_map().begin(), itMEnd = parm_map().end();
        itM != itMEnd;
        ++itM){
      cout<<"key: "<<(*itM).first<<", type: "<<(*itM).second->type().name()<<", value: ";
      (*itM).second->write_repn<C,R,Z>(cout);
      cout<<"\n";
    }
  }
  else
    cout<<"  object pointer is NULL";  
  cout<<endl;
}


template <class T>
inline void extensible_parameters_base<T>::swap(extensible_parameters_base& other)
{ std::swap(pobject_, other.pobject_); }

template <class T>
void extensible_parameters_base<T>::copy(const extensible_parameters_base& other)
{
  if (pobject_ != NULL){
    delete pobject_;
    pobject_ = NULL;
  }
  if (other.pobject_ != NULL)
    pobject_ = other.pobject_->clone();   
}

template <class T>
inline extensible_parameters_base<T>& extensible_parameters_base<T>::operator=(const extensible_parameters_base& other)
{
  copy(other);
  return *this;
}


template <class T>
bool extensible_parameters_base<T>::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{
  bool status(true);
  if (pobject_ != NULL){
    delete pobject_;
    pobject_ = NULL;
  }

  status = (status && simple_object_base::readBinaryVirtual(fp, pobject_));
  if (status && !enforce_usage())
    throw std::string("extensible_parameters_base<T>::readBinary: object as read does not have correct structure: \n"
                      "   (<name>, <coefficients list>, ...)");
  return status;     
}

template <class T>
bool extensible_parameters_base<T>::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string)      
{
  bool status(true);
  status = (status && simple_object_base::writeBinaryVirtual(fp, pobject_));
  return status;
}

template <class T>
size_t extensible_parameters_base<T>::binarySize(void)const throw(std::string)
{
  size_t val(0);
    
  val += simple_object_base::binarySizeVirtual(pobject_);
  return val;
}


template <class T>
template <class C, class R, class Z>
void extensible_parameters_base<T>::read(std::istream& is) throw(std::string)
{
  std::string src;
  scanner_util::load_group_from_stream(src, is, '(', ')');
  if (is){
    if (pobject_ != NULL){
      delete pobject_;
      pobject_ = NULL;
    }
    pobject_ = simple_object_base::read_repn_virtual<C,R,Z>(src);  
    if (!enforce_usage())
      throw std::string("extensible_parameters_base<T>::read: object as read does not have correct structure: \n"
                      "   (<name>, <coefficients list>, ...)");
  }
  // non-exception status returned via stream flags.
}

template <class T>
template <class C, class R, class Z>
void extensible_parameters_base<T>::write(std::ostream& os)const throw(std::string)
{ 
  assert(pobject_ != NULL);
  pobject_->write_repn<C,R,Z>(os);
  // non-exception status returned via stream flags.
}


template <class T>
extensible_parameters_base<T>::~extensible_parameters_base(void)
{
  if (pobject_ != NULL){
    delete pobject_;    
    pobject_ = NULL;
  } 
}

template <class T>
extensible_parameters_base<T>::extensible_parameters_base(bool initialize)
  : pobject_(NULL)
{
  if (initialize)
    init(); 
}

template <class T>
extensible_parameters_base<T>::extensible_parameters_base(const extensible_parameters_base& other)
  : pobject_(NULL)
{
  // note: "init()" _not_ required prior to copy.
  copy(other);
}

template <class T>
extensible_parameters_base<T>::extensible_parameters_base(const std::string& name, const std::vector<T>* pcoeff)
  : pobject_(NULL)
{
  init(&name, pcoeff);
}

template <class T>
extensible_parameters_base<T>::extensible_parameters_base(const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
  : pobject_(NULL)
{ 
  clone_from(pobject, transfer_ownership);
}


// ------------------------------- command-line options implementation: -----------------------------------------------------
template <class T>
inline void options_map<T>::init(bool initialize)
{
  assert(pobject_ == NULL);
  if (initialize)
    pobject_ = new simple_object<object_map>();
  set_cache_current(false);
}

template <class T>
template <class C, class R, class Z>
void options_map<T>::init(const std::string& src, const std::string* ppreamble, const std::string* ppostscript) throw(std::string)
{
  assert(pobject_ == NULL);
  pobject_ = simple_object_base::template read_repn_virtual<C,R,Z>(src, ppreamble, ppostscript); 
  if (!enforce_usage())
    throw std::string("options_map<T>::init: object as initialized does not have correct structure: \n"
                      "   <object map>");  
  set_cache_current(false);
  update_cache();
}

template <class T>
template <class C, class R, class Z>
void options_map<T>::init(std::istream& src, 
          std::istream* ppreamble, 
          std::istream* ppostscript) throw(std::string)
{
  assert(pobject_ == NULL);
  pobject_ = simple_object_base::template read_repn_virtual<C,R,Z>(src, ppreamble, ppostscript);   
  if (!enforce_usage())
    throw std::string("options_map<T>::init: object as initialized does not have correct structure: \n"
                      "   <object map>");  
  set_cache_current(false);
  update_cache();
}          

template <class T>
inline bool options_map<T>::enforce_usage(void)const
{
  bool test(true);
  if (pobject_ == NULL || !pobject_->is<object_map>())
    test = false;
  return test;   
}

// named-parm map:
template <class T>
inline simple_object_base::object_map& options_map<T>::parm_map(void) 
{ return (pobject_->as<object_map>()); }

template <class T>
inline const simple_object_base::object_map& options_map<T>::parm_map(void)const 
{ return (pobject_->as<object_map>()); }


template <class T>
void options_map<T>::dump_parm_map(void)const
{
  using std::cout;
  using std::endl;
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  cout<<"dump of parm map: \n";
  if (pobject_ != NULL){
    for(object_map::const_iterator itM = parm_map().begin(), itMEnd = parm_map().end();
        itM != itMEnd;
        ++itM){
      cout<<"key: "<<(*itM).first<<", type: "<<(*itM).second->type().name()<<", value: ";
      (*itM).second->write_repn<C,R,Z>(cout);
      cout<<"\n";
    }
  }
  else
    cout<<"  object pointer is NULL";  
  cout<<endl;
}


template <class T>
inline bool options_map<T>::cache_current(void)const
{ return cache_current_; }

template <class T>
inline void options_map<T>::set_cache_current(bool flag)
{ cache_current_ = flag; }

template <class T>
void options_map<T>::update_cache(bool derived, bool to_cache) throw(std::string)
{ 
  // base_class: nothing to do except update the flag:
  if (!derived)
    cache_current_ = true; 
}

    
template <class T>
inline const simple_object_base* options_map<T>::object_pointer(void)const
{ return pobject_; }

// "clone_from" method clones object pointer from simple_object_base*, using its "clone" method,
//    while enforcing structure of this base class (see usage comment above):
// (transfer_ownership => assign object pointer: use with care)
template <class T>
void options_map<T>::clone_from(const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
{ 
  if (pobject_ != NULL)
    delete pobject_;
    
  if (!transfer_ownership)
    pobject_ = pobject->clone();
  else
    pobject_ = const_cast<simple_object_base*>(pobject); 
     
  if (!enforce_usage())
    throw std::string("options_map<T>::clone_from: object as cloned does not have correct structure: \n"
                      "   <object map>"); 
  
  set_cache_current(false);
  update_cache();
}

template <class T>
inline options_map<T>* options_map<T>::clone(void)const throw(std::string)
{ 
  options_map<T> *val = new options_map(*this); 
  // note: copy constructor updates cache
  return val;
}

// test if a parm exists:
template <class T>
inline bool options_map<T>::has_named_parm(const std::string& key)const
{ return simple_object_base::has_named_parm(parm_map(), key); }

template <class T>
inline const std::type_info& options_map<T>::named_parm_type(const std::string& key)const throw(std::string) 
{ return simple_object_base::named_parm_type(parm_map(), key); }

// single_entity types U (see simple_object_traits<T>):
// (note: use of this form with list_entity or object_entity types throws exception)
template <class T>
template <class U>
inline const U& options_map<T>::get_named_parm(const std::string& key)const throw(std::string)
{ return simple_object_base::template get_named_parm<U>(parm_map(), key); }

// container types U (see simple_object_traits<T>):
template <class T>
template <class U>
inline void options_map<T>::get_named_parm(const std::string& key, U& val)const throw(std::string)
{ simple_object_base::template get_named_parm<U>(parm_map(), key, val); }

template <class T>
template <class U>
inline void options_map<T>::set_named_parm(const std::string& key, const U& val, bool from_cache) throw(std::string)
{ 
  simple_object_base::template set_named_parm<U>(parm_map(), key, val); 
  if (!from_cache){
    set_cache_current(false);
    update_cache();
  }
}

// object pointers themselves:
template <class T>
inline const simple_object_base* options_map<T>::get_named_object(const std::string& key)const throw(std::string)
{ return simple_object_base::get_named_object(parm_map(), key); }

template <class T>
inline void options_map<T>::set_named_object(const std::string& key, const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
{ 
  simple_object_base::set_named_object(parm_map(), key, pobject, transfer_ownership); 
  set_cache_current(false);
  update_cache();  
}


template <class T>
inline void options_map<T>::swap(options_map& other)
{ 
  std::swap(pobject_, other.pobject_); 
  other.set_cache_current(false);
  other.update_cache(); 
  set_cache_current(false);
  update_cache(); 
}


/*
 * derived flag: set to true if called _within_ derived-class copy method:
 * (important:
 *    -- for purposes of parameter caching alone, derived-class does _not_ need to implement a copy method
 *       (it just needs to implement "update_cache") 
 *    -- _if_ implemented, however, derived-class "copy" method needs to first call this base class "copy" with derived _true_, 
 *          and then call its own "invalidate_cache", and "update_cache" (i.e. cache update should happen _last_)
 *    -- derived class "copy constructor" should just call "copy" _without_ setting derived, and further,. it should _not_ construct the base class
 *         in its CTOR using the base-class copy constructor 
 *         (which would work, but would generate multiple calls to base_class "copy" and "update_cache").
 */
template <class T>
void options_map<T>::copy(const options_map& other, bool derived) throw(std::string)
{
  if (pobject_ != NULL){
    delete pobject_;
    pobject_ = NULL;
  }
  if (other.pobject_ != NULL)
    pobject_ = other.pobject_->clone();   

  if (!derived){
    // this clause (or its analogy) should only execute at at level of most-derived class:
    set_cache_current(false);
    update_cache();
  } 
}

template <class T>
inline options_map<T>& options_map<T>::operator=(const options_map& other) throw(std::string)
{
  copy(other);
  return *this;
}

/*
 * usage note: parameter cache does not need to participate in binary I/O;
 *   it is updated at end of successful "readBinary"
 *
 */
template <class T>
inline bool options_map<T>::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{
  bool status(true);
  if (pobject_ != NULL){
    delete pobject_;
    pobject_ = NULL;
  }

  status = (status && simple_object_base::readBinaryVirtual(fp, pobject_));
  if (status){
    if (!enforce_usage())
    throw std::string("options_map<T>::readBinary: object as read does not have correct structure: \n"
                      "   <object map>");    
    set_cache_current(false);
    update_cache();   
  }
  
  return status;     
}

template <class T>
inline bool options_map<T>::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string) 
{
  bool status(true);
  status = (status && simple_object_base::writeBinaryVirtual(fp, pobject_));
  return status;
}
     
template <class T>
inline size_t options_map<T>::binarySize(void)const throw(std::string)
{
  size_t val(0);
  
  val += simple_object_base::binarySizeVirtual(pobject_);
  return val;
}

// for the following methods, break-out number types as explicit parameters
// (this keeps module independent from namespace TMatrix, otherwise, I need
//    TMatrix::numberTraits<T> to obtain the types)
template <class T>
template <class C, class R, class Z>
void options_map<T>::read(std::istream& is) throw(std::string)
{
  std::string src;
  scanner_util::load_group_from_stream(src, is, '{', '}');
  if (is){
    if (pobject_ != NULL){
      delete pobject_;
      pobject_ = NULL;
    }
    pobject_ = simple_object_base::read_repn_virtual<C,R,Z>(src);  

    if (!enforce_usage())
      throw std::string("options_map<T>::read: object as read does not have correct structure: \n"
                        "   <object map>");

    set_cache_current(false);
    update_cache();
  }
  // non-exception error status returned via stream flags.
}


template <class T>
template <class C, class R, class Z>   
void options_map<T>::write(std::ostream& os)const throw(std::string)
{ 
  if (pobject_ == NULL)
    os<<"NULL";
  else
    pobject_->write_repn<C,R,Z>(os);
  // non-exception status returned via stream flags.
}

template <class T>
options_map<T>::~options_map(void)
{
  if (pobject_ != NULL){
    delete pobject_;    
    pobject_ = NULL;
  } 
}


template <class T>
inline options_map<T>::options_map(bool initialize)
  : pobject_(NULL), cache_current_(false)
{ init(initialize); }

template <class T>
inline options_map<T>::options_map(const options_map& other) throw(std::string)
  : pobject_(NULL), cache_current_(false)
{
  // notes: 
  // (1) "init()" _not_ required prior to copy.
  // (2) copy calls _virtual_ "update_cache": this will work correctly provided any class
  //   derived from this one, calls its _own_ update_cache in its copy / copy-constructor,
  //   in which case to also use this base-class "copy" it should call base_class::copy(other, true)
  //   to indicate this usage.
  copy(other);
}

#if 1
template <class T>
template <class C, class R, class Z>
inline options_map<T> options_map<T>::parse_options(const std::string& src, const std::string* ppreamble, const std::string* ppostscript) throw(std::string)
{
  options_map<T> val;
  val.init<C,R,Z>(src, ppreamble, ppostscript);
  return val;
}

template <class T>
template <class C, class R, class Z>
inline options_map<T> options_map<T>::parse_options(
            std::istream& src, 
            std::istream* ppreamble, 
            std::istream* ppostscript) throw(std::string)
{
  options_map<T> val;
  val.init<C,R,Z>(src, ppreamble, ppostscript);
  return val;
}

template <class T>
template <class C, class R, class Z>
inline options_map<T> options_map<T>::parse_options(
            int argc, const char *argv[], 
            std::istream* ppreamble, 
            std::istream* ppostscript) throw(std::string)
{
  options_map<T> val;
  std::string src_;
  scanner_util::argument_string(src_, argc, argv);
  if (src_.empty())
    src_ = "{}";  // init to empty map    
  std::istringstream src(src_);
  
  val.init<C,R,Z>(src, ppreamble, ppostscript);
  return val;
}
#else
template <class T>
template <class C, class R, class Z>
inline options_map<T>::options_map(const std::string& src, const std::string* ppreamble, const std::string* ppostscript) throw(std::string)
  : pobject_(NULL), cache_current_(false)
{
  init<C,R,Z>(src, ppreamble, ppostscript);
}

template <class T>
template <class C, class R, class Z>
inline options_map<T>::options_map(std::istream& src, 
            std::istream* ppreamble, 
            std::istream* ppostscript) throw(std::string)
  : pobject_(NULL), cache_current_(false)
{
  init<C,R,Z>(src, ppreamble, ppostscript);
}

template <class T>
template <class C, class R, class Z>
inline options_map<T>::options_map(
            int argc, const char *argv[], 
            std::istream* ppreamble, 
            std::istream* ppostscript) throw(std::string)
  : pobject_(NULL), cache_current_(false)
{
  std::string src_;
  scanner_util::argument_string(src_, argc, argv);
  if (src_.empty())
    src_ = "{}";  // init to empty map    
  std::istringstream src(src_);
  
  init<C,R,Z>(src, ppreamble, ppostscript);
}
#endif

template <class T>
options_map<T>::options_map(const simple_object_base* object, bool transfer_ownership) throw(std::string)
  : pobject_(NULL), cache_current_(false)
{ 
  clone_from(pobject_, transfer_ownership);
  if (!enforce_usage())
    throw std::string("options_map<T>::options_map: object as cloned does not have correct structure: \n"
                    "   <object map>");
}


} // namespace python_util

#endif // __python_util_template__h
