// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util.cpp,v $

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


#if 0
// ================================================================
#define EXCLUDE_TEMPLATE_BODIES
// ================================================================
#endif

#include <portability.h>

#include <pthread.h>

#if defined(__USE_MPI)
  #include <mpi.h>
#endif

#if defined(__USE_OMP) 
  #include <omp.h>
#endif

#include <Python.h>
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
// #include <numpy/ufuncobject.h>
// #include <numpy/__multiarray_api.h>
// #include <numpy/__ufunc_api.h>

#include <assert.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <complex>
#include <string>
#include <list>

#include <functional>
#if defined(__PGI)
  // PGI compiler uses STLport standard template library; extensions are in "std" namespace (see _STL_EXT_NAMESPACE_ macro in "portability.h"):
  #include <algorithm>
  #include <numeric>
  // #include <slist>
  #include <hash_map>
#else
  #include <algorithm>
  #include <numeric>
  #include <ext/algorithm>  
  #include <ext/numeric>
  // #include <ext/slist>
  #include <ext/hash_map>
#endif

#include <typeinfo>
#include <tr1/type_traits>

#include <commUtil.h>
using commUtil::abstractCommHandle;
using commUtil::readBinary;
using commUtil::writeBinary;

#include <numericalConstants.h>

// --------- for gmm::binarySize: -------------
#include <gmm/gmm_kernel.h>
#include <gmm/gmm_dense_lu.h>

#include <ghostIterator.h>
#include <gmm_ext.h>
// --------------------------------------------


#include <ntuple.h>
using linalg::ntuple;


#include <scanner_util.h>

#if 0
// ================================================================
#undef EXCLUDE_TEMPLATE_BODIES
// ------- python_util template bodies: this source module: -------
// ================================================================
#endif

#include <python_util.h>

#if 0
// ================================================================
#define EXCLUDE_TEMPLATE_BODIES
// ================================================================
#endif

namespace python_util{

// -------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------- EXPLICIT template and template function instantiations: ---------------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------

#if defined(__USE_MERE)
template bool simple_object_base::is<mere::C>(void)const;
template bool simple_object_base::is<mere::R>(void)const;
template bool simple_object_base::is<mere::Z>(void)const;
#endif
template bool simple_object_base::is<std::complex<double> >(void)const;
template bool simple_object_base::is<double>(void)const;
template bool simple_object_base::is<long>(void)const;

template bool simple_object_base::is<bool>(void)const;
template bool simple_object_base::is<std::string>(void)const;

template bool simple_object_base::is<simple_object_base::object_map>(void)const;
template bool simple_object_base::is<simple_object_base::object_list>(void)const;
// -------------------------------------------------------------

#if 0 
// --------- form for template instantiation: ------------
template class Foo<int>;

// --------- form for function template instantiation: ---
// (may need a forward declaration with "extern" to prevent instantiation at point of first use)
template ostream& operator <<(ostream&, const Foo<int>&);


// --------- list of types for instantiation: -----------
mere::C
mere::R
mere::Z

std::complex<double>
double
long

bool
std::string

std::vector<std::complex<double> >
std::list<std::complex<double> >
ntuple<std::complex<double>,3>
ntuple<std::complex<double>,2>
ntuple<std::complex<double>,1>

std::vector<double>
std::list<double>
ntuple<double,3>
ntuple<double,2>
ntuple<double,1>

simple_object_base
simple_object_base::object_map
simple_object_base::object_list

generic_object< ... above types ... >


// ------------------------------------------------------
#endif

// -------------------------------------------------------------------------------------------------------------------------------------------------------
// ------------------------------------------- end: EXPLICIT template and template function instantiations: ----------------------------------------------
// -------------------------------------------------------------------------------------------------------------------------------------------------------


const char* python_error::what(void)const throw()
{ return what_.c_str(); }

python_error::~python_error(void) throw()
{ }

python_error::python_error(void) throw()
  : what_("unspecified python_error")
{ }

python_error& python_error::operator=(const python_error& other) throw()
{ 
  what_ = other.what_; 
  return *this;
}

python_error::python_error(const python_error& other) throw()
{ operator=(other); }

python_error::python_error(const std::string& msg) throw()
  : what_(msg)
{
  // append any error information from the python C/API error methods:
  what_ += " ";
  what_ += implementation_module::python_error_string();
}


simple_object_base::RECOGNIZED_TYPE simple_object_base::enum_from_type_info(const std::type_info& type_)
{
  /*
   * IMPLEMENTATION NOTE: new gnu C++ rtti ABI _requires_ either comparison by name,
   *   or _excruciatingly_ careful linking restrictions for this section to work correctly
   *   with cross-module objects (e.g. with shared-library files).  For the moment, I simply use
   *   comparison by-name.
   *
   */
   
  RECOGNIZED_TYPE val(NONE_TYPE);
  
  if (0 == strcmp(type_.name(),  typeid(void).name()))
    val = NONE_TYPE;
  #if defined(__USE_MERE) // _allow_ simultaneous utilization of mere types and double-precision types.
  else                    //   === this needs re-working... ===
  if (0 == strcmp(type_.name(), typeid(mere::C).name()))
    val = COMPLEX_TYPE;
  else
  if (0 == strcmp(type_.name(),  typeid(mere::R).name()))
    val = REAL_TYPE;  
  else
  if (0 == strcmp(type_.name(),  typeid(mere::Z).name()))
    val = INTEGER_TYPE;  
  #endif  
  else
  if (0 == strcmp(type_.name(),  typeid(std::complex<double>).name()))
    val = COMPLEX_TYPE;
  else
  if (0 == strcmp(type_.name(),  typeid(double).name()))
    val = REAL_TYPE;  
  else
  if (0 == strcmp(type_.name(),  typeid(long).name()))
    val = INTEGER_TYPE;    
  else
  if (0 == strcmp(type_.name(),  typeid(bool).name()))
    val = BOOL_TYPE;  
  else
  if (0 == strcmp(type_.name(),  typeid(std::string).name()))
    val = STRING_TYPE;  
  else
  if (0 == strcmp(type_.name(),  typeid(object_list).name()))
    val = OBJECT_LIST_TYPE;  
  else
  if (0 == strcmp(type_.name(),  typeid(object_map).name()))
    val = OBJECT_MAP_TYPE;  
  else
    throw std::string("simple_object_base::enum_from_type_info: unrecognized type_info: ") + type_.name();
    
  return val;
}

const std::type_info& simple_object_base::type_info_from_enum(RECOGNIZED_TYPE etype_)
{
  const std::type_info* pval(&typeid(void));
  switch (etype_){
    case NONE_TYPE:
      pval = &typeid(void);
    break;
    #if defined(__USE_MERE) // _allow_ simultaneous utilization of mere types and double-precision types.
    case COMPLEX_TYPE:      //   here this is a bit trickier: default to arbitrary precision, if it is in-use
      pval = &typeid(mere::C);
    break;
    case REAL_TYPE:
      pval = &typeid(mere::R);
    break;
    case INTEGER_TYPE:
      pval = &typeid(mere::Z);
    break;
    #else
    case COMPLEX_TYPE:
      pval = &typeid(std::complex<double>);
    break;
    case REAL_TYPE:
      pval = &typeid(double);
    break;
    case INTEGER_TYPE:
      pval = &typeid(long);
    break;    
    #endif
    case BOOL_TYPE:
      pval = &typeid(bool);
    break;
    case STRING_TYPE:
      pval = &typeid(std::string);
    break;
    case OBJECT_LIST_TYPE:
      pval = &typeid(object_list);
    break;
    case OBJECT_MAP_TYPE:
      pval = &typeid(object_map);
    break;
    default:
      throw std::string("simple_object_base::type_info_from_enum: unrecognized type enum");
    // break;    
  }
  
  return *pval;  
}


simple_object_base::generic_ptr_base* simple_object_base::generic_ptr_base::clone(void)const
{ throw std::string("simple_object_base::generic_ptr_base::clone: generic pointer has undefined type"); }

simple_object_base::generic_ptr_base::~generic_ptr_base(void)
{ }


// ----------------------------- simple_object_base: static variables: -------------------------------------------------
simple_object_base::extract_func simple_object_base::default_extract_func_ = simple_object_base::NULL_extract_func;
simple_object_base::insert_func simple_object_base::default_insert_func_ = simple_object_base::NULL_insert_func;
// ---------------------------------------------------------------------------------------------------------------------

simple_object_base* simple_object_base::NULL_extract_func(const PyObject* src) throw(std::runtime_error)
{ throw std::runtime_error("extract from python object to simple_object_base*: default extract method uninitialized"); }
    
PyObject* simple_object_base::NULL_insert_func(const simple_object_base* src)  throw(std::runtime_error)
{ throw std::runtime_error("insert from simple_object_base* to python object: default insert method uninitialized"); }

        
void simple_object_base::free_pointers_(object_list& l)
{
  for(object_list::iterator itL = l.begin(), itLEnd = l.end();
      itL != itLEnd;
      ++itL)
    if (*itL != NULL) // allow uninitialized object_list pointers (?).   
      delete *itL;
  l.clear();  
}

void simple_object_base::free_pointers_(object_map& m)
{
  for(object_map::iterator itM = m.begin(), itMEnd = m.end();
      itM != itMEnd;
      ++itM)
    if ((*itM).second != NULL)
      delete (*itM).second; // allow uninitialized object_map pointers (?).
  m.clear();  
}


bool simple_object_base::readBinary_(commUtil::abstractCommHandle *fp, object_list& l)
{
   using commUtil::readBinary;
   bool status = true;
   
   free_pointers_(l);
   
   // read size
   size_t size_;
   status = (status && readBinary(fp, size_));
   l.resize(size_,NULL); 

   if (0 < size_)
     for (object_list::iterator itL = l.begin(), itLEnd = l.end();     
          status && (itL != itLEnd);          
          ++itL){           
       status = (status && readBinaryVirtual(fp, *itL));
     }
   return status;
}

bool simple_object_base::writeBinary_(commUtil::abstractCommHandle *fp, const object_list& l)
{
  using commUtil::writeBinary;
  bool status(true);
  
  size_t size_(l.size()); // explicitely control output-type
  status = (status && writeBinary(fp, size_));  
   if (0 < size_)
     for (object_list::const_iterator itL = l.begin(), itLEnd = l.end();     
          status && (itL != itLEnd);          
          ++itL)
       status = (status && writeBinaryVirtual(fp, *itL));

   return status;
}

      
size_t simple_object_base::binarySize_(const object_list& l) throw(std::string)
{
  size_t val(0);
  size_t size_(l.size()); 
  
  val += sizeof(size_t);
  if (0 < size_)
    for (object_list::const_iterator itL = l.begin(), itLEnd = l.end();     
         itL != itLEnd;          
         ++itL)
      val += binarySizeVirtual(*itL);

  return val;
}


bool simple_object_base::readBinary_(commUtil::abstractCommHandle *fp, object_map& m)
{
  using commUtil::readBinary;
  bool status = true;

  free_pointers_(m);

  // read size
  size_t size_;
  status = (status && readBinary(fp, size_));

  std::string key;
  if (0 < size_)
    for (size_t n = 0; status && (n < size_); ++n){
      simple_object_base *val(NULL); // will be cloned when read (init to NULL, otherwise object deleted by "readBinaryVirtual"!)
      status = (status && readBinary(fp, key)); 
      status = (status && readBinaryVirtual(fp, val));
      if (status)
        m[key] = val;
    }     

  return status;
}

bool simple_object_base::writeBinary_(commUtil::abstractCommHandle *fp, const object_map& m)
{
  using commUtil::writeBinary;
  bool status = true;

  // write size (explicitely control output type):
  size_t size_(m.size());
  status = (status && writeBinary(fp, size_));

  if (0 < size_)
    for (object_map::const_iterator itM = m.begin(), itMEnd = m.end();
         status && (itM != itMEnd);
         ++itM){
      status = (status && writeBinary(fp, (*itM).first)); 
      status = (status && writeBinaryVirtual(fp, (*itM).second));
    }     

  return status;
}

      
size_t simple_object_base::binarySize_(const object_map& m) throw(std::string)
{
  throw std::string("simple_object_base::binarySize_: recheck this implementation: \n"
    "  any utilization must take account of the fact that binarySize(std::string) is not constant!");
#if 0
  size_t val(0);
  size_t size_(m.size()); 
  val += sizeof(size_t);
  if (0 < size_)
    for (object_map::const_iterator itM = m.begin(), itMEnd = m.end();     
         itM != itMEnd;          
         ++itM){
      val += binarySize_((*itM).first); // ummm this doesn't work -- need constant key size, I think!
      val += binarySizeVirtual((*itM).second);
    }
  return val;
#endif
}

    
// ---------------------------------------- arbitrary named parameter utility methods: -------------------------------------------

const std::type_info& simple_object_base::named_parm_type(const object_map& map, const std::string& key) throw(std::string) 
{ 
  const std::type_info *ptype(NULL);
  object_map::const_iterator itP = map.find(key);
  if (itP != map.end())
    ptype = &((*itP).second->type());
  else
    throw std::string("named_parm_type: no entry found for specified key: ") + key;   
  return *ptype;
}

// object pointers themselves:
const simple_object_base* simple_object_base::get_named_object(const object_map& map, const std::string& key) throw(std::string)
{
  // note on form: dispatcher-type partial specializations only allowed at namespace scope.
  object_map::const_iterator itP = map.find(key);
  const simple_object_base *pobject(NULL);
  if (itP != map.end())
    pobject = (*itP).second;
  else  
    throw std::string("get_named_object: no entry found for specified key: ") + key;
  
  return pobject; 
}

void simple_object_base::set_named_object(object_map& map, const std::string& key, const simple_object_base* pobject, bool transfer_ownership) throw(std::string)
{
  // Notes on form: 
  //   -- dispatcher-type partial specializations only allowed at namespace scope;
  //   -- one cannot assume that a map-entry pointer is initialized to NULL;
  //   -- no attempt is made to re-use allocation.
  
  std::pair<object_map::iterator, bool> test = map.insert(std::pair<std::string, simple_object_base*>(key,NULL));
  object_map::iterator &itP(test.first);
  if (!test.second && ((*itP).second != NULL))
    delete (*itP).second;
  
  if (transfer_ownership)
    (*itP).second = const_cast<simple_object_base*>(pobject);
  else     
    (*itP).second = pobject->clone();   
}

// ---------------------------------------- end: arbitrary named parameter utility methods: ------------------------------------



// methods to allow binary read and write from pointer to base-class:
bool simple_object_base::readBinaryVirtual(commUtil::abstractCommHandle *fp, simple_object_base*& pobject) throw(std::string)
{
  using commUtil::readBinary;
  bool status(true);
  
  long etype_(static_cast<long>(NONE_TYPE));
  status = (status && readBinary(fp, etype_));

  if (status){
    #if 0
    if (pobject != NULL){
      delete pobject;
      pobject = NULL;
    }
    #else
    assert(pobject == NULL);  // require object allocation management _outside_ of this method!
    #endif
    switch (static_cast<RECOGNIZED_TYPE>(etype_)){
      case NONE_TYPE:
        pobject = NULL;
        #if 1
        throw std::string("simple_object_base::readBinaryVirtual: read of NONE_TYPE object");
        #endif
      // break;
      #if defined(__USE_MERE)
      case COMPLEX_TYPE:
        // use default (zero size) constructors to avoid immediate re-allocation in size > 1 (i.e. array) case:
        pobject = new simple_object<mere::C>(); // (mere::C::zero());
      break;
      case REAL_TYPE:
        pobject = new simple_object<mere::R>(); // (mere::R::zero());
      break;
      case INTEGER_TYPE:
        pobject = new simple_object<mere::Z>(); // (mere::Z::zero());
      break;
      #else
      case COMPLEX_TYPE:
        pobject = new simple_object<std::complex<double> >; // (std::complex<double>(0.0, 0.0));
      break;
      case REAL_TYPE:
        pobject = new simple_object<double>(); // (0.0);
      break;
      case INTEGER_TYPE:
        pobject = new simple_object<long>(); // (0L);
      break;    
      #endif
      case BOOL_TYPE:
        pobject = new simple_object<bool>(); // (false);
      break;
      case STRING_TYPE:
        pobject = new simple_object<std::string>();
      break;
      case OBJECT_LIST_TYPE:
        pobject = new simple_object<object_list>();
      break;
      case OBJECT_MAP_TYPE:
        pobject = new simple_object<object_map>();
      break;
      default:
        throw std::string("simple_object_base::type_info_from_enum: unrecognized type enum");
      // break;    
    }
    // call the appropriate virtual readBinary:
    status = (status && pobject->readBinary(fp));    
  }
  
  return status;  
}

bool simple_object_base::writeBinaryVirtual(commUtil::abstractCommHandle *fp, const simple_object_base* pobject) throw(std::string)
{
  using commUtil::writeBinary;
  bool status(true);
  
  status = (status && writeBinary(fp, static_cast<long>(enum_from_type_info(pobject->type()))));
  status = (status && pobject->writeBinary(fp)); // virtual writeBinary(...)
  return status;  
}

size_t simple_object_base::binarySizeVirtual(const simple_object_base* pobject) throw(std::string)
{
  size_t val(0);
  
  val = sizeof(long); // type enum
  val += pobject->binarySize(); // virtual binarySize()
  return val;
}

bool simple_object_base::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{ 
  #if 0
  throw std::string("simple_object_base::readBinary: _pure_(?!) simple_object_base with ") 
    + (pv_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer"); 
  #else
  throw std::string("simple_object_base::readBinary: _pure_(?!) simple_object_base with ") 
    + (pg_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer");   
  #endif
}

bool simple_object_base::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string)      
{ 
  #if 0
  throw std::string("simple_object_base::writeBinary: _pure_(?!) simple_object_base with ") 
    + (pv_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer"); 
  #else
  throw std::string("simple_object_base::writeBinary: _pure_(?!) simple_object_base with ") 
    + (pg_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer");   
  #endif
}

size_t simple_object_base::binarySize(void)const throw(std::string)
{ 
  #if 0
  throw std::string("simple_object_base::binarySize: _pure_(?!) simple_object_base with ") 
    + (pv_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer"); 
  #else
  throw std::string("simple_object_base::binarySize: _pure_(?!) simple_object_base with ") 
    + (pg_ == NULL?"NULL derived-class pointer":"non-NULL derived-class pointer");   
  #endif
}

const std::type_info& simple_object_base::type(void)const
{ return typeid(void); }

size_t simple_object_base::size(void)const
{ return 1; }
    
/**
 * @brief Data ownership state for the simple_object.
 */
bool simple_object_base::own_data(void)const
{ return true; }

/**
 * @brief Reset data-ownership flag associated with the simple_object.
 *   - At present, this is allowed for all single-entity simple_object 
 *    (i.e. not simple_object<object_list> or simple_object<object_map>).
 *   - Calling release_data for non single-entity objects, or when the 
 *     ownership flag is not set is an error.
 *   .
 */
void simple_object_base::release_data(void)const throw(std::runtime_error)
{
  throw std::runtime_error("simple_object_base::release_data: _pure_(?!) simple_object_base");
}

simple_object_base* simple_object_base::clone(void)const
{
  return new simple_object_base(pg_->clone()); // default is "void" object 
}

simple_object_base::~simple_object_base(void)
{
  assert(NULL != pg_);
  delete pg_; 
}

simple_object_base::simple_object_base(generic_ptr_base* pg)
  : pg_(pg)
{ }  


// ------------------------- simple_object<object_list>: -----------------------------------------

// note: the read/writeBinary methods are "virtual": they cannot be inlined:

bool simple_object<simple_object_base::object_list>::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{ return simple_object_base::readBinary_(fp, as_<object_list>()); }

bool simple_object<simple_object_base::object_list>::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string)      
{ return simple_object_base::writeBinary_(fp, as_<object_list>()); }

size_t simple_object<simple_object_base::object_list>::binarySize(void)const throw(std::string)
{ return simple_object_base::binarySize_(as_<object_list>()); }


const std::type_info& simple_object<simple_object_base::object_list>::type(void)const
{
  return typeid(simple_object_base::object_list);
}

size_t simple_object<simple_object_base::object_list>::size(void)const
{
  return 1;
}

/**
 * @brief Data ownership state for the simple_object.
 */
bool simple_object<simple_object_base::object_list>::own_data(void)const
{ return true; }

/**
 * @brief Reset data-ownership flag associated with the simple_object.
 *   - At present, this is allowed for all single-entity simple_object 
 *    (i.e. not simple_object<object_list> or simple_object<object_map>).
 *   - Calling release_data for non single-entity objects, or when the 
 *     ownership flag is not set is an error.
 *   .
 */
void simple_object<simple_object_base::object_list>::release_data(void)const throw(std::runtime_error)
{
  throw std::runtime_error("simple_object<simple_object_base::object_list>::release_data: \n"
    "  only supported for single-entity simple_object types");
}

#if 0
/**
 * @brief Work-around for the GNU C++ RTTI implementation problem that typeid(T) in different modules are not comparable.
 *    (i.e. test simple_object<T>::static_type == type() rather than typeid(T) == type() for cross-module comparison).
 */
const std::type_info& simple_object<simple_object_base::object_list>::static_type(void) 
{
  return typeid(simple_object_base::object_list);
}
#endif

simple_object_base* simple_object<simple_object_base::object_list>::clone(void)const
{
  return new simple_object<simple_object_base::object_list>(*this);
}

simple_object<simple_object_base::object_list>& 
  simple_object<simple_object_base::object_list>::operator=(const simple_object& other)
{
  // maintain ownership of pointer objects:
  assign_with_clone_(other.as_<object_list>());
  return *this;
}

simple_object<simple_object_base::object_list>& simple_object<simple_object_base::object_list>::operator=(const object_list& other)
{
  // maintain ownership of pointer objects:
  assign_with_clone_(other);
  return *this;
}

simple_object<simple_object_base::object_list>::~simple_object(void)
{
  // enforce usage (see note at simple_object_base::size() declaration):
  assert(size() == 1);
   
  // additional clean-up for lists and maps of simple_object_base pointer:
  free_pointers_();    

  delete ptr_<object_list>(); 
  ptr_<object_list>() = NULL;
}

simple_object<simple_object_base::object_list>::simple_object(void)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_list>())
{
  ptr_<object_list>() = new object_list(); 
}

simple_object<simple_object_base::object_list>::simple_object(const simple_object& other)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_list>())
{
  // this format allows pointer cloning (in "operator="):
  ptr_<object_list>() = new object_list();
  operator=(other); 
}

simple_object<simple_object_base::object_list>::simple_object(const object_list& l)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_list>())
{
  // this format allows pointer cloning (in "operator="):
  ptr_<object_list>() = new object_list();
  operator=(l); 
}

// ------------------------- simple_object<object_map>: -----------------------------------------

// note: the read/writeBinary methods are "virtual": they cannot be inlined:

bool simple_object<simple_object_base::object_map>::readBinary(commUtil::abstractCommHandle *fp) throw(std::string)
{ return simple_object_base::readBinary_(fp, as_<object_map>()); }

bool simple_object<simple_object_base::object_map>::writeBinary(commUtil::abstractCommHandle *fp)const throw(std::string)      
{ return simple_object_base::writeBinary_(fp, as_<object_map>()); }

size_t simple_object<simple_object_base::object_map>::binarySize(void)const throw(std::string)
{ return simple_object_base::binarySize_(as_<object_map>()); }

const std::type_info& simple_object<simple_object_base::object_map>::type(void)const
{
  return typeid(simple_object_base::object_map);
}

size_t simple_object<simple_object_base::object_map>::size(void)const
{
  return 1;
}

/**
 * @brief Data ownership state for the simple_object.
 */
bool simple_object<simple_object_base::object_map>::own_data(void)const
{ return true; }

/**
 * @brief Reset data-ownership flag associated with the simple_object.
 *   - At present, this is allowed for all single-entity simple_object 
 *    (i.e. not simple_object<object_list> or simple_object<object_map>).
 *   - Calling release_data for non single-entity objects, or when the 
 *     ownership flag is not set is an error.
 *   .
 */
void simple_object<simple_object_base::object_map>::release_data(void)const throw(std::runtime_error)
{
  throw std::runtime_error("simple_object<simple_object_base::object_map>::release_data: \n"
    "  only supported for single-entity simple_object types");
}

#if 0
/**
 * @brief Work-around for the GNU C++ RTTI implementation problem that typeid(T) in different modules are not comparable.
 *    (i.e. test simple_object<T>::static_type == type() rather than typeid(T) == type() for cross-module comparison).
 */
const std::type_info& simple_object<simple_object_base::object_map>::static_type(void) 
{
  return typeid(simple_object_base::object_map);
}
#endif

simple_object_base* simple_object<simple_object_base::object_map>::clone(void)const
{
  return new simple_object(*this);
}

simple_object<simple_object_base::object_map>& simple_object<simple_object_base::object_map>::operator=(const simple_object& other)
{
  // maintain ownership of pointer objects:
  assign_with_clone_(other.as_<object_map>());
  return *this;
}

simple_object<simple_object_base::object_map>& simple_object<simple_object_base::object_map>::operator=(const object_map& other)
{
  // maintain ownership of pointer objects:
  assign_with_clone_(other);
  return *this;
}


simple_object<simple_object_base::object_map>::~simple_object(void)
{ 
  // enforce usage (see note at simple_object_base::size() declaration):
  assert(size() == 1);
  
  // additional clean-up for lists and maps of simple_object_base pointer:
  free_pointers_();    
  delete ptr_<object_map>(); 
  ptr_<object_map>() = NULL;
}

simple_object<simple_object_base::object_map>::simple_object(void)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_map>())
{
  ptr_<object_map>() = new object_map(); 
}

simple_object<simple_object_base::object_map>::simple_object(const simple_object& other)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_map>())
{
  // this format allows pointer cloning (in "operator="):
  ptr_<object_map>() = new object_map();
  operator=(other); 
}

simple_object<simple_object_base::object_map>::simple_object(const object_map& m)
  : simple_object_base(new simple_object_base::generic_ptr<simple_object_base::object_map>())
{
  // this format allows pointer cloning (in "operator="):
  ptr_<object_map>() = new object_map();
  operator=(m); 
}


} // namespace python_util


