#if !defined(__python_util_inline__h)
#define __python_util_inline__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util_inline.h,v $

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

//! simple_object is a RANK=0 object (though not necessarily a number object):
inline bool simple_object_base::is_scalar(void)const
{ return (size() == 1) && !is<object_list>(); }

/**
  * @brief test if simple_object is one of: RANK>0 object (with N_components = 0), empty list, empty map.
  */
inline bool simple_object_base::is_empty(void)const
{ 
  return 
    (size() == 0) 
    || (is<object_list>() && as_<object_list>().size() == 0) 
    || (is<object_map>() && as_<object_map>().size() == 0);
}

/**
 * @brief Set default method for extraction from python object.
 */
inline void simple_object_base::set_default_extractor(extract_func func)
{ default_extract_func_ = func; }

/**
 * @brief Convert the python object to a simple_object_base* using the default extractor: returns a new pointer.
 *    The "PyObject* style" extract signature (and meaning) is used here (read/write array references version).
 */
inline simple_object_base* simple_object_base::extract(const PyObject* src)
{ return default_extract_func_(src); }

/**
 * @brief Set default method for insertion to python object.
 */
inline void simple_object_base::set_default_inserter(insert_func func)
{ default_insert_func_ = func; }

/**
 * @brief Convert the simple_object_base* to a python object using the default inserter: returns a new reference.
 *    The "PyObject* style" insert signature (and meaning) is used here.
 *    Warning: any python arrays created to reference RANK>0 number data in the source simple_object 
 *    will have read/write access; data ownership of these objects is transferred to the python object.
 */
inline PyObject* simple_object_base::insert(const simple_object_base* src)
{ return default_insert_func_(src); }

// ------------------------- generic object methods: -------------------------------------------------

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
template <class C, class R, class Z>
template <class T>
inline bool generic_object<C,R,Z>::is(void)const
{ return ptr_->is<T>(); }

template <class C, class R, class Z>
template <class T>
inline const T& generic_object<C,R,Z>::as(void)const
{ return ptr_->as<T>(); }

template <class C, class R, class Z>
template <class T>
inline T& generic_object<C,R,Z>::as(void)
{ return ptr_->as<T>(); }
#endif
 
template <class C, class R, class Z>
inline generic_object<C,R,Z>::operator simple_object_base*(void)const
{ return ptr_; }

template <class C, class R, class Z>
inline generic_object<C,R,Z>::operator simple_object_base*(void)
{ return ptr_; }

template <class C, class R, class Z>
inline simple_object_base* generic_object<C,R,Z>::ptr(void)const
{ return ptr_; }

template <class C, class R, class Z>
inline simple_object_base* generic_object<C,R,Z>::ptr(void)
{ return ptr_; }


template <class C, class R, class Z>
inline bool generic_object<C,R,Z>::readBinary(commUtil::abstractCommHandle *fp)
{ 
  if (NULL != ptr_){
    // all re-allocation handled outside of "readBinaryVirtual": 
    delete ptr_;
    ptr_ = NULL;
  }
  return simple_object_base::readBinaryVirtual(fp, ptr_); 
}

template <class C, class R, class Z>
inline bool generic_object<C,R,Z>::writeBinary(commUtil::abstractCommHandle *fp)const      
{ return simple_object_base::writeBinaryVirtual(fp, ptr_); }

template <class C, class R, class Z>
inline size_t generic_object<C,R,Z>::binarySize(void)const
{ return simple_object_base::binarySizeVirtual(ptr_); }


// ------------------- conversion to and from PyObject*: ----------------------------

/**
 * @brief Convert a generic_object to a python object: returns a new reference.
 */
template <class C, class R, class Z>
inline PyObject* generic_object<C,R,Z>::insert(void)const
{
  return implementation_module::insert<C,R,Z>(ptr_);
}



/**
 * @brief Convert python object to a generic_object.
 */
template <class C, class R, class Z>
inline void generic_object<C,R,Z>::extract(const PyObject* src)
{
  if (NULL != ptr_){
    delete ptr_;
    ptr_ = NULL;
  }  
  ptr_ = implementation_module::extract<C,R,Z>(src);
}  

// ----------------------------------------------------------------------------------

template <class C, class R, class Z>
inline generic_object<C,R,Z>& generic_object<C,R,Z>::operator=(const simple_object_base* other)
{
  clone_from(other);
  return *this; 
}

template <class C, class R, class Z>
inline void generic_object<C,R,Z>::clone_from(const simple_object_base* other, bool transfer_ownership)
{
  if (NULL != ptr_){
    delete ptr_;
    ptr_ = NULL;
  }
  
  if (transfer_ownership)
    ptr_ = const_cast<simple_object_base*>(other);
  else  
  if (NULL != other)
    ptr_ = other->clone();  
}

template <class C, class R, class Z>
inline generic_object<C,R,Z>& generic_object<C,R,Z>::operator=(const generic_object& other)
{ return operator=(other.ptr()); }

template <class C, class R, class Z>
inline generic_object<C,R,Z>::~generic_object(void)
{
  if (NULL != ptr_){
    delete ptr_; 
    ptr_ = NULL;
  }
}

template <class C, class R, class Z>
inline generic_object<C,R,Z>::generic_object(const generic_object& other)
  : ptr_(NULL)
{ operator=(other); }

template <class C, class R, class Z>
inline generic_object<C,R,Z>::generic_object(const simple_object_base *p, bool transfer_ownership)
  : ptr_(NULL)
{
  if (transfer_ownership)
    ptr_ = const_cast<simple_object_base*>(p);
  else  
  if (NULL != p)
    ptr_ = p->clone(); 
}

template <class C, class R, class Z>
inline generic_object<C,R,Z>::generic_object(void)
  : ptr_(NULL)
{ }    

// ------------------------- end: generic object methods: --------------------------------------------


// test if a parm exists:
inline bool simple_object_base::has_named_parm(const object_map& map, const std::string& key)
{
  return (map.find(key) != map.end()); 
}


inline void simple_object<simple_object_base::object_list>::free_pointers_(void)
{
  simple_object_base::free_pointers_(as_<object_list>()); 
}


inline void simple_object<simple_object_base::object_list>::clear(void)
{
  free_pointers_(); 
}

inline void simple_object<simple_object_base::object_map>::free_pointers_(void)
{
  simple_object_base::free_pointers_(as_<object_map>());  
}

inline void simple_object<simple_object_base::object_map>::clear(void)
{
  free_pointers_();  
}

inline void simple_object<simple_object_base::object_list>::assign_with_clone_(const object_list& other)
{
  // clone the object_list:
  free_pointers_();
  object_list &l(as_<object_list>());  
  l.resize(other.size(),NULL);
  object_list::const_iterator itOther = other.begin();
  for(object_list::iterator itL = l.begin(), itLEnd = l.end();
      itL != itLEnd;
      ++itL, ++itOther)
    *itL = (*itOther)->clone();    
}

inline void simple_object<simple_object_base::object_map>::assign_with_clone_(const object_map& other)
{
  // clone the object-map's values:
  free_pointers_();
  object_map &m(as_<object_map>());
  for(object_map::const_iterator itOther = other.begin(), itOtherEnd = other.end();
      itOther != itOtherEnd;
      ++itOther)
    m[(*itOther).first] = (*itOther).second->clone();    
}


} // namespace python_util

#endif // __python_util_inline__h
