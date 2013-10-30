#if !defined(__python_util_module__h)
#define __python_util_module__h
//
// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util_module.h,v $

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


//! Modular separation of specifics of python_util implementation.
/**
 * Utilities to support I/O methods in template class "python_list". 
 * Modularization of interface allows isolation of implementation using ``boost::python'' classes. 
 * boost::python utilization occurs only at source-file level 
 *  (this is intended to eventually function as a stand-alone library module, 
 *   and/or possibly to be replace by a separate custom parser of restricted python syntax).
 *
 * Implementation note: due to its wide-use Python C/API type "PyObject*" _will_ be visible outside this module;
 *   participating in "insert" (to the PyObject*) and "extract" (from the PyObject*) methods.
 *
 */

// use LHS = RHS syntax convention for assignment-like methods:

/* 
 * naming convention: 
 *   "insert" (to the PyObject*)
 *   "extract" (from the PyObject*)
 *
 *   "insert_simple_object" (to a std::string)
 *   (could also have named this "repn")
 *
 *   "extract_simple_object" (from a std::string)
 *  
 */
 
namespace implementation_module{

using python_util::python_error;

// "PyObject*" related methods (interpreter initialization occurs _outside_ these methods):

/**
 * @brief Allow transfer-of-ownership of "array" new allocated data to python
 *   (via registration of associated deallocation callback-object as the base-object of the python array).
 */
template <class T>
class python_deallocator{
  public:
    // place this first, assuming that there are be alignment requirements:
    PyObject_HEAD
    
  private:
  
    static PyObject *class_type_object_; // PyTypeObject for the entire class
    T *ptr_;

    static void dealloc_(python_deallocator<T> *self);

    /**
     * @brief unique python name for this deallocation class.
     */
    static const char* class_name(void);
    
    //! private constructor (presently not used anywhere; all initialization is in "register_deallocator"):
    python_deallocator(T* ptr);
      
  public:

    /**
     * @brief Initialize and register the python type associated with this deallocator class.
     *   usage (equivalent to explicit template instantiation): 
     *     in _one_ source-module only:
     *       PyObject *python_deallocator<T>::class_type_object_ = python_deallocator<T>::registered_type_object();
     *     other source modules (and in-fact in the header), must ensure:
     *       extern template class python_deallocator<T>;
     */
    static PyObject* registered_type_object(void) throw(python_error);

    /**
     * @brief Implicitely create a python_deallocator and attach to the owner python-object (which must be an array):
     *   @param[in] owner   owner python-object (which must be an array)
     *     when this object is destroyed the dealloc_ method will be called
     *   @param[in] ptr  pointer to allocated data controlled by this deallocator
     */
    static void attach_deallocator(PyObject *owner, T* ptr) throw(python_error);

}; // class python_deallocator

// =========== EXPLICITELY DISSALLOW MULTIPLE INSTANTIATIONS: =======================

// --- these are the only supported types, other types should throw exception: ---
extern template class python_deallocator<std::complex<double> >;
extern template class python_deallocator<double >;
extern template class python_deallocator<long >;

extern template class python_deallocator<linalg::ntuple<std::complex<double>,1> >;
extern template class python_deallocator<linalg::ntuple<double,1> >;
extern template class python_deallocator<linalg::ntuple<long,1> >;
extern template class python_deallocator<linalg::ntuple<size_t,1> >;

extern template class python_deallocator<linalg::ntuple<std::complex<double>,2> >;
extern template class python_deallocator<linalg::ntuple<double,2> >;
extern template class python_deallocator<linalg::ntuple<long,2> >;
extern template class python_deallocator<linalg::ntuple<size_t,2> >;

extern template class python_deallocator<linalg::ntuple<std::complex<double>,3> >;
extern template class python_deallocator<linalg::ntuple<double,3> >;
extern template class python_deallocator<linalg::ntuple<long,3> >;
extern template class python_deallocator<linalg::ntuple<size_t,3> >;

// ==================================================================================

/**
 * @brief Check the type of a python object.
 *   - for non-pointer types, this checks whether the python-object is directly convertible to the C type;
 *   - for simple pointer-types, this checks direct in-memory correspondence of the python-object type to the C type;
 *     (note: for RANK > 0 types T: call type_check using the linalg::tensor_traits<T>::scalar_type* ).
 *   - for simple_object_base* this checks whether the type of the python-object is in the subset of types
 *     representable as a simple_object.
 *   . 
 */
template <class T>
bool type_check(const PyObject* src); 


/**
 * @brief return the python/numpy data-type enum corresponding to type T.
 */
template <class T>
int datatype_enum(void);


template <class T> struct python_array_flags;
/**
 * @brief the python array initialization-flags corresponding to type T (T: mutable type)
 */
template <class T>
struct python_array_flags{
  static const int flags = NPY_CARRAY;
};

/**
 * @brief the python array initialization-flags corresponding to type T (T: immutable type)
 */
template <class T>
struct python_array_flags<const T>{
  static const int flags = NPY_CARRAY_RO;
};

/**
 * @brief: Convert type T to python object.
 *   @param[in]  t: value to be converted
 *   @return  new reference to the python object
 */
template <class T>
PyObject* insert(const T& t) throw(python_error);

/**
 * @brief: Convert a pointer to type T to python/numpy array object.
 *   A direct in-memory correspondence must exist between type T and python array data-type (see "type_check").
 *   @param[in]  ptr  value to be converted
 *   @param[in]  len  number of elements
 *   @param[in] transfer_ownership  
 *     - true: register a deallocation method with the python object (using "delete[]");
 *     - false: copy the object
 *     .
 *   @return  new reference to the python object
 *   Usage note: shape specification of multidimensional arrays should be completed outside of this method.
 */
template <class T>
PyObject* insert(T* ptr, size_t len, bool transfer_ownership=false) throw(python_error, std::runtime_error);

/**
 * @brief Convert an ntuple<T,DIM> to a python object as a python array.
 */
template <class T,size_t DIM>
PyObject* insert(const linalg::ntuple<T,DIM>& src) throw(python_error, std::runtime_error);

/**
 * @brief: Convert simple_object_base* to python object, using the specified number types (const version).
 *   WARNING: Any python arrays referencing RANK>0 number objects in to the incoming simple_object will have read-write access:
 *     data ownership will be transferred to the python object.
 *   @tparam C: complex number type
 *   @tparam R: real number type
 *   @tparam Z: integer number type
 *
 *   @param[in]  src: value to be converted
 *   @return  new reference to the python object
 */
template <class C, class R, class Z>
PyObject* insert(const simple_object_base* src) throw(python_error);

#if 0 
// -------------- not finished yet: ------------------------------------
/**
 * @brief: Convert STL or gmm-interfaced container type to python object.
 *   @param[in]  u: value to be converted
 *   @return  new reference to the python object
 *   converts to python/numpy array type
 */
template <class U>
PyObject* insert(const U& u) throw(python_error);
#endif

/**
 * @brief: Convert python object to type T.
 *   @param[in]  object: python object to be converted
 *   @return converted value
 */
template <class T>
T extract(const PyObject* object) throw(python_error);

/**
 * @brief: Convert a python/numpy array object to a pointer to type T.
 *   A direct in-memory correspondence must exist between the python array data-type and type T (see "type_check").
 *   @param[out]  ptr  destination pointer
 *   @param[out]  len  number of elements
 *   @param[in]  own_data  specify data ownership for destination pointer;
 *     - true: copy the object (there is no straighforward way to transfer ownership from a python object),
 *         this implies allocation with "new[]";
 *     - false: reference the object: WARNING: in this case destination pointer must not be additionally deallocated.
 *     .
 *   Usage note: shape specification for multidimensional arrays must be obtained outside of this method.
 */
template <class T>
void extract(T*& ptr, size_t& len, const PyObject* src, bool own_data) throw(python_error, std::runtime_error);


/**
 * @brief: Convert a python/numpy array object to a pointer to type T: const version.
 */
template <class T>
inline void extract(const T*& ptr, size_t& len, const PyObject* src, bool own_data) throw(python_error, std::runtime_error)
{ extract<T>(const_cast<T*&>(ptr), len, src, own_data); }

/**
 * @brief Convert a python array object to an ntuple<T,DIM>.
 */
template <class T,size_t DIM>
void extract(linalg::ntuple<T,DIM>& dest, const PyObject* src) throw(python_error, std::runtime_error);

/**
 * @brief Convert a python object to  a simple_object, using the specified number classes.
 *    Any RANK>0 number data in the created simple_object 
 *      referencing arrays in the python object will be read-only,
 *      this is enforced at the non-const simple_object_base::ptr<T>() method. 
 */
template <class C, class R, class Z>
simple_object_base* extract(const PyObject* object) throw(python_error);

#if 0 
// -------------- not finished yet: ------------------------------------
/**
 * @brief: Convert python object to STL or gmm-interfaced container type.
 *   @param[in]   object python-object to be converted
 *   @param[out]  u  destination container
 *   Converts from either python/numpy array type or from sequence or sequence of sequences.
 *   It is an error if all sub-objects in the sequences are not of type convertible to linalg_traits<U>::value_type.
 */
template <class U>
void extract(const PyObject* object, U& u) throw(python_error);
#endif

template <class C, class R, class Z>
void extract_simple_object(simple_object_base *&dest, const std::string& src,
                           const std::string* ppreamble = NULL, const std::string* ppostscript = NULL); 

template <class C, class R, class Z>
void insert_simple_object(std::string& dest, const simple_object_base *src); 


template <class C, class R, class Z>
void read_pickle(simple_object_base *&dest, const std::string& filename);

template <class C, class R, class Z>
void write_pickle(const std::string& filename, const simple_object_base *src);


std::string python_error_string(void) throw();

void unload_interpreter(void);

// ------------ full and partial template specializations: ---------------------------

template < >
bool type_check<simple_object_base*>(const PyObject* src); 

template < >
bool type_check<bool>(const PyObject* src); 

template < >
bool type_check<std::string>(const PyObject* src); 

template < >
bool type_check<size_t>(const PyObject* src); 

/**
 * @brief Convert from python object to simple_object_base*, using the default extractor for class simple_object_base: returns new pointer.
 *    See: "simple_object_base::set_default_extractor" method.
 */
template < >
simple_object_base* extract<simple_object_base*>(const PyObject* src) throw(python_error);

template < >
bool extract<bool>(const PyObject* src) throw(python_error);

template < >
std::string extract<std::string>(const PyObject* src) throw(python_error);

#if 0 
// -------------- not finished yet: ------------------------------------
/**
 * @brief: Convert python object to vector of ntuple.
 *   @param[in]   object python-object to be converted
 *   @param[out]  v  destination vector of ntuple<T,DIM>
 *   Converts from either python/numpy array type or from sequence of sequences.
 *   It is an error if all sub-objects in the sequences are not of type convertible to type T.
 */
template <class T,size_t DIM>
void extract(const PyObject* object, std::vector<ntuple<T,DIM> >& v) throw(python_error);

/**
 * @brief: Convert python object to dense_vector_ref.
 *   @param[in]   object python-object to be converted
 *   @param[out]  v  destination dense-vector-ref
 *    - converts from python/numpy array type only.
 *    - python array must have value-type with direct in-memory equivalence 
        to type linalg_traits<dense_vector_ref<PT> >::value_type;
 *    - actual shape of array is ignored, except for total length,
 *      excepting when PT is a *ntuple<T,DIM>
 *      in which case the last dimension of shape must be DIM.
 *    .
 */
template <class PT>
void extract(const PyObject* object, gmm::dense_vector_ref<PT>& v) throw(python_error);
template <class T,size_t DIM>
void extract(const PyObject* object, gmm::dense_vector_ref<linalg::ntuple<T,DIM>*>& v) throw(python_error);

/**
 * @brief: Convert python object to dense_vector_ref with associated shape information.
 *   @param[in]   object python-object to be converted
 *   @param[out]  v  destination dense-vector-ref
 *   @param[out]  shape  dimension limits
 *    - converts from python/numpy array type only.
 *    - python array must have value-type with direct in-memory equivalence 
        to type linalg_traits<dense_vector_ref<PT> >::value_type;
 *    - when PT is a *ntuple<T,DIM'>
 *      and the last dimension of shape is DIM' (here an additional "DIM" from that of the template parameter), 
 *      the returned shape will be the actual array shape[:-1].
 *    .
 */
template <class PT,size_t DIM>
void extract(const PyObject* object, gmm::dense_vector_ref<PT>& v, ntuple<size_t,DIM>& shape) throw(python_error);
template <class T,size_t DIM1,size_t DIM2>
void extract(const PyObject* object, gmm::dense_vector_ref<ntuple<T,DIM2>*>& v, ntuple<size_t,DIM1>& shape) throw(python_error);

/**
 * @brief: Convert python object to dense_matrix_ref.
 *   @param[in]   object python-object to be converted
 *   @param[out]  m  destination dense-matrix-ref
 *    - converts from python/numpy array type only.
 *    - python array must have value-type with direct in-memory equivalence 
        to type linalg_traits<dense_matrix_ref<PT> >::value_type;
 *    - actual shape of array must have length 2, excepting when PT is a *ntuple<T,DIM>
 *      in which case the shape has length 3,and the last dimension must have value DIM.
 *      (note: for conversion to/from N-dimensional arrays use the "dense_vector_ref" form, 
 *       and pass the shape explicitely as a separate parameter)
 *    .
 */
template <class PT>
void extract(const PyObject* object, gmm::dense_matrix_ref<PT>& m) throw(python_error);
template <class T,size_t DIM>
void extract(const PyObject* object, gmm::dense_matrix_ref<ntuple<T,DIM>*>& m) throw(python_error);
 
// ------------- end: not finished yet: ------------------------------------
#endif

/**
 * @brief Convert from simple_object_base* to python object, using the default inserter for class simple_object_base: returns new reference.
 *    See: "simple_object_base::set_default_inserter" method.
 *    Warning: any arrays in the python object referencing RANK>0 number objects in the source simple_object
 *      will have read/write access, provided the source object also had such access;
 *      data-ownership will be transferred to the python object if possible.
 */
template < >
PyObject* insert<simple_object_base* const>(simple_object_base* const& src) throw(python_error);

template < >
PyObject* insert<bool>(const bool& flag) throw(python_error);

template < >
PyObject* insert<std::string>(const std::string& s) throw(python_error);


#if 0 
// -------------- not finished yet: ------------------------------------

/**
 * @brief: Convert from vector of ntuple to python object.
 *   @param[in]   v  vector to be converted
 *   @return  equivalent python object
 *     - converts to python/numpy array type.
 */
template <class T,size_t DIM>
PyObject* insert(const std::vector<ntuple<T,DIM> >& v) throw(python_error);

/**
 * @brief: Convert dense_vector_ref to python object.
 *   @param[in]  v  dense-vector-ref to be converted.
 *    - converts to python/numpy array type.
 *    - array will be 1D, except in when PT is a *ntuple<T,DIM>
 *      in which case a 2D array will be created with last dimension DIM.
 *    - ownership of data will be transferred to python object, 
 *      which will handle deallocation using standard python methods
 *      (i.e. any _explicit_ deallocation of the referenced object 
 *       after using this method is a serious error).
 *    .
 */
template <class PT>
PyObject* insert(const gmm::dense_vector_ref<PT>& v) throw(python_error);
template <class T,size_t DIM>
PyObject* insert(const gmm::dense_vector_ref<ntuple<T,DIM>*>& v) throw(python_error);

/**
 * @brief: Convert dense_vector_ref to python N-dimensional array of requested shape.
 *   @param[in]  v  dense-vector-ref to be converted.
 *   @param[in] shape  dimension limits
 *    - converts to python/numpy array type.
 *    - when PT is *ntuple<T,DIM'>
 *      the specified shape will be appended with the last dimension DIM' (here DIM' is a separate "DIM" from that of the tempate parameter).
 *    - ownership of data will be transferred to python object, 
 *      which will handle deallocation using standard python methods
 *      (i.e. any _explicit_ deallocation of the referenced object 
 *       after using this method is a serious error).
 *    .
 */
template <class PT,size_t DIM>
PyObject* insert(const gmm::dense_vector_ref<PT>& v, ntuple<size_t,DIM>& shape) throw(python_error);
template <class T,size_t DIM1,size_t DIM2>
PyObject* insert(const gmm::dense_vector_ref<ntuple<T,DIM2>*>& v, ntuple<size_t,DIM1>& shape) throw(python_error);

/**
 * @brief: Convert dense_matrix_ref to python object.
 *   @param[in]  m  dense-matrix-ref to be converted
 *    - converts to python/numpy array type.
 *    - array will be 2D, except in when PT is a *ntuple<T,DIM>
 *      in which case a 3D array will be created with last dimension DIM.
 *    - ownership of data will be transferred to python object, 
 *      which will handle deallocation using standard python methods
 *      (i.e. any _explicit_ deallocation of the referenced object 
 *       after using this method is a serious error).
 *    .
 */
template <class PT>
PyObject* insert(const gmm::dense_matrix_ref<PT>& m) throw(python_error);
template <class T,size_t DIM>
PyObject* insert(const gmm::dense_matrix_ref<ntuple<T,DIM>*>& m) throw(python_error);

// -------------- end: not finished yet: ------------------------------------
#endif

//      ------------- double precision specializations: -----------------------

template < >
bool type_check<std::complex<double> >(const PyObject* src);
 
template < >
bool type_check<std::complex<double>*>(const PyObject* src); 
 
template < >
bool type_check<const std::complex<double>*>(const PyObject* src); 


template < >
bool type_check<double>(const PyObject* src);
 
template < >
bool type_check<double*>(const PyObject* src); 

template < >
bool type_check<const double*>(const PyObject* src); 


template < >
bool type_check<long>(const PyObject* src);
 
template < >
bool type_check<long*>(const PyObject* src); 
 
template < >
bool type_check<const long*>(const PyObject* src); 

 
template < >
bool type_check<size_t*>(const PyObject* src); 
 
template < >
bool type_check<const size_t*>(const PyObject* src); 


template < >
int datatype_enum<std::complex<double> >(void);

template < >
int datatype_enum<double>(void);

template < >
int datatype_enum<long>(void);

// WARNING: size_t converts to NPY_LONG, and this assumes in-memory compatibility:
template < >
int datatype_enum<size_t>(void);


template < >
PyObject* insert<std::complex<double>, double, long>(const simple_object_base* src) throw(python_error);

template < >
PyObject* insert<std::complex<double> >(const std::complex<double>& c) throw(python_error);
template < >
PyObject* insert<double>(const double& r) throw(python_error);
template < >
PyObject* insert<long>(const long& z) throw(python_error);


template < >
simple_object_base* extract<std::complex<double>, double, long>(const PyObject* object) throw(python_error);


template < >
std::complex<double> extract<std::complex<double> >(const PyObject* object) throw(python_error);

template < >
double extract<double>(const PyObject* object) throw(python_error);

template < >
long extract<long>(const PyObject* object) throw(python_error);


template <>
void extract_simple_object<std::complex<double>, double, long>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript); 

template <>
void insert_simple_object<std::complex<double>, double, long>
  (std::string& dest, const simple_object_base *src); 


//      ------------- arbitrary precision specializations: -----------------------

#if defined(__USE_MERE)
/*
 * notes: 
 *   - mere classes can co-exist with double, so these definitions must work independently;
 *   - important: at present "bignum" or its python equivalent is not used on the python side;
 *       these conversions go to and from python "complex", "float" and "long" types.
 *   .
 */


template < >
bool type_check<mere::C>(const PyObject* src); 

template < >
bool type_check<mere::C*>(const PyObject* src); 

template < >
bool type_check<mere::R>(const PyObject* src);
 
template < >
bool type_check<mere::R*>(const PyObject* src); 

template < >
bool type_check<mere::Z>(const PyObject* src);
 
template < >
bool type_check<mere::Z*>(const PyObject* src); 



template < >
PyObject* insert<mere::C, mere::R, mere::Z>(const simple_object_base* src) throw(python_error);

template < >
PyObject* insert<mere::C>(const mere::C& c) throw(python_error);

template < >
PyObject* insert<mere::R>(const mere::R& r) throw(python_error);

template < >
PyObject* insert<mere::Z>(const mere::Z& z) throw(python_error);


template < >
simple_object_base* extract<mere::C, mere::R, mere::Z>(const PyObject* object) throw(python_error);


template <>
mere::C extract<mere::C>(const PyObject* object) throw(python_error);

template <>
mere::R extract<mere::R>(const PyObject* object) throw(python_error);

template <>
mere::Z extract<mere::Z>(const PyObject* object) throw(python_error);

template <>
void extract_simple_object<mere::C, mere::R, mere::Z>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript); 

template <>
void insert_simple_object<mere::C, mere::R, mere::Z>
  (std::string& dest, const simple_object_base *src); 

#endif

}  // namespace implementation_module

}  // namespace python_util

#if !defined(EXCLUDE_TEMPLATE_BODIES)
#include <python_util_module_template.h>
#endif

#endif // __python_util_module__h
