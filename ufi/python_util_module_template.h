#if !defined(__python_util_module_template__h)
#define __python_util_module_template__h
//
// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util_module_template.h,v $

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
namespace implementation_module{


template <class T>
void python_deallocator<T>::dealloc_(python_deallocator<T> *self)
{
  // I allow this: NULL pointer corresponds to size == 0 array object:
  #if 0
  assert(NULL != self->ptr_); // NULL indicates a usage error.
  #endif
  
  if (NULL != self->ptr_){
    delete[] self->ptr_; // always use "array" new and delete
    self->ptr_ = NULL;
  }
     
  self->ob_type->tp_free((PyObject *)self);
}

/**
 * @brief unique python name for this deallocation class.
 */
template <class T>
const char* python_deallocator<T>::class_name(void)
{ return typeid(python_deallocator<T>).name(); }

//! private constructor (presently not used anywhere; all initialization is in "register_deallocator"):
template <class T>
python_deallocator<T>::python_deallocator(T* ptr)
  : ptr_(ptr)
{ }

/**
* @brief Initialize and register the python type associated with this deallocator class.
*   usage (equivalent to explicit template instantiation): 
*     in _one_ source-module only:
*       PyObject *python_deallocator<T>::class_type_object_ = python_deallocator<T>::registered_type_object();
*     other source modules (and in-fact in the header), must ensure:
*       extern template class python_deallocator<T>;
*/
template <class T>
PyObject* python_deallocator<T>::registered_type_object(void) throw(python_error)
{
  // create a stack copy and then copy to the heap; this allows
  //   initializer syntax to be used to intialize the C/API object:

  PyTypeObject class_type = {
    PyObject_HEAD_INIT(NULL)
    0,                                          /*ob_size*/
    const_cast<char*>(class_name()),                       /*tp_name*/
    sizeof(python_deallocator<T>),              /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)python_deallocator<T>::dealloc_,            /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    const_cast<char*>("python_deallocator<T> deallocation object"),           /* tp_doc */
  };  
  class_type.tp_new = PyType_GenericNew;
  
  // note: PyObject_HEAD_INIT initializes the object's ob_refcnt to 1, so it is not
  //   necessary to increment the count prior to the copy to the on-heap object.

  // allocate an object on the heap with matching contents:
  PyTypeObject* rval = PyObject_New(PyTypeObject, &PyType_Type);
  if (NULL == rval)
    throw python_error("python_deallocator<T>::registered_type_object: ");

  memcpy(reinterpret_cast<void*>(rval), reinterpret_cast<void*>(&class_type), sizeof(PyTypeObject));
  
	#if 0 // this seems to be a problem if called on the _static_ (then _followed_ by the memcpy):
  if (PyType_Ready(&class_type) < 0)
    throw python_error("python_deallocator<T>::registered_type_object: ");
	#else
  if (PyType_Ready(rval) < 0)
    throw python_error("python_deallocator<T>::registered_type_object: ");	
	#endif

	return reinterpret_cast<PyObject*>(rval);
}

/**
* @brief Implicitely create a python_deallocator and attach to the owner python-object (which must be an array):
*   @param[in] owner   owner python-object (which must be an array)
*     when this object is destroyed the dealloc_ method will be called
*   @param[in] ptr  pointer to allocated data controlled by this deallocator
*/
template <class T>
void python_deallocator<T>::attach_deallocator(PyObject *owner, T* ptr) throw(python_error)
{
  assert(!(PyArray_FLAGS(owner)&NPY_OWNDATA));
  
  PyObject *newobj = reinterpret_cast<PyObject*>(PyObject_New(python_deallocator<T>, reinterpret_cast<PyTypeObject*>(class_type_object_)));
  if (newobj == NULL)
    throw python_error("python_deallocator<T>::attach_deallocator: ");

  // initialize the unique attributes:  
  reinterpret_cast<python_deallocator<T> *>(newobj)->ptr_ = ptr;

  // assign as the base-class of the owning array 
  //  (assume that this is done with a ref-count = 1, newly-created array):
  assert(NULL == PyArray_BASE(owner));
  PyArray_BASE(owner) = newobj;
}


/**
 * @brief Check the type of a python object.
 *   - for non-pointer types, this checks whether the python-object is directly convertible to the C type;
 *   - for simple pointer-types, this checks in-memory compatibility of the python-object type to the C type;
 *   - for simple_object_base* this checks whether the type of the python-object is in the subset of types
 *     representable as a simple_object.
 *   . 
 */
template <class T>
bool type_check(const PyObject* src) 
{
  /* throw unexpected */
  throw std::runtime_error("python_util::implementation_module:: no \"type_check\" is available for specified type;\n"
                           "  if a direct in-memory mapping does not exist (e.g. mere types) transfer to/from python\n"
                           "  should be implemented using non-array sequence types only");  
}

/**
 * @brief return the python/numpy data-type enum corresponding to type T.
 */
template <class T>
int datatype_enum(void)
{ 
  /* throw unexpected */
  throw std::runtime_error("python_util::implementation_module:: no \"datatype_enum\" is available for specified type;\n"
                           "  if a direct in-memory mapping does not exist (e.g. mere types) transfer to/from python\n"
                           "  should be implemented using non-array sequence types only"); 
}

/**
 * @brief: Convert type T to python object.
 *   @param[in]  t: value to be converted
 *   @return  new reference to the python object
 */
template <class T>
PyObject* insert(const T& t) throw(python_error)
{
  return t.insert();
}

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
PyObject* insert(T* ptr, size_t len, bool transfer_ownership) throw(python_error, std::runtime_error)
{
  // the value-type (may have RANK>0):
  // (this will also be used in the following for the non-const type)
  typedef typename std::iterator_traits<T*>::value_type value_type;
  
  // the non-const pointer type:
  typedef value_type* PTR; // need the non-const pointer type
  
  // the associated number type (RANK=0):
  typedef typename linalg::tensor_traits<value_type>::scalar_type S;


  PyObject *dest(NULL);
  
  // number of dimensions for destination array:
  size_t DIM(1); 
  npy_int dims[2] = {len,0}; // worst-case number of dims (this implementation limited to RANK <= 1). 
  if (linalg::tensor_traits<value_type>::rank > 0){
    if (linalg::tensor_traits<value_type>::rank == 1){
      // transfer ntuple<T,1> as 1-dimension:
      if (linalg::tensor_traits<value_type>::N_components > 1){
        DIM += 1;
        dims[DIM-1] = linalg::tensor_traits<value_type>::N_components;
      }
    }  
    else
      throw std::runtime_error("python_util::implementation_module::insert for RANK > 1 objects not implemented");      
  }

  #if !defined(EMBEDDED_PYTHON) 
  if (!transfer_ownership){
    // default should be to create array as C-style contiguous:
    dest = PyArray_SimpleNew(DIM, reinterpret_cast<npy_intp*>(&dims), datatype_enum<S>());
    if (NULL == dest)
      throw python_error("python_util::implementation_module::insert: ");
      
    std::copy(ptr, ptr+len, reinterpret_cast<PTR>(PyArray_DATA(dest)));
  }
  else{
    // create python array with reference to incoming data pointer (if T is "const", this will be a read-only array):
    dest = PyArray_New(&PyArray_Type, DIM, reinterpret_cast<npy_intp*>(&dims), datatype_enum<S>(),
                       NULL, reinterpret_cast<void*>(const_cast<PTR>(ptr)), 0, python_array_flags<T>::flags, NULL);
    
    // dest = PyArray_SimpleNewFromData(DIM, reinterpret_cast<npy_intp*>(&dims), datatype_enum<S>(), ptr);
    if (NULL == dest)
      throw python_error("python_util::implementation_module::insert: ");
      
    // attach the deallocation method (use "value_type", as "T" may be const):
    python_deallocator<value_type>::attach_deallocator(dest, const_cast<PTR>(ptr));  
  }
  #else
  // ------------- embedded python: multiarray.so module load currently not working: -------------------
  // create corresponding python object as tuple:
  dest = PyTuple_New(dims[0]);
  if (NULL == dest)
    throw python_error("python_util::implementation_module::insert: ");

  if (DIM==1){
    int n(0);
    for(T *it = ptr, *itEnd = ptr + len;
        it != itEnd;
        ++it, ++n){
      PyObject *item_ = insert(*it);

      // "PyTuple_SetItem" steals ref:
      if (PyTuple_SetItem(dest, n, item_)){
        Py_DECREF(dest);
        Py_DECREF(item_);
        throw python_error("python_util::implementation_module::insert: "); 
      }      
    }
  }
  else{
    assert(DIM==2); // at present, only rank==0, 1 are supported:
    int n(0);
    for(T *it = ptr, *itEnd = ptr + len;
        it != itEnd;
        ++it, ++n){
      PyObject *item_ = insert(it, dims[1], transfer_ownership);

      // "PyTuple_SetItem" steals ref:
      if (PyTuple_SetItem(dest, n, item_)){
        Py_DECREF(dest);
        Py_DECREF(item_);
        throw python_error("python_util::implementation_module::insert: "); 
      }      
    }
  }  
  #endif

  return dest;
}

/**
 * @brief Convert an ntuple<T,DIM> to a python object as a python array.
 *   WARNING: python value_type (see "datatype_enum<T>") must be in-memory compatible to C type; no conversion occurs.
 */
template <class T,size_t DIM>
PyObject* insert(const linalg::ntuple<T,DIM>& src) throw(python_error, std::runtime_error)
{ return insert(src.begin(), DIM, false); }

/**
 * @brief: Convert python object to type T.
 *   @param[in]  object: python object to be converted
 *   @return converted value
 */
template <class T>
T extract(const PyObject* src) throw(python_error)
{
  T rval;
  rval.extract(src);
  return rval;
}

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
// implementation note: "own_data" changed from optional "own_data=true" because of ambiguity between "false" and NULL. 
template <class T>
void extract(T*& ptr, size_t& len, const PyObject* src_, bool own_data) throw(python_error, std::runtime_error)
{
  // python/C api does not use "const:
  PyObject *src(const_cast<PyObject*>(src_));
  
  assert((NULL == ptr) && (NULL != src)); // usage-error trap.

  // distinguish between type "T" and its non-const type:  
  typedef typename std::iterator_traits<T*>::value_type value_type;  // non-const type (possibly RANK>0)
  typedef value_type *PTR;                                           // non-const pointer
  typedef typename linalg::tensor_traits<value_type>::scalar_type S; // non-const RANK=0 number type
  
  // at this point (with respect to the simple_object implementation), non-owned data-references
  //   are read-only, and so it suffices to call type_check for the lesser requirement, which is that
  //   the python array be at least read-only:
  if (type_check<const S*>(src)){
  
    if (linalg::tensor_traits<value_type>::rank > 0){
      if (linalg::tensor_traits<value_type>::rank == 1){
        if (!(((linalg::tensor_traits<value_type>::N_components == 1) && (static_cast<size_t>(PyArray_NDIM(src)) == 1))
              || (linalg::tensor_traits<value_type>::N_components == static_cast<size_t>(PyArray_DIM(src,PyArray_NDIM(src)-1)))))
          throw std::runtime_error("python_util::implementation_module::extract: RANK=1: \n"
            "  array and non-scalar must either both be 1D, or last dimension of python array must match N_components of non-scalar");   
      }
      else
        throw std::runtime_error("python_util::implementation_module::extract for RANK > 1 objects not implemented");
    }
  
    #if !defined(EMBEDDED_PYTHON)
    if (own_data){
      const PTR data(reinterpret_cast<PTR>(PyArray_DATA(src)));
      if (PyErr_Occurred())
        throw python_error("python_util::implementation_module::extract: ");
        
      len = static_cast<size_t>(PyArray_SIZE(src));
      if (PyErr_Occurred())
        throw python_error("python_util::implementation_module::extract: ");
        
      if (len % linalg::tensor_traits<value_type>::N_components != 0)
        throw std::runtime_error("python_util::implementation_module::extract: RANK > 0 destination:\n"
          "  number of components in destination incompatible with number of elements in python array");  
      len /= linalg::tensor_traits<value_type>::N_components;
      
      ptr = new value_type[len];
      
      std::copy(data, data+len, ptr);
    }
    else{
      ptr = reinterpret_cast<T*>(PyArray_DATA(src));
      if (PyErr_Occurred())
        throw python_error("python_util::implementation_module::extract: ");
        
      len = static_cast<size_t>(PyArray_SIZE(src));
      if (PyErr_Occurred())
        throw python_error("python_util::implementation_module::extract: ");
      len /= linalg::tensor_traits<value_type>::N_components;
    }  
    #else
    // import of "multiarray.so" currently not working, no option but to copy from the python sequence:
    assert(own_data); 

    if (PySequence_Check(src)){  
        int size_(PySequence_Size(src));
        if (-1 == size_)
          throw python_error("python_util::implementation_module::extract: ");
        ptr = new value_type[size_];
        
        if (linalg::tensor_traits<value_type>::N_components == 1){
          int n(0);
          for(value_type *it = ptr, *itEnd = ptr + size_;
              it != itEnd;
              ++it, ++n){
            // "PySequence_GetItem" returns new ref:
            PyObject *item_ = PySequence_GetItem(src, n);
            if (NULL == item_)
              throw python_error("python_util::implementation_module::extract: ");
            *it = extract<value_type>(item_);
            Py_DECREF(item_);
          }    
        }
        else{
          int n(0);
          for(value_type *it = ptr, *itEnd = ptr + size_;
              it != itEnd;
              ++it, ++n){
            // "PySequence_GetItem" returns new ref:
            PyObject *item_ = PySequence_GetItem(src, n);
            if (NULL == item_)
              throw python_error("python_util::implementation_module::extract: ");
            
            if (PySequence_Check(item_)){
              if (static_cast<size_t>(PySequence_Size(item_)) != linalg::tensor_traits<value_type>::N_components){
                Py_DECREF(item_);
                throw python_error("python_util::implementation_module::extract: sub-item size does not match RANK 1 size");              
              }
              
              int S_n(0);
              for(S *S_it = reinterpret_cast<S*>(it), *S_itEnd = S_it + linalg::tensor_traits<value_type>::N_components;
                  S_it != S_itEnd;
                  ++S_it, ++S_n){

                // "PySequence_GetItem" returns new ref:
                PyObject *S_item_ = PySequence_GetItem(item_, n);
                if (NULL == S_item_)
                  throw python_error("python_util::implementation_module::extract: ");
                *S_it = extract<S>(S_item_);
                Py_DECREF(S_item_);
              }                
            }
            else{
              Py_DECREF(item_);
              throw python_error("python_util::implementation_module::extract: sub-item is not a sequence");
            }
                        
            Py_DECREF(item_);
          }            
        }
    }
    else 
      throw std::string("python_util::implementation_module::extract: python object is not a sequence"); 
    #endif
    
  }
  else{
    #if 1
    // *** DEBUG ***
    std::cerr<<"array-check:   "<<PyArray_Check(src)<<"\n"
             <<"  ISCARRAY_RO: "<<PyArray_ISCARRAY_RO(src)<<"\n"
             <<"  ISFARRAY_RO: "<<PyArray_ISFARRAY_RO(src)<<"\n"
             <<"  ISCARRAY:    "<<PyArray_ISCARRAY(src)<<"\n"
             <<"  ISFARRAY:    "<<PyArray_ISFARRAY(src)<<"\n"
             <<"  datatype_enum:       "<<datatype_enum<double>()<<"\n"
             <<"  array-datatype-enum: "<<PyArray_TYPE(src)<<"\n"<<std::endl;
    #endif
    throw std::runtime_error("python_util::implementation_module::extract: python object has no direct in-memory mapping to target type");
  }
}

/**
 * @brief Convert a python array object to an ntuple<T,DIM>.
 *    - Data types must have in-memory compatibility: no conversion occurs (see "datatype_enum<T>");
 *    - array must be either 1D, or have 1 x DIM, or DIM x 1 dimensions.
 *    .
 */
template <class T,size_t DIM>
void extract(linalg::ntuple<T,DIM>& dest, const PyObject* src) throw(python_error, std::runtime_error)
{
  if (type_check<const T*>(src) 
      && (DIM == PyArray_SIZE(src))
      && ( (PyArray_NDIM(src) == 1) 
           || ((PyArray_NDIM(src) == 2) && (PyArray_DIM(src,0) == 1)) 
           || ((PyArray_NDIM(src) == 2) && (PyArray_DIM(src,1) == 1)))){
    const T *data(reinterpret_cast<const T*>(PyArray_DATA(src)));
    if (PyErr_Occurred())
      throw python_error("python_util::implementation_module::extract: ");

    std::copy(data, data+DIM, dest.begin());
  }
  else
    throw std::runtime_error("python_util::implementation_module::extract: python object has no direct in-memory mapping to ntuple<T,DIM>");
}

} // namespace implementation_module
} // namespace python_util


// additional STL / gmm-interfaced insert and extract methods:
#include "python_util_module_gmm_template.h"


#endif // __python_util_module_template__h
