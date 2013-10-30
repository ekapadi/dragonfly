// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util_module.cpp,v $

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


#include <portability.h>

#include <pthread.h>

#if defined(__USE_MPI)
  #include <mpi.h>
#endif

#if defined(__USE_OMP) 
  #include <omp.h>
#endif

#if 1 
#include <Python.h>
#if !defined(EMBEDDED_PYTHON)
#define NO_IMPORT_ARRAY
#else
// turned _off_ for the moment ("multiarray.so" import problems):
// #define PY_ARRAY_UNIQUE_SYMBOL PyArray_API
#endif

#include <numpy/arrayobject.h>
// #include <numpy/ufuncobject.h>
// #include <numpy/__multiarray_api.h>
// #include <numpy/__ufunc_api.h>
#endif 

#include <assert.h>

#include <stdexcept>

#include <iostream>
#include <sstream>
#include <iomanip>
using std::cout;
using std::endl;

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

#if 0
#include <gmm/gmm_kernel.h>
#include <gmm/gmm_dense_lu.h>
#include <vector_view.h>
#include <matrix_view.h>
#endif

#include <ntuple.h>
using linalg::ntuple;

#if 0
using linalg::ntuple_interval;
#endif

#include <scanner_util.h>
#include <python_util.h>


#if 0
// =========================================== BOOST::PYTHON version ===================================================================================

#include <boost/python.hpp>


/*
http://stackoverflow.com/questions/1008343/boost-python-and-python-exceptions

The proper way to do this under Boost.Python is not to use the C API directly, but install an exception translator: register_exception_translator<RuntimeException>(my_runtime_exception_translator); void my_runtime_exception_translator(RuntimeException const& ex) { PyErr_SetString(PyExc_RuntimeError, ex.toString().c_str()); }
 */

using namespace boost::python;

namespace python_util{
namespace implementation_module{

std::string python_error_string(void) throw()
{
  // assume interpreter is initialized:
  std::string msg;
  try{
    if (NULL != PyErr_Occurred()){

      // save (and clear) the error indicator:
      PyObject *type_save_(NULL), *value_save_(NULL), *traceback_save_(NULL);

      // "PyErr_Fetch" returns new refs (any of which may be NULL, but not all):
      PyErr_Fetch(&type_save_, &value_save_, &traceback_save_);

      // "PyImport_AddModule" returns borrowed ref:
      // "PyModule_GetDict" returns borrowed ref:

      // borrowed refs to global and local namespace dict:
      PyObject 
        *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
        *locals_ = PyEval_GetLocals();

      PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

      // re-direct stderr to string-stream:
      rval_ = PyRun_String("import sys; import StringIO; stderr_save = sys.stderr; sys.stderr = StringIO.StringIO()", Py_file_input, globals_, locals_);
      if (NULL == rval_){
        Py_XDECREF(type_save_); 
        Py_XDECREF(value_save_);
        Py_XDECREF(traceback_save_);
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      }
      Py_DECREF(rval_);

      // restore the error indicator (consumes the refs):
      PyErr_Restore(type_save_, value_save_, traceback_save_);

      // print the error to the redirected stderr, re-clear the error indicator:
      PyErr_Print();

      rval_ = PyRun_String("msg_ = sys.stderr.getvalue()", Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      Py_DECREF(rval_);

      // "PyDict_GetItemString" returns borrowed ref:
      PyObject* msg_ = PyDict_GetItemString(globals_, "msg_");
      if (NULL == msg_) 
        throw std::runtime_error("python_error_string: PyDict_GetItemString error return");// throw unexpected...  
      msg = extract<std::string>(msg_);

      // restore stderr to its original stream:
      rval_ = PyRun_String("sys.stderr.close(); sys.stderr = stderr_save", Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      Py_DECREF(rval_); // usually "None"
    }
    else
      msg = "python error-state NOT set";
  }
  catch(std::string& new_error){
    msg = "python_error_string: unable to extract error information from python C/API:\n  "
    msg += new_error;
  }
  
  return msg;
}


void unload_interpreter(void){
  if (Py_IsInitialized())
    Py_Finalize();
}


// file-scope only prototypes:
void extract_simple_object_(simple_object_base *&dest, const object &src); 
void insert_simple_object_(object& dest, const simple_object_base *src);

// ------------ explicit template specializations: ---------------------------

// -------------- _always_ define double-precision versions: -----------------
template <>
void extract_simple_object<std::complex<double>, double, long>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // evaluate the string[s] in the global namespace:
    if (ppreamble != NULL)
      exec(ppreamble->c_str(), global, global);
    std::string exec_str(
      "try:\n"
      "  object = " + src + "\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    exec(exec_str.c_str(), global, global);
    if (ppostscript != NULL)
      exec(ppostscript->c_str(), global, global);    
    object extractable = global["object"];
    #if 0
    // *** DEBUG ***
    exec("obj_str = repr(object)",global,global);
    std::string stmp = extract<std::string>(global["obj_str"]);
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    extract_simple_object_(dest, extractable);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("python_util::extract_simple_object:");
  }
  catch (...)
  {
    cout<<"unknown python exception: error string:"<<endl;
    cout<<"    "<<python_error_string()<<endl;
  }    
} 
    
void extract_simple_object_(simple_object_base *&dest, const object &src) 
{
  using boost::python::extract;
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  

  PyObject *src_(src.ptr());
  if (PyBool_Check(src_)){
    dest = new simple_object<bool>(extract<bool>(src));
  }    
  else
  if (PyNumber_Check(src_)){
    if (PyComplex_Check(src_))
      dest = new simple_object<C>(extract<C>(src));
    else
    if (PyFloat_Check(src_))
      dest = new simple_object<R>(extract<R>(src));
    else
    if (PyInt_Check(src_))
      dest = new simple_object<Z>(extract<Z>(src));
  }
  else
  if (PySequence_Check(src_)){  
    if (PyString_Check(src_))
      dest = new simple_object<std::string>(extract<std::string>(src));
    else
    if (PyTuple_Check(src_) || PyList_Check(src_)){
      dest = new simple_object<LIST>();
      for(size_t n = 0; n < static_cast<size_t>(len(src)); ++n){
        simple_object_base *dest_;
        extract_simple_object_(dest_, src[n]);
        dest->as<LIST>().push_back(dest_);
      }  
    }  
  }
  else
  if (PyDict_Check(src_)){
    // extract dictionary with _string_ key and arbitrary value:
    dest = new simple_object<MAP>();
    std::string key;
    simple_object_base *val;
    object dict_items(src.attr("items")());
    for(size_t nitem = 0; nitem < static_cast<size_t>(len(dict_items)); ++nitem){
      object item(dict_items[nitem]);
      key = extract<std::string>(item[0]);
      extract_simple_object_(val, item[1]);
      dest->as<MAP>()[key] = val;
    }    
  }
  else{
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));
    global["bad_object"] = src;
    exec("obj_str = str(bad_object)",global,global);
    std::string obj_str = extract<std::string>(global["obj_str"]);
    throw python_error("extract_simple_object: object type not in { <number>, <sequence>, <dictionary> }: ") + obj_str;
  }
} 


template <>
void insert_simple_object<std::complex<double>, double, long>
  (std::string& dest, const simple_object_base *src)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // convert simple_object_base to boost::python object:
    object dest_;
    insert_simple_object_(dest_, src);

    // get the string-representation and convert to a std::string:
    global["object"] = dest_;
    exec("obj_str = repr(object)",global,global);
    dest = extract<std::string>(global["obj_str"]);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("insert_simple_object:");
  }  
}  

void insert_simple_object_(object& dest, const simple_object_base *src)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  

  //  PyObject *src_(src.ptr());
  if (src->is<bool>())
    dest = object(src->as<bool>());
  else
  if (src->is<C>())
    dest = object(src->as<C>());
  else
  if (src->is<R>())
    dest = object(src->as<R>());
  else
  if (src->is<Z>())
    dest = object(src->as<Z>());
  else
  if (src->is<std::string>())
    dest = object((src->as<std::string>()).c_str());
  else
  if (src->is<LIST>()){
    object dest_ = list(), obj;
    const LIST &l(src->as<LIST>());
    for(LIST::const_iterator itl = l.begin(), itlEnd = l.end();
        itl != itlEnd;
        ++itl){
      insert_simple_object_(obj, *itl);
      dest_.attr("append")(obj);   
    }
    dest = tuple(dest_);  // repn as tuple rather than a list    
  }
  else
  if (src->is<MAP>()){
    dest = dict();
    std::string key;
    object val;
    const MAP &m(src->as<MAP>());
    for(MAP::const_iterator itm = m.begin(), itmEnd = m.end();
        itm != itmEnd;
        ++itm){
      key = (*itm).first;  
      insert_simple_object_(val, (*itm).second);
      dest[key] = val;   
    }   
  }
  else  
    throw python_error("insert_simple_object: object type not in { std::complex<double>, double, long, std::string, object_list, object_map }");
} 


template < >
void read_pickle<std::complex<double>, double, long>
  (simple_object_base *&dest, const std::string& filename)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // evaluate the string in the global namespace:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  infile = open('" + filename + "','rb')\n"
      "  object = PK.load(infile)\n"
      "  infile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    exec(exec_str.c_str(), global, global);
    object extractable = global["object"];
    #if 0
    // *** DEBUG ***
    exec("obj_str = repr(object)",global,global);
    std::string stmp = extract<std::string>(global["obj_str"]);
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    extract_simple_object_(dest, extractable);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("python_util::read_pickle:");
  }
  catch (...)
  {
    cout<<"python_util::read_pickle:  unknown python exception: error string:"<<endl;
    cout<<"    "<<python_error_string()<<endl;
  } 
}

template < >
void write_pickle<std::complex<double>, double, long>
  (const std::string& filename, const simple_object_base *src)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // convert simple_object_base to boost::python object:
    object src_;
    insert_simple_object_(src_, src);

    global["object"] = src_;
    
    // get the string-representation and convert to a std::string:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  outfile = open('" + filename + "','wb')\n"
      "  PK.dump(object, outfile, PK.HIGHEST_PROTOCOL)\n"
      "  outfile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");                  
    exec(exec_str.c_str(),global,global);
    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("python_util::write_pickle:");
  }  
}


#if defined(__USE_MERE) && 0 // ------------------ temporarily turned-off --------------------------------

// -------------------------- arbitrary-precision versions: ----------------------------------------------
// ---------------- note: these definitions may co-exist with the double-precision definitions -----------

template <>
void extract_simple_object<mere::C, mere::R, mere::Z>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // evaluate the string[s] in the global namespace:
    if (ppreamble != NULL)
      exec(ppreamble->c_str(), global, global);
    exec((std::string("object = ") + src).c_str(), global, global);
    if (ppostscript != NULL)
      exec(ppostscript->c_str(), global, global);    
    object extractable = global["object"];
    #if 0
    // *** DEBUG ***
    exec("obj_str = repr(object)",global,global);
    std::string stmp = extract<std::string>(global["obj_str"]);
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    extract_simple_object_(dest, extractable);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("extract_simple_object:");
  }  
} 

void extract_simple_object_(simple_object_base *&dest, const object &src) 
{
  using boost::python::extract;
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  

  PyObject *src_(src.ptr());
  if (PyBool_Check(src_)){
    dest = new simple_object<bool>(extract<bool>(src)); 
  }   
  else
  if (PyNumber_Check(src_)){
    if (PyComplex_Check(src_))
      dest = new simple_object<C>(extract<C>(src));
    else
    if (PyFloat_Check(src_))
      dest = new simple_object<R>(extract<R>(src));
    else
    if (PyInt_Check(src_))
      dest = new simple_object<Z>(extract<Z>(src));
  }
  else
  if (PySequence_Check(src_)){  
    if (PyString_Check(src_))
      dest = new simple_object<std::string>(extract<std::string>(src));
    else
    if (PyTuple_Check(src_) || PyList_Check(src_)){
      dest = new simple_object<LIST>();
      for(size_t n = 0; n < static_cast<size_t>(len(src)); ++n){
        simple_object_base *dest_;
        extract_simple_object_(dest_, src[n]);
        dest->as<LIST>().push_back(dest_);
      }  
    }  
  }
  else
  if (PyDict_Check(src_)){
    // extract dictionary with _string_ key and arbitrary value:
    dest = new simple_object<MAP>();
    std::string key;
    simple_object_base *val;
    object dict_items(src.attr("items")());
    for(size_t nitem = 0; nitem < static_cast<size_t>(len(dict_items)); ++nitem){
      object item(dict_items[nitem]);
      key = extract<std::string>(item[0]);
      extract_simple_object_(val, item[1]);
      dest->as<MAP>()[key] = val;
    }    
  }
  else{
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));
    global["bad_object"] = src;
    exec("obj_str = str(bad_object)",global,global);
    std::string obj_str = extract<std::string>(global["obj_str"]);
    throw python_error("extract_simple_object: object type not in { <number>, bool, <sequence>, <dictionary> }: ") + obj_str;
  }
} 


template <>
void insert_simple_object<mere::C, mere::R, mere::Z>
  (std::string& dest, const simple_object_base *src)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // convert simple_object_base to boost::python object:
    object dest_;
    insert_simple_object_(dest_, src);

    // get the string-representation and convert to a std::string:
    global["object"] = dest_;
    exec("obj_str = repr(object)",global,global);
    dest = extract<std::string>(global["obj_str"]);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("insert_simple_object:");
  }  
}  

void insert_simple_object_(object& dest, const simple_object_base *src)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  

  //  PyObject *src_(src.ptr());
  if (src->is<bool>())
    dest = object(src->as<bool>());
  else
  if (src->is<C>())
    dest = object(src->as<C>());
  else
  if (src->is<R>())
    dest = object(src->as<R>());
  else
  if (src->is<Z>())
    dest = object(src->as<Z>());
  else
  if (src->is<std::string>())
    dest = object((src->as<std::string>()).c_str());
  else
  if (src->is<LIST>()){
    object dest_ = list(), obj;
    const LIST &l(src->as<LIST>());
    for(LIST::const_iterator itl = l.begin(), itlEnd = l.end();
        itl != itlEnd;
        ++itl){
      insert_simple_object_(obj, *itl);
      dest_.attr("append")(obj);   
    }
    dest = tuple(dest_);  // repn as tuple rather than a list    
  }
  else
  if (src->is<MAP>()){
    dest = dict();
    std::string key;
    object val;
    const MAP &m(src->as<MAP>());
    for(MAP::const_iterator itm = m.begin(), itmEnd = m.end();
        itm != itmEnd;
        ++itm){
      key = (*itm).first;  
      insert_simple_object_(val, (*itm).second);
      dest[key] = val;   
    }   
  }
  else  
    throw python_error("insert_simple_object: object type not in { std::complex<double>, double, long, bool, std::string, object_list, object_map }");
}



template < >
void read_pickle<mere::C, mere::R, mere::Z>
  (simple_object_base *&dest, const std::string& filename)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // evaluate the string in the global namespace:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  infile = open('" + filename + "','rb')\n"
      "  object = PK.load(infile)\n"
      "  infile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    exec(exec_str.c_str(), global, global);
    object extractable = global["object"];
    #if 0
    // *** DEBUG ***
    exec("obj_str = repr(object)",global,global);
    std::string stmp = extract<std::string>(global["obj_str"]);
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    extract_simple_object_(dest, extractable);

    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("python_util::read_pickle:");
  }
  catch (...)
  {
    cout<<"python_util::read_pickle:  unknown python exception: error string:"<<endl;
    cout<<"    "<<python_error_string()<<endl;
  } 
}

template < >
void write_pickle<mere::C, mere::R, mere::Z>
  (const std::string& filename, const simple_object_base *src)
{
  using boost::python::extract;
  try{
    if (!Py_IsInitialized())
      Py_Initialize();
    // get global namespace;
    object main_module = import(const_cast<char*>("__main__"));
    object global(main_module.attr(const_cast<char*>("__dict__")));

    // convert simple_object_base to boost::python object:
    object src_;
    insert_simple_object_(src_, src);

    global["object"] = src_;

    // get the string-representation and convert to a std::string:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  outfile = open('" + filename + "','wb')\n"
      "  PK.dump(object, outfile, PK.HIGHEST_PROTOCOL)\n"
      "  outfile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");            
    exec(exec_str.c_str(),global,global);
    #if 0 // try this later 
    Py_Finalize();
    #endif
  }
  catch(error_already_set const &)
  {
    throw python_error("write_pickle:");
  }  
}

 
#endif

} // namespace implementation_module

} // namespace python_util

// =========================================== end: BOOST::PYTHON version ==============================================================================
#else
// ====================================== python C/API only version ====================================================================================

namespace python_util{
namespace implementation_module{

std::string python_error_string(void) throw()
{
  // assume interpreter is initialized:
  std::string msg;
  
  try{
    if (NULL != PyErr_Occurred()){

      // save (and clear) the error indicator:
      PyObject *type_save_(NULL), *value_save_(NULL), *traceback_save_(NULL);

      // "PyErr_Fetch" returns new refs (any of which may be NULL, but not all):
      PyErr_Fetch(&type_save_, &value_save_, &traceback_save_);

      // "PyImport_AddModule" returns borrowed ref:
      // "PyModule_GetDict" returns borrowed ref:

      // borrowed refs to global and local namespace dict:
      PyObject 
        *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
        *locals_ = PyEval_GetLocals();

      PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

      // re-direct stderr to string-stream:
      rval_ = PyRun_String("import sys; import StringIO; stderr_save = sys.stderr; sys.stderr = StringIO.StringIO()", Py_file_input, globals_, locals_);
      if (NULL == rval_){
        Py_XDECREF(type_save_); 
        Py_XDECREF(value_save_);
        Py_XDECREF(traceback_save_);
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      }
      Py_DECREF(rval_);

      // restore the error indicator (consumes the refs):
      PyErr_Restore(type_save_, value_save_, traceback_save_);

      // print the error to the redirected stderr, re-clear the error indicator:
      PyErr_Print();

      rval_ = PyRun_String("msg_ = sys.stderr.getvalue()", Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      Py_DECREF(rval_);

      // "PyDict_GetItemString" returns borrowed ref:
      PyObject* msg_ = PyDict_GetItemString(globals_, "msg_");
      if (NULL == msg_) 
        throw std::runtime_error("python_error_string: PyDict_GetItemString error return"); // throw unexpected...  
      msg = extract<std::string>(msg_);

      // restore stderr to its original stream:
      rval_ = PyRun_String("sys.stderr.close(); sys.stderr = stderr_save", Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::runtime_error("python_error_string: PyRun_String error return"); // throw unexpected...
      Py_DECREF(rval_); // usually "None"
    }
    else
      msg = "python error-state NOT set";
  }
  catch(std::string& new_error){
    msg = "python_error_string: unable to extract error information from python C/API:\n  ";
    msg += new_error;
  }
  
  return msg;
}

void unload_interpreter(void){
  if (Py_IsInitialized())
    Py_Finalize();
}


// ------------ explicit template specializations: ---------------------------

// -------------- _always_ define double-precision versions: -----------------
template <>
void extract_simple_object<std::complex<double>, double, long>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  
  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

    // evaluate the string[s] in the global namespace:
    if (ppreamble != NULL){
      rval_ = PyRun_String(ppreamble->c_str(), Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
      Py_DECREF(rval_);
    }
    std::string exec_str(
      "try:\n"
      "  object = " + src + "\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); // usually "None"


    if (ppostscript != NULL){
      rval_ = PyRun_String(ppostscript->c_str(), Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
      Py_DECREF(rval_);
    }  

    // "PyDict_GetItemString" returns borrowed ref:
    PyObject* extractable_ = PyDict_GetItemString(globals_, "object");
    #if 0
    // *** DEBUG ***
    rval_ = PyRun_String("obj_str = repr(object)", Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_);
    std::string stmp = extract<std::string>(PyDict_GetItemString(globals_, "obj_str");
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    dest = extract<C,R,Z>(extractable_);
    
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
} 

template <>
void insert_simple_object<std::complex<double>, double, long>
  (std::string& dest, const simple_object_base *src)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  
  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

    // convert simple_object_base to python object:
    PyObject* dest_ = insert<C,R,Z>(src);

    // get the string-representation and convert to a std::string:
    if (PyDict_SetItemString(globals_, "object", dest_)){
      Py_DECREF(dest_);
      throw std::string("insert_simple_object<C,R,Z>: PyDict_SetItemString error return");
    }
    
    rval_ = PyRun_String("obj_str = repr(object)", Py_file_input, globals_, locals_);
    if (NULL == rval_){
      Py_DECREF(dest_);
      throw std::string("insert_simple_object<C,R,Z>: PyRun_String error return");
    }
    Py_DECREF(rval_);

    // borrowed ref to "obj_str":
    PyObject* obj_str_ = PyDict_GetItemString(globals_, "obj_str");
    if (NULL == obj_str_)
      throw std::string("insert_simple_object<C,R,Z>: PyDict_GetItemString error return");

    dest = extract<std::string>(obj_str_);
    
    Py_DECREF(dest_);
   
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}  

template < >
void read_pickle<std::complex<double>, double, long>
  (simple_object_base *&dest, const std::string& filename)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  
  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.
  
    // evaluate the string in the global namespace:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  infile = open('" + filename + "','rb')\n"
      "  object = PK.load(infile)\n"
      "  infile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("read_pickle<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); // usually "None"

    // "PyDict_GetItemString" returns borrowed ref:
    PyObject* extractable = PyDict_GetItemString(globals_, "object");
    if (NULL == extractable)
      throw std::string("read_pickle<C,R,Z>: PyDict_GetItemString error return");
 
    dest = extract<C,R,Z>(extractable);

    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}

template < >
void write_pickle<std::complex<double>, double, long>
  (const std::string& filename, const simple_object_base *src)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;

  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

    // convert simple_object_base* to python object:
    PyObject* src_ = insert<C,R,Z>(src);

    if (PyDict_SetItemString(globals_, "object", src_)){
      Py_DECREF(src_);
      throw std::string("write_pickle<C,R,Z>: PyDict_SetItemString error return");
    }
    
    // get the string-representation and convert to a std::string:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  outfile = open('" + filename + "','wb')\n"
      "  PK.dump(object, outfile, PK.HIGHEST_PROTOCOL)\n"
      "  outfile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");                      
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("write_pickle<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); 
    
    Py_DECREF(src_);
    
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}


#if defined(__USE_MERE)

// -------------------------- arbitrary-precision versions: ----------------------------------------------
// ---------------- note: these definitions may co-exist with the double-precision definitions -----------

template <>
void extract_simple_object<mere::C, mere::R, mere::Z>
  (simple_object_base *&dest, const std::string& src,
   const std::string* ppreamble, const std::string* ppostscript)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  
  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.
    
    // evaluate the string[s] in the global namespace:
    if (ppreamble != NULL){
      rval_ = PyRun_String(ppreamble->c_str(), Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
      Py_DECREF(rval_);
    }
    std::string exec_str(
      "try:\n"
      "  object = " + src + "\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); // usually "None"


    if (ppostscript != NULL){
      rval_ = PyRun_String(ppostscript->c_str(), Py_file_input, globals_, locals_);
      if (NULL == rval_)
        throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
      Py_DECREF(rval_);
    }  

    // "PyDict_GetItemString" returns borrowed ref:
    PyObject* extractable_ = PyDict_GetItemString(globals_, "object");
    #if 0
    // *** DEBUG ***
    rval_ = PyRun_String("obj_str = repr(object)", Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("extract_simple_object<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_);
    std::string stmp = extract<std::string>(PyDict_GetItemString(globals_, "obj_str");
    cout<<"object in embedded interpreter: "<<stmp<<endl;
    #endif

    dest = extract<C,R,Z>(extractable_);
    
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
} 

template <>
void insert_simple_object<mere::C, mere::R, mere::Z>
  (std::string& dest, const simple_object_base *src)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  
  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

    // convert simple_object_base to python object:
    PyObject* dest_ = insert<C,R,Z>(src);

    // get the string-representation and convert to a std::string:
    if (PyDict_SetItemString(globals_, "object", dest_)){
      Py_DECREF(dest_);
      throw std::string("insert_simple_object<C,R,Z>: PyDict_SetItemString error return");
    }
    
    rval_ = PyRun_String("obj_str = repr(object)", Py_file_input, globals_, locals_);
    if (NULL == rval_){
      Py_DECREF(dest_);
      throw std::string("insert_simple_object<C,R,Z>: PyRun_String error return");
    }
    Py_DECREF(rval_);

    // borrowed ref to "obj_str":
    PyObject* obj_str_ = PyDict_GetItemString(globals_, "obj_str")
    if (NULL == obj_str_)
      throw std::string("insert_simple_object<C,R,Z>: PyDict_GetItemString error return");

    dest = extract<std::string>(obj_str_);
    
    Py_DECREF(dest_);
   
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}  

template < >
void read_pickle<mere::C, mere::R, mere::Z>
  (simple_object_base *&dest, const std::string& filename)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;

  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.
  
    // evaluate the string in the global namespace:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  infile = open('" + filename + "','rb')\n"
      "  object = PK.load(infile)\n"
      "  infile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("read_pickle<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); // usually "None"

    // "PyDict_GetItemString" returns borrowed ref:
    PyObject* extractable = PyDict_GetItemString(globals_, "object");
    if (NULL == extractable)
      throw std::string("read_pickle<C,R,Z>: PyDict_GetItemString error return");
 
    dest = extract<C,R,Z>(extractable);

    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}

template < >
void write_pickle<mere::C, mere::R, mere::Z>
  (const std::string& filename, const simple_object_base *src)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;

  try{
    if (!Py_IsInitialized()){
      Py_Initialize();
      #if defined(EMBEDDED_PYTHON) && 0
      import_array();
      #endif
    }

    // "PyImport_AddModule" returns borrowed ref:
    // "PyModule_GetDict" returns borrowed ref:
    
    // borrowed refs to global and local namespace dict:
    PyObject 
      *globals_ = PyModule_GetDict(PyImport_AddModule(const_cast<char*>("__main__"))),
      *locals_ = NULL;

    PyObject* rval_(NULL); // return values are usually "None", but I assume the refs are owned.

    // convert simple_object_base* to python object:
    PyObject* src_ = insert<C,R,Z>(src);

    if (PyDict_SetItemString(globals_, "object", src_)){
      Py_DECREF(src_);
      throw std::string("write_pickle<C,R,Z>: PyDict_SetItemString error return");
    }
    
    // get the string-representation and convert to a std::string:
    std::string exec_str(
      "import sys\n"
      "try:\n"
      "  import pickle as PK\n"
      "  outfile = open('" + filename + "','wb')\n"
      "  PK.dump(object, outfile, PK.HIGHEST_PROTOCOL)\n"
      "  outfile.close()\n"
      "except:\n"
      "  print 'embedded python exception:', sys.exc_info()\n");                      
    rval_ = PyRun_String(exec_str.c_str(), Py_file_input, globals_, locals_);
    if (NULL == rval_)
      throw std::string("write_pickle<C,R,Z>: PyRun_String error return");
    Py_DECREF(rval_); 
    
    Py_DECREF(src_);
    
    #if 1
    Py_Finalize();
    #endif
  }
  catch(std::string& msg)
  {
    #if 0
    PyErr_Print();
    Py_Finalize();
    throw msg;
    #else
    msg += ":\n" + python_error_string();
    Py_Finalize();
    throw msg;
    #endif
  }   
}


#endif

} // namespace implementation_module

} // namespace python_util

// ====================================== end: python C/API only version ===============================================================================
#endif


// ============================= non-exclusive python C/API functions: =================================================================================
namespace python_util{
namespace implementation_module{


// =========== EXPLICIT TEMPLATE INSTANTIATIONS: =============================================

// --- these are the only supported types, other types should throw exception: ---
template class python_deallocator<std::complex<double> >;
template class python_deallocator<double >;
template class python_deallocator<long >;

template class python_deallocator<ntuple<std::complex<double>,1> >;
template class python_deallocator<ntuple<double,1> >;
template class python_deallocator<ntuple<long,1> >;
template class python_deallocator<ntuple<size_t,1> >;

template class python_deallocator<ntuple<std::complex<double>,2> >;
template class python_deallocator<ntuple<double,2> >;
template class python_deallocator<ntuple<long,2> >;
template class python_deallocator<ntuple<size_t,2> >;

template class python_deallocator<ntuple<std::complex<double>,3> >;
template class python_deallocator<ntuple<double,3> >;
template class python_deallocator<ntuple<long,3> >;
template class python_deallocator<ntuple<size_t,3> >;

// --- associated static variables (see comments at definition of "python_deallocator<T>"): ---

template< >
PyObject* python_deallocator<std::complex<double> >::class_type_object_ 
  = python_deallocator<std::complex<double> >::registered_type_object();
template< >
PyObject* python_deallocator<double>::class_type_object_ 
  = python_deallocator<double>::registered_type_object();
template< >
PyObject* python_deallocator<long>::class_type_object_ 
  = python_deallocator<long>::registered_type_object();

template< >
PyObject* python_deallocator<ntuple<std::complex<double>,1> >::class_type_object_ 
  = python_deallocator<ntuple<std::complex<double>,1> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<double,1> >::class_type_object_ 
  = python_deallocator<ntuple<double,1> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<long,1> >::class_type_object_ 
  = python_deallocator<ntuple<long,1> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<size_t,1> >::class_type_object_ 
  = python_deallocator<ntuple<size_t,1> >::registered_type_object();

template< >
PyObject* python_deallocator<ntuple<std::complex<double>,2> >::class_type_object_ 
  = python_deallocator<ntuple<std::complex<double>,2> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<double,2> >::class_type_object_ 
  = python_deallocator<ntuple<double,2> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<long,2> >::class_type_object_ 
  = python_deallocator<ntuple<long,2> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<size_t,2> >::class_type_object_ 
  = python_deallocator<ntuple<size_t,2> >::registered_type_object();

template< >
PyObject* python_deallocator<ntuple<std::complex<double>,3> >::class_type_object_ 
  = python_deallocator<ntuple<std::complex<double>,3> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<double,3> >::class_type_object_ 
  = python_deallocator<ntuple<double,3> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<long,3> >::class_type_object_ 
  = python_deallocator<ntuple<long,3> >::registered_type_object();
template< >
PyObject* python_deallocator<ntuple<size_t,3> >::class_type_object_ 
  = python_deallocator<ntuple<size_t,3> >::registered_type_object();

// ==================================================================================



// ----------------- general specializations: ---------------------------

template < >
bool type_check<simple_object_base*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  if (PyNumber_Check(src)
      || PyBool_Check(src)
      || PyString_Check(src)
      || PySequence_Check(src)
      || PyMapping_Check(src))
    return true;
    
  return false;
}

template < >
bool type_check<bool>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // no conversion to bool: require actual bool type:
  return PyBool_Check(src); 
}

template < >
bool type_check<std::string>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  return PyString_Check(src);
}


template < >
bool type_check<size_t>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // here we do not require in-memory compatibility (that check is at "type_check<size_t*>"):
  return PyLong_Check(src) || PyInt_Check(src);
}

template < >
int datatype_enum<std::complex<double> >(void)
{ return NPY_CDOUBLE; }

template < >
int datatype_enum<double>(void)
{ return NPY_DOUBLE; }

template < >
int datatype_enum<long>(void)
{ return NPY_LONG; }

// WARNING: size_t converts to NPY_LONG, and this assumes in-memory compatibility of these types:
template < >
int datatype_enum<size_t>(void)
{ 
#if !defined(__x86_64__)
  // note: on 32-bit platforms, NPY_LONG is _not_ equivalent to python long type (which apparently is 64-bit):
  return NPY_LONG; 
#else
  return NPY_LONG;
#endif
}

/**
 * @brief Convert from python object to simple_object_base*, using the default extractor for class simple_object_base: returns new pointer.
 *   See: "simple_object_base::set_default_extractor" method.
 *   Any RANK>0 number objects created in the simple-object that reference array data in the python object
 *     will have read-only access (enforced at non-const simple_object_base::ptr<T>(void) method).
 */
template < >
simple_object_base* extract<simple_object_base*>(const PyObject* src) throw(python_error)
{ return simple_object_base::extract(src); }

template < >
bool extract<bool>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  bool dest(false);
  
  if (PyBool_Check(src))
    dest = (src == Py_True? true: false);
  else
    throw python_error("extract<bool>: PyObject* is not a boolean");
  
  return dest;
}


template < >
std::string extract<std::string>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  std::string dest;
  
  if (PyString_Check(src))
    dest = PyString_AsString(const_cast<PyObject*>(src));  
  else
    throw python_error("extract<std::string>: PyObject* is not a string");
  
  return dest;
}


/**
 * @brief Convert from simple_object_base* to python object, using the default inserter for class simple_object_base: returns new reference.
 *    WARNING: read/write references are created to any RANK>0 number data in the source simple_object;
 *      data ownership is transferred to the python object, if possible.
 *    See: "simple_object_base::set_default_inserter" method.
 */
template < >
PyObject* insert<simple_object_base* const>(simple_object_base* const& src) throw(python_error)
{ return simple_object_base::insert(src); }

template < >
PyObject* insert<bool>(const bool& flag) throw(python_error)
{
  PyObject* dest(NULL);
  
  dest = (flag? Py_True: Py_False);
  Py_INCREF(dest);
  
  return dest;
}


template < >
PyObject* insert<std::string>(const std::string& s) throw(python_error)
{
  PyObject* dest(NULL);
  
  dest = PyString_FromString(s.c_str());
  if (NULL == dest)
    throw python_error("insert<std::string>: PyString_FromString error return");
    
  return dest;    
}

// --------------------------- double-precision specializations: ---------------------------------


template < >
bool type_check<std::complex<double> >(const PyObject* src_)
{ 
  PyObject *src = const_cast<PyObject*>(src_);
  return PyComplex_Check(src); 
}
 
template < >
bool type_check<std::complex<double>*>(const PyObject* src_)
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_CDOUBLE:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY(src) || PyArray_ISFARRAY(src))
      && (datatype_enum<std::complex<double> >() == PyArray_TYPE(src)))
    return true;
  
  return false;
} 
 
template < >
bool type_check<const std::complex<double>*>(const PyObject* src_)
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_CDOUBLE:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY_RO(src) || PyArray_ISFARRAY_RO(src))
      && (datatype_enum<std::complex<double> >() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<double>(const PyObject* src_)
{ 
  PyObject *src = const_cast<PyObject*>(src_);
  return PyFloat_Check(src); 
}
 
template < >
bool type_check<double*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_DOUBLE:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY(src) || PyArray_ISFARRAY(src))
      && (datatype_enum<double>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<const double*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_DOUBLE:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY_RO(src) || PyArray_ISFARRAY_RO(src))
      && (datatype_enum<double>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<long>(const PyObject* src_)
{ 
  PyObject *src = const_cast<PyObject*>(src_);
  return PyLong_Check(src) || PyInt_Check(src); 
}
 
template < >
bool type_check<long*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_LONG:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY(src) || PyArray_ISFARRAY(src))
      && (datatype_enum<long>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<const long*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_LONG:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY_RO(src) || PyArray_ISFARRAY_RO(src))
      && (datatype_enum<long>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<size_t*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_LONG:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY(src) || PyArray_ISFARRAY(src))
      && (datatype_enum<size_t>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}
 
template < >
bool type_check<const size_t*>(const PyObject* src_) 
{
  PyObject *src = const_cast<PyObject*>(src_);
  // supported conversion is from numpy array of type NPY_LONG:
  //   minimum requirement for meaningful conversion is that array is contiguous, aligned, in machine-byte order;
  //   it _also_ may or may-not be writeable (RO (read-only) is the lesser requirement).
  if (PyArray_Check(src)
      && (PyArray_ISCARRAY_RO(src) || PyArray_ISFARRAY_RO(src))
      && (datatype_enum<size_t>() == PyArray_TYPE(src)))
    return true;
  
  return false;
}

template < >
PyObject* insert<std::complex<double>, double, long>(const simple_object_base* src) throw(python_error)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  
  PyObject *dest(NULL);
  
  // Implementation note: "switch" statement by "src->type()" (or address thereof) does not work across shared library boundaries
  //   (without some really difficult work to ensure the type_info objects are not duplicated...).
  
  if (NULL == src){
    // allow "None" to correspond to NULL simple object: not sure if this is a good idea...
    dest = Py_None;
    Py_INCREF(dest);  
    std::cerr<<"WARNING: insert<C,R,Z>: converting NULL simple_object_base* to python object \"None\""<<std::endl;
  }
  else
  if (src->is<bool>()){
    assert(src->size() == 1); // RANK = 1 bool entities not supported...
    dest = insert(src->as<bool>());
  }  
  else
  if (src->is<C>()){
    if (src->size() == 1)    
      dest = insert(src->as<C>()); 
    else{
      // Implementation note: 
      //   RANK = 1 number entity: 
      //   Default practice is to transfer any data ownership to python object, 
      //     which can always be accomplished if the data is actually owned by the simple_object.
      //   Other forms of data transfer (i.e. which involve copy) can still be accomplished
      //     through the usage of the single-entity (RANK = 1) insert signature explicitly: 
      //       "PyObject *dest = insert<T>(T*, len, false);".
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<C>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }   
  else
  if (src->is<R>()){
    if (src->size() == 1)    
      dest = insert(src->as<R>()); 
    else{
      // RANK = 1 number entity:
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<R>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }    
  else
  if (src->is<Z>()){
    if (src->size() == 1)    
      dest = insert(src->as<Z>()); 
    else{
      // RANK = 1 number entity:
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<Z>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }     
  else
  if (src->is<std::string>()){
    assert(src->size() == 1); // RANK = 1 std::string entities not supported...
    dest = insert(src->as<std::string>());
  }  
  else
  if (src->is<LIST>()){
    // create corresponding python object as tuple:
    dest = PyTuple_New(src->as<LIST>().size());
    if (NULL == dest)
      throw python_error("insert<C,R,Z>: PyTuple_New error return");
    
    const LIST &l(src->as<LIST>());
    size_t n(0);
    for(LIST::const_iterator itl = l.begin(), itlEnd = l.end();
        itl != itlEnd;
        ++itl, ++n){
      PyObject *item_ = insert<C,R,Z>(*itl);

      // "PyTuple_SetItem" steals ref:
      if (PyTuple_SetItem(dest, n, item_)){
        Py_DECREF(dest);
        Py_DECREF(item_);
        throw python_error("insert<C,R,Z>: PyTuple_SetItem error return"); 
      }  
    }
  }
  else
  if (src->is<MAP>()){
    dest = PyDict_New();
    if (NULL == dest)
      throw python_error("insert<C,R,Z>: PyDict_New error return");
      
    const MAP &m(src->as<MAP>());
    for(MAP::const_iterator itm = m.begin(), itmEnd = m.end();
        itm != itmEnd;
        ++itm){
      const std::string &key_((*itm).first); 
      
      PyObject *val_ = insert<C,R,Z>((*itm).second);
       
      // PyDict_SetItemString borrows ref:
      if (PyDict_SetItemString(dest, key_.c_str(), val_)){
        Py_DECREF(dest);
        Py_DECREF(val_);
        throw python_error("insert<C,R,Z>: PyDict_SetItemStr error return"); 
      }
      Py_DECREF(val_);      
    }   
  }
  else  
    throw std::runtime_error("insert<simple_object_base*>: object type not in \n    { std::complex<double>, double, long, NULL, bool, std::string, object_list, object_map }");

  return dest;
}


template < >
PyObject* insert<std::complex<double> >(const std::complex<double>& c) throw(python_error)
{
  PyObject *dest(NULL);
  
  dest = PyComplex_FromDoubles(real(c), imag(c));
  if (NULL == dest)
    throw python_error("insert<std::complex<double> >: PyComplex_FromCComplex error return");
    
  return dest;   
}

template < >
PyObject* insert<double>(const double& r) throw(python_error)
{
  PyObject *dest(NULL);
  
  dest = PyFloat_FromDouble(r);
  if (NULL == dest)
    throw python_error("insert<double>: PyFloat_FromDouble error return");
    
  return dest;   
}


template < >
PyObject* insert<long>(const long& z) throw(python_error)
{
  PyObject *dest(NULL);
  
  dest = PyLong_FromLong(z);
  if (NULL == dest)
    throw python_error("insert<long>: PyLong_FromLong error return");
    
  return dest;   
}

//! Convert a python object to a simple_object (array references are read-only).
template < >
simple_object_base* extract<std::complex<double>, double, long>(const PyObject* src_) throw(python_error)
{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;

  #if 1
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  #endif
  
  simple_object_base *dest(NULL);

  if (NULL == src)
    throw std::runtime_error("extract<C,R,Z>: NULL PyObject*"); 

  if (Py_None == src){
    std::cerr<<"WARNING: extract<C,R,Z>: converting python object \"None\" to NULL"<<std::endl;
  }
  else
  if (PyBool_Check(src)){
    dest = new simple_object<bool>(extract<bool>(src));
  }  
  else
  if (PyNumber_Check(src)){
    #if !defined(EMBEDDED_PYTHON) // turned _off_ for the moment (problems with "multiarray.so" import):
    // array moved here (array supports number-protocol => PyNumber_Check(<array object>)):
    if (PyArray_Check(src)){
      if (type_check<const C*>(src)){
        const C *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<C>(ptr, len, false);
      }
      else
      if (type_check<const R*>(src)){
        const R *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<R>(ptr, len, false);
      }      
      else
      if (type_check<const Z*>(src)){
        const Z *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<Z>(ptr, len, false);
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
        throw std::runtime_error("extract<simple_object_base*>: either python array cannot be referenced as in-memory mapping,\n"
                                 "  or its data-type is not supported");  
      }                              
    }    
    else
    #endif
    if (PyComplex_Check(src))
      dest = new simple_object<C>(extract<C>(src));
    else
    if (PyFloat_Check(src))
      dest = new simple_object<R>(extract<R>(src));
    else
    if (PyLong_Check(src) || PyInt_Check(src))
      dest = new simple_object<Z>(extract<Z>(src));
  }
  else
  if (PySequence_Check(src)){  
    if (PyString_Check(src)){
      dest = new simple_object<std::string>(extract<std::string>(src));
    }  
    else{
      dest = new simple_object<LIST>();
      int size_(PySequence_Size(src));
      if (-1 == size_)
        throw python_error("extract<simple_object_base*>: PySequence_Size error return");

      for(size_t n = 0; n < static_cast<size_t>(size_); ++n){        
        // "PySequence_GetItem" returns new ref:
        PyObject *item_ = PySequence_GetItem(src, n);
        if (NULL == item_)
          throw python_error("extract<simple_object_base*>: PySequence_GetItem error return");
        dest->as<LIST>().push_back(extract<C,R,Z>(item_));
        Py_DECREF(item_);
      }  
    }  
  }
  else
  if (PyMapping_Check(src)){
    // extract object-map with _string_ key and arbitrary value:
    dest = new simple_object<MAP>();
    std::string key;
    simple_object_base *val;

    // "PyMapping_Items" returns new ref:
    PyObject *items_ = PyMapping_Items(src);
    if (NULL == items_)
      throw python_error("extract<simple_object_base*>: PyMapping_Items error return");    
    
    int size_ = PyMapping_Size(src);
    if (-1 == size_){
      Py_DECREF(items_);
      throw python_error("extract<simple_object_base*>: PyMapping_Size error return");
    }
          
    for(size_t n = 0; n < static_cast<size_t>(size_); ++n){
      
      // "PyList_GetItem" returns borrowed ref:
      PyObject *item_ = PyList_GetItem(items_, n);
      if (NULL == item_){
        Py_DECREF(items_);
        throw python_error("extract<simple_object_base*>: PyList_GetItem error return");
      }
      
      // "PyTuple_GetItem" returns borrowed ref:
      PyObject *key_ = PyTuple_GetItem(item_,0);
      PyObject *val_ = PyTuple_GetItem(item_,1);
      if (NULL == key_ || NULL == val_){
        Py_DECREF(items_);
        throw python_error("extract<simple_object_base*>: PyTuple_GetItem error return");
      }

      key = extract<std::string>(key_);
      val = extract<C,R,Z>(val_);
      dest->as<MAP>()[key] = val;
    }
    
    Py_DECREF(items_);    
  }
  else{
    // obtain a string representation of the bad-object, to add to exception string:
    
    PyObject *format_ = PyString_FromString("%s");
    if (NULL == format_)
      throw python_error("extract<simple_object_base*>: PyString_FromString error return");
    
    PyObject *arg_ = PyTuple_New(1);
    if (NULL == arg_){
      Py_DECREF(format_);
      throw python_error("extract<simple_object_base*>: PyTuple_New error return");
    }
    
    Py_INCREF(src);
    // "PyTuple_SetItem" steals ref:
    if (PyTuple_SetItem(arg_, 0, src)){
      Py_DECREF(format_);
      Py_DECREF(arg_);
      Py_DECREF(src); 
      throw python_error("extract<simple_object_base*>: PyTuple_SetItem error return");    
    }

    PyObject* str_ = PyString_Format(format_, arg_);
    if (NULL == str_){
      Py_DECREF(format_);
      Py_DECREF(arg_); // at this point: arg_ owns its ref to src
      throw python_error("extract<simple_object_base*>: PyString_Format error return");         
    }
    
    const char* obj_str_ = PyString_AsString(str_);
    if (NULL == obj_str_){
      Py_DECREF(format_);
      Py_DECREF(arg_);
      Py_DECREF(str_);
      throw python_error("extract<simple_object_base*>: PyString_AsString error return");               
    }
    
    // copy to the output std::string before releasing the ref to "str_"    
    std::string msg("extract<simple_object_base*>: object type not in { <number>, bool, <sequence>, <mapping> }: ");
    msg += obj_str_;
    
    Py_DECREF(format_);
    Py_DECREF(arg_);
    Py_DECREF(str_);    
    
    throw msg;        
  }

  if (NULL == dest)
    throw std::runtime_error("extract<simple_object_base*>: object extracted as NULL pointer");
    
  return dest;
} 


template < >
std::complex<double> extract<std::complex<double> >(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));

  Py_complex c_;

  if (PyComplex_Check(src))
    c_ = PyComplex_AsCComplex(src);
  else
    throw python_error("extract<std::complex<double> >: PyObject* is not a complex number");

  return std::complex<double>(c_.real, c_.imag);
}

template < >
double extract<double>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));

  double r;

  if (PyFloat_Check(src))
    r = PyFloat_AsDouble(src);
  else
    throw python_error("extract<double>: PyObject* is not a real number");
  
  return r;
}

template < >
long extract<long>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  long z;

  if (PyLong_Check(src))
    z = PyLong_AsLong(src);
  else
  if (PyInt_Check(src)){
    // also allow python integer type:
    PyObject *long_ = PyNumber_Long(src);
    if (NULL == long_)
      throw python_error("extract<std::string>: PyNumber_Long error return");
    z = PyLong_AsLong(long_);

    Py_DECREF(long_);
  }  
  else  
    throw python_error("extract<long>: PyObject* is not an integer");

  return z;
}

#if defined(__USE_MERE)
// ================ definitions of arbitrary-precision versions must co-exist with double-precision definitions ===========================

/*
 * notes: 
 *   - mere classes can co-exist with double, so these definitions must work independently;
 *   - important: at present "bignum" or its python equivalent is not used on the python side;
 *       these conversions go to and from python "complex", "float" and "long" types.
 *   .
 */


template < >
bool type_check<mere::C>(const PyObject* src)
{ return PyComplex_Check(src); } 

template < >
bool type_check<mere::C*>(const PyObject* src) 
{
  // Implementation note: 
  //   this could simply return false, but the exception provides more information;
  //   generally calling this method would indicate a usage-error of some kind.
  throw std::runtime_error("python_util::implementation_module::type_check: \n"
    "  no direct in-memory mapping is available between type mere::C and python object type.");
  /* throw unexpected */
}

template < >
bool type_check<mere::R>(const PyObject* src)
{ return PyFloat_Check(src); } 
 
template < >
bool type_check<mere::R*>(const PyObject* src) 
{
  // Implementation note: 
  //   this could simply return false, but the exception provides more information;
  //   generally calling this method would indicate a usage-error of some kind.
  throw std::runtime_error("python_util::implementation_module::type_check: \n"
    "  no direct in-memory mapping is available between type mere::R and python object type.");
  /* throw unexpected */
}

template < >
bool type_check<mere::Z>(const PyObject* src)
{ return PyLong_Check(src) || PyInt_Check(src); } 
 
template < >
bool type_check<mere::Z*>(const PyObject* src) 
{
  // Implementation note: 
  //   this could simply return false, but the exception provides more information;
  //   generally calling this method would indicate a usage-error of some kind.
  throw std::runtime_error("python_util::implementation_module::type_check: \n"
    "  no direct in-memory mapping is available between type mere::Z and python object type.");
  /* throw unexpected */
}


template < >
PyObject* insert<mere::C, mere::R, mere::Z>(const simple_object_base* src) throw(python_error)
{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
  
  PyObject *dest(NULL);

  if (NULL == src){
    // allow "None" to correspond to NULL simple object: not sure if this is a good idea...
    dest = Py_None;
    Py_INCREF(dest);  
    std::cerr<<"WARNING: insert<C,R,Z>: converting NULL simple_object_base* to python object \"None\""<<std::endl;
  }
  else
  if (src->is<bool>()){
    assert(src->size() == 1); // RANK = 1 bool entities not supported...
    dest = insert(src->as<bool>());
  }  
  else
  if (src->is<C>()){
    if (src->size() == 1)    
      dest = insert(src->as<C>()); 
    else{
      // Implementation note: 
      //   RANK = 1 number entity: 
      //   Default practice is to transfer any data ownership to python object, 
      //     which can always be accomplished if the data is actually owned by the simple_object.
      //   Other forms of data transfer (i.e. which involve copy) can still be accomplished
      //     throw the usage of the single-entity (RANK = 1) insert signature explicitly: 
      //       "PyObject *dest = insert<T>(T*, len, false);".
      
      // I'll leave these sections in, for the moment, 
      //   but not having an arbitrary-precision representation on the python-side,
      //   with direct in-memory mapping to C types, these insert methods will throw an exception
      //   (at either of "type_check" or "datatype_enum").
      
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<C>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }   
  else
  if (src->is<R>()){
    if (src->size() == 1)    
      dest = insert(src->as<R>()); 
    else{
      // RANK = 1 number entity:
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<R>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }    
  else
  if (src->is<Z>()){
    if (src->size() == 1)    
      dest = insert(src->as<Z>()); 
    else{
      // RANK = 1 number entity:
      bool transfer_data(src->own_data());
      dest = insert(src->ptr<Z>(), src->size(), transfer_data);
      if (transfer_data)
        src->release_data();      
    }
  }     
  else
  if (src->is<std::string>()){
    assert(src->size() == 1); // RANK = 1 std::string entities not supported...
    dest = insert(src->as<std::string>());
  }  
  else
  if (src->is<LIST>()){
    // create corresponding python object as tuple:
    dest = PyTuple_New(src->as<LIST>().size());
    if (NULL == dest)
      throw python_error("insert<C,R,Z>: PyTuple_New error return");
    
    const LIST &l(src->as<LIST>());
    size_t n(0);
    for(LIST::const_iterator itl = l.begin(), itlEnd = l.end();
        itl != itlEnd;
        ++itl, ++n){
      PyObject item_ = insert<C,R,Z>(*itl);

      // "PyTuple_SetItem" steals ref:
      if (PyTuple_SetItem(dest, n, item_)){
        Py_DECREF(dest);
        Py_DECREF(item_);
        throw python_error("insert<C,R,Z>: PyTuple_SetItem error return"); 
      }  
    }
  }
  else
  if (src->is<MAP>()){
    dest = PyDict_New();
    if (NULL == dest)
      throw python_error("insert<C,R,Z>: PyDict_New error return");
      
    const MAP &m(src->as<MAP>());
    for(MAP::const_iterator itm = m.begin(), itmEnd = m.end();
        itm != itmEnd;
        ++itm){
      const std::string &key_((*itm).first); 
      
      PyObject *val_ = insert<C,R,Z>((*itm).second);
       
      // PyDict_SetItemString borrows ref: 
      if (PyDict_SetItemString(dest, key_.c_str(), val_)){
        Py_DECREF(dest);
        Py_DECREF(val_);
        throw python_error("insert<C,R,Z>: PyDict_SetItemStr error return"); 
      }
      Py_DECREF(val_);      
    }   
  }
  else  
    throw python_error("insert<simple_object_base*>: object type not in \n    { std::complex<double>, double, long, NULL, bool, std::string, object_list, object_map }");

  return dest;
}

template < >
PyObject* insert<mere::C>(const mere::C& c) throw(python_error)
{
  PyObject *dest(NULL);
  
  // use identical forms for all C,R,Z:
  Py_complex c_;
  mere::conv(c_.real, real(c));
  mere::conv(c_.imag, imag(c));
  
  dest = PyComplex_FromCComplex(c_);
  if (NULL == dest)
    throw python_error("insert<std::complex<double> >: PyComplex_FromCComplex error return");
    
  return dest;   
}

template < >
PyObject* insert<mere::R>(const mere::R& r) throw(python_error)
{
  PyObject *dest(NULL);
  
  // use identical forms for all C,R,Z:
  double r_;
  mere::conv(r_, r);
  
  dest = PyFloat_FromDouble(r_);
  if (NULL == dest)
    throw python_error("insert<double>: PyFloat_FromDouble error return");
    
  return dest;   
}


template < >
PyObject* insert<mere::Z>(const mere::Z& z) throw(python_error)
{
  PyObject *dest(NULL);
  
  // use identical forms for all C,R,Z:
  long z_;
  mere::conv(z_, z);
  
  dest = PyLong_FromLong(z_);
  if (NULL == dest)
    throw python_error("insert<long>: PyLong_FromLong error return");
    
  return dest;   
}


//! Convert a python object to a simple_object (array references are read/write).
template < >
simple_object_base* extract<mere::C, mere::R, mere::Z>(const PyObject* src_) throw(python_error)
{
  using boost::python::extract;
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef simple_object_base::object_list LIST;
  typedef simple_object_base::object_map  MAP;
    
  #if 1  
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  #endif
  
  simple_object_base *dest(NULL);

  if (NULL == src)
    throw python_error("extract<C,R,Z>: NULL PyObject*"); 

  if (Py_None == src){
    std::cerr<<"WARNING: extract<C,R,Z>: converting python object \"None\" to NULL"<<std::endl;
  }
  else
  if (PyBool_Check(src)){
    dest = new simple_object<bool>(extract<bool>(src)); 
  } 
  else
  if (PyNumber_Check(src)){      
    #if !defined(EMBEDDED_PYTHON) // turned _off_ for the moment (problems with "multiarray.so" import):
    // array moved here (array supports number-protocol => PyNumber_Check(<array object>)):
    if (PyArray_Check(src)){
      // implementation note: at present, no direct in-memory mapping is available for mere types from python types;
      //   however, this section will be left in for informational purposes; 
      //   type_check will throw an exception if this conversion is ever requested.
      if (type_check<const C*>(src)){
        const C *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<C>(ptr, len, false);
      }
      else
      if (type_check<const R*>(src)){
        const R *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<R>(ptr, len, false);
      }      
      else
      if (type_check<const Z*>(src)){
        const Z *ptr(NULL);
        size_t len(0);
        // create simple_object as a reference to the python array data:
        extract(ptr, len, src, false);
        dest = new simple_object<Z>(ptr, len, false);
      }
      else
        throw std::runtime_error("extract<simple_object_base*>: either python array cannot be referenced as in-memory mapping,\n"
                                 "  or its data-type is not supported");     
    }
    else
    #endif
    if (PyComplex_Check(src))
      dest = new simple_object<C>(extract<C>(src));
    else
    if (PyFloat_Check(src))
      dest = new simple_object<R>(extract<R>(src));
    else
    if (PyLong_Check(src) || PyInt_Check(src))
      dest = new simple_object<Z>(extract<Z>(src));
  }
  else
  if (PySequence_Check(src)){  
    if (PyString_Check(src)){
      dest = new simple_object<std::string>(extract<std::string>(src));
    }
    else{
      dest = new simple_object<LIST>();
      int size_(PySequence_Size(src));
      if (-1 == size_)
        throw python_error("extract<simple_object_base*>: PySequence_Size error return");
      for(size_t n = 0; n < static_cast<size_t>(size_); ++n){
        
        // "PySequence_GetItem" returns new ref:
        PyObject *item_ = PySequence_GetItem(src, n);
        if (NULL == item_)
          throw python_error("extract<simple_object_base*>: PySequence_GetItem error return");
        dest->as<LIST>().push_back(extract<C,R,Z>(item_));
        Py_DECREF(item_);
      }  
    }  
  }
  else
  if (PyMapping_Check(src)){
    // extract object-map with _string_ key and arbitrary value:
    dest = new simple_object<MAP>();
    std::string key;
    simple_object_base *val;

    // "PyMapping_Items" returns new ref:
    PyObject *items_ = PyMapping_Items(src);
    if (NULL == items_)
      throw python_error("extract<simple_object_base*>: PyMapping_Items error return");    
    
    int size_ = PyMapping_Size(src);
    if (-1 == size_){
      Py_DECREF(items_);
      throw python_error("extract<simple_object_base*>: PyMapping_Size error return");
    }
          
    for(size_t n = 0; n < static_cast<size_t>(size_); ++n){

      // "PyList_GetItem" returns borrowed ref:
      PyObject *item_ = PyList_GetItem(items_, n);
      if (NULL == item_){
        Py_DECREF(items_);
        throw python_error("extract<simple_object_base*>: PyList_GetItem error return");
      }
      
      // "PyTuple_GetItem" returns borrowed ref:
      PyObject *key_ = PyTuple_GetItem(item_,0);
      PyObject *val_ = PyTuple_GetItem(item_,1);
      if (NULL == key_ || NULL == val_){
        Py_DECREF(items_);
        throw python_error("extract<simple_object_base*>: PyTuple_GetItem error return");
      }

      key = extract<std::string>(key_);
      val = extract<C,R,Z>(val_);
      dest->as<MAP>()[key] = val;
    }
    
    Py_DECREF(items_);    
  }
  else{
    // obtain a string representation of the bad-object, to add to exception string:
    
    PyObject *format_ = PyString_FromString("%s");
    if (NULL == format_)
      throw python_error("extract<simple_object_base*>: PyString_FromString error return");
    
    PyObject *arg_ = PyTuple_New(1);
    if (NULL == arg_){
      Py_DECREF(format_);
      throw python_error("extract<simple_object_base*>: PyTuple_New error return");
    }
    
    Py_INCREF(src);
    // "PyTuple_SetItem" steals ref:
    if (PyTuple_SetItem(arg_, 0, src)){
      Py_DECREF(format_);
      Py_DECREF(arg_);
      Py_DECREF(src_); 
      throw python_error("extract<simple_object_base*>: PyTuple_SetItem error return");    
    }

    PyObject* str_ = PyString_Format(format_, arg_);
    if (NULL == str_){
      Py_DECREF(format_);
      Py_DECREF(arg_); // at this point: arg_ owns its ref to src
      throw python_error("extract<simple_object_base*>: PyString_Format error return");         
    }
    
    const char* obj_str_ = PyString_AsString(str_);
    if (NULL == obj_str_){
      Py_DECREF(format_);
      Py_DECREF(arg_);
      Py_DECREF(str_);
      throw python_error("extract<simple_object_base*>: PyString_AsString error return");               
    }
    
    // copy to the output std::string before releasing the ref to "str_"    
    std::string msg("extract<simple_object_base*>: object type not in { <number>, bool, <sequence>, <mapping> }: ") + obj_str_;
    
    Py_DECREF(format_);
    Py_DECREF(arg_);
    Py_DECREF(str_);    
    
    throw msg;        
  }

  if (NULL == dest)
    throw std::runtime_error("extract<simple_object_base*>: object extracted as NULL pointer (=> any call to clone will _segfault_!)");

  return dest;
} 


template < >
mere::C extract<mere::C>(const PyObject* src_) throw(python_error)
{  
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  mere::C c;

  // use identical forms for all C,R,Z:  
  if (PyComplex_Check(src)){
    Py_complex c_ = PyComplex_AsCComplex(src);
    mere::conv(real(c), c_.real);
    mere::conv(imag(c), c_.imag);    
  }
  else
    throw python_error("extract<std::complex<mere::C> >: PyObject* is not a complex number");

  return c;
}

template < >
double extract<mere::R>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  mere::R r;

  // use identical forms for all C,R,Z:  
  if (PyFloat_Check(src)){
    double r_ = PyFloat_AsDouble(src);
    mere::conv(r, r_);
  }
  else
    throw python_error("extract<mere::R>: PyObject* is not a real number");

  return r;
}

template < >
long extract<mere::Z>(const PyObject* src_) throw(python_error)
{
  // python C/API does not use "const":
  PyObject *src(const_cast<PyObject*>(src_));
  
  mere::Z z;

  // use identical forms for all C,R,Z:  
  if (PyLong_Check(src)){
    long z_ = PyLong_AsLong(src);
    mere::conv(z, z_);
  }
  else
  if (PyInt_Check(src)){
    // also allow python integer type:
    PyObject *long_ = PyNumber_Long(src);
    if (NULL == long_)
      throw python_error("extract<std::string>: PyNumber_Long error return");
    long z_ = PyLong_AsLong(long_);
    conv(z, z_);
    Py_DECREF(long_);
  }  
  else  
    throw python_error("extract<mere::Z>: PyObject* is not an integer");

  return r;
}

#endif


} // namespace implementation_module
} // namespace python_util
// ============================= end: non-exclusive python C/API functions =============================================================================
