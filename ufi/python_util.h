#if !defined(__python_util__h)
#define __python_util__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/python_util.h,v $

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


//! Methods to parse a restricted python ``list'' syntax, for use in parameter and command-line options I/O.
/**
 * Classes and methods to allow selected interface to an embedded python interpreter.  The specific objective of these classes
 * is to allow ease of processing of complex input and output parameter sets.
 */


#if !defined(_TYPEINFO_H) && !defined(_TYPEINFO)
  #error #include <typeinfo> prior to this file!
#endif
 

#if defined(__USE_MERE)
  // warning: this may screw-up the general header-file sequence:
  #include <mere.h>
#endif 
 
#if !defined(__commUtil__h)
// commUtil.h:
namespace commUtil{
  class abstractCommHandle; // forward declaration
}

using commUtil::abstractCommHandle;  // forward declaration
#endif 

#if !defined(Py_PYTHON_H)
  // Python.h:
  struct PyObject; // forward decl
#endif 
  
namespace python_util{


/**
 * @brief Exception class for methods using python C/API (derived from std::exception).
 */
class python_error: public std::exception{
  private:
    std::string what_;
  public:
    virtual const char* what(void)const noexcept (true);
    
    python_error& operator=(const python_error& other);    
    
    virtual ~python_error(void);
    
    python_error(void);
    python_error(const python_error& other);
    python_error(const std::string& msg);
}; // class python_error 


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
PyObject*  wrapped_external_functor(PyObject *self, PyObject *args);
 

class simple_object_base{
  public:
    // note: changed from std::list to std::vector to allow the convenience of subscripting.
    //  (after initial implementation it turned out that "insert" was never necessary)
    typedef std::vector<simple_object_base*> object_list;
    
    typedef std::unordered_map<std::string, simple_object_base*> object_map;
    
    /**
     * @brief Extract function typedef: convert a python object to a simple_object_base*: returns new pointer.
     */
    typedef simple_object_base* (*extract_func)(const PyObject* src);

    /**
     * @brief Type of insert function: convert simple_object_base* to a python object: returns new reference.
     *   Warning: any RANK>0 number objects in simple_object will be referenced (in python) as read/write python arrays;
     *     data ownership will be transfered to the python object.
     */
    typedef PyObject* (*insert_func)(const simple_object_base* src);
 
  protected:
  
    // (Note: _implicitely_ this enforces that the _only_ possible 
    //    simple_object<T> are the enumerated types in "RECOGNIZED_TYPE";
    //    possibly a better design-scheme would be to do this _explicitely_?)
    // (Further note: this explicit design starts below at namespace scope with "simple_object_traits"...) 
    
    // enum of object-types for use by "writeBinaryVirtual" / "readBinaryVirtual" methods.
    // (note: fix enum values to assist in implementation of associated python methods.)
    enum RECOGNIZED_TYPE {NONE_TYPE=0, COMPLEX_TYPE=1, REAL_TYPE=2, INTEGER_TYPE=3, BOOL_TYPE=4, STRING_TYPE=5, OBJECT_LIST_TYPE=6, OBJECT_MAP_TYPE=7};
    
    static RECOGNIZED_TYPE enum_from_type_info(const std::type_info& type_);
    static const std::type_info&  type_info_from_enum(RECOGNIZED_TYPE etype_);

    // facilitate strict aliasing (technique in analogy to boost::any):
    class generic_ptr_base{
      public:
        
        virtual generic_ptr_base* clone(void)const;
        
        virtual ~generic_ptr_base(void);
    };
    
    /*
     * IMPLEMENTATION NOTE: apart from its generic nature, this class should behave 
     *   _exactly_ like any other pointer class; which is why "ptr" is public.  
     *   Therefore, the following functionality is avoided:
     *     - ownership tracking;
     *     - deep copy on assignment (or clone);
     *     - automatic deletion on destruction.
     *     .
     */
    template <class T>
    class generic_ptr: public generic_ptr_base{
      public:
        T *ptr;
        
        virtual generic_ptr_base* clone(void)const
        { return new generic_ptr(*this); }
        
        virtual ~generic_ptr(void)
        {
          #if 0
          if (NULL != ptr){
            delete ptr;
            ptr = NULL;
          }
          #else
          assert(NULL == ptr);
          #endif
        }
        
        generic_ptr(void)
          : ptr(NULL)
        { }
        
        generic_ptr(const generic_ptr& other)
          #if 0
          : ptr(other.ptr != NULL? new T(*other.ptr): NULL)
          #else
          : ptr(other.ptr)
          #endif
        { }  
    };
    
  private:
    #if 0
    void * pv_;
    #else
    // strict aliasing:
    generic_ptr_base *pg_;
    #endif

    static extract_func default_extract_func_;
    static insert_func default_insert_func_;
        
  protected:

    
    template <class T>
    inline const T& as_(void)const;
 
    template <class T>
    inline T& as_(void);

    template <class T>
    inline T*& ptr_(void);

    template <class T>
    inline const T* ptr_(void)const;
    
    // for use by simple_object<T>:

    static void free_pointers_(object_list& l);    
    static void free_pointers_(object_map& m);


    static simple_object_base* NULL_extract_func(const PyObject* src);
    
    static PyObject* NULL_insert_func(const simple_object_base* src);
    
    static bool writeBinary_(std::ostream &out, const object_list& l);      
    static bool readBinary_(std::istream &in, object_list& l);
    static size_t binarySize_(const object_list& l);

    static bool writeBinary_(std::ostream &out, const object_map& m);      
    static bool readBinary_(std::istream &in, object_map& m);
    static size_t binarySize_(const object_map& m);
              
  public:

    template <class T>
    inline bool is(void)const;

    /**
      * @brief test if simple_object is a RANK=0 object (though not necessarily a number object).
      */
    inline bool is_scalar(void)const;

    /**
      * @brief test if simple_object is one of: RANK>0 object (with N_components = 0), empty list, empty map.
      */
    inline bool is_empty(void)const;

    
    template <class T>
    inline const T& as(void)const;

    template <class T>
    inline T& as(void);

    /**
     * @brief Access object pointer as pointer to T.
     *   - throws an exception if object is not a simple_object<T>.
     *   .
     */
    template <class T>
    inline const T* ptr(void)const;

    /**
     * @brief Access object pointer as pointer to T.
     *   - throws an exception if object is not a simple_object<T>.
     *   - throws an exception if access is to RANK>0 number data
     *     _not_ owned by the object.
     *   .
     */
    template <class T>
    inline T* ptr(void);
    
    // ----------------------- object_map utility methods: ----------------------------------------

    // test if a parm exists in the map:
    static bool has_named_parm(const object_map& map, const std::string& key);

    static const std::type_info& named_parm_type(const object_map& map, const std::string& key); 

    // single_entity types U (see simple_object_traits<T>):
    // (note: use of this form with list_entity or object_entity types throws exception)
    template <class U>
    static const U& get_named_parm(const object_map& map, const std::string& key);

    template <class U>
    static void set_named_parm(object_map& map, const std::string& key, const U& val);

    // container types U (see simple_object_traits<T>):
    template <class U>
    static void get_named_parm(const object_map& map, const std::string& key, U& val);

    // object pointers themselves:
    static const simple_object_base* get_named_object(const object_map& map, const std::string& key);

    static void set_named_object(object_map& map, const std::string& key, 
                                 const simple_object_base* pobject, bool transfer_ownership=false);

    // ----------------------- conversion to and from python C/API PyObject*: ---------------------    
    
    /*
     * Implementation note:
     *   Extract and insert methods require information about the number types.
     *   As I do not want to make this base class a template class, I will use the following
     *     function-pointer based technique instead.
     *   Note that the default (i.e. not explicitely initialized) extractor and inserter throw exceptions.
     */
    
    /**
     * @brief Set default method for extraction from python object.
     */
    inline static void set_default_extractor(extract_func func);
    
    /**
     * @brief Convert the python object to a simple_object_base* using the default extractor: returns a new pointer.
     *    The "PyObject* style" extract signature (and meaning) is used here.
     */
    inline static simple_object_base* extract(const PyObject* src);
    
    /**
     * @brief Set default method for insertion to python object.
     *    Warning: any references to RANK>0 number data in the created python object will have read/write access:
     *    data ownership will be transferred to the python object.
     */
    inline static void set_default_inserter(insert_func func);
         
    /**
     * @brief Convert the simple_object_base* to a python object using the default inserter: returns a new reference.
     *    The "PyObject* style" insert signature (and meaning) is used here.
     */
    inline static PyObject* insert(const simple_object_base* src);
     
    
    // ----------------------- I/O methods (primary justification for this class) -----------------
        
    // initialize a simple_object_base pointer from a stream representation as python syntax 
    // (Allows arbitrary python code, including preamble and postscript. Note: name of python object is "object".):
    template <class C, class R, class Z>
    static simple_object_base* read_repn_virtual(const std::string& src, const std::string* ppreamble = NULL, const std::string* ppostscript = NULL);

    template <class C, class R, class Z>
    static simple_object_base* read_repn_virtual(std::istream& src, 
                                                 std::istream* ppreamble=NULL, 
                                                 std::istream* ppostscript=NULL);
    
    // output the object to a stream representation (as python syntax):
    template <class C, class R, class Z>
    void write_repn(std::ostream &dest)const;    
    
    // output the object to "pickle" binary representation:
    template <class C, class R, class Z>
    void write_pickle(const std::string& filename)const;    
 
    // construct object from "pickle" binary representation:
    template <class C, class R, class Z>
    static simple_object_base* read_pickle_virtual(const std::string& filename);    
       

    // output the object to a stream representation (in syntax appropriate for display to end-user):
    template <class C, class R, class Z>
    void write(std::ostream &dest)const;        

    // methods to allow binary read and write from pointer to base-class:
    static bool writeBinaryVirtual(std::ostream &out, const simple_object_base* pobject);
    static bool readBinaryVirtual(std::istream &in, simple_object_base*& pobject);
    static size_t binarySizeVirtual(const simple_object_base* pobject);

    virtual bool writeBinary(std::ostream &out) const;      
    virtual bool readBinary(std::istream &in);
    virtual size_t binarySize(void)const;

    // --------------------------------------------------------------------------------------------
    
       
    virtual const std::type_info& type(void)const;
    
    /**
     * @brief Allow contiguous storage arrays (type method will return type of value-type).
     *   IMPLEMENTATION NOTES:
     *     - allowing size > 1 will complicate extraction for target list types, but these _target_ types in general
     *         are single value-type container types, so it shouldn't be too difficult to implement;
     *     - this method shall return 1 for simple_object<object_map>, simple_object<object_list> specializations,
     *       which will also _hide_ the associated (*T, len, own_data=true) constructor signature;
     *     - all allocation management is at simple_object.
     */
    virtual size_t size(void)const;
    
    
    /**
     * @brief Data ownership state for the simple_object.
     */
    virtual bool own_data(void)const;

    /**
     * @brief Reset data-ownership flag associated with the simple_object.
     *   - At present, this is allowed for all single-entity simple_object 
     *    (i.e. not simple_object<object_list> or simple_object<object_map>).
     *   - Calling release_data for non single-entity objects, or when the 
     *     ownership flag is not set is an error.
     *   .
     */
    virtual void release_data(void)const;
   
    virtual simple_object_base* clone(void)const;

    virtual ~simple_object_base(void);

  protected:
    
    // protected constructor: transfer ownership of the generic_ptr<T>:    
    simple_object_base(generic_ptr_base *pg); 
};

/**
 * @brief A wrapper-class for virtual-base-class pointers so that read/writeBinary semantics can be used directly.
 *
 *   For the moment, this class holds its number system information in its template parameters; probably there is a better way to do this,
 *     but for "extraction" operations from string and PyObject* types, this information is required.
 *
 *  Usage: generic_object owns its underlying simple_object_base* => "clone()" occurs at assignment and deletion occurs at destruction. 
 *
 */
template <class C, class R, class Z>
class generic_object{
  private:
  
  simple_object_base *ptr_;
 
  public:

#if !defined(__INTEL_COMPILER) && !defined(__PGI)
  template <class T>
  inline bool is(void)const;

  template <class T>
  inline const T& as(void)const;

  template <class T>
  inline T& as(void);
#else
  // INTEL C++ compiler and PGI C++ compiler
  // template functions of template classes must have definitions at point of primary declaration:
  template <class T>
  inline bool is(void)const
  { return ptr_->is<T>(); }

  template <class T>
  inline const T& as(void)const
  { return ptr_->as<T>(); }

  template <class T>
  inline T& as(void)
  { return ptr_->as<T>(); }
#endif
 
  inline operator simple_object_base*(void) const;

  inline operator simple_object_base*(void);
  
  inline simple_object_base* ptr(void)const;
  
  inline simple_object_base* ptr(void);
  
  
  inline bool writeBinary(std::ostream &out) const;      

  inline bool readBinary(std::istream &in);
    
  inline size_t binarySize(void) const;


  // ------------------- conversion to and from PyObject*: ----------------------------
  
  /**
   * @brief Convert a generic_object to a python object: returns a new reference.
   */
  inline PyObject* insert(void)const;  
  
  /**
   * @brief Convert a python_object to a generic_object.
   */
  inline void extract(const PyObject* src);
  
  // ----------------------------------------------------------------------------------
  
  inline generic_object& operator=(const simple_object_base* other);
  
  inline void clone_from(const simple_object_base* other, bool transfer_ownership=false);
  
  inline generic_object& operator=(const generic_object& other);
  
  inline ~generic_object(void);
  
  inline generic_object(const generic_object& other);
 
  inline generic_object(const simple_object_base *p, bool transfer_ownership=false);
  
  inline generic_object(void);  
};
// ---------------------------------------------------------------------------------------------------------------------------------------------



template <class T>
class simple_object: public simple_object_base{
  private:
    
    // "size_" and "own_data_" attributes should not exist for "object_list" and "object_map" specializations:
    size_t size_;
    
    mutable bool own_data_;
    
  protected:
  
    // allocation control methods:
    void free_pointers_(void);
    void resize_(size_t new_size, bool new_own_data);
    
  public:
  
    typedef simple_object_base base_class;
    typedef base_class::object_list object_list;
    typedef base_class::object_map object_map;

    virtual bool writeBinary(std::ostream &out) const;      
    virtual bool readBinary(std::istream &in);
    virtual size_t binarySize(void)const;

    virtual const std::type_info& type(void)const;

    /**
     * @brief Allow contiguous storage arrays (type method will return type of value-type).
     *   IMPLEMENTATION NOTES:
     *     - allowing size > 1 will complicate extraction for target list types, but these _target_ types in general
     *         are single value-type container types, so it shouldn't be too difficult to implement;
     *     - this method shall return 1 for simple_object<object_map>, simple_object<object_list> specializations,
     *       which will also _hide_ the associated (*T, len, own_data=true) constructor signature;
     *     - all allocation management is at simple_object.
     */
     virtual size_t size(void)const;

    
    /**
     * @brief Data ownership state for the simple_object.
     */
    virtual bool own_data(void)const;

    /**
     * @brief Reset data-ownership flag associated with the simple_object.
     *   - At present, this is allowed for all single-entity simple_object 
     *    (i.e. not simple_object<object_list> or simple_object<object_map>).
     *   - Calling release_data for non single-entity objects, or when the 
     *     ownership flag is not set is an error.
     *   .
     */
    virtual void release_data(void)const;
 
    #if 0
    /**
     * @brief Work-around for the GNU C++ RTTI implementation problem that typeid(T) in different modules are not comparable.
     *    (i.e. test simple_object<T>::static_type == type() rather than typeid(T) == type() for cross-module comparison).
     *    ------------ THIS DOESN'T SOLVE THE PROBLEM! --------------
     */
    static const std::type_info& static_type(void); 
    #endif
    
    virtual simple_object_base* clone(void)const;
    
    simple_object& operator=(const simple_object& other);
    
    simple_object& operator=(const T& t);
    
    virtual ~simple_object(void);


    simple_object(void);
    
    simple_object(const simple_object<T>& other);
    
    simple_object(const T& t);
    
    /**
     * @brief Construct simple object from contiguous in-memory data.
     *   @param[in] ptr  data pointer
     *   @param[in] len  number of data elements
     *   @param[in] own_data  true => allocate, otherwise reference the data
     *   @param[in] copy_data  copy from the in-memory data
     */
    simple_object(const T* ptr, size_t len, bool own_data=true, bool copy_data=true);    
};

// ---------------------------------- explicit specialization for object_list: -----------------------------
template <>
class simple_object<simple_object_base::object_list>: public simple_object_base{
  private:

    void free_pointers_(void);

    void assign_with_clone_(const object_list& other);
  
  protected:
    
  public:
  
    typedef simple_object_base base_class;
    typedef base_class::object_list object_list;
    typedef base_class::object_map object_map;

    virtual bool writeBinary(std::ostream &out) const;      
    virtual bool readBinary(std::istream &in);
    virtual size_t binarySize(void) const;

    virtual const std::type_info& type(void) const;

    /**
     * @brief Allow contiguous storage arrays (type method will return type of value-type).
     *   IMPLEMENTATION NOTES:
     *     - allowing size > 1 will complicate extraction for target list types, but these _target_ types in general
     *         are single value-type container types, so it shouldn't be too difficult to implement;
     *     - this method shall return 1 for simple_object<object_map>, simple_object<object_list> specializations,
     *       which will also _hide_ the associated (*T, len, own_data=true) constructor signature;
     *     - all allocation management is at simple_object.
     */
    virtual size_t size(void)const;

    
    /**
     * @brief Data ownership state for the simple_object.
     */
    virtual bool own_data(void)const;
    
    /**
     * @brief Reset data-ownership flag associated with the simple_object.
     *   - At present, this is allowed for all single-entity simple_object 
     *    (i.e. not simple_object<object_list> or simple_object<object_map>).
     *   - Calling release_data for non single-entity objects, or when the 
     *     ownership flag is not set is an error.
     *   .
     */
    virtual void release_data(void)const;
 
    #if 0
    /**
     * @brief Work-around for the GNU C++ RTTI implementation problem that typeid(T) in different modules are not comparable.
     *    (i.e. test simple_object<T>::static_type == type() rather than typeid(T) == type() for cross-module comparison).
     */
    static const std::type_info& static_type(void); 
    #endif
    
    virtual simple_object_base* clone(void)const;
    
    simple_object& operator=(const simple_object& other);
    
    simple_object& operator=(const object_list& l);

    inline void clear(void);
    
    virtual ~simple_object(void);

    simple_object(void);
    
    simple_object(const simple_object& other);
    
    simple_object(const object_list& l);
    
    // U is any STL-style container class:
    template <class U>
    simple_object(const U& u); 
};

// ---------------------------------- explicit specialization for object_map: -----------------------------
template <>
class simple_object<simple_object_base::object_map>: public simple_object_base{
  private:

    void free_pointers_(void);

    void assign_with_clone_(const object_map& other);
    
  protected:
    
  public:
  
    typedef simple_object_base base_class;
    typedef base_class::object_list object_list;
    typedef base_class::object_map object_map;

    virtual bool writeBinary(std::ostream &out) const;      
    virtual bool readBinary(std::istream &in);
    virtual size_t binarySize(void) const;

    virtual const std::type_info& type(void) const;

    /**
     * @brief Allow contiguous storage arrays (type method will return type of value-type).
     *   IMPLEMENTATION NOTES:
     *     - allowing size > 1 will complicate extraction for target list types, but these _target_ types in general
     *         are single value-type container types, so it shouldn't be too difficult to implement;
     *     - this method shall return 1 for simple_object<object_map>, simple_object<object_list> specializations,
     *       which will also _hide_ the associated (*T, len, own_data=true) constructor signature;
     *     - all allocation management is at simple_object.
     */
     virtual size_t size(void)const;

    
    /**
     * @brief Data ownership state for the simple_object.
     */
    virtual bool own_data(void)const;

    /**
     * @brief Reset data-ownership flag associated with the simple_object.
     *   - At present, this is allowed for all single-entity simple_object 
     *    (i.e. not simple_object<object_list> or simple_object<object_map>).
     *   - Calling release_data for non single-entity objects, or when the 
     *     ownership flag is not set is an error.
     *   .
     */
    virtual void release_data(void)const;
    
    #if 0
    /**
     * @brief Work-around for the GNU C++ RTTI implementation problem that typeid(T) in different modules are not comparable.
     *    (i.e. test simple_object<T>::static_type == type() rather than typeid(T) == type() for cross-module comparison).
     */
    static const std::type_info& static_type(void); 
    #endif
    
    virtual simple_object_base* clone(void)const;
    
    simple_object& operator=(const simple_object& other);
    
    simple_object& operator=(const object_map& m);

    inline void clear(void);
    
    virtual ~simple_object(void);

    simple_object(void);
    
    simple_object(const simple_object& other);
    
    simple_object(const object_map& m);
    
    // M is any STL-style map class:
    template <class M>
    simple_object(const M& m); 
};

// ------------------------------------------ extensible_parameters_base: -------------------------------------------

/*
 * a base-class for parameters structures that require arbitrary extensibility:
 *   examples: 
 *     material_info (potential extensions for use with mixture formula permittivities)
 *     particle_info (potential extensions for use with arbitrary regions)
 *
 *   usage:
 *    associated simple_object_base* always is an object_list* at outer level (see "init" and "enforce_usage" methods)
 *    std::string <name>: allows determination of what is in LIST[2:]
 *    std::vector<T> <coefficients>: always present as LIST[1]
 *      (the assumption here is that at least a parameter list of simple number will always be required, 
 *         and further, this allows backwards compatibility with the existing usages).
 *    usage validation failures throw exceptions.
 *   warning: there are unexpected requirements for complex-number usage.
 *     Unless it is somehow _known_ that every eventual usage of the sub-classes will be
 *       with only real numbers, this class should be instantiated with a complex number-type. 
 *     (for example: material_info use with permittivityFunctor, 
 *            in the general case, requires complex parameters.
 *         particle_info use with regionSelector<T>, where regionSelector<T>::scale is complex (normal case for metals)
 *           and regionSelector<T> methods dealing with helical basis are complex).
 *     This does not mean that sub-classes shouldn't where appropriate _cast_ the coefficient parameters
 *       to real numbers (or ensure that they actually are real numbers), just
 *       that limiting the sub-class to real-number type may cause problems with
 *         eventual extensions at some later time.
 */
 
template <class T>
class extensible_parameters_base{
  public:
  
    typedef simple_object_base::object_list object_list;
    typedef simple_object_base::object_map object_map;
    
  private:
    simple_object_base *pobject_;
    
    // set default configuration:
    void init(const std::string* pname = NULL, const std::vector<T>* pcoeff = NULL);
    
    bool enforce_usage(void)const;
  
  protected:
  
    // named-parm map:
    object_map& parm_map(void); 
    const object_map& parm_map(void)const; 
      
    void dump_parm_map(void)const;
      
  public:
        
    const simple_object_base* object_pointer(void)const;
    
    // "clone_from" method clones object pointer from simple_object_base*, using its "clone" method,
    //    while enforcing structure of this base class (see usage comment above):
    // (transfer_ownership => assign object pointer: use with care)
    virtual void clone_from(const simple_object_base* pobject, bool transfer_ownership=false);
    
    virtual extensible_parameters_base* clone(void)const;
     
     
    std::string& name(void);
    const std::string& name(void)const;
    
    // U is any iterable container type of value_type convertable to T via direct assignment T = U::value_type:
    template <class U>
    void set_coeff(const U& src);

    // U is any iterable container type of value_type convertable to T via direct assignment U::value_type = T:
    template <class U>
    void get_coeff(U& dest)const;

    
    // arbitrary named parameters:
    
    // test if a parm exists:
    bool has_named_parm(const std::string& key)const;
 
    const std::type_info& named_parm_type(const std::string& key)const; 
       
    // single_entity types U (see simple_object_traits<T>):
    // (note: use of this form with list_entity or object_entity types throws exception)
    template <class U>
    const U& get_named_parm(const std::string& key)const;
   
    // container types U (see simple_object_traits<T>):
    template <class U>
    void get_named_parm(const std::string& key, U& val)const;

    template <class U>
    void set_named_parm(const std::string& key, const U& val);

    // object pointers themselves:
    const simple_object_base* get_named_object(const std::string& key)const;
    
    void set_named_object(const std::string& key, const simple_object_base* pobject, bool transfer_ownership=false);
    

    void swap(extensible_parameters_base& other);

    void copy(const extensible_parameters_base& other);
    extensible_parameters_base& operator=(const extensible_parameters_base& other);


    bool writeBinary(std::ostream &out) const;      
    bool readBinary(std::istream &in);
    size_t binarySize(void) const;

    // for the following methods, break-out number types as explicit parameters
    // (this keeps module independent from namespace TMatrix, otherwise, I need
    //    TMatrix::numberTraits<T> to obtain the types)
    template <class C, class R, class Z>
    void read(std::istream& is);
    template <class C, class R, class Z>   
    void write(std::ostream& os)const;

    virtual ~extensible_parameters_base(void);
    
    extensible_parameters_base(bool initialize=false);
    extensible_parameters_base(const extensible_parameters_base& other);
    extensible_parameters_base(const std::string& name, const std::vector<T>* pcoeff = NULL);
    extensible_parameters_base(const simple_object_base* object, bool transfer_ownership=false);
}; 


// ---------------------------------- command-line options implementation: ---------------------------------- 
// (note: with the addition of the parameter caching methods, this class will also serve as the base class for
//    a functor parameters class (e.g. where the given functor is the primary object of a test shell))
template <class T>
class options_map{
  public:
  
    typedef simple_object_base::object_list object_list;
    typedef simple_object_base::object_map object_map;
    
  private:
    simple_object_base *pobject_;

    bool cache_current_;
    
    
    void init(bool initialize=false);
    
    template <class C, class R, class Z>
    void init(const std::string& src, const std::string* ppreamble = NULL, const std::string* ppostscript = NULL);
    
    template <class C, class R, class Z>
    void init(std::istream& src, 
              std::istream* ppreamble=NULL, 
              std::istream* ppostscript=NULL);

    bool enforce_usage(void)const;
  
  protected:
  
    // named-parm map:
    object_map& parm_map(void); 
    const object_map& parm_map(void)const; 
      
    void dump_parm_map(void)const;
    
    
    void set_cache_current(bool flag=true);
          
  public:

    bool cache_current(void)const;

    virtual void update_cache(bool derived=false, bool to_cache=true);
        
    const simple_object_base* object_pointer(void)const;
    
    // "clone_from" method clones object pointer from simple_object_base*, using its "clone" method,
    //    while enforcing structure of this base class (see usage comment above):
    // (transfer_ownership => assign object pointer: use with care)
    virtual void clone_from(const simple_object_base* pobject, bool transfer_ownership=false);
    
    virtual options_map* clone(void)const;
    
    // test if a parm exists:
    bool has_named_parm(const std::string& key)const;
 
    const std::type_info& named_parm_type(const std::string& key)const; 
       
    // single_entity types U (see simple_object_traits<T>):
    // (note: use of this form with list_entity or object_entity types throws exception)
    template <class U>
    const U& get_named_parm(const std::string& key)const;
   
    // container types U (see simple_object_traits<T>):
    template <class U>
    void get_named_parm(const std::string& key, U& val)const;

    template <class U>
    void set_named_parm(const std::string& key, const U& val, bool from_cache=false);

    // object pointers themselves:
    const simple_object_base* get_named_object(const std::string& key)const;
    
    // note: assumes that object attributes are _not_ cached => no "from_cache" flag is provided for this method:
    void set_named_object(const std::string& key, const simple_object_base* pobject, bool transfer_ownership=false);
    

    void swap(options_map& other);

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
    void copy(const options_map& other, bool derived=false);
    options_map& operator=(const options_map& other);


    bool writeBinary(std::ostream &out) const;      
    bool readBinary(std::istream &in);
    size_t binarySize(void) const;

    // for the following methods, break-out number types as explicit parameters
    // (this keeps module independent from namespace TMatrix, otherwise, I need
    //    TMatrix::numberTraits<T> to obtain the types)
    template <class C, class R, class Z>
    void read(std::istream& is);
    
    template <class C, class R, class Z>   
    void write(std::ostream& os)const;


    virtual ~options_map(void);

    
    options_map(bool initialize=false);
    
    options_map(const options_map& other);

    #if 1
    // ---------- factory methods: -----------------
    template <class C, class R, class Z>
    static options_map parse_options(const std::string& src, const std::string* ppreamble = NULL, const std::string* ppostscript = NULL);
    
    template <class C, class R, class Z>
    static options_map parse_options(std::istream& src, 
                std::istream* ppreamble=NULL, 
                std::istream* ppostscript=NULL);
                
    template <class C, class R, class Z>
    static options_map parse_options(int argc, const char *argv[], 
                std::istream* ppreamble=NULL, 
                std::istream* ppostscript=NULL);
    // ---------------------------------------------
    #else
        
    template <class C, class R, class Z>
    options_map(const std::string& src, const std::string* ppreamble = NULL, const std::string* ppostscript = NULL);
    
    template <class C, class R, class Z>
    options_map(std::istream& src, 
                std::istream* ppreamble=NULL, 
                std::istream* ppostscript=NULL);
                
    template <class C, class R, class Z>
    options_map(int argc, const char *argv[], 
                std::istream* ppreamble=NULL, 
                std::istream* ppostscript=NULL);
    #endif
                
    options_map(const simple_object_base* object, bool transfer_ownership=false);
}; // class options_map


// ---------------------------------------- dispatcher classes: ----------------------------------
// traits class and support-types for expression templates to allow common types to be converted to and from simple_object:
// (note, this also explicitely defines recognized types)

struct abstract_true{
  abstract_true(void){}
};
struct abstract_false{
  abstract_false(void){}
};
struct abstract_null_entity{
  abstract_null_entity(void){}
};
struct single_entity{
  single_entity(void){}
};
struct list_entity{
  list_entity(void){}
};
struct map_entity{
  map_entity(void){}
};
struct complex_number_entity{
  complex_number_entity(void){}
};
struct real_number_entity{
  real_number_entity(void){}
};
struct integral_number_entity{
  integral_number_entity(void){}
};

// default:
template <class T>
struct simple_object_traits{
  typedef abstract_null_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

// explicit instantiations:

#if defined(__ntuple__h)
template <class T, size_t DIM>
struct simple_object_traits<linalg::ntuple<T,DIM> >{
  typedef list_entity entity_type;
  static inline bool is_single(void) { return false; }
};
#endif

#if defined(__USE_MERE)
template <>
struct simple_object_traits<mere::C>{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef single_entity entity_type;
  typedef complex_number_entity number_type;
  static inline bool is_single(void) { return true; }
};
template <>
struct simple_object_traits<mere::R>{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef single_entity entity_type;
  typedef real_number_entity number_type;
  static inline bool is_single(void) { return true; }
};
template <>
struct simple_object_traits<mere::Z>{
  typedef mere::C C;
  typedef mere::R R;
  typedef mere::Z Z;
  typedef single_entity entity_type;
  typedef integral_number_entity number_type;
  static inline bool is_single(void) { return true; }
};
#endif

// -------------- always define for double-precision types: ---------------------
template <>
struct simple_object_traits<std::complex<double> >{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef single_entity entity_type;
  typedef complex_number_entity number_type;
  static inline bool is_single(void) { return true; }
};
template <>
struct simple_object_traits<double>{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef single_entity entity_type;
  typedef real_number_entity number_type;
  static inline bool is_single(void) { return true; }
};
template <>
struct simple_object_traits<long>{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef single_entity entity_type;
  typedef integral_number_entity number_type;
  static inline bool is_single(void) { return true; }
};

// define size_t as an integral_number_entity, but only in a limited manner:
// (I'm not sure if this is a good idea...)
template <>
struct simple_object_traits<size_t>{
  typedef std::complex<double> C;
  typedef double R;
  typedef long Z;
  typedef single_entity entity_type;
  typedef integral_number_entity number_type;
  static inline bool is_single(void) { return true; }
};

#if 0 // -------------- omit: directly specialize the "new_simple_object" and "extract" for these types: ---------------------

// allow extraction of simple_object_base* itself:
template <>
struct simple_object_traits<simple_object_base*>{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

// ... and equivalently, of generic_object:
template <>
template <class C, class R, class Z>
struct simple_object_traits<generic_object<C,R,Z> >{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

// allow extraction and insertion of options_map<T> to and from simple_object_base*:
template <class T>
struct simple_object_traits<options_map<T> >{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

// allow extraction and insertion of extensible_parameters_base from and to simple_object_base*:
template <>
struct simple_object_traits<extensible_parameters_base>{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

#endif // -----------------end: omit ------------------------------------------------------------------------------------------

template <>
struct simple_object_traits<bool>{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

template <>
struct simple_object_traits<std::string>{
  typedef single_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return true; }
};

// std::vector, std::list, ntuple (traits defined in its ntuple.h, but in this namespace)
template <class T>
struct simple_object_traits<std::vector<T> >{
  typedef list_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return false; }
};

template <class T>
struct simple_object_traits<std::list<T> >{
  typedef list_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return false; }
};

// std::unordered_map with std::string as key (should this be generalized to allow other types as keys?):
template <class T>
struct simple_object_traits<std::unordered_map<std::string, T> >{
  typedef map_entity entity_type;
  typedef abstract_null_entity number_type;
  static inline bool is_single(void) { return false; }
};




// type-specific extraction from simple_object_base* :
template <class T>
void extract(T& t, const simple_object_base* pobject);

// -------- specializations of extract: -----------------------------------------------------------------------
template < >
inline void extract<PyObject*>(PyObject*& dest, const simple_object_base* pobject);

template < >
inline void extract<simple_object_base*>(simple_object_base*& dest, const simple_object_base* pobject);

template <class C, class R, class Z>
inline void extract(generic_object<C,R,Z>& dest, const simple_object_base* pobject);

template <class T>
inline void extract(extensible_parameters_base<T>& dest, const simple_object_base* pobject);

template <class T>
inline void extract(options_map<T>& dest, const simple_object_base* pobject);

// ------------------------------------------------------------------------------------------------------------

template <class T, class ET>
void dispatch_extract(T& t, const simple_object_base *pobject, ET);

template <class T>
void dispatch_extract(T& t, const simple_object_base *pobject, single_entity);

template <class U>
void dispatch_extract(U& u, const simple_object_base *pobject, list_entity);

template <class M>
void dispatch_extract(M& m, const simple_object_base *pobject, map_entity);

template <class T>
void extract_single(T& t, const simple_object_base *pobject);

template <class T, class NT>
inline void dispatch_extract_single(T& dest, const simple_object_base *pobject, NT);

template <class T>
inline void dispatch_extract_single(T& dest, const simple_object_base *pobject, abstract_null_entity);

// default: extract container from RANK=1 number object
//   or alternatively from object-list format: 
//   NT is simple_object_traits<typename U::value_type>::number_type
template <class U, class NT>
void extract_list(U& u, const simple_object_base *pobject, NT);

// special case: extract non-number container from object-list (or alternative format for number container):
template <class U>
void extract_list(U& u, const simple_object_base *pobject, abstract_null_entity);

template <class M>
void extract_map(M& m, const simple_object_base *pobject);


// generic allocation of simple_object_base*
template <class T>
simple_object_base* new_simple_object(const T& t);

// ------------- specializations of "new_simple_object" _bypassing_ dispatcher : --------------------------------

// simple_object_base* itself (using clone()):
inline simple_object_base* new_simple_object(const simple_object_base* t);

// generic_object (using clone()):
template <class C, class R, class Z>
inline simple_object_base* new_simple_object(const generic_object<C,R,Z>& t);

// extensible_parameters_base (using clone()):
template <class T>
inline simple_object_base* new_simple_object(const extensible_parameters_base<T>& t);

// options_map<T> (using clone()):
template <class T>
inline simple_object_base* new_simple_object(const options_map<T>& t);

// --------------------------------------------------------------------------------------------------------------

template <class T, class ET>
simple_object_base* new_object_dispatch(const T& t, ET);

template <class T>
simple_object_base* new_object_dispatch(const T& t, single_entity);

template <class U>
simple_object_base* new_object_dispatch(const U& u, list_entity);

template <class M>
simple_object_base* new_object_dispatch(const M& m, map_entity);

template <class T>
simple_object_base* new_single(const T& t);

// default: create a RANK=1 number object: 
//  "NT" is 
//     simple_object_traits<
//       typename linalg::tensor_traits<
//         typename gmm::linalg_traits<U>::value_type>::scalar_type>::number_type
template <class U, class NT>
simple_object_base* new_list(const U& u, NT);

// special case: not a number-type: create an object-list:
template <class U>
inline simple_object_base* new_list(const U& src, abstract_null_entity);


template <class M>
simple_object_base* new_map(const M& m);


// promotion to more-general number classes:
template <class T>
void extract_number(T& t, const simple_object_base* pobject);

template <class T, class NT>
void dispatch_extract_number(T& t, const simple_object_base* pobject, NT);

template <class T>
void dispatch_extract_number(T& t, const simple_object_base* pobject, complex_number_entity);

template <class T>
void dispatch_extract_number(T& t, const simple_object_base* pobject, real_number_entity);

template <class T>
void dispatch_extract_number(T& t, const simple_object_base* pobject, integral_number_entity);

template <class T>
void extract_complex_number(T& t, const simple_object_base* pobject);

template <class T>
void extract_real_number(T& t, const simple_object_base* pobject);

template <class T>
void extract_integral_number(T& t, const simple_object_base* pobject);

template < >
void extract_integral_number<size_t>(size_t& t, const simple_object_base* pobject);

// ---------------------------------------- end: dispatcher classes: ----------------------------------



} // namespace python_util

#include <python_util_module.h>

// "python_util_inline.h" requires "python_util_module.h"
#include <python_util_inline.h>

#if !defined(EXCLUDE_TEMPLATE_BODIES)
#include <python_util_template.h>
#endif

#endif // __python_util__h
