=================================
*ekapadi/dragonfly/ufi*
=================================

The universal function interface (UFI) codes were developed to provide an interface technique to allow the ``numericalFunctor``-derived
classes in the dissertation code-base (see `Optical Scattering From Nanoparticle Aggregates`__), to be quickly interfaced from a scripting language such as Python.  In this manner, the specific initialization details of the large number of numerical simulations required could be conveniently laid out in Python, while taking complete advantage of the highly-optimized C++ code-base running on the MPI cluster.

.. __: http://repositories.tdl.org/tdl-ir/handle/2152/ETD-UT-2010-12-2247


The primary requirement for the UFI design stems from the fact that many of the ``numericalFunctor``-derived classes require instantiation of a *parameters* structure containing *default*-state and *initialization*-state values prior to the actual onset of the computation itself.  The details of these parameters structures are generally quite complex and depend on the specific ``numericalFunctor`` variant used. For example, parameters structures can contain anything from sub-parameter range limits or enumeration values selecting between multiple algorithm-variants, to complete geometry-list information (e.g., in the case of an *aggregate T-matrix* or a *Monte-Carlo* calculation).  

It was found that alternative interfacing techniques, such as Simplified Wrapper and Interface Generator (SWIG), were not well suited to quickly providing interfaces to these complex parameters structures.  In cases where such SWIG-based interfaces would only be used a few times, the time investment in interface development was prohibitive.  In contrast, use of the UFI design merely requires the instantiation of a *dictionary* corresponding to the functor's parameters in Python,  and then implementation of the parameters-initialization and functor-evaluation sequence using the analogous ``std::map``-based container from the C++ side.

Although it was convenient to use an embedded Python interpreter to implement certain aspects of the UFI classes themselves, these classes provide a general method to interface between C++ and *many* other programming languages.


Source files (main ones):
-------------------------

  In order to be consistent with the design goals of the entire T-matrix codebase, 
  the UFI codes are fully-templated with respect to *real* and *complex* number classes.   
  This templating thereby optionally facilitates arbitrary-precision numerical calculation.
  Further, these classes support the ``abstractCommHandle`` class, which provides a unified interface to disk-resident, memory-resident, and MPI inter-process data transfer.
  
  python_util.h, python_util_template.h:
    These header files define the ``simple_object`` template variant, which serves as a type-safe wrapper for <real number>, <complex number>, <integer>, <bool>, <string>, <object list> and <object map>, where the latter are lists and maps of ``simple_object`` itself.  Additional utilities are provided to allow streaming of ``simple_object`` to and from a representation as Python syntax, or as an efficient binary representation.  The ``extensible_parameters`` class then builds on the ``simple_object`` implementation to provide a *generic* function-parameters implementation, and the ``option_map`` class is an analogous *generic* command-line options implementation.
  
  python_util_module.h, python_util_module_template.h, python_util_module.cpp:
    These source files *modularize* the dependence of the ``simple_object`` implementation on the embedded Python interpreter, 
    thereby allowing alternative embedding techniques to be used.  (During initial development work, it was found that ``boost::python`` was *unstable* with respect to
    unloading, and therefore an alternative hard-coded implementation was completed.) 
    Additionally, utility methods are implemented here to allow efficient *by-reference* transfer between C++ and Python/numpy array objects. 
  
