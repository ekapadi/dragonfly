==================================
*ekapadi/dragonfly/mjp*:
==================================

 .. image:: MJP_2D_1.png
   :height: 500px
   :width: 500 px
   :scale: 50 %
   :alt: 2D MJP example
   :align: right

The *dragonfly/mjp* repository contains codes implementing an event-driven dynamics simulation of maximally-jammed particle packing (MJP).


These codes were initially developed to provide high statistical-quality
*random* aggregate instantiations for use in computational electrodynamic simulations.  The algorithm used
is based on pseudo-code and discussion from the reference::

    Boris D. Lubachevsky and Frank H. Stillinger,
    "Geometric Properties of Random Disk Packings",
    Journal of Statistical Physics, Vol. 60, Nos. 5/6, Pp. 561-583, (1990).

A video representation of the *evolution* of a single algorithm run (for the *3D* case) is viewable here: `MJP instantiation`__ [1]_.

.. __: file:./mjp_simulation.mpg

In order to generate particle-packing instantiations within a specified boundary region, this *hypercube* implementation of the MJP algorithm requires that an additional region-constraining *filter* step be applied to the simulation end-point particle distribution after the simulation's completion.  Clearly this is only a stop-gap solution and a more accurate calculation of volume-constrained particle packing requires extension of the hypercube version of the algorithm.  Such an extension would be straightforward, and would for example allow particle-packing evolution within a volume defined by a closed triangular-mesh. (Such an extension is planned for the near future, however please contact Ekapadi if you have an *urgent* need for such an extension.)

.. [1] For some browser installations, it may be required to download the video file prior to viewing.


MJP algorithm, Python version:
===================================

The original MJP *reference* algorithm was developed in Python.  In order to obtain familiarity with the
algorithm or for cases where only a few MJP simulations are required, the use of the Python version is recommended [2]_.

Source files
------------

  mjp.py:
    MJP ``system`` class and supporting ``event``, ``sphere``, ``state``, ``cell_array``, and ``cell`` classes and associated utility methods.
    (At present, ``multi-threaded`` and ``parallel_python multi-threaded`` versions of these MJP classes are for experimental use only, and should not be used for production simulation.)
    
  vector_linalg.py:
    Miscellaneous linear algebra utilities such as cartesian to spherical co-ordinate transformation and Euler rotation and its inverse.

  parallel_util.py:
    Python implementation of selected C++ ``parallelUtil`` classes.
    
    
MJP algorithm, C++ version:
================================

When a large number of MJP simulations are required, the use of the C++ version of the algorithm is recommended.  This version is about 50 times faster than the Python version.  In addition, the C++ version is multi-threaded and can efficiently utilize multi-core CPU configurations [2]_.

Source files (main ones):
-------------------------

  In order to be consistent with the design goals of the entire T-matrix codebase, 
  the C++ MJP codes are fully-templated with respect to *real* and *complex* number classes.   
  This templating thereby optionally facilitates arbitrary-precision numerical calculation.
  Further, these classes support the ``abstractCommHandle`` class, which provides a unified interface to disk-resident, memory-resident, and MPI inter-process data transfer.
  
  particle_packing.h:
    This header file includes the declaration of the abstract ``system`` class with nested ``state``, ``particle``, ``cell``, and ``event`` classes, as well as the ``mjp_system`` specialization of the ``system`` class to spherical particles. Additional utility routines are provided which allow all structures to be *persistant* both to disk, and with respect to inter-process communication (e.g., for use by MPI).  This latter utility support is consistent with that used by the entire ``TMatrix`` namespace, although the ``particle_packing`` namespace has been deliberately designed to function in a stand-alone manner (i.e., without any dependence on the ``TMatrix`` namespace itself).
    
  particle_packing_template.h:
    Template class method definitions for classes defined in ``particle_packing.h``.
  
  commUtil.h:    
    Template class method definitions supporting a unified low-level binary interface to 
    disk-resident, memory-resident, and MPI inter-process data transfer.

  parallelUtil.h 
    Support classes and methods for object-oriented multi-thread and multi-process coding.

  ntuple.h:
    N-dimensional fully-templated numerical tuple implementation.
  
  linalgUtil.h:
    Miscellaneous linear algebra utility methods.  
  
  ND_matrix.h:
    N-dimensional matrix and associated interpolation classes.
       
External library dependancies:
------------------------------
  In addition to depending on the *Universal Function Interface* (ufi) provided in the present *github* repository at *ekapadi/ufi*, 
  the C++ MJP implementation is also dependent on the *Generic Matrix Methods* (``gmm++``) template library, which is available from `GetFEM++`__ 

.. __: http://download.gna.org/getfem/html/homepage/gmm.html

.. [2] Usage examples for the Python and C++ versions of the MJP codes will be added to the repository as soon as possible.


=====================================================
Concurrent *Programming Languages* benchmark project:
=====================================================

Through adjustment of its initialization parameter values, the computational *concurrency* requirements of the MJP algorithm can be made to range from *embarassingly parallel* computation (i.e., CPU-bound and requiring only rare inter-process communication,) to completely inter-process I/O-bound computation.  As such, it provides an excellent *benchmark* algorithm for comparing the effectiveness of concurrency solutions for multiple computer languages.
To serve this end, development work is in progress to provide reference versions of the MJP algorithm in several additional languages, namely *OCaml*, *Haskell*, and *Clojure*.  As soon as these implementations are available they will be added to this repository, and the associated benchmark results will be published at the http://www.ekapadi.com website. 

