===============================================
*ekapadi/dragonfly* ( |qing1ting2| ):
===============================================

The *dragonfly*  repository is Ekapadi's public target repository for:
  #. Computational electrodynamics code from Kort Travis' dissertation project `Optical Scattering From Nanoparticle Aggregates`__ [#]_.
  #. Small stand-alone open source projects and experiments.
  #. Other open-source projects sponsored by Ekapadi LLC.
  
.. __: http://repositories.tdl.org/tdl-ir/handle/2152/ETD-UT-2010-12-2247

.. [#] As soon as is practical, the complete code-base associated with the dissertation
  will be released into the public domain.  Release specifics depend on the publication scheduling of forthcoming articles.
  
.. |qing1ting2| unicode:: 0x8713 .. dragonfly

Table of Contents
=================

`ekapadi/dragonfly/mjp`__ :
-----------------------------

.. __: https://github.com/ekapadi/dragonfly/tree/master/mjp

The *dragonfly/mjp* directory contains Python and C++ source code implementing an event-driven dynamics simulation of maximally-jammed particle packing (MJP).  Although the statistical physics of such closely-packed structures is interesting in its own right, these structures serve as important *input* material for other numerical modeling efforts ranging from the simulation of advanced solid rocket-fuel compositions to the simulation of nanoparticle aggregate structures as used by the present code-base. 

`ekapadi/dragonfly/ufi`__ :
-----------------------------

.. __: https://github.com/ekapadi/dragonfly/tree/master/ufi

The *dragonfly/ufi* directory contains source code implementing a universal function interface (UFI).  The UFI was developed as an interface technique to allow the ``numericalFunctor``-derived
classes in the dissertation code-base, to be quickly interfaced from a scripting language such as Python.

`ekapadi/dragonfly/base`__ :
------------------------------

.. __: https://github.com/ekapadi/dragonfly/tree/master/base

The *dragonfly/base* directory contains source code shared among the various ``dragonfly`` sub-projects.
