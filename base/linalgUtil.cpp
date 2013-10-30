// $Source: /usr/data0/leipzig_work/tmat_cvs/src/linalgUtil.cpp,v $

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

#if defined(__USE_PTHREAD)
  #include <pthread.h>
#endif

#if defined(__USE_MPI)
  #include <mpi.h>
#endif
  
#if defined(__USE_OMP) 
  #include <omp.h>
#endif

#if defined(__NO_THREAD__)
  #define __THREAD_ATTRIBUTE__
#else
  #define __THREAD_ATTRIBUTE__ __thread
#endif


#include <assert.h>
#include <errno.h>
// #include <typeinfo.h>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <locale>
#include <cmath>
#include <string>
#include <complex>
#include <vector>
#include <list>

#include <functional>
#include <algorithm>
#include <numeric>
#if defined(__PGI)
  // TACC port: see discussion in "portability.h" at definition of _STL_EXT_NAMESPACE_ macro:
  // #include <slist>
  #include <hash_map>
#else
  #include <ext/algorithm>
  #include <ext/numeric>
  // #include <ext/slist>
  #include <ext/hash_map>
#endif

#include <sstream>
#include <typeinfo>
#include <tr1/type_traits>

// #define __verbose__ 0

#include <commUtil.h>
using commUtil::abstractCommHandle;
using commUtil::fileHandle;
using commUtil::processCommHandle;
using commUtil::open;
using commUtil::read;
using commUtil::write;
using commUtil::close;
using commUtil::readBinary;
using commUtil::writeBinary;

#include <parallelUtil.h>
#include <statusUtil.h>

#if 0
  #include <scanner_util.h>
  #include <python_util.h>
#endif

#include <md5sum.h>


#include <tmp_src/ProcessFileParms.h>
#include <tmp_src/Status.h>

#include <numericalConstants.h>
#if 0
  #include <numericalFunctor.h>
  #include <md5sum_template_forward.h>
#endif

#include <gmm/gmm_kernel.h>
#include <gmm/gmm_dense_lu.h>

#include <ghostIterator.h>
#include <gmm_ext.h>

#include <ntuple.h>
using linalg::ntuple;

#include <linalgUtil.h>

#if 0 // -- delta --
  #include <numericalFunctor.h>
  #include <md5sum_template_forward.h>
#endif

#if 0
	#include <multipoleSextuple.h>
	#include <besselFunctor.h>
	#include <stirlingsAppxFunctor.h>
	#include <gauntFunctor.h>
	#include <wigner3JFunctor.h>
	#include <legendreFunctor.h>
	#include <jacobiFunctor.h>
	#include <vectorSphericalHarmonicFunctor.h>
#endif

#if 0 // -- delta --
#include <eulerRotationFunctor.h>
#include <factorialFunctor.h>
#endif

using number::zero;
using number::one;
using number::conj;

using std::distance;
using std::abs;
using std::conj;

#if 0 // -- delta --
#include <regionSelector.h>
#endif
#include <ghostIterator.h>

using std::cout;
using std::endl;
using std::setw;
using std::setprecision;

namespace linalg{

// linear index from N-dimensional row-major indices in shape (i.e. dimension limits):
//   here I assume indices are "size_t"
size_t row_major_index(const std::vector<size_t>& indices, const std::vector<size_t>& shape)
{
  size_t n(0), ndim(0);
  for(std::vector<size_t>::const_iterator itN = indices.begin(), itNEnd = indices.end();
      itN != itNEnd;
      ++itN, ++ndim){
    n += (*itN) * product(gmm::sub_vector(shape, gmm::sub_interval(ndim + 1, (gmm::vect_size(shape) - ndim - 1) )));
  }
  return n;      
}


// N-dimensional row-major indices (constrained by shape) from linear index:  
void inverse_row_major_index(size_t n, const std::vector<size_t>& shape, std::vector<size_t>& indices)
{
  indices.clear();
  indices.resize(shape.size(),0);
  size_t n_rem(n);
  for(size_t ndim = 0; ndim < shape.size(); ++ndim){
    // number of elements per increment at this ndim:
    size_t N_elts(product(gmm::sub_vector(shape, gmm::sub_interval(ndim + 1, (gmm::vect_size(shape) - ndim - 1) ))));
    indices[ndim] = n_rem/N_elts;  // truncation
    n_rem -= indices[ndim]*N_elts; 
  }
}

// increment a specified dimension of a set of indices subject to an interval-list constraint:
// (false return value => increment went out-of-range;
//    specifically: out-of-range return leaves highest dimension index at end-of-range, 
//    others at start-of-range positions)
bool increment_row_major_indices(const std::vector<size_t>& src, size_t ndim, 
                                 const std::vector<std::pair<size_t,size_t> >& intervals, 
                                 std::vector<size_t>& dest)
{
  bool status(true);
  if (&src != &dest)
    dest = src;
  dest[ndim] += 1;
  if (dest[ndim] >= intervals[ndim].second){
    // all indices except highest dimension revert to start-of-range:
    // (highest dimension will be "carry" marker)
    if (ndim > 0)
      dest[ndim] = intervals[ndim].first;
    if (!((ndim > 0) && increment_row_major_indices(dest, ndim-1, intervals, dest)))
      status = false; // out-of-range
  }
  return status;  
}

} // namespace linalg
