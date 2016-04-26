
// $Source: $

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

#include <assert.h>
#include <stdexcept>

#include <cstddef>
using std::size_t;

#include <type_traits>
#include <algorithm>
#include <limits>
#include <complex>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;
using std::setw;

#include <gmm/gmm.h>

#include "cross_platform.h"
#include "numberTraits.h"
#if 1
// most of the following are required in order to move "ntuple" to "namespace linalg":
using number::epsilon;
using number::zero;
using number::one;
using number::pi;
using number::integer;
using number::ratio;
using number::mod;
using number::conv;
using number::numberTraits;
#endif
#include "factorial.h"


#include "commUtil.h"
using commUtil::abstractCommHandle;

#include "hash_combine.h"
#include "ntuple.h"
#include "ND_matrix.h"

namespace linalg{




} // namespace linalg
