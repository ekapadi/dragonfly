"""
  Python classes for event-driven dynamics simulation to generate maximally-jammed particle packings.
  
  Based on pseudo-code and discussion from the reference::
  
    Boris D. Lubachevsky and Frank H. Stillinger,
    "Geometric Properties of Random Disk Packings",
    Journal of Statistical Physics, Vol. 60, Nos. 5/6, Pp. 561-583, (1990)
  
  $Source: /usr/data0/leipzig_work/tmat_cvs/src/mjp.py,v $
"""

# /* *********************************************************************** */
# /*                                                                         */
# /* Copyright (C) 2010  Kort Travis                                         */
# /*                                                                         */
# /*                                                                         */
# /* This program is free software; you can redistribute it and/or modify    */
# /* it under the terms of the GNU Lesser General Public License as          */
# /* published by the Free Software Foundation; version 2.1 of the License.  */
# /*                                                                         */
# /* This program is distributed in the hope that it will be useful,         */
# /* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
# /* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
# /* GNU Lesser General Public License for more details.                     */
# /*                                                                         */
# /* You should have received a copy of the GNU Lesser General Public        */
# /* License along with this program; if not, write to the Free Software     */
# /* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
# /* USA.                                                                    */
# /*                                                                         */
# /* *********************************************************************** */

import sys
from numpy import *
import scipy as SP
import pylab as PL
import vector_linalg as LA
import parallel_util as PU
import threading
import pdb 


# define parallel python methods _only_ if the module is installed:
global pp_loaded
try:
  import pp as PP
  pp_loaded = True
except ImportError:
  pp_loaded = False
  
_debug_print_ = False


# ----------------------- OPTIMIZATION (pylab versions > 15% initial profile) ---------
def norm(l1):
  ' 2-norm of list elements '
  val = 0.0
  for x in l1:
    val += x * x
  return sqrt(val)  


def dot(l1,l2):
  ' dot product between two lists or tuples '
  # assert len(l1) == len(l2)
  val = 0.0
  for ndim in range(len(l1)):
    val += l1[ndim]*l2[ndim]
  return val  
  
def dist(l1,l2):
  ' distance between points (list or n-tuple parms) ' 
  # assert len(l1) == len(l2)
  val = 0.0
  tmp = 0.0
  for ndim in range(len(l1)):
    tmp = (l1[ndim] - l2[ndim])
    val += tmp * tmp
  return sqrt(val) 
  
# ---------------------- > 10% initial profile ---------------
def product(l1):
  """
  product of list entries 
    : where product([]) \defas 1.0
  """
  val = 1.0
  for x in l1:
    val *= x
  return val 
    
# --------------------------------------------------------------------------------------
  
def row_major_index(indices, shape):
  """Linear index from index-tuple constrained by shape.
  """
  n = 0
  for ndim,n_ in enumerate(indices):
    n += n_ * int(product(shape[ndim+1::]))
  return n  
  
def inverse_row_major_index(n, shape, shift=0):
  """Index-tuple from linear index constrained by shape.
  
  <index list> = inverse_row_major_index(<raw linear index>, <shape>)
    shift => starting offset for index intervals
  """
  indices = zeros((len(shape),),int)
  n_rem = n
  for ndim in range(len(shape)):
    # number of elements per increment at this ndim:
    N_elts = product(shape[ndim+1::])
    x = int(n_rem/N_elts)
    indices[ndim] = int(x + shift)
    n_rem -= x * N_elts
  return indices 

def row_major_index_wrap(indices0, indices1, shape):
  """Tests wrap-around of any sub-index between indices0 and indices1 index tuples.
  
    return values
      test: True if any dimension wraps around
      wrap_dim: +-1 or 0 for each dimension, corresponding to wrap direction
  """
  test = False
  wrap_dim = array([0.0 for n in range(len(indices0))], double)
  for ndim in range(len(indices0)):
    if (indices0[ndim] == shape[ndim]-1 and indices1[ndim] == 0):
      test = True
      wrap_dim[ndim] = 1.0
    elif (indices0[ndim] == 0 and indices1[ndim] == shape[ndim]-1):
      test = True
      wrap_dim[ndim] = -1.0
  return test, wrap_dim  
  
class cell(object):
  """ 
    Local sub-hypercube boundary cell for particle-in-cell simulation.
     
    attributes:
      particle_: <list of particles in the cell>
      neighbor_: <neighboring cells (where offset (0, 0, ...) is current cell)>
      positive_neighbor_: <neighboring cells with >= 0 offset indices (where offset (0, 0, ...) is current cell)>      
      index_: linear index of the cell in its parent cell_array
      wrap_: True => corresponding _positive_ neighbor is across periodic B.C.
      wrap_dim_: +=1 or 0 for each dimension if wrap_
  """
  __slots__ = ['particle_', 'neighbor_', 'positive_neighbor_', 'index_', 'indices_', 'wrap_', 'wrap_dim_']
  
  @staticmethod 
  def in_cube(r, p, cube, use_size=False):
    """
    test = in_cube(p, cube)
    input parms:
      r: <particle radius>
      p: <position n-tuple>
      cube: n-dimensional cube as n-tuple of interval (e.g. ((min_x, max_x), (min_y, max_y), ...) ) 
    """
    if use_size:
      r_ = r
    else:
      r_ = 0.0
    # ------------ OPTIMIZATION (> 5% first profile) ---------------------------
    # return all([(interval[0] <= (p[ndim]-r_)) and ((p[ndim]+r_) <= interval[1]) for ndim, interval  in enumerate(cube)])
    test = True
    for ndim in range(len(cube)):
      interval = cube[ndim]
      if (interval[0] > (p[ndim]-r_)) or ((p[ndim]+r_) > interval[1]):
        test = False
        break
    return test

  
  @staticmethod 
  def N_face(N_dim):
    return N_dim*2
    
  @staticmethod 
  def N_neighbor(N_dim):
    return 3**N_dim - 1

  @staticmethod 
  def N_positive_neighbor(N_dim):
    return 2**N_dim - 1
  
  def neighbor(self, index_, N_dim):
    """ neighboring cell from relative ND index """
    # central offset (cell.N_neighbor(N_dim)/2) is cell_ itself (uninitialized)
    return self.neighbor_[row_major_index(index_, [3 for ndim in range(N_dim)])] 

  def face_neighbor(self, dimension, face, N_dim):
    """ neighboring cell from dimension and face index """
    # relative row-major index: for each dimension -face is "0", cell-center is "1", and +face is "2":
    #   (note: new conditional syntax not in python version 2.4)
    x = [((ndim != dimension) and 1) or ((face == 1) and 2) or 0  for ndim in range(N_dim)]    
    return self.neighbor(x, N_dim)

  def transfer_particle(self, event_, N_dim):
    ' transfer particle to another cell ' 
    self.particle_.remove(event_.particle0_)
    event_.cell_.face_neighbor(event_.dimension_, event_.face_, N_dim).particle_.append(event_.particle0_)
    
  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(cell, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self, cell_number, cell_indices):
    self.particle_ = []
    self.neighbor_ = []           # For uniformity, "neighbor_", "positive_neighbor_", "wrap_" and "wrap_dim_" now initialized via "list.append";
    self.positive_neighbor_ = []  #   note that for convenience during usage, after initialization each will also have a position for the cell itself.
    self.index_ = cell_number
    self.indices_ = cell_indices
    self.wrap_ = []
    self.wrap_dim_ = [] 
     
  def __del__(self):
    for cell_ in self.neighbor_:
      cell_ = []# re-init to [], removing circular ref to self
    for cell_ in self.positive_neighbor_:
      cell_ = []
  
  def __getstate__(self):
    return (self.particle_, self.neighbor_, self.positive_neighbor_, self.index_, self.indices_, self.wrap_, self.wrap_dim_)
    
  def __setstate__(self,state):
    self.particle_ = state[0] 
    self.neighbor_ = state[1]
    self.positive_neighbor_ = state[2]
    self.index_ = state[3]
    self.indices_ = state[4] 
    self.wrap_ = state[5]
    self.wrap_dim_ = state[6]

    
class cell_array(object):
  """ N-dimensional boundary-cell array, represents hypercube.
  """
  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(cell_array, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self, N, N_dim):
    # adjust number of cells so that
    # hypercube side length is an integer:
    
    self.N_dim_ = N_dim
    self.N1_ = int(ceil(double(N)**(1.0/double(N_dim))))
    self.N_ = int(self.N1_**N_dim)
    self.shape = [self.N1_  for ndim in range(self.N_dim_) ]
  
    # initialize cell array:
    self.data_ = [cell(cell_number, inverse_row_major_index(cell_number, self.shape)) for cell_number in range(self.N_)]
    for cell_ in self.data_:
      self.__init_neighbors(cell_)
     
  # note: default "__get_state__" and "__set_state__" should work (i.e. using self.__dict__)

  def __init_neighbors(self, cell_):
    """ Initialize cell neighbor-lists.
     
      neighbor_: _all_ adjacent cells;
      positive_neighbor_: group of adjacent cells "owned" by cell 
        (i.e. which don't overlap with analogous groups from adjacent cells))
    """
    nx = cell_.indices_
    offset_shape = [3 for ndim in range(self.N_dim_)]
    # print 'init neighbors for cell %d' % cell_number
    for neighbor in range(cell.N_neighbor(self.N_dim_)+1):
      # note: initialize neighbor-index (0,0,...) to self
      
      nx_offset = inverse_row_major_index(neighbor, offset_shape, shift=-1)
      
      # application of modulus => periodic boundary conditions for the global hypercube:
      nx_neighbor = [(nx[ndim] + nx_offset[ndim]) % self.shape[ndim]  for ndim in range(self.N_dim_)] 
      cell_.neighbor_.append(self.data_[row_major_index(nx_neighbor, self.shape)])

      # "positive" neighbors: any neighbor in the wrapped 1/2-space corresponding to the 1st dimension:
      if (nx_offset[0] >= 0):
        cell_.positive_neighbor_.append(self.data_[row_major_index(nx_neighbor, self.shape)])
        wrap, wrap_dim = row_major_index_wrap(nx, nx_neighbor, self.shape)
        cell_.wrap_.append(wrap)
        cell_.wrap_dim_.append(wrap_dim)

        
  def cell_edge_length(self, L):
    """Length of boundary-cell edge.
    """
    return L/double(self.N1_)
    
  def cell_center(self, cell_, L):
    """Center position of cell at ND-index.
    
      <position n-tuple> = cell_center(cell_, L)
        input parms:
          cell_: cell for which to find center
          L: length of global boundary hypercube edge
    """
    nx = cell_.indices_
    p = array([((double(nx[ndim]) + 1.0/2.0)/double(self.shape[ndim]))*L for ndim in range(self.N_dim_)])
    return p
  
  def cell_cube(self, cell_, L):
    """Cell interior intervals at ND-index.
     
      <interval n-tuple> = cell_cube(index_, L)
        input parms:
          cell_: cell for which to find boundary cube
          L: length of global boundary hypercube edge
    """
    nx = cell_.indices_
    cube = [((double(nx[ndim])/double(self.shape[ndim]))*L, (double(nx[ndim] + 1)/double(self.shape[ndim]))*L) for ndim in range(self.N_dim_)]
    return cube
 
class state(object):
  """Particle kinemetic state.
  """
  __slots__ = ['position_', 'velocity_', 'time_'] # underscore: direct attribute access
  @staticmethod 
  def __new__(cls):
    result =  super(state, cls).__new__(cls)
    result.position_ = []
    result.velocity_ = []
    result.time_ = []
    return result
  
  def __init__(self, position, velocity, time):
    self.position_ = position
    self.velocity_ = velocity
    self.time_ = time
  
  def __getstate__(self):
    return (self.position_, self.velocity_, self.time_)

  def __setstate__(self, state_):
    self.__init__(*state_)
  
class sphere(object):
  """Sphere physical parameters and kinematic state.
  
     Note on object-oriented analysis:
       "sphere" is a specialization of "particle". However, for this reference implementation
       there is no abstract "particle" base-class.  
  """
  __slots__ = ['radius_', 'dr_', 'state_', 'current_state_', 'index_']  # index_ used as id

  def __deepcopy__(self,memo={}):
    from copy import deepcopy
    result = type(self).__new__(sphere)
    memo[id(self)] = result
    result.__init__(self.radius_, self.dr_, deepcopy(self.current_state().position_, memo), deepcopy(self.current_state().velocity_, memo), id(result))
    return result   

  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(sphere, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self, radius, dr, position, velocity, sphere_number):
    self.radius_ = radius
    self.dr_ = dr
    self.state_ = [state.__new__(state), state.__new__(state)] # old/new state double-buffer
    self.current_state_ = 0
    self.current_state().position_ = position
    self.current_state().velocity_ = velocity
    self.current_state().time_ = 0.0
    self.index_ = sphere_number
    
    
  def __getstate__(self):
    return (self.radius_, self.dr_, self.state_, self.current_state_, self.index_)
    
  def __setstate__(self,state):
    self.radius_ = state[0]
    self.dr_ = state[1]
    self.state_ = state[2]
    self.current_state_ = state[3]
    self.index_ = state[4]
    
  def current_state(self):
    return self.state_[self.current_state_]
    
  def rotate_state_buffer(self):
    self.current_state_ = (self.current_state_ + 1) % 2
      
  def collision_time(self, other):
    """ test, t, jam = collision_time(other) """    
    # from pylab import dist
    s1 = self.current_state() 
    s2 = other.current_state()   
    test = None
    dt = None
    jam = False

    if _debug_print_:    
      # *** DEBUG ***
      print 'collision time (%d -> %d): ' % (self.index_, other.index_)
      print '  dist, r1, r2: %s, %s, %s' % (dist(s1.position_, s2.position_), self.radius_, other.radius_)

    r12_0 = dist(s1.position_, s2.position_) - (self.radius_ + other.radius_)
    if (abs(r12_0) < 1.0*finfo(double).eps):
      r12_0 = 0.0
      
    if (r12_0 >= 0.0):

      n12 = (s2.position_ - s1.position_)
      n12 /= norm(n12)
      # center of mass system:
      v_12 = dot(s2.velocity_ - s1.velocity_, n12) - (self.dr_ + other.dr_)

      # solve the equation:
      #   t: 
      #     r12(t) = |(s1.position_ + s1.velocity_*t) - (s2.position_ + s2.velocity_*t)| -  ((r1 + dr1*t) + (r2 + dr2*t))
      #            = 0
      #     \rightarrow
      #   r12(0) = (v_12:i - dr1 - dr2)*t   # i: intrinsic (i.e. without radial expansion)
      #          = v_{1 2}*t

      if (r12_0 > 0.0):
        # three cases: (1) not moving, (2) moving towards each other, and (3) moving apart:
        test, dt = ((abs(v_12) < finfo(double).eps) and (False, inf)) \
                      or ((v_12 < 0.0) and (True, r12_0/abs(v_12))) \
                      or ((v_12 > 0.0) and (False, inf))
        assert dt >= 0.0
      elif ((self.dr_ + other.dr_) - abs(dot(s2.velocity_ - s1.velocity_, n12)) > 10.0*finfo(double).eps):
        jam = True
        test, dt = True, 0.0
        if _debug_print_:
          # *** DEBUG ***
          print 'GLANCING COLLISION'
      else:
        test, dt = False, inf
        
    else:
      test = False
      dt = 0.0
      
      # Eventually, this should be "unreachable code" and appears to occur rarely due to some type of round-off error. 
      #   Regardless, it represents an ERROR condition: 
      print 'WARNING: OVERLAP(%d, %d): | %s - %s | = %s < %s + %s = %s' % \
        (self.index_, other.index_, s1.position_, s2.position_, dist(s1.position_, s2.position_), self.radius_, other.radius_, self.radius_ + other.radius_)
    
    # assert dist(s1.position_, s2.position_) >= (self.radius_ + other.radius_)    

    if _debug_print_:
      # *** DEBUG ***
      """
      if (isnan(v_12)):
        pdb.set_trace()
      """
      print 'collision time (%d -> %d), dt: %s, %s:' % (self.index_, other.index_, s1.time_ + dt, dt)
      if (r12_0 >= 0.0):
        print '  r12_0, v12, v12_eff: %s, %s, %s' % (r12_0, dot(s2.velocity_ - s1.velocity_, n12), v_12)
      else:
        print '  OVERLAP: r12_0: %s' % r12_0

    return test, s1.time_ + dt, jam    

    
  def collide(self, other, wrap_position):
    """Implement sphere-pair collision.
     
      collide(other, wrap_position)
      input parms:
        other: sphere to collide with 
        wrap_position: other.current_state().position, accounting for periodic B.C. (virtual extension)
      note: assumes positions and times are at correct pre-collision values and only modifies velocity vectors  
    """
    s1 = self.current_state() 
    s2 = other.current_state()  
       
    n12 = (wrap_position - s1.position_)
    n12 /= norm(n12)
    # separate out components perpendicular to collision direction (these will be unmodified):
    v1_p = n12.copy()
    v1_p *= dot(s1.velocity_, n12)
    v1_s = s1.velocity_.copy()
    v1_s -= v1_p
    v2_p = n12.copy()
    v2_p *= dot(s2.velocity_, n12)
    v2_s = s2.velocity_.copy()
    v2_s -= v2_p
    
    # retain perpendicular component:
    v1 = v1_s
    v2 = v2_s
    
    # explicit kinetic energy increment for radial growth:
    v1 -= n12*other.dr_
    v2 += n12*self.dr_
    
    # exchange perpendicular components
    v1 += v2_p     
    v2 += v1_p

    # transfer information to new states:
    self.rotate_state_buffer();
    other.rotate_state_buffer();
   
    self.current_state().__init__(s1.position_, v1, s1.time_)
    other.current_state().__init__(s2.position_, v2, s2.time_) # here _actual_ position is used, _not_ wrap_position!

  def jam(self, other, wrap_position):
    """Implement sphere-pair jamming collision.
       
      jam(other, wrap_position)
      input parms:
        other: sphere to collide with 
        wrap_position: other.current_state().position, accounting for periodic B.C. (virtual extension)
      notes: jam is defined as a glancing collision with ((dr1 + dr2) > (v2 - v1)\cdot n_{1 2} and intersurface distance < epsilon;
        assumes positions and times are at correct pre-collision values and only modifies velocity vectors  
    """
    s1 = self.current_state() 
    s2 = other.current_state()  
       
    n12 = (wrap_position - s1.position_)
    n12 /= norm(n12)
    # separate out components perpendicular to collision direction (these will be unmodified):
    v1_p = n12.copy()
    v1_p *= dot(s1.velocity_, n12)
    v1_s = s1.velocity_.copy()
    v1_s -= v1_p
    v2_p = n12.copy()
    v2_p *= dot(s2.velocity_, n12)
    v2_s = s2.velocity_.copy()
    v2_s -= v2_p
    
    # retain perpendicular component:
    v1 = v1_s
    v2 = v2_s
    
    # explicit kinetic energy increment for radial growth:
    v1 -= n12*other.dr_
    v2 += n12*self.dr_
    
    # zero relative perpendicular component
    # (this case deals with continuous growth-induced collisions, in glancing case)
    v_p = v1_p 
    v_p += v2_p
    v_p /= 2.0
    v1 += v_p
    v2 += v_p
    
    if _debug_print_:
      # *** DEBUG ***
      print 'GLANCE, v1_p, v2_p, v12_eff: %s, %s, %s' % (-n12*other.dr_ + v_p, n12*self.dr_ + v_p, dot((v2-v1),n12) - (self.dr_ + other.dr_))
      
    # transfer information to new states:
    self.rotate_state_buffer();
    other.rotate_state_buffer();
   
    self.current_state().__init__(s1.position_, v1, s1.time_)
    other.current_state().__init__(s2.position_, v2, s2.time_) # here _actual_ position is used, _not_ wrap_position!


  def move(self, t, L):
    """ modify position and radius using current velocities """
    s1 = self.current_state()
    dt = t - s1.time_
    if dt > finfo(double).eps: # dt == 0.0 for any collision partner that has already been moved (rare event)
      s1.position_ += s1.velocity_*dt
      s1.position_ %= L # apply periodic B.C.
      s1.time_ = t
      self.radius_ += self.dr_*dt
    else:
      print 'move: dt < epsilon'     
              
  def cell_exit_time(self, cube):
    """Evaluate sphere cell-exit time to adjacent cells.
    
    test, t, dimension, face = cell_exit_time(cube)
    input parms:  
      cube: n-dimensional cube as n-tuple of interval (e.g. ((min_x, max_x), (min_y, max_y), ...) )
    returns:
      test: intersection exists (i.e. t is finite)
      t: entry time into adjacent cell (entry time => intersection time + epsilon)
      dimension: dimension index (i.e. 1st dimension has index 0)
      face: face index (as interval index)      
    """   
    # assume: 
    #   p, v as 1D-array or n-tuple (i.e. _not_ single row matrices)
    #   len(p) == len(v) == len(cube)
    #   all(interval[0] <= interval[1] for interval in enumerate(cube))

    # only deal with cube-exit case:
    if _debug_print_:
      # *** DEBUG ***
      if not cell.in_cube(self.radius_, self.current_state().position_, cube):
        print '%s not in cube %s' % (self.current_state().position_, cube)
    assert cell.in_cube(self.radius_, self.current_state().position_, cube)  

    test = None
    dt = finfo(double).max
    dimension = None
    face = None
    
    s1 = self.current_state()
    
    # for cell usage, exit event is determined by center-location,
    #   sphere radius is ignored:
    for ndim in range(0,len(cube)):
      x = s1.position_[ndim]
      v = s1.velocity_[ndim]
      interval = cube[ndim]
      face_ = None
      test_ = None
      if abs(v) < finfo(double).eps:
        test_ = False
      elif (v > 0.0):
        test_ = True
        dt_ = (interval[1] - x)/v + finfo(double).eps
        face_ = 1
      else:
        test_ = True
        dt_ = (x - interval[0])/abs(v) + finfo(double).eps
        face_ = 0

      if test_ and (dt_ < dt):
        test = True
        dt = dt_
        dimension = ndim
        face = face_
    
    assert (s1.time_ + dt) > 0.0
    return test, s1.time_ + dt, dimension, face
 


class event(object):
  """Record information about a specific dynamic event. 
  """
  
  # either multiple spheres, or cell-exit information will be valid (not both):
  __slots__ = ['time_', 'cell_', 'particle0_', 'particle1_', 'wrap_position_', 'dimension_', 'face_', 'jam_']  
  
  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(event, cls).__new__(cls, *args, **kwargs)

  """
  @staticmethod
  def new():
    ' work-around for parallel python global name problems '
    return event.__new__(event)
  """
  
  def __init__(self, time, cell0, particle0, particle1=None, wrap_position=None, dimension=None, face=None, jam=None):
    self.time_ = time
    self.cell_ = cell0
    self.particle0_ = particle0
    self.particle1_ = particle1
    self.wrap_position_ = wrap_position
    self.dimension_ = dimension
    self.face_ = face
    self.jam_ = jam

  def __getstate__(self):
    return (self.time_, self.cell_, self.particle0_, self.particle1_, self.wrap_position_, self.dimension_, self.face_, self.jam_)
      
  def __setstate__(self, state):
    self.__init__(*state)
        
  def is_cell_exit(self):
    """ is cell-exit event """
    return self.particle1_ == None
    
  def is_collision(self):
    return self.particle1_ != None

  def debug_print(self):
    if self.is_cell_exit():
      print 'cell exit event at time %s:' % self.time_
      print '  particle %d --> (%d, %d), velocity %s' % \
        (self.particle0_.index_, self.dimension_, self.face_, self.particle0_.current_state().velocity_)
    elif self.is_collision():
      print 'collision event at time %s:' % self.time_
      print '  particles %d, %d; velocities %s, %s' % \
        (self.particle0_.index_, self.particle1_.index_, self.particle0_.current_state().velocity_, self.particle1_.current_state().velocity_)
    else:
      print 'unknown event'
    if (self.particle0_.current_state().time_ == self.time_):
      print 'ZERO TIME STEP'
    if (self.jam_):
      print 'GLANCING COLLISION'
    
class system(object):
  """System of spheres.
  """
  @staticmethod 
  def __new__(cls, *args, **kwargs):
    return super(system, cls).__new__(cls, *args, **kwargs)
  
  def __init__(self, N, L, N_dim, v, dr):
    """
     system(N, L, N_dim, v, dr)
     input parms:
       N: number of spheres
       L: length of global hypercube boundary edge
       N_dim: dimensionality of hypercube boundary 
       v: initial velocity magnitude
       dr: initial radial growth velocity 
       ... other init variants using kwargs ...    
    """
    import warnings
    
    self.__N = N
    self.__L = L
    self.__N_dim = N_dim
    # init cell array, n-dimensional list of spheres:
    self.__cells = cell_array(N, N_dim)
    self.__init_state(dr, v)

    if (dr >= v):
      warnings.warn('radial growth velocity >= thermal velocity',category=UserWarning)
      
    """
    # *** DEBUG ***
    for cell_ in self.__cells.data_:
      print 'cell %d, cube: %s' % (cell_.index_, self.__cells.cell_cube(cell_, self.__L)) 
    """


  # default "__getstate__" and "__setstate__" should work (i.e. using self.__dict__)
  
      
  def __init_state(self, dr, v):
    """ place spheres at initial velocity and radius """  
    # for now, just place a sphere at the center of the first "N" cells:
    sphere_number = 0
    for cell_ in self.__cells.data_:
      if (cell_.index_ >= self.__N):
        break
      nv = (SP.rand(self.__N_dim)*2.0 - 1.0)
      nv /= dot(nv, nv)
      cell_.particle_.append(sphere(0.0, dr, self.__cells.cell_center(cell_, self.__L), nv*v, sphere_number))
      sphere_number += 1 

  def step(self):
    """Evolve system through next event.
    
      <completion test> = step()
    """
    
    jam_list, event_ = self.__next_event() 
    self.__step(jam_list, event_)
    
    return self.complete()

    
  def step_mt(self):
    """Evolve system through next event (multi-threaded version).
     
      <completion test> = step_mt()
    """

    event_results = PU.threaded_block(target=self.__next_event_mt)
    
    # condense threaded results:
    jam_list = []
    event_ = None
    for er in event_results:
      jam_list.extend(er[0])
      if event_:
        test_event = er[1]
        if test_event.time_ < event_.time_:
          event_ = test_event
      else:
        event_ = er[1]
    
    # note: "__step_mt" is sub-threaded:    
    self.__step_mt(jam_list, event_)
    
    return self.complete()

    
  def step_pp(self):
    """Evolve system through next event (parallel-python multi-threaded version).
      
      <completion test> = step_pp()
    """
    import parallel_util, mjp # use fully-qualified classes for parallel python
    
    pp_server = PP.Server(ncpus=PU.worker_thread.NUM_THREADS())
    
    # event_results = PU.pp_threaded_block(pp_server, target=self.__next_event_mt, depfuncs=(mjp.event, mjp.event.__new__), modules=('parallel_util','mjp'), globals=globals())
    
    globals_ = globals().copy()
    globals_.update(locals())

    #  depfuncs=(mjp.event, mjp.event.__new__)
    event_results = PU.pp_threaded_block(pp_server, target=self._system__next_event_mt, depfuncs=(), modules=('parallel_util','mjp'), globals=globals_)
    
    # condense threaded results:
    jam_list = []
    event_ = None
    for er in event_results:
      jam_list.extend(er[0])
      if event_:
        test_event = er[1]
        if test_event.time_ < event_.time_:
          event_ = test_event
      else:
        event_ = er[1]
    
    # note: "__step_pp" is sub-threaded:    
    self.__step_pp(pp_server, jam_list, event_)
    
    del pp_server
    
    return self.complete()        


  def show(self, label=True):
    """.2D snap-shot of system.
    """
    if self.__N_dim > 2:
      raise 'show: dimensions > 2 not supported'

    hfig = PL.figure(); PL.hold(1)
    for cell_ in self.__cells.data_:
      for particle_ in cell_.particle_:
        system.show_particle(particle_, hfig, label)
    
    PL.hold(0)
    PL.draw()
    return hfig    

  @staticmethod
  def show_particle(p, hfig, label=True):
    """Show particle as circle on figure.
    """
    pos = p.current_state().position_
    r = p.radius_
    id = p.index_
    nv = p.current_state().velocity_
    nv /= sqrt(dot(nv,nv))

    if _debug_print_:
      # *** DEBUG ***
      print 'pos, r, nv_%d: %s, %s, %s' % (id, pos, r, nv)
    
    PL.figure(hfig.number)
    phi = linspace(0.0, 2.0*pi, 50)
    x = (cos(phi)*r + pos[0]).tolist(); y = (sin(phi)*r + pos[1]).tolist();
    PL.plot(x,y)
    if label:
      PL.text(pos[0], pos[1], str(id), fontsize=20)
    ha = PL.Arrow(pos[0], pos[1], nv[0]*r, nv[1]*r, width=0.05)
    haxes = PL.gca()
    haxes.add_patch(ha)
    # PL.draw()

  def radius_interval(self):
    """Calculate interval containing particle radii.
    """
    result = [finfo(double).max, 0.0]
    for cell_ in self.__cells.data_:
      for particle_ in cell_.particle_:
        result[0] = min(particle_.radius_, result[0])
        result[1] = max(particle_.radius_, result[1])
    return result
              
  def complete(self):
    """Simulation is complete when particle diameter >= maximum diameter.
    """
    import warnings
    test = False
    ri = self.radius_interval()
    L1 = self.__cells.cell_edge_length(self.__L)
    if 2.0*ri[1] > L1:
      print 'Warning: maximum particle diameter %s exceeds cell dimension %s' % (2.0*ri[1], L1)
      test = True
    return test

    
  def __next_event(self):
    """Generate list of zero-time events and soonest finite-time event.
      
      jam_list, next = __next_event()
    """

    # list of zero-time, jam-adjustment events:
    jam_list = [] 
    
    # soonest finite-time event:
    next = event.__new__(event)
    
    time = finfo(double).max
    # iterate cell array:
    for cell0 in self.__cells.data_:
      """
      # *** DEBUG ***
      print 'scanning cell: %d' % cell0.index_; sys.stdout.flush()
      """
      
      cube = self.__cells.cell_cube(cell0, self.__L)
      for particle0 in cell0.particle_:
        test, time_, dimension, face = particle0.cell_exit_time(cube)
        if (test and (time_ < time)):
          time = time_
          next.__init__(time, cell0, particle0, dimension=dimension, face=face)

        # ----------- OPTIMIZATION: look only at 1/2 of neighbors for each cell: ------------        
        for nn, cell1 in enumerate(cell0.positive_neighbor_):
          """
          # *** DEBUG ***
          print 'scanning neighbor %d: cell %d' % (nn, cell1.index_); sys.stdout.flush()
          """
          # note: cell1 now may be cell0 itself
                  
          for particle1 in cell1.particle_:
            if not (particle1 is particle0):
              wrap_particle1 = self.__wrap_particle(particle1, cell0.wrap_[nn], cell0.wrap_dim_[nn])
              test, time_, jam = particle0.collision_time(wrap_particle1)
              if test:
                if jam:
                  jam_list.append(event(time_, cell0, particle0, particle1, wrap_particle1.current_state().position_, jam=True))
                elif (time_ < time):
                  time = time_
                  next.__init__(time, cell0, particle0, particle1, wrap_particle1.current_state().position_, jam=False)

    if (time < finfo(double).max):
      if _debug_print_:
        # *** DEBUG ***
        print 'next: ',
        next.debug_print()
      else:
        pass    
    else:
      print 'ALL EVENTS ARE JAMMED'
   
    return jam_list, next    

    
  def __next_event_mt(self, kwargs_=None):
    """Generate list of zero-time events and soonest finite-time event (multi-threaded version).
    
      jam_list, next = __next_event_mt(**kwargs)
        **kwargs: presently accepts "NUM_THREADS" and "thread_offset" (required for use with parallel python)
    """

    # ---------------------- parallel python work-around: ----------------------------
    _debug_print_ = False

    import parallel_util as PU
    import mjp
    import scipy as SP
    import numpy as NP

    # parallel python work-around:
    if kwargs_:
      kwargs = kwargs_
    else:
      kwargs = {}
      
    # ---------------------------------------------------------------------------------  
    # list of zero-time, jam-adjustment events:
    jam_list = [] 
              
    # soonest finite-time event:
    next = mjp.event.__new__(mjp.event)
    
    time = NP.finfo(SP.double).max
    # iterate cell array:
    for cell0 in PU.threaded_iterator(self._system__cells.data_, **kwargs): 
      """
      # *** DEBUG ***
      print 'scanning cell: %d' % cell0.index_; sys.stdout.flush()
      """
      
      cube = self._system__cells.cell_cube(cell0, self._system__L)
      for particle0 in cell0.particle_:
        test, time_, dimension, face = particle0.cell_exit_time(cube)
        if (test and (time_ < time)):
          time = time_
          next.__init__(time, cell0, particle0, dimension=dimension, face=face)

        # ----------- OPTIMIZATION: look only at 1/2 of neighbors for each cell: ------------        
        for nn, cell1 in enumerate(cell0.positive_neighbor_):
          """
          # *** DEBUG ***
          print 'scanning neighbor %d: cell %d' % (nn, cell1.index_); sys.stdout.flush()
          """
          # note: cell1 now may be cell0 itself
                  
          for particle1 in cell1.particle_:
            if not (particle1 is particle0):
              wrap_particle1 = self._system__wrap_particle(particle1, cell0.wrap_[nn], cell0.wrap_dim_[nn])
              test, time_, jam = particle0.collision_time(wrap_particle1)
              if test:
                if jam:
                  jam_list.append(event(time_, cell0, particle0, particle1, wrap_particle1.current_state().position_, jam=True))
                elif (time_ < time):
                  time = time_
                  next.__init__(time, cell0, particle0, particle1, wrap_particle1.current_state().position_, jam=False)

    if (time < NP.finfo(SP.double).max):
      if _debug_print_:
        # *** DEBUG ***
        print 'next: ',
        next.debug_print()
      else:
        pass    
    else:
      print 'ALL EVENTS ARE JAMMED'
   
    return jam_list, next    
    
    
  def __wrap_particle(self, particle, wrap, wrap_dim):
    """Project particle state across periodic B.C. 
    """
    from copy import deepcopy
    if not wrap:
      val = particle
    else:
      val = deepcopy(particle)
      
      # *** DEBUG ***
      val.index_ = -particle.index_
      
      # image position across periodic B.C. (i.e. extend virtually):
      val.current_state().position_ += wrap_dim * self.__L
            
    return val  
    
            
  def __step(self, jam_list, event_):
    """Process jam-list and next event.
    """
        
    for jam_event in jam_list:
      jam_event.particle0_.jam(jam_event.particle1_, jam_event.wrap_position_)
    
    # move all particles to kinematic coordinates of finite-time event:
    for cell0 in self.__cells.data_:
      # iterate particles in cell:
      for np, particle0 in enumerate(cell0.particle_):
        """
        if particle0 is event_.particle0_:
          if _debug_print_:
            # *** DEBUG ***
            print '%d: ** event particle **' % particle0.index_
        else:
          if _debug_print_:
            # *** DEBUG ***
            print '%d:    vanilla move' % particle0.index_          
        """
        particle0.move(event_.time_, self.__L)
       
    # actualize event:
    if event_.is_cell_exit(): 
      if _debug_print_:
        # *** DEBUG ***
        print 'cell exit(%d, %d): exit cell %d, enter cell %d at position %s' % \
          (event_.dimension_, event_.face_, event_.cell_.index_, event_.cell_.face_neighbor(event_.dimension_, event_.face_, self.__N_dim).index_,
             event_.particle0_.current_state().position_)

      # transfer particle0 to new cell
      event_.cell_.transfer_particle(event_, self.__N_dim)
    elif event_.is_collision():
      if _debug_print_:
        # *** DEBUG ***
        print 'collision'
 
      event_.particle0_.collide(event_.particle1_, event_.wrap_position_)


  def __step_mt(self, jam_list, event_):
    """Process jam-list and next event (multi-threaded version).
    """
    
    
    PU.threaded_block(target=self.process_jam_events_, args=(jam_list,))
    
    # move all particles to kinematic coordinates of finite-time event:
    PU.threaded_block(target=self.move_particles_, args=(event_,))
       
    # actualize event:
    if event_.is_cell_exit(): 
      if _debug_print_:
        # *** DEBUG ***
        print 'cell exit(%d, %d): exit cell %d, enter cell %d at position %s' % \
          (event_.dimension_, event_.face_, event_.cell_.index_, event_.cell_.face_neighbor(event_.dimension_, event_.face_, self.__N_dim).index_,
             event_.particle0_.current_state().position_)

      # transfer particle0 to new cell
      event_.cell_.transfer_particle(event_, self.__N_dim)
    elif event_.is_collision():
      if _debug_print_:
        # *** DEBUG ***
        print 'collision'
 
      event_.particle0_.collide(event_.particle1_, event_.wrap_position_)


  def __step_pp(self, pp_server, jam_list, event_):
    """Process jam-list and next event (parallel_python multi-threaded version).
    """
    
    globals_ = globals().copy()
    globals_.update(locals()) 
       
    pdb.set_trace()   
    PU.pp_threaded_block(pp_server, target=self.process_jam_events_, args=(jam_list,), modules=('parallel_util','mjp'), globals=globals_)
    
    # move all particles to kinematic coordinates of finite-time event:
    PU.pp_threaded_block(pp_server, target=self.move_particles_, args=(event_,), modules=('parallel_util','mjp'), globals=globals_)
       
    # actualize event:
    if event_.is_cell_exit(): 
      if _debug_print_:
        # *** DEBUG ***
        print 'cell exit(%d, %d): exit cell %d, enter cell %d at position %s' % \
          (event_.dimension_, event_.face_, event_.cell_.index_, event_.cell_.face_neighbor(event_.dimension_, event_.face_, self.__N_dim).index_,
             event_.particle0_.current_state().position_)

      # transfer particle0 to new cell
      event_.cell_.transfer_particle(event_, self.__N_dim)
    elif event_.is_collision():
      if _debug_print_:
        # *** DEBUG ***
        print 'collision'
 
      event_.particle0_.collide(event_.particle1_, event_.wrap_position_)


  def process_jam_events_(self, jam_list, kwargs_=None):
    """Process zero-time jam events.
    """

    # ---------------------- parallel python work-around: ----------------------------
    _debug_print_ = False

    import parallel_util as PU
    import mjp
    import scipy as SP
    import numpy as NP

    # parallel python work-around:
    if kwargs_:
      kwargs = kwargs_
    else:
      kwargs = {}
      
    # ---------------------------------------------------------------------------------  
    
    for jam_event in PU.threaded_iterator(jam_list, **kwargs):
      jam_event.particle0_.jam(jam_event.particle1_, jam_event.wrap_position_)
               
  def move_particles_(self, event_, kwargs_=None):
    """Move all particles to kinematic coordinates of next finite-time event.
    
    move_particles_(event_)
    """

    # ---------------------- parallel python work-around: ----------------------------
    _debug_print_ = False

    import parallel_util as PU
    import mjp
    import scipy as SP
    import numpy as NP

    # parallel python work-around:
    if kwargs_:
      kwargs = kwargs_
    else:
      kwargs = {}
      
    # ---------------------------------------------------------------------------------  

    for cell0 in PU.threaded_iterator(self._system__cells.data_, **kwargs):
      # iterate particles in cell:
      for np, particle0 in enumerate(cell0.particle_):
        particle0.move(event_.time_, self._system__L)    
    
