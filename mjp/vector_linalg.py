"""
PYTHON vector linear algebra utilities.  -- parameters assume numpy arrays or matrices (i.e. _not_ ntuples)
  $Source: /usr/data0/leipzig_work/tmat_cvs/src/vector_linalg.py,v $
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


from numpy import *

def min_non_zero(v, eps = finfo(double).eps):
  """ minimum non-zero element """
  result = finfo(double).max
  for d in v:
    if (abs(d) > eps) and (d < result):
      result = d
  return result
 
def radius(p):
  """ radius from cartesian point """
  return sqrt(sum(p**2,0))

def transverse_radius(p):
  """ transverse radius from cartesian point """
  return sqrt(sum(p[0:2]**2,0))


def cart2sphere_(x, y, z):
  """
  r, t, p = cart2sphere_(x, y, z)
     spherical coordinate point from cartesian 
  """  
  
  # spherical coordinate point <-- (r, <zenith angle>, <azimuth angle>)
  r = sqrt(x**2 + y**2 + z**2).real; t = 0.0; p = 0.0 

  # special case flags:
  zx = (x == 0.0); zy = (y == 0.0); zz = (z == 0.0) 

  if (zx and zy and zz): # r==0
    r = 0.0
    t = 0.0
    p = 0.0 # define this case.
  else:
    if (zx and zy): # rho==0
      p = 0.0 # define this case
    else:
      p = math.atan2(y, x) # domain: \phi \in [0, 2*\pi]
      if (0.0 > p):
        p += 2.0*pi
        if (p >= 2.0*pi):
          p = 0.0 # this catches case: <x or y> < 0 && < x or y within eps of zero >

    t = math.acos(z/r)  

  return r, t, p  
  

def sphere2cart_(r, t, p):
  """
  x, y, z = sphere2cart_(r, t, p)
    cartesian coordinate point from spherical
  """

  # cartesian coordinates <-- (x, y, z)
  x = 0.0; y = 0.0; z = 0.0

  rho = r*sin(t)
  x = rho*cos(p)
  y = rho*sin(p)
  z = r*cos(t)

  return x, y, z
  
  
def cart2sphere(cp):
  """
  sp = cart2sphere(cp)
     spherical coordinate point[s] from cartesian 
  """
  if (len(cp.shape)==1) or (cp.shape[0]==1):
    cp_ = asmatrix(cp) # unify indexing as (1,3) matrix
    # define only for single 3D points:
    if (cp_.shape[1] != 3):
      raise 'cart2sphere: parm must be 3D point'

     # return value of same type as input parm:
    sp = zeros_like(cp)  # sp <-- (r, <zenith angle>, <azimuth angle>)
    sp_ = asmatrix(sp) # unify indexing as (1,3) matrix
    sp_[0,0], sp_[0,1], sp_[0,2] = cart2sphere_(cp_[0,0], cp_[0,1], cp_[0,2])
  else:
    # multiple points:
    sp = zeros_like(cp)
    for np    in    range(cp.shape[0]):
      sp[np,:] = cart2sphere(cp[np,:])
        
  return sp
  
  
def sphere2cart(sp):
  """
  cp = sphere2cart(sp)
    cartesian coordinate point[s] from spherical
  """

  if (len(sp.shape)==1) or (sp.shape[0]==1):
    sp_ = asmatrix(sp) # unify indexing as (1,3) matrix
    # define only for single 3D points:
    if (sp_.shape[1] != 3):
      raise 'sphere2cart: parm must be 3D point'

     # return value of same type as input parm:
    cp = zeros_like(sp) # cp <-- (x, y, z)
    cp_ = asmatrix(cp) # unify indexing as (1,3) matrix
    cp_[0,0], cp_[0,1], cp_[0,2] = sphere2cart_(sp_[0,0], sp_[0,1], sp_[0,2])
  else:
    # multiple points:
    cp = zeros_like(sp)
    for np    in    range(sp.shape[0]):
      cp[np,:] = sphere2cart(sp[np,:])

  return cp
  
  
"""  
#if 1 // moved from numericalConstants_template.h because of methods using std::vector:
template <class T>
inline void cart2sphere(const T& x, const T& y, const T& z,
                   T& r, T& theta, T& phi)
{
 typedef typename numberTraits<T>::magnitudeType R;
 
 // allow in-place:
 R r_(real(sqrt( x*x + y*y + z*z ))), theta_(zero<R>()), phi_(zero<R>());

 // special case flags:
 bool zx( zero<T>()==x ), zy( zero<T>()==y ), zz( zero<T>()==z );
 if (zx&&zy&&zz){ // r==0
   theta_ = zero<R>();
   phi_ = zero<R>(); // define this case.
 }
 else{ // rho==0
  if (zx&&zy){
    phi_ = zero<R>();
  }
  else{ 
    phi_ = real( atan2(y, x) ); // domain: \phi \in [0, 2*\pi]
    if (zero<R>() > phi_){
			const R pi2(integer<R>(2)*pi<R>());
			phi_ += pi2;
      if (phi_ >= pi2)    // this catches case: <x or y> < 0 && < x or y within eps of zero >
			  phi_ = zero<R>();  
		}
  }
  theta_ = real( acos(z/r_) );
 }
 
 // transfer to return values:
 r = r_;
 theta = theta_;
 phi = phi_; 
}

template <class T>
inline void sphere2cart(const T& r_, const T& theta_, const T& phi_,
                   T& x, T& y, T& z)
{
 // allow in-place
 T r(r_), theta(theta_), phi(phi_);
 
 T rho( r*sin(theta) );
 x = rho*cos(phi);
 y = rho*sin(phi);
 z = r*cos(theta);
}
#endif
"""


def cartVect2sphere_(x, y, z, v_x, v_y, v_z):
  """
  r, t, p, v_r, v_t, v_p = cartVect2sphere_(x, y, z, v_x, v_y, v_z)
    spherical coordinate position and vector components from cartesian 
  """
  sphere_pos, sphere_vec = cartVect2sphere(matrix([x,y,z]), matrix([v_x, v_y, v_z]))
  return sphere_pos[0,0], sphere_pos[0,1], sphere_pos[0,2], \
         sphere_vec[0,0], sphere_vec[0,1], sphere_vec[0,2]
         

def sphereVect2cart_(r, t, p, v_r, v_t, v_p):
  """
  x, y, z, v_x, v_y, v_z = sphereVect2cart_(r, t, p, v_r, v_t, v_p)
    cartesian coordinate position and vector components from spherical 
  """
  cart_pos, cart_vec = sphereVect2cart(matrix([r, t, p]), matrix([v_r, v_t, v_p]))
  return cart_pos[0,0], cart_pos[0,1], cart_pos[0,2], \
         cart_vec[0,0], cart_vec[0,1], cart_vec[0,2]

  
def cartVect2sphere(cart_pos, cart_vec):
  """
  sphere_pos, sphere_vec = cartVect2sphere(cart_pos, cart_vec)
    spherical coordinate position and vector components from cartesian 
  """
  sphere_pos = cart2sphere(cart_pos)
  sphere_pos_ = asmatrix(sphere_pos) # unify indexing as (1,3) matrix
  r = sphere_pos_[0,0]; t = sphere_pos_[0,1]; p = sphere_pos_[0,2]  
  st = sin(t); ct = cos(t)
  sp = sin(p); cp = cos(p)
  m_cart2sphere = matrix([[st*cp, st*sp, ct], 
                         [ct*cp, ct*sp, -st], 
                         [-sp, cp, 0.0]])
  sphere_vec = m_cart2sphere*matrix(cart_vec,copy=False).T
  # match return type to parm type:
  if isinstance(cart_vec,matrix):
    sphere_vec_ = sphere_vec.T
  else:
    sphere_vec_ = asarray(sphere_vec.T).reshape(cart_vec.shape) # reshape allows mono-shape input arrays  
  return sphere_pos, sphere_vec_

def sphereVect2cart(sphere_pos, sphere_vec):
  """
  cart_pos, cart_vec = sphereVect2cart(sphere_pos, sphere_vec)
    cartesian coordinate position and vector components from spherical 
  """
  cart_pos = sphere2cart(sphere_pos)
  sphere_pos_ = asmatrix(sphere_pos) # unify indexing as (1,3) matrix
  r = sphere_pos_[0,0]; t = sphere_pos_[0,1]; p = sphere_pos_[0,2]  
  st = sin(t); ct = cos(t)
  sp = sin(p); cp = cos(p)
  m_sphere2cart = matrix([[st*cp, ct*cp, -sp],
                          [st*sp, ct*sp, cp],
                          [ct, -st, 0.0]]) 
  cart_vec = m_sphere2cart*matrix(sphere_vec,copy=False).T
  # match return type to parm type:
  if isinstance(sphere_vec,matrix):
    cart_vec_ = cart_vec.T
  else:
    cart_vec_ = asarray(cart_vec.T).reshape(sphere_vec.shape) # reshape allows mono-shape input arrays  
  return cart_pos, cart_vec_

   
"""
// for 3D vectors:
template <class T>
inline void cartVect2sphere(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV)
{
  typedef typename numberTraits<T>::magnitudeType R;
  
  cart2sphere(vSrcP,vDestP);

  // see note at "cart2sphere" w.r.t. complex coordinates; at this point the angular coordinates will be real-valued:  
  const R st(sin(real(vDestP[1]))), ct(cos(real(vDestP[1]))),
	  sp(sin(real(vDestP[2]))), cp(cos(real(vDestP[2])));
	
	gmm::dense_matrix<T> aCartToSphere(3,3);  // could be of type "R", but leave as "T" to simplify "mult" below
	/*
  aCartToSphere = [ st*cp, st*sp,  ct ;...
                    ct*cp, ct*sp,  -st;...
                    -sp,   cp,     0           ];
	*/
	aCartToSphere(0,0) = st*cp; aCartToSphere(0,1) = st*sp; aCartToSphere(0,2) = ct;
	aCartToSphere(1,0) = ct*cp; aCartToSphere(1,1) = ct*sp; aCartToSphere(1,2) = -st;
	aCartToSphere(2,0) = -sp;   aCartToSphere(2,1) = cp;    aCartToSphere(2,2) = zero<T>();
	
	if (vDestV.size() != 3) // allows in-place (via gmm, in-place...)
	  vDestV.resize(3,zero<T>());
	gmm::mult(aCartToSphere, vSrcV, vDestV);									
}	
                 
template <class T>
inline void sphereVect2cart(const std::vector<T>& vSrcP, const std::vector<T>& vSrcV, 
                            std::vector<T>& vDestP, std::vector<T>& vDestV)    
{
  sphere2cart(vSrcP,vDestP);

  // see note at "sphere2cart" w.r.t. complex coordinates; here the angular coordinates should be real-valued:  

  const R st(sin(real(vSrcP[1]))), ct(cos(real(vSrcP[1]))),
	  sp(sin(real(vSrcP[2]))), cp(cos(real(vSrcP[2])));
	
	gmm::dense_matrix<T> aSphereToCart(3,3); // could be of type "R", but leave as "T" to simplify "mult" below
	/*
    aSphereToCart = [ sin(t)*cos(p), cos(t)*cos(p),       -sin(p);...
                     sin(t)*sin(p), cos(t)*sin(p),        cos(p);...
                     cos(t),       -sin(t),             0];	
  */
	aSphereToCart(0,0) = st*cp; aSphereToCart(0,1) = ct*cp; aSphereToCart(0,2) = -sp;
	aSphereToCart(1,0) = st*sp; aSphereToCart(1,1) = ct*sp; aSphereToCart(1,2) = cp;
	aSphereToCart(2,0) = ct;    aSphereToCart(2,1) = -st;   aSphereToCart(2,2) = zero<T>();
	
	if (vDestV.size() != 3) // allows in-place (via gmm, in-place...)
	  vDestV.resize(3,zero<T>());

	gmm::mult(aSphereToCart, vSrcV, vDestV);	
}	
		
"""

def cartMesh_sphereVect2cart(aXX, aYY, aZZ, avRR, avTT, avPP):
  """
  avXX, avYY, avZZ = cartMesh_sphereVect2cart(aXX, aYY, aZZ, avRR, avTT, avPP)
    cartesian vector components from spherical components on a cartesian mesh
  """
  avXX = zeros_like(avRR); avYY = zeros_like(avRR); avZZ = zeros_like(avRR);
  for nx in ndindex(aXX.shape):
    x,y,z = aXX[nx], aYY[nx], aZZ[nx]
    r,t,p = cart2sphere_(x,y,z)
    x,y,z,avXX[nx],avYY[nx],avZZ[nx] = sphereVect2cart_(r,t,p,avRR[nx],avTT[nx],avPP[nx])
    
  return avXX, avYY, avZZ
  
  
def cartMesh_cartVect2sphere(aXX, aYY, aZZ, avXX, avYY, avZZ):
  """
  avRR, avTT, avPP = cartMesh_cartVect2sphere(aXX, aYY, aZZ, avXX, avYY, avZZ)
    spherical vector components from cartesian components on a cartesian mesh
  """
  avRR = zeros_like(avXX); avTT = zeros_like(avXX); avPP = zeros_like(avXX);
  for nx in ndindex(aXX.shape):
    x,y,z = aXX[nx], aYY[nx], aZZ[nx]
    r,t,p,avRR[nx],avTT[nx],avPP[nx] = cartVect2sphere_(x,y,z,avXX[nx],avYY[nx],avZZ[nx])
    
  return avRR, avTT, avPP
    

def apply_rotation(aRot, p):
  """ apply rotation matrix while preserving form of input point array/matrix to destination """
  if (len(p.shape)==1) or (p.shape[0]==1):
    pp_ = (aRot*matrix(p,copy=False).T).T               
  else:
    # vectorize multiple points:
    # (note: no breakout here to allow reuse of aRot)
    pp_ = asmatrix(zeros(p.shape,p.dtype))
    for np in range(p.shape[0]):
      pp_[np,:] = (aRot*matrix(p[np,:],copy=False).T).T
    # match return type to parm type:

  if isinstance(p,matrix):
    pp = pp_
  else:
    pp = asarray(pp_).reshape(p.shape) # reshape allows mono-shape input array

  return pp    


def euler_rotation_matrix(alpha, beta, gamma):
  """
    generate Euler rotation matrix
    to rotate cartesian coordinates following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
  """
  ca = cos(alpha); sa = sin(alpha)
  cb = cos(beta); sb = sin(beta)
  cg = cos(gamma); sg = sin(gamma)
  aRot = matrix([[cg*cb*ca - sg*sa, cg*cb*sa + sg*ca, -cg*sb],
                 [-sg*cb*ca - cg*sa, -sg*cb*sa + cg*ca, sg*sb], # 19.08.2011: found transcription-error from Arfken, sign of this term

                 [sb*ca, sb*sa, cb]])
  return aRot

 
def euler_rotation(p, alpha, beta, gamma):
  """
  <cartesian position>' = euler_rotation(<cartesian position>, alpha, beta, gamma)
    rotate cartesian coordinates following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
  """
  aRot =  euler_rotation_matrix(alpha, beta, gamma)
  return apply_rotation(aRot, p)    


def euler_rotation_(x, y, z, alpha, beta, gamma):
  """
  x', y', z' = euler_rotation_(x, y, z, alpha, beta, gamma)
    rotate cartesian coordinates following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
  """
  pp = euler_rotation(matrix([x, y, z]), alpha, beta, gamma)
  return pp[0,0], pp[0,1], pp[0,2]


def inverse_euler_rotation_matrix(alpha, beta, gamma):
  """
  generate inverse Euler rotation matrix
    to produce the inverse of the cartesian coordinate rotation following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
    (for the _forward_ transformation)
  """
  ca = cos(alpha); sa = sin(alpha)
  cb = cos(beta); sb = sin(beta)
  cg = cos(gamma); sg = sin(gamma)
  aInvRot = matrix([[cg*cb*ca - sg*sa, sg*cb*ca - cg*sa, sb*ca], # 19.08.2011: see "euler_rotation_matrix" 
                                                                 #   re changes due to  transcription-error from Arfken (not yet verified...)
                 [cg*cb*sa + sg*ca, -sg*cb*sa + cg*ca, sb*sa],
                 [cg*sb, sg*sb, cb]])
  return aInvRot
  
  
def inverse_euler_rotation(p, alpha, beta, gamma):
  """
  <cartesian position>' = inverse_euler_rotation(<cartesian position>, alpha, beta, gamma)
    inverse of the cartesian coordinate rotation following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
    (for the _forward_ transformation)
  """
  aInvRot =  inverse_euler_rotation_matrix(alpha, beta, gamma)
  return apply_rotation(aInvRot, p)
 

def inverse_euler_rotation_(x, y, z, alpha, beta, gamma):
  """
  x', y', z' = inverse_euler_rotation_(x, y, z, alpha, beta, gamma)
    inverse of the cartesian coordinate rotation following Z(alpha), X'(beta), Z''(gamma) rotation convention
    positive angles produce +helical (counter-clockwise) rotation of axes or -helical rotation of coordinates
    (for the _forward_ transformation)
  """
  pp = inverse_euler_rotation(matrix([x, y, z]), alpha, beta, gamma)
  return pp[0,0], pp[0,1], pp[0,2]
  
