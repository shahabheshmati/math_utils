#!/usr/bin/env python

import roslib
import rospy
from numpy import zeros, array, asarray, mean
from numpy import sum as npsum

from math import sqrt, atan2, sin, pi, asin, cos
from geometry_msgs.msg import Vector3, Quaternion

# ***** General functions *****
def anglewrap(x):
	y=x
	if x > pi:
		y=x-2.0*pi
	elif x<-pi:
		y=x+2.0*pi

	return y	


# sign function
def sign(x): 
	if x==0:
		return 0
	elif x>0:
		return 1
	else:
		return -1

# saturation
def saturation(val,minval,maxval):
	return max(minval,min(val,maxval))

# ***** Quaternions - Euler - Rotation Matrix *****
#Quaternion norm
def quat_norm (q):
	return sqrt(q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z)

#Quaternion normalize
def quat_normalize (q):
	qn=Quaternion()  	
  	qnorm = quat_norm(q)
  	qn.w = q.w/qnorm
 	qn.x = q.x/qnorm
  	qn.y = q.y/qnorm
  	qn.z = q.z/qnorm

	return qn

#Quaternion multiply
def quat_multiply (q,r):
	qc=Quaternion()
	qc.w = q.w*r.w - q.x*r.x - q.y*r.y - q.z*r.z
	qc.x = q.w*r.x + r.w*q.x + q.y*r.z - q.z*r.y
	qc.y = q.w*r.y + r.w*q.y + q.z*r.x - q.x*r.z
	qc.z = q.w*r.z + r.w*q.z + q.x*r.y - q.y*r.x

	return quat_normalize(qc)

#Quaternion inverse (conj)
def quat_inverse (q):
	q_inv=Quaternion() 
	norm = quat_norm(q)
  	q_inv.w = +q.w / norm
  	q_inv.x = -q.x / norm
  	q_inv.y = -q.y / norm
  	q_inv.z = -q.z / norm

	return q_inv

#Quaternion divide
def quat_divide (q,r):
	return quat_multiply(quat_inverse(r),q)

#Rotate a quaternion by axis angles - used also in quaternion differential equation
def quat_rotation (axis):
	qr=Quaternion()  
	qr.w=1.0
	wnorm = vector3_norm(axis)

	if (wnorm>=1e-7):
		sinhwnorm = sin(0.5*wnorm)
		qr.w = cos(0.5*wnorm)
		qr.x = (axis.x/wnorm)*sinhwnorm
		qr.y = (axis.y/wnorm)*sinhwnorm
		qr.z = (axis.z/wnorm)*sinhwnorm

	return qr

#Rotate a quaternion by a vector and qet the rotated vector
def quat_vector3_rotate (q,vec):
	vecq = Quaternion(vec.x,vec.y,vec.z,0.0)

	qb = quat_multiply(quat_multiply(q,vecq),q_quat_inverse(q))

	return Vector3(qb.x,qb.y,qb.z)

#Quaternion normalize and auto-complete w
def quat_normal_wcomplete (q):
	q.w = sqrt (1.0 - (q.x*q.x + q.y*q.y + q.z*q.z))

#Quaternion duality - get equivalent with positive w
def quat_equiv_wpos_get (q):
	if q.w < 0.0:
		q.w = -q.w
		q.x = -q.x
		q.y = -q.y
		q.z = -q.z

#Quaternion duality - get equivalent with positive w
def quat_equiv_wneg_get (q):
	if q.w > 0.0:
		q.w = -q.w
		q.x = -q.x
		q.y = -q.y
		q.z = -q.z

#Quaternion to Rotation Matrix - (Reb)
def quat2rotmtx (quat):
	rotmtx=zeros((3,3))
  	rotmtx[0][0] = quat.w*quat.w + quat.x*quat.x - quat.y*quat.y - quat.z*quat.z
  	rotmtx[0][1] = 2.0*(quat.x*quat.y - quat.w*quat.z)
  	rotmtx[0][2] = 2.0*(quat.x*quat.z + quat.w*quat.y)
  	rotmtx[1][0] = 2.0*(quat.x*quat.y + quat.w*quat.z)
  	rotmtx[1][1] = quat.w*quat.w - quat.x*quat.x + quat.y*quat.y - quat.z*quat.z
  	rotmtx[1][2] = 2.0*(quat.y*quat.z - quat.w*quat.x)
  	rotmtx[2][0] = 2.0*(quat.x*quat.z - quat.w*quat.y)
  	rotmtx[2][1] = 2.0*(quat.y*quat.z + quat.w*quat.x)
  	rotmtx[2][2] = quat.w*quat.w - quat.x*quat.x - quat.y*quat.y + quat.z*quat.z
	
	return rotmtx

#Rotation Matrix to Quaternion - (Reb)
def rotmtx2quat (rotmtx):
	return eulerZYX2quat(rotmtx2eulerZYX(rotmtx))

'''def rotmtx2quat (rotmtx):
	q = Quaternion()
	q.w = 1.0

	tr = rotmtx[0][0] + rotmtx[1][1] + rotmtx[2][2]

	if tr>0:
		sqtrp1 = sqrt(tr + 1.0);
		sqtrp1x2 = 2.0*sqtrp1;

		q.w = 0.5*sqtrp1;
		q.x = (rotmtx[2][1] - rotmtx[1][2])/sqtrp1x2
		q.y = (rotmtx[0][2] - rotmtx[2][0])/sqtrp1x2
		q.z = (rotmtx[1][0] - rotmtx[0][1])/sqtrp1x2

		return q
	else:
		d0=rotmtx[0][0]
		d1=rotmtx[1][1]
		d2=rotmtx[2][2]

		if ((d1 > d0) and (d1 > d2)):
			sqdip1 = sqrt(d1 - d0 - d2 + 1.0 )
			q.y = 0.5*sqdip1

			if abs(sqdip1)>1e-7:
				sqdip1 = 0.5/sqdip1

			q.w = (rotmtx[0][2] - rotmtx[2][0])*sqdip1
			q.x = (rotmtx[1][0] + rotmtx[0][1])*sqdip1
			q.z = (rotmtx[2][1] + rotmtx[1][2])*sqdip1
			
			return q

		elif (d2 > d0):
			#max value at R(3,3)
			sqdip1 = sqrt(d2 - d0 - d1 + 1.0 )

			q.z = 0.5*sqdip1;

			if abs(sqdip1)>1e-7:
				sqdip1 = 0.5/sqdip1

			q.w = (rotmtx[1][0] - rotmtx[0][1])*sqdip1
			q.x = (rotmtx[0][2] + rotmtx[2][0])*sqdip1
			q.y = (rotmtx[2][1] + rotmtx[1][2])*sqdip1

			return q
			
		else:
			# max value at R(1,1)
			sqdip1 = sqrt(d0 - d1 - d2 + 1.0 )

			q.x = 0.5*sqdip1

			if abs(sqdip1) > 1e-7:
				sqdip1 = 0.5/sqdip1

			q.w = (rotmtx[2][1] - rotmtx[1][2])*sqdip1
			q.y = (rotmtx[1][0] + rotmtx[0][1])*sqdip1
			q.z = (rotmtx[0][2] + rotmtx[2][0])*sqdip1

			return q'''

#Quaternion to Euler-ZYX
def quat2eulerZYX (q):
    euler = Vector3()
    tol = quat2eulerZYX.tolerance
    
    qww, qxx, qyy, qzz = q.w*q.w, q.x*q.x, q.y*q.y, q.z*q.z
    qwx, qxy, qyz, qxz= q.w*q.x, q.x*q.y, q.y*q.z, q.x*q.z
    qwy, qwz = q.w*q.y, q.w*q.z

    test = -2.0 * (qxz - qwy)
    if test > +tol:
        euler.x = atan2 (-2.0*(qyz-qwx), qww-qxx+qyy-qzz)
        euler.y = +0.5 * pi
        euler.z = 0.0

	return euler

    elif test < -tol:
        euler.x = atan2 (-2.0*(qyz-qwx), qww-qxx+qyy-qzz)
        euler.y = -0.5 * pi
        euler.z = tol

	return euler

    else:
        euler.x = atan2 (2.0*(qyz+qwx), qww-qxx-qyy+qzz)
        euler.y = asin (test)
        euler.z = atan2 (2.0*(qxy+qwz), qww+qxx-qyy-qzz)

	return euler
quat2eulerZYX.tolerance=0.99999

#Euler-ZYX to Quaternion
def eulerZYX2quat(euler):
	q = Quaternion()
	cpsi = cos (0.5 * euler.z)
	spsi = sin (0.5 * euler.z)
  	ctheta = cos (0.5 * euler.y)
	stheta = sin (0.5 * euler.y)
  	cphi = cos (0.5 * euler.x)
	sphi = sin (0.5 * euler.x)
  	q.w = cphi*ctheta*cpsi + sphi*stheta*spsi
  	q.x = sphi*ctheta*cpsi - cphi*stheta*spsi
  	q.y = cphi*stheta*cpsi + sphi*ctheta*spsi
  	q.z = cphi*ctheta*spsi - sphi*stheta*cpsi

	return quat_normalize(q)

#Euler-ZYX to Rotation Matrix - (Reb)
def eulerZYX2rotmtx(euler):
	rotmtx=zeros((3,3))
  	cpsi   = cos(euler.z)
	spsi   = sin(euler.z)
  	ctheta = cos(euler.y)
	stheta = sin(euler.y)
  	cphi   = cos(euler.x)
	sphi   = sin(euler.x)
  	#Calculate rotation matrix
  	rotmtx[0][0] = cpsi * ctheta
  	rotmtx[0][1] = sphi * cpsi * stheta - cphi * spsi
  	rotmtx[0][2] = cphi * cpsi * stheta + sphi * spsi
  	rotmtx[1][0] = spsi * ctheta
  	rotmtx[1][1] = sphi * spsi * stheta + cphi * cpsi
  	rotmtx[1][2] = cphi * spsi * stheta - sphi * cpsi
  	rotmtx[2][0] = -stheta
  	rotmtx[2][1] = sphi * ctheta
  	rotmtx[2][2] = cphi * ctheta

	return rotmtx

#Rotation Matrix to Euler-ZYX - (Reb)
def rotmtx2eulerZYX (mtx):

	tolerance = rotmtx2eulerZYX.tolerance
	test = -mtx [2][0]
	
	if test > +tolerance:
        	phi = atan2 (-mtx[1][2], mtx[1][1])
        	theta = +0.5 * pi
        	psi = 0.0
    	elif test < -tolerance:
        	phi = atan2 (-mtx[1][2], mtx[1][1])
        	theta = -0.5 * pi
        	psi = 0.0
    	else:
        	phi = atan2 (mtx[2][1], mtx[2][2])
        	theta = asin (-mtx[2][0])
        	psi = atan2 (mtx[1][0], mtx[0][0])

    	return Vector3(phi, theta, psi)
rotmtx2eulerZYX.tolerance = 0.99999

# ***** Matrix Manipualtion *****
#Create diag block matrix by N nxn matrices
def block_diag (*arrs):
	"""Create a new diagonal matrix from the provided arrays.

	Parameters
	----------
	a, b, c, ... : ndarray
	Input arrays.

	Returns
	-------
	D : ndarray
	Array with a, b, c, ... on the diagonal.
	"""
	arrs = [asarray(a) for a in arrs]
	shapes = array([a.shape for a in arrs])
	out = zeros(npsum(shapes, axis=0))

	r, c = 0, 0
	for i, (rr, cc) in enumerate(shapes):
		out[r:r + rr, c:c + cc] = arrs[i]
		r += rr
		c += cc
	return out

#***** Vector 3 msgs Manipulation *****
# From Vector3 to np.array
def Vector32Array(vector):
	return array([vector.x,vector.y,vector.z],'float').reshape(3,1)

# From np.array to Vector3
def Array2Vector3 (array):
	return Vector3(array[0],array[1],array[2])

#Norm of Vector3
def vector3_norm (vec):
	return sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z)

#Euclidean Distance of two vectors (Vector3)
def vector3_distance (veca,vecb):
	vec = Vector3()
	vec.x = veca.x - vecb.x
	vec.y = veca.y - vecb.y
	vec.z = veca.z - vecb.z
	return vector3_norm(vec)

#Cross product of two vectors (Vector3)
def vector3_cross (a,b):
	c = Vector3()	
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;
	
	return c


class FilterMA():
	def __init__(self,n):
		self.n=n
		self.points = zeros((n,1))
		self.reset()

	def feed(self,p):
		for x in range(0,self.n-1):
			self.points[x] = self.points[x+1]

		self.points[self.n-1] = p

		return mean(self.points)

	def reset(self):
		for x in range(0,self.n):
			self.points[x] = 0.0
