#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFRAP.

#PyFRAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

"""PyFRAP module for simple geometric operations, such as

	* Angle computation.
	* Normal computation.
	
"""
#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#String
import string

#PyFRAP modules
from .pyfrp_term_module import *
from . import pyfrp_misc_module

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def getAngle(vec1,vec2):
	
	r"""Returns angle between two vectors in radians.
	
	Angle is calculated via
	
	.. math:: \phi = \frac{v_1 \dot v_2}{|v_1| |v_2|}
	
	.. note:: Checks for numerical errors and corrects them if necessary.
	
	Args:
		vec1 (numpy.ndarray): Vector 1.
		vec2 (numpy.ndarray): Vector 2.
	
	Returns: 
		float: Angle.
	
	"""

	x=np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
	
	if x>1.0:
		if x>1.+1E-4:
			printError("getAngle: x out of bounds.")
		x=1.0
	
	if x==1.0:
		return 0.0
	if x==-1.0:
		return np.pi
	
	a=np.arccos(x)
	
	if a<0:
		return getAngle(vec2,vec1)

	return a

def flipCoordinate(x,destAxis,origAxis='x',debug=False):
	
	"""Transforms coodinate from one axis to another by
	rolling the coordinates, e.g. clockwise turning the 
	point.
	
	``destAxis`` and ``origAxis`` are given as one of 
	``x,y,z``. 
	
	Args:
		x (numpy.ndarray): Coordinate to turn.
		destAxis (str): Destination axis.
		
	Keyword Args:	
		origAxis (str): Original axis.
		debug (bool): Print debugging output.
	
	Returns:
		numpy.ndarray: Transformed coordinate.
	
	"""
	
		
	# Calculate differences between axis
	try:
		axisDiff=abs(string.lowercase.index(destAxis)-string.lowercase.index(origAxis))
	except AttributeError:
		axisDiff=abs(string.ascii_lowercase.index(destAxis)-string.ascii_lowercase.index(origAxis))
		
		
		
	# Roll
	xnew=np.roll(x,axisDiff)
	
	# Print debugging messages
	if debug:
		print("Transforming coordinate " , x, " from axis ", origAxis, " to axis ", destAxis , ".")
		print("axisDiff = ", axisDiff)
		print("xnew = ", xnew)
	
	return xnew 

def newellsMethod(vertices):
	
	"""Computes normal using Newell's method.
		
	Adapted from http://stackoverflow.com/questions/39001642/calculating-surface-normal-in-python-using-newells-method.
	
	Vertices should be given as 
	
	>>> vertices=[[x1,y1,z1],[x2,y2,z2],...]
	
	Args:
		vertices (list): List of vertex coordinates.
	
	Returns:
		numpy.ndarray: Normal vector to surface.
	"""
	
	#Newell's method
	n = [0.0, 0.0, 0.0]

	for i, x in enumerate(vertices):
		
		x2 = vertices[(i+1) % len(vertices)]
		
		nBefore=list(n)
		
		n[0] += (x[1] - x2[1]) * (x[2] + x2[2])
		n[1] += (x[2] - x2[2]) * (x[0] - x2[0])
		n[2] += (x[0] - x2[0]) * (x[1] - x2[1])
		
		
		#This is for closed loops that somehow return a zero normal vector.
		#Don't know if this is a good workaround
		if sum(nBefore)!=0 and sum(n)==0:
			n=list(nBefore)
		
	normalised = [i/sum(n) for i in n]
	
	return normalised

def computeNormal(vertices,method="cross"):
	
	"""Computes normal.

	Vertices should be given as 
	
	>>> vertices=[[x1,y1,z1],[x2,y2,z2],...]
	
	Currently there are two methods available:
		
		* ``cross``, see also :py:func:`normalByCross`.
		* ``newells``, see also :py:func:`newells`.
	
	If method is unknown, will fall back to ``cross``.
	
	Args:
		vertices (list): List of vertex coordinates.
	
	Keyword Args:
		method (str): Method of normal computation.
	
	Returns:
		numpy.ndarray: Normal vector.
	"""
	
	if method.lower() not in ["cross","newells"]:
		printWarning("Unknown method " + method + " for normal computation. Will use method=cross." )
		method=cross
		
	if method.lower()=="cross":
		vec1=np.array(vertices[1])-np.array(vertices[0])
		vec2=np.array(vertices[2])-np.array(vertices[0])

		return normalByCross(vec1,vec2)
		
	elif method.lower()=="newells":
		return newellsMethod(vertices)
		
def normalByCross(vec1,vec2):
	
	r"""Returns normalised normal vectors of plane spanned by two vectors.
	
	Normal vector is computed by:
	
	.. math:: \mathbf{n} = \frac{\mathbf{v_1} \times \mathbf{v_2}}{|\mathbf{v_1} \times \mathbf{v_2}|}
	
	.. note:: Will return zero vector if ``vec1`` and ``vec2`` are colinear.
	
	Args:
		vec1 (numpy.ndarray): Vector 1.
		vec2 (numpy.ndarray): Vector 2.
	
	Returns:
		numpy.ndarray: Normal vector.
	"""
	
	if checkColinear(vec1,vec2):
		
		printWarning("Can't compute normal of vectors, they seem to be colinear. Returning zero.")
		return np.zeros(np.shape(vec1))
		
	return np.cross(vec1,vec2)/np.linalg.norm(np.cross(vec1,vec2))
	
def checkColinear(vec1,vec2):
	
	"""Returns True if two vectors are colinear.
	
	Args:
		vec1 (numpy.ndarray): Vector 1.
		vec2 (numpy.ndarray): Vector 2.
	
	Returns:
		bool: True if colinear.
	
	"""
	
	return sum(np.cross(vec1,vec2))==0
 
def getRotMatrix(n1,n2):

	"""Builds rotation matrix for the rotation of a vector n2 
	onto n1.
	
	Taken from http://stackoverflow.com/questions/9423621/3d-rotations-of-a-plane .
	
	Args:
		n1 (numpy.ndarray): Vector to be rotated to.
		n2 (numpy.ndarray): Vector to rotate.
		
	Returns:
		numpy.ndarray: Rotation matrix

	"""
	
	# Rotation axis
	rotAxis=np.cross(n2,n1)
	rotAxis=rotAxis/np.linalg.norm(rotAxis)
	
	#Rotation angle
	c=np.dot(n2,n1)/(np.linalg.norm(n2)*np.linalg.norm(n1))
	s=np.sqrt(1-c**2)
	C=1-c
	
	x=rotAxis[0]
	y=rotAxis[1]
	z=rotAxis[2]
	
	# Build rotation matrix
	rmat = np.array([[x*x*C+c,x*y*C-z*s,x*z*C+y*s],[y*x*C+z*s,y*y*C+c,y*z*C-x*s],[z*x*C-y*s,z*y*C+x*s,z*z*C+c]])
	
	return rmat.T

def decodeEuclideanBase(d):
	
	"""Decodes a euclidean base vector given as a literal.
	
	Example:
	
	>>> decodeEuclideanBase('z')
	>>> array([ 0.,  0.,  1.])
	
	Args:
		d (str): Direction ("x"/"y"/"z")
		
	Returns:
		numpy.ndarray: Base vector.
	
	"""
	
	if d=='x':
		return np.array([1.,0.,0.])
	elif d=='y':
		return np.array([0.,1.,0.])
	elif d=='z':
		return np.array([0.,0.,1.])
	else:
		printError("Unknown direction "+d)
		raise ValueError
		


	
