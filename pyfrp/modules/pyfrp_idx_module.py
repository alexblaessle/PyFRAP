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

"""Indexing module for PyFRAP toolbox. Mainly contains functions that help finding either

	* image indices
	* mesh indices 
	* extended indices
	
for all types of ROIs, such as

	* circular :py:class:`pyfrp.subclasses.pyfrp_ROI.radialROI`
	* rectangular :py:class:`pyfrp.subclasses.pyfrp_ROI.rectangleROI`
	* quadratic :py:class:`pyfrp.subclasses.pyfrp_ROI.squarelROI`
	* polygonial :py:class:`pyfrp.subclasses.pyfrp_ROI.polyROI`
	
Also provides functions to handle indices in case of quadrant reduction and a powerful suite
of *check* functions that help to figure out if a list of coordinates is inside a ROI., using 
``numpy.where``.
	
"""	

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#numpy/Scipy
import numpy as np

#Plotting
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Misc
import os, os.path
import sys

#Only import scipy functions if environment is not RTD
if os.environ.get('READTHEDOCS', None) != 'True':
	from scipy.interpolate.interpnd import _ndim_coords_from_arrays
	from scipy.spatial import cKDTree
	from scipy.spatial import Delaunay
	
#PyFRAP modules
from . import pyfrp_misc_module as pyfrp_misc
from . import pyfrp_plot_module as pyfrp_plt

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def getCircleIdxImg(center,radius,res,debug=False):
		
	"""Returns all indices of image that lie within given circle.
	
	Args:
		center (numpy.ndarray): Center of circle.
		radius (float): Radius of circle.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
			
			* ind_circ_x (list): List of indices inside circle in x-direction.
			* ind_circ_y (list): List of indices inside circle in y-direction.

	"""
		
	#Empty index vectors
	ind_circ_x=[]
	ind_circ_y=[]
	
	if debug:
		ind_slice_debug=np.zeros((res,res))
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in circle
			if checkInsideCircle(i+1,j+1,center,radius):
				ind_circ_x.append(i)
				ind_circ_y.append(j)
				if debug:
					ind_circ_debug[i,j]=1
	
	if debug:
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Circle"],sup="getCircleIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)
		
	return ind_circ_x,ind_circ_y

def getRectangleIdxImg(offset,sidelengthX,sidelengthY,res,debug=False):
	
	"""Returns all indices of image that lie within given rectangle.
	
	.. note:: Offset is set to be bottom left corner.
	
	Args:
		offset (numpy.ndarray): Offset of rectangle.
		sidelengthX (float): Sidelength in x-directiion.
		sidelengthY (float): Sidelength in y-directiion.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
			
			* indX (list): List of indices inside rectangle in x-direction.
			* indY (list): List of indices inside rectangle in y-direction.

	"""
	
	indX=[]
	indY=[]
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in square
			if checkInsideRectangle(i+1,j+1,offset,sidelengthX,sidelengthY):
				indX.append(i)
				indY.append(j)
				
				if debug:
					indDebug[i,j]=1
					
	if debug:
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Rectangle"],sup="getRectangleIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)			
		
	return indX,indY	


def getSquareIdxImg(offset,sidelength,res,debug=False):
	
	"""Returns all indices of image that lie within given square.
	
	.. note:: Offset is set to be bottom left corner.
	
	Args:
		offset (numpy.ndarray): Offset of square.
		sidelengtX (float): Sidelength.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
			
			* indX (list): List of indices inside square in x-direction.
			* indY (list): List of indices inside square in y-direction.

	"""
	
	
	indX=[]
	indY=[]
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in square
			if checkInsideSquare(i+1,j+1,offset,sidelength):
				indX.append(i)
				indY.append(j)
				
				if debug:
					indDebug[i,j]=1
					
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["Square"],sup="getSquareIdxImg debugging output")
		axes[0].contourf(ind_circ_debug)			
		
	return indX,indY

def getAllIdxImg(res,debug=False):
	
	"""Returns all indices of image.
	
	Args:
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
			
			* indX (list): List of indices inside image in x-direction.
			* indY (list): List of indices inside image in y-direction.

	"""
	
	#Empty index vectors
	indX=[]
	indY=[]
	
	if debug:
		indDebug=zeros((res,res))
	
	#Go through all pixels
	for i in range(int(res)):
		for j in range(int(res)):
			
			#Check if in circle
			indX.append(i)
			indY.append(j)
			if debug:
				indDebug[i,j]=1
				
	if debug:
		
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,1],titles=["slice"],sup="getAllIdxImg debugging output")
		axes[0].contourf(indDebug)
		
	return indX,indY

def getPolyIdxImg(corners,res,debug=False):
	
	"""Returns all indices of image that lie within given polygon.
	
	Args:
		corners (list): List of (x,y)-coordinates of corners.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
			
			* indX (list): List of indices inside polygon in x-direction.
			* indY (list): List of indices inside polygon in y-direction.

	"""
	
	#Convert to np array if necessary
	corners=np.asarray(corners)
	
	#Define polygonial patch
	poly = ptc.Polygon(corners,edgecolor='r',facecolor=(1,0,0,.2),)
	
	#Create grid
	x_int=np.arange(1,res+1,1)
	y_int=np.arange(1,res+1,1)
	g = np.meshgrid(x_int, y_int)
	
	#Zip them into coordinate tuples
	coords = list(zip(*(c.flat for c in g)))
	
	
	
	#Check which point is inside
	pts = np.vstack([p for p in coords if poly.contains_point(p, radius=0)])
	

	
	indX,indY= pts[:,0],pts[:,1]
	
	return indX,indY
		
def getCircleIdxMesh(center,radius,mesh,zmin="-inf",zmax="inf",debug=False):
	
	"""Returns all indices of mesh that lie within given circle and between ``zmin`` and
	``zmax``.
	
	Args:
		center (numpy.ndarray): Center of circle.
		radius (float): Radius of circle.
		mesh (pyfrp.subclasses.pyfrp_mesh.mesh): Mesh.
		
	Keyword Args:
		zmin (float): Minimal z-coordinate.
		zmax (float): Maximal z-coordinate.
		debug (bool): Print debugging messages.
		
	Returns:
		list: List of mesh indices inside circle. 
		
	"""
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.getCellCenters()
	
	#Convert into complex numbers
	c=np.array([np.complex(xc,yc) for xc,yc in zip(x,y)])
	centerC=np.complex(center[0],center[1])
	
	#Get indices in Circle
	indCircle=np.where(np.abs(c-centerC)<radius)[0]

	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	

	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indCircle)
	
	return indFinal


def getSliceIdxMesh(z,zmin,zmax,debug=False):
	
	"""Returns all indices of mesh that lie within given slice between ``zmin`` and
	``zmax``.
	
	Args:
		z (float): z-coordinates of mesh.
		zmin (float): Minimal z-coordinate.
		zmax (float): Maximal z-coordinate.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		list: List of mesh indices inside slice. 
		
	"""
	
	indSlice=np.where((z<zmax) & (z > zmin))[0]
	return indSlice

def getRectangleIdxMesh(sidelengthX,sidelengthY,offset,mesh,zmin="-inf",zmax="inf",debug=False):
	
	"""Returns all indices of mesh that lie within given rectangle and between ``zmin`` and
	``zmax``.
	
	.. note:: Offset is set to be bottom left corner.
	
	Args:
		offset (numpy.ndarray): Offset of rectangle.
		sidelengthX (float): Sidelength in x-directiion.
		sidelengthY (float): Sidelength in y-directiion.
		mesh (pyfrp.subclasses.pyfrp_mesh.mesh): Mesh.
		
	Keyword Args:
		zmin (float): Minimal z-coordinate.
		zmax (float): Maximal z-coordinate.
		debug (bool): Print debugging messages.
		
	Returns:
		list: List of mesh indices inside rectangle. 
		
	"""
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.getCellCenters()
	
	#Getting indices
	indSquare=np.where((offset[0]<x) & (x<offset[0]+sidelengthX) & (offset[1]<y) & (y<offset[1]+sidelengthY))[0]
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indSquare)
	
	return indFinal

def getSquareIdxMesh(sidelength,offset,mesh,zmin="-inf",zmax="inf",debug=False):
	
	"""Returns all indices of mesh that lie within given square and between ``zmin`` and
	``zmax``.
	
	.. note:: Offset is set to be bottom left corner.
	
	Args:
		offset (numpy.ndarray): Offset of square.
		sidelength (float): Sidelength.
		mesh (pyfrp.subclasses.pyfrp_mesh.mesh): Mesh.
		
	Keyword Args:
		zmin (float): Minimal z-coordinate.
		zmax (float): Maximal z-coordinate.
		debug (bool): Print debugging messages.
		
	Returns:
		list: List of mesh indices inside square. 
		
	"""
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.getCellCenters()
	
	#Getting indices
	indSquare=np.where((offset[0]<x) & (x<offset[0]+sidelength) & (offset[1]<y) & (y<offset[1]+sidelength))[0]
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indSquare)
	
	return indFinal

def getPolyIdxMesh(corners,mesh,zmin="-inf",zmax="inf",debug=False):
	
	"""Returns all indices of mesh that lie within given polygon and between ``zmin`` and
	``zmax``.
	
	Args:
		corners (list): List of (x,y)-coordinates of corners.
		mesh (pyfrp.subclasses.pyfrp_mesh.mesh): Mesh.
		
	Keyword Args:
		zmin (float): Minimal z-coordinate.
		zmax (float): Maximal z-coordinate.
		debug (bool): Print debugging messages.
		
	Returns:
		list: List of mesh indices inside polygon. 
		
	"""
	
	#Checking that zmin/zmax are converted into numpy floats
	zmin=pyfrp_misc.translateNPFloat(zmin)
	zmax=pyfrp_misc.translateNPFloat(zmax)
	
	#Grabbing cellCenters of mesh
	x,y,z=mesh.getCellCenters()
	
	coords=list(zip(np.array(x),np.array(y)))
	
	#Bookkeeping list
	#indPoly=[]
	
	
	#coords = list(zip(*(c.flat for c in g)))
	#Define polygonial patch
	poly = ptc.Polygon(corners,edgecolor='r',facecolor=(1,0,0,.2),)
	
	#Check which point is inside
	indPoly = [i for i,p in enumerate(coords) if poly.contains_point(p, radius=0)]
	
	#Loop through coordinates and check if inside
	#for i in range(len(x)):
		#if checkInsidePoly(x[i],y[i],corners):
			#indPoly.append(i)
	
	#Get indices in Slice
	indSlice=getSliceIdxMesh(z,zmin,zmax)
	
	#Get matches indices
	indFinal=pyfrp_misc.matchVals(indSlice,indPoly)
		
	return indFinal	
		
def checkInsideCircle(x,y,center,radius):
	
	"""Checks if coordinate (x,y) is in circle with given radius and center.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	   
	Args:
		x (numpy.ndarray): Array of x-coordinates.
		y (numpy.ndarray): Array of y-coordinates.
		center (numpy.ndarray): Center of circle.
		radius (float): Radius of circle.
		
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	return np.sqrt((x-center[0])**2+(y-center[1])**2)<radius

def checkInsideSquare(x,y,offset,sidelength):
	
	"""Checks if coordinate (x,y) is in square with given offset and sidelength.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	   
	.. note:: Offset is set to be bottom left corner.   
	   
	Args:
		x (numpy.ndarray): Array of x-coordinates.
		y (numpy.ndarray): Array of y-coordinates.
		offset (numpy.ndarray): Offset of square.
		sidelength (float): Sidelength.
		
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	return (x<=offset[0]+sidelength) * (offset[0]<=x) * (y<=offset[1]+sidelength) * (offset[1]<=y)
		
def checkInsideRectangle(x,y,offset,sidelengthX,sidelengthY):
	
	"""Checks if coordinate (x,y) is in rectangle with given offset and sidelength.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	   
	.. note:: Offset is set to be bottom left corner.   
	   
	Args:
		x (numpy.ndarray): Array of x-coordinates.
		y (numpy.ndarray): Array of y-coordinates.
		offset (numpy.ndarray): Offset of rectangle.
		sidelengthX (float): Sidelength in x-direction.
		sidelengthY (float): Sidelength in y-direction.
		
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	return (x<=offset[0]+sidelengthX) * (offset[0]<=x) * (y<=offset[1]+sidelengthY) * (offset[1]<=y)
	
def checkInsideImg(x,y,res,offset=[0,0]):
	
	"""Checks if coordinate (x,y) is inside image.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	      
	Args:
		x (numpy.ndarray): Array of x-coordinates.
		y (numpy.ndarray): Array of y-coordinates.
		res (int): Resolution of image (e.g. 512).
		
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	return (x>offset[0]) * (x<offset[0]+res)*(y>offset[1]) * (y<offset[1]+res)

def checkInsidePolyVec(x,y,poly):	
	
	"""Checks if coordinate (x,y) is inside polyogn, checks first if vector or just value.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	      
	Args:
		poly (list): List of (x,y)-coordinates of corners.
		x (numpy.ndarray): Array of x-coordinates.
		y (numpy.ndarray): Array of y-coordinates.
			
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	try:
		len(x)
	except TypeError:
		return checkInsidePoly(x,y,poly)
	
	vec=[]
	for i in range(len(x)):
		vec.append(checkInsidePoly(x[i],y[i],poly))
	     
	return vec

def checkInsidePolyOld(x,y,poly):
	
	"""Checks if coordinate (x,y) is inside polyogn.
	
	Adapted from http://www.ariel.com.au/a/python-point-int-poly.html.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	      
	Args:
		poly (list): List of (x,y)-coordinates of corners.
		x (float): x-coordinate.
		y (float): y-coordinate.
			
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	n = len(poly)
	inside =False
		
	
	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y
	
	return inside

def checkInsidePoly(x,y,poly):
	
	"""Checks if coordinate (x,y) is inside polyogn.
	
	Adapted from http://www.ariel.com.au/a/python-point-int-poly.html.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	      
	Args:
		poly (list): List of (x,y)-coordinates of corners.
		x (float): x-coordinate.
		y (float): y-coordinate.
			
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	try:
		x[0]
		wasVec=True
	except:
		x=[x]
		y=[y]
		wasVec=False	
		
	coords=list(zip(np.array(x),np.array(y)))
	
	#Define polygonial patch
	poly = ptc.Polygon(poly,edgecolor='r',facecolor=(1,0,0,.2),)
	
	#Check which point is inside
	ret= [poly.contains_point(p, radius=0) for i,p in enumerate(coords)]
	
	if not wasVec:
		ret=ret[0]
		
	return ret
	
def checkQuad(x,y,res):
	
	"""Checks if coordinate (x,y) is inside first quadrant.
	
	.. note:: If ``x`` and ``y`` are ``float``, will return ``bool``, otherwise
	   ``numpy.ndarray`` of booleans.
	      
	Args:
		x (float): x-coordinate.
		y (float): y-coordinate.
		res (int): Resolution of image (e.g. 512).
		
	Returns:
		bool: True if inside, otherwise False.
			
	"""
	
	return (res/2.-1<=x) & (res/2.-1<=y)
	
def checkSquareSize(ind_sq_x,ind_sq_y,sidelength):
	
	"""Checks if square described by indices has right size.
	      
	Args:
		ind_sq_x (list): Image indices in x-direction.
		ind_sq_y (list): Image indices y x-direction.
		sidelength (float): Sidelength of square.
		
	Returns:
		bool: True if equal size, otherwise False.
			
	"""

	if floor(sidelength)**2==len(ind_sq_x) and floor(sidelength)**2==len(ind_sq_x) and floor(sidelength)**2==len(unique(ind_sq_x))*len(unique(ind_sq_y)):
		return True
	else:
		return False
		
	
def checkSquareCenteredFromInd(ind_sq_x,ind_sq_y,res):
	
	"""Checks if square described by indices is cenntered in image.
	      
	Args:
		ind_sq_x (list): Image indices in x-direction.
		ind_sq_y (list): Image indices y x-direction.
		res (int): Resolution of image (e.g. 512).
		
	Returns:
		bool: True if centered, otherwise False.
			
	"""
	
	#first bracket should be zero and second too
	return not bool((min(ind_sq_x)-(res-1-max(ind_sq_x)))+(min(ind_sq_y)-(res-1-max(ind_sq_y))))

def checkSquCentered(offset,sidelength,res):
	
	"""Checks if square is centered in image.
	
	.. note:: Need a correction by .5 because there is a difference between pixels and coordinates.
	      
	Args:
		offset (numpy.ndarray): Offset of square.
		sidelength (float): Sidelength.
		res (int): Resolution of image (e.g. 512).
		
	Returns:
		bool: True if centered, otherwise False.
			
	"""
	
	return not bool(((res-sidelength)-2*(offset[0]-0.5)) and ((res-sidelength)-2*(offset[1]-0.5)))


def idx2QuadImg(indX,indY,res,debug=False):
	
	"""Reduces indices found for whole domain to first quadrant.
	
	.. note:: Need a correction by .5 because there is a difference between pixels and coordinates.
	      
	Args:
		indX (list): List of indices in x-direction.
		indY (list): List of indices in y-direction.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* indXQuad (list): List of reduced indices.
			* indYQuad (list): List of reduced indices.
			
	"""
	
	#Convert to np array
	indX=np.asarray(indX)
	indY=np.asarray(indY)
	
	#Find indices of first quadrant
	inds=where((res/2.-1<=indX) & (res/2.-1<=indY))[0]
	
	#Assign new indices
	indXQuad=indX[inds]
	indYQuad=indY[inds]
	
	#Debugging plot if necessary
	if debug:
	
		#Create figure
		fig,axes = pyfrp_plt.make_subplot([1,3],titles=["Original indices","Flipped indices","Original and flipped indices"],sup="ind2quad debugging output")
		
		#Plot original and flipped indices
		axes[0].plot(indX,indY,'y*')
		axes[1].plot(indXQuad,indYQuad,'r*')
		axes[2].plot(indX,indY,'y*')
		axes[2].plot(indXQuad,indYQuad,'r*')
		
		#Draw
		plt.draw()
	
	return indXQuad, indYQuad

def regions2quad(inds,res,debug=False):
	
	"""Reduces indices of regions specified in ``inds`` found for whole domain to first quadrant.
	      
	Args:
		inds (list): List of indices-list doubles.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		list: List of reduced indices-list doubles.
			
	"""
	
	
	inds_quad=[]
	for ind in inds:
		ind_quad=ind2quad(ind[0],ind[1],res,debug=debug)
		inds_quad.append(ind_quad)	
	
	return inds_quad

def ind2mask(vals,ind_x,ind_y,val):
	
	"""Converts indices lists into mask.
	      
	Args:
		vals (numpy.ndarray): Array on which to be masked.
		ind_x (list): Indices in x-direction.
		ind_y (list): Indices in y-direction.
		val (float): Value that is assigned to pixels in indices-lists.
	
	Returns:
		numpy.ndarray: Masked array.
			
	"""
	
	vals[ind_x,ind_y]=val
	
	return vals

def mask2ind(mask,res):
	
	"""Converts mask into indices list.
	      
	Args:
		mask (numpy.ndarray): Mask array.
		res (int): Resolution of image (e.g. 512).
	
	Returns:
		tuple: Tuple containing:
		
			* indX_new (list): List of x-indices.
			* indY_new (list): List of y-indices.
			
	"""
	
	#To bool
	mask=mask.astype(bool)

	#idx grid
	x=np.arange(res)
	y=np.arange(res)	
	X,Y=np.meshgrid(x,y)

	#Slice idx grid
	idxX_new=X[mask].flatten().astype(int)
	idxY_new=Y[mask].flatten().astype(int)
	
	return idxX_new, idxY_new

def getExtendedPixelsSquare(offset,sidelength,res,debug=False):
	
	"""Finds theoretical pixels that could be filled up with rim
	concentration for a square.
	      
	Args:
		offset (numpy.ndarray): Offset of square.
		sidelength (float): Sidelength.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* indX (list): List of x-indices.
			* indY (list): List of y-indices.
			
	"""
	
	indX=[]
	indY=[]
	
	x=np.arange(np.floor(offset[0]),np.ceil(offset[0]+sidelength))
	y=np.arange(np.floor(offset[1]),np.ceil(offset[1]+sidelength))
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res):
				indX.append(i)
				indY.append(j)
				
	return indX,indY



def getExtendedPixelsRectangle(offset,sidelengthX,sidelengthY,res,debug=False):
	
	"""Finds theoretical pixels that could be filled up with rim
	concentration for a rectangle.
	      
	Args:
		offset (numpy.ndarray): Offset of rectangle.
		sidelengthX (float): Sidelength in x-direction.
		sidelengthY (float): Sidelength in y-direction.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* indX (list): List of x-indices.
			* indY (list): List of y-indices.
			
	"""
	
	x=np.arange(np.floor(offset[1]),np.ceil(offset[0]+sidelengthX))
	y=np.arange(np.floor(offset[1]),np.ceil(offset[1]+sidelengthY))
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkInsideRectangle(i,j,offset,sidelengthX,sidelengthY):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

def getExtendedPixelsCircle(center,radius,res,debug=False):
	
	"""Finds theoretical pixels that could be filled up with rim
	concentration for a circle.
	      
	Args:
		center (numpy.ndarray): Center of circle.
		radius (float): Radius of circle.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* indX (list): List of x-indices.
			* indY (list): List of y-indices.
			
	"""
	
	x=np.arange(np.ceil(center[0]-radius),np.floor(center[0]+radius))
	y=np.arange(np.ceil(center[1]-radius),np.floor(center[1]+radius))
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkInsideCircle(i,j,center,radius):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

def getExtendedPixelsPolygon(corners,res,debug=False):
	
	"""Finds theoretical pixels that could be filled up with rim
	concentration for a polygon.
	      
	Args:
		corners (list): List of (x,y)-coordinates of corners.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* indX (list): List of x-indices.
			* indY (list): List of y-indices.
			
	"""
	
	cornersNP=np.array(corners)
	
	xmax=cornersNP[:,0].max()
	xmin=cornersNP[:,0].min()
	ymax=cornersNP[:,1].max()
	ymin=cornersNP[:,1].min()
	
	x=np.arange(np.floor(xmin),np.ceil(xmax))
	y=np.arange(np.floor(ymin),np.ceil(ymax))
	
	
	indX=[]
	indY=[]
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res) and checkPolygon(i,j,corners):
				indX.append(i)
				indY.append(j)
				
	return indX,indY

def getCommonExtendedPixels(ROIs,res,debug=False,procedures=None):
	
	"""Finds theoretical pixels that could be filled up with rim concentration 
	for a list of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` ROIs. 
	
	The ``procedures`` input is only necessary if ROI is a :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI`,
	combining multiple ROIs via addition/substraction.
	
	Args:
		ROIs (list): List of ROIs.
		res (int): Resolution of image (e.g. 512).
	
	Keyword Args:
		debug (bool): Print debugging messages.
		procedures (list): List of addition/substraction procedures.
		
	Returns:
		tuple: Tuple containing:
		
			* indX (list): List of x-indices.
			* indY (list): List of y-indices.
			
	"""
	
	xExtend,yExtend=getCommonXYExtend(ROIs,debug=debug)
	
	x=np.arange(np.floor(xExtend[0]),np.ceil(xExtend[1]))
	y=np.arange(np.floor(yExtend[0]),np.ceil(yExtend[1]))
	
	indX=[]
	indY=[]
	
	if procedures==None:
		procedures=np.ones(np.shape(ROIs))
	
	for i in x:
		for j in y:
			if not checkInsideImg(i,j,res):
				b=True
				for k,r in enumerate(ROIs):
					if 1+procedures[k]:
						b=b and r.checkXYInside(i,j)
					else:					
						b=b and not r.checkXYInside(i,j)
						
				if b:
					indX.append(i)
					indY.append(j)
	
	return indX,indY

def getCommonXYExtend(ROIs,debug=False):
	
	"""Finds common x-y-extend of a list of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` ROIs..
	
	Args:
		ROIs (list): List of ROIs.
		
	Returns:
		tuple: Tuple containing:
		
			* xExtend (list): ``[minx,maxx]``.
			* yExtend (list): ``[miny,maxy]``.
			
	"""
	
	xExtends,yExtends=[],[]
	for r in ROIs:
		xExtend,yExtend=r.computeXYExtend()
		xExtends.append(xExtend)
		yExtends.append(yExtend)
	
	xExtends=np.array(xExtends)
	yExtends=np.array(yExtends)
	
	xExtend=[xExtends[:,0].min(),xExtends[:,1].max()]
	yExtend=[yExtends[:,0].min(),yExtends[:,1].max()]
	
	if debug:
		print("xExtend = ", xExtend)
		print("yExtend = ", yExtend)
	
	return xExtend, yExtend

def remRepeatedImgIdxs(idxX,idxY,debug=False):
	
	"""Remove repeated indices tupels from index lists for images.
	
	Args:
		idxX (list): List of x-indices.
		idxY (list): List of y-indices.
	
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
		
			* idxX (list): Filtred list of x-indices.
			* idxY (list): Filtered list of y-indices.
	
	"""
	
	idx=list(zip(idxX,idxY))
	idx=pyfrp_misc.remRepeatsList(idx)
	idxX,idxY=pyfrp_misc.unzipLists(idx)
	return idxX,idxY

def maskMeshByDistance(x,y,d,grid):
	
	"""Filters all (x,y) coordinates that are more than d 
	in meshgrid given some actual coordinates (x,y).
	
	Args:
		x (numpy.ndarray): x-coordinates.
		y (numpy.ndarray): y-coordinates.
		d (float): Maximum distance.
		grid (numpy.ndarray): Numpy meshgrid.
	
	Returns:
		idxs (list): List of booleans.
			
	"""
	
	#Convert x/y into useful array
	xy=np.vstack((x,y))
	
	#Compute distances to nearest neighbors
	tree = cKDTree(xy.T)
	xi = _ndim_coords_from_arrays(tuple(grid))
	dists, indexes = tree.query(xi)
	
	#Get indices of distances that are further apart than d
	idxs = (dists > d)
	
	return idxs,

def nearestNeighbour3D(xi,yi,zi,x,y,z,k=1,minD=None):
	
	"""Finds k nearest neighbour to points.
	
	Uses http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.query.html#scipy.spatial.cKDTree.query .
	
	Example:
	
	>>> from pyfrp.modules import pyfrp_idx_module
	>>> import numpy as np
	>>> import matplotlib.pyplot as plt

	>>> N=50

	>>> x=np.random.random(N)
	>>> y=np.random.random(N)
	>>> z=np.random.random(N)

	>>> xi=[0.,1.]
	>>> yi=[0.,1.]
	>>> zi=[0.,1.]

	>>> idx,dist=pyfrp_idx_module.nearestNeighbour3D(xi,yi,zi,x,y,z,k=1)

	>>> fig=plt.figure()
	>>> ax=fig.add_subplot(111,projection='3d')
	>>> ax.scatter(x,y,z,color='k')

	>>> ax.scatter(x[idx[0]],y[idx[0]],z[idx[0]],color='b',s=50)
	>>> ax.scatter(x[idx[1]],y[idx[1]],z[idx[1]],color='r',s=50)
	>>> ax.scatter(xi,yi,zi,color='g',s=50)
	
	>>> plt.show()
	
	.. image:: ../imgs/pyfrp_idx_module/neighbours3D.png
	
	
	.. note:: To avoid that the a point is nearest neighbour to itself, one can choose 
	   minD=0 to avoid that points with distances ``0`` are returned. Note, if this could
	   be the possibility, one might have at least ``k=2`` to avoid an empty return.
	
	Args:
		xi (numpy.ndarray): x-coordinates of points to find neighbours to.
		yi (numpy.ndarray): x-coordinates of points to find neighbours to.
		zi (numpy.ndarray): x-coordinates of points to find neighbours to.
		x (numpy.ndarray): x-coordinates of possible neighbours.
		y (numpy.ndarray): x-coordinates of possible neighbours.
		z (numpy.ndarray): x-coordinates of possible neighbours.
	
	Keyword Args:
		k (int): Number of neighbours to find per point.
		minD (float): Minimum distance nearest neighbour must be away.
		
	Returns:
		tuple: Tuple containing:
		
			* indexes (list): Indices of k -closest neighbours per point.
			* dists (list): Distances of k -closest neighbours per point.

	"""
	
	#Convert x/y into useful array
	xyz=np.vstack((x,y,z))
	xyzi=np.vstack((xi,yi,zi))
	
	#Compute distances to nearest neighbors
	tree = cKDTree(xyz.T)
	dists, indexes = tree.query(xyzi.T,k=k)
	
	if minD!=None:
		indexes=indexes[dists>minD]
		dists=dists[dists>minD]
		
	return indexes, dists
	
def triangulatePoly(coords,addPoints=False,iterations=2,debug=False):
	
	"""Triangulates a polygon with given coords using a Delaunay triangulation.
	
	Uses http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html#scipy.spatial.Delaunay to 
	calculate triangulation, then filters the triangles actually lying within the polygon.
	
	Args:
		coords (list): List of (x,y)-coordinates of corners.
	
	Keyword Args:
		addPoints (bool): Allow incremental addition of points.
		iterations (int): Number of iterations of additional point adding.
		debug (boo): Print debugging messages.
		
	Returns:
		tuple: Tuple containing:
		
			* triFinal (list): List of found triangles.
			* coordsTri (list): List of vertex coordinates.
	"""
	
	#Bookkeeping list
	triFinal=[]
		
	#Triangulate
	tri=Delaunay(coords,incremental=addPoints)
	
	#Backup original coordinates
	coordsOrg=list(coords)
	
	if debug:
		print("Found ", len(tri.simplices.copy()), "triangles in initial call.")

	#Incrementally refine triangulation
	if addPoints:
		for i in range(iterations):
			
			mids=[]
			for j in range(len(tri.simplices)):
				mid=getCenterOfMass(coords[tri.simplices.copy()[j]])
				mids.append(mid)
				
			coords=np.asarray(list(coords)+mids)
			tri.add_points(mids,restart=True)
				
		if debug:
			print("Found ", len(tri.simplices.copy()), "triangles after iterations.")
		
			
	#Remember assigment of points by traingulation function
	coordsTri=tri.points
	
	midsIn=[]
	midsOut=[]
	triOutIdx=[]
	triInIdx=[]
	
	for i in range(len(tri.simplices)):
	
		# Get COM of triangle
		mid=getCenterOfMass(coordsTri[tri.simplices.copy()[i]])
		
		#Check if triangle is inside original polygon
		if checkInsidePolyVec(mid[0],mid[1],coordsOrg):
			triFinal.append(tri.simplices.copy()[i])
			midsIn.append(mid)
			triInIdx.append(i)		
		else:
			triOutIdx.append(i)
			midsOut.append(mid)
	
	if debug:
		print("Removed ", len(tri.simplices.copy())-len(triFinal), "triangles through COM criteria.")	
		print("Returning ", len(triFinal), "triangles.")	
		
		
	return triFinal,coordsTri
	
def getCenterOfMass(xs,axis=0,masses=None):
	
	r"""Computes center of mass of a given set of points.
	
	.. note:: If ``masses==None``, then all points are assigned :math:`m_i=1`.
	
	Center of mass is computed by:
	
	.. math:: C=\frac{1}{M}\sum\limits_i m_i (x_i)^T
	
	where 
	
	.. math:: M = \sum\limits_i m_i
	
	Args:
		xs (numpy.ndarray): Coordinates.
		
	Keyword Args:
		masses (numpy.ndarray): List of masses.
	
	Returns:
		numpy.ndarray: Center of mass.
	
	"""
	
	if masses==None:
		masses=np.ones(xs.shape[0])
	
	return (masses*xs.T).T.sum(axis=axis)/float(xs.shape[axis])

			



				
				
	
	
	
