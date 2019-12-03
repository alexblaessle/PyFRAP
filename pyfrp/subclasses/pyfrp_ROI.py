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

#ROI module for PyFRAP toolbox, including following ROI objects:

#(1)  ROI
#(2)  radialROI
#(3)  sliceROI
#(4)  radialSliceROI
#(5)  squareROI
#(6)  squareSliceROI
#(7)  rectangleROI
#(8)  rectangleSliceROI
#(9)  polyROI
#(10)  polySliceROI
#(11)  customROI

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module 
from pyfrp.modules import pyfrp_idx_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_img_module
from pyfrp.modules import pyfrp_integration_module
from pyfrp.modules import pyfrp_fit_module
from pyfrp.modules import pyfrp_gmsh_geometry
from pyfrp.modules import pyfrp_gmsh_module
from pyfrp.modules import pyfrp_openscad_module

from pyfrp.modules.pyfrp_term_module import *

#Plotting
import matplotlib.pyplot as plt
import matplotlib.patches as ptc

#Time 
import time

#Copy
import copy

#OS
import os
import shutil

#Solid/Opescad
import solid
import solid.utils

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main ROI class

class ROI(object):

	def __init__(self,embryo,name,Id,zmin='-inf',zmax='inf',color='b'):
		
		#Name/Id
		self.name=name
		self.Id=Id
		self.embryo=embryo
		self.color=color
		
		#zExtend
		self.zmin=pyfrp_misc_module.translateNPFloat(zmin)
		self.zmax=pyfrp_misc_module.translateNPFloat(zmax)
		
		#Idxs from data analysis/simulation
		self.imgIdxX=[]
		self.imgIdxY=[]
		self.extImgIdxX=[]
		self.extImgIdxY=[]
		self.meshIdx=[]
		
		#Mask
		self.imgMask=None
		self.extMask=None
		
		#Number of extended pixels
		self.numExt=0 # None
		
		#Result Dataseries
		self.dataVec=[]
		self.simVec=[]
		self.dataVecPinned=[]
		self.simVecPinned=[]
		
		#Rim concentration
		self.useForRim=False
	
	def setId(self,Id):
		
		"""Sets Id of ROI.
		
		Args:
			Id (int): New Id.
		
		Returns:
			int: New Id.
		
		"""
		
		self.Id=Id
		return self.Id
	
	def setName(self,n):
		
		"""Sets name of ROI.
		
		Args:
			n (str): New name.
			
		Returns:
			str: New name.
		
		"""
		
		self.name=n
		return self.name
	
	def setZExtend(self,zmin,zmax):
		
		"""Sets extend in z-direction.
		
		Args:
			zmin (float): Minimum z-coordinate.
			zmax (float): Maximum z-coordinate.
		
		Returns:
			list: New z-extend given by ``[zmin,zmax]``.
		
		"""
		
		self.zmin=zmin
		self.zmax=zmax
		
		return [self.zmin,self.zmax]
	
	def getZExtend(self):
		
		"""Returns extend in z-direction.
	
		Returns:
			list: Z-extend given by ``[zmin,zmax]``.
		
		"""
		
		return [self.zmin,self.zmax]
	
	def getRealZExend(self):
		
		r"""Returns real extend in z-direction.
	
		Real extend returns
		
		.. math:: z_{\mathrm{min}}=\mathrm{max} (z_{\mathrm{min,ROI}},z_{\mathrm{min,geometry}})
		
		and 
		
		.. math:: z_{\mathrm{max}}=\mathrm{min} (z_{\mathrm{max,ROI}},z_{\mathrm{max,geometry}})
		
		Returns:
			list: Z-extend given by ``[zmin,zmax]``.
		
		"""
		
		zminGeo,zmaxGeo=self.embryo.geometry.getZExtend()
	
		zmin=max([zminGeo,self.zmin])
		zmax=min([zmaxGeo,self.zmax])
		
		return [zmin,zmax]
	
	def getOpenscadZExtend(self):
		
		"""Returns extend in z-direction suitable for rendering 
		the ROI via openscad.
		
		If either ``zmin`` or ``zmax`` is infinity, then uses 
		:py:func:getRealZExend to return more meaningful extend.
		
		Returns:
			list: Z-extend given by ``[zmin,zmax]``.
		
		"""
		
		if np.inf==abs(self.zmin):
			zminReal,zmaxReal=self.getRealZExend()
			zmin=zminReal
		else:
			zmin=self.zmin
		if np.inf==abs(self.zmax):
			zminReal,zmaxReal=self.getRealZExend()
			zmax=zmaxReal
		else:
			zmax=self.zmax
		
		return zmin,zmax
	
	def getId(self):
		
		"""Returns Id of ROI.
		
		Returns:
			int: Id.
			
		"""	
		
		return self.Id
	
	def getName(self):
		
		"""Returns name of ROI.
		
		Returns:
			str: Current name.
		"""	
		
		return self.name
	
	def getImgIdx(self):
		
		"""Returns image indices of ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		return self.imgIdxX,self.imgIdxY
	
	def getExtImgIdx(self):
		
		"""Returns extended image indices of ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* extImgIdxX (list): Extended image indices in x-direction.
				* extImgIdxY (list): Extended image indices in y-direction.
				
		"""
		
		return self.extImgIdxX,self.extImgIdxY
	
	def getMeshIdx(self):
		
		"""Returns mesh indices of ROI.
		
		Returns:
			list: Mesh indices of ROI.
				
		"""
		
		return self.meshIdx
	
	def getMeshDensity(self):
		
		r"""Returns average mesh density inside ROI.
		
		Mesh density is defined by 
		
		.. math:: \rho=N/V,
		
		where :math:`N` is the number of mesh nodes inside ROI 
		and :math:`V` is the volume of ROI, see also 
		:py:func:`getVolume`.
		
		Returns:
			float: Mesh density.
				
		"""
		
		volume=self.getVolume()
		return len(self.meshIdx)/volume
		
	def getVolume(self):
		
		r"""Returns volume of ROI.
		
		Since ROIs only behave linearly in z-direction, volume
		is given by
		
		.. math:: V = A * h,
		
		where :math:`h` is ROI height (see :py:func:`getROIHeight`) and 
		:math:`A` is ROI area (see :py:func:`getArea`).
		
		Returns:
			float: ROI volume.
		
		"""
		
		area=self.getArea()
		height=self.getROIHeight()
		
		return area*height
		
	def getROIHeight(self):
		
		"""Returns height of ROI.
		
		Returns:
			float: Height of ROI.
		
		"""
	
		if np.isfinite(self.zmax):
			zMax=self.zmax
		else:
			dump,zMax=self.getMeshIdxZExtend()
		
		if np.isfinite(self.zmin):
			zMin=self.zmin
		else:
			zMin,dump=self.getMeshIdxZExtend()
			
		return abs(zMax-zMin)
	
	def getArea(self):
		
		"""Returns area of ROI.
		
		Area is computed as area covered by
		``imgMask + extMask`` 
		
		Returns:
			float: Area of ROI.
		
		"""
		
		if self.imgMask==None:
			self.computeImgMask()
		if self.extMask==None:
			self.computeExtMask()
			if self.extMask==None:
				return self.imgMask.sum()
		
		return self.imgMask.sum()+self.extMask.sum()
		
	def getMeshIdxZExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in z-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum z-coordinate.
				* (float): Maximum z-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh.mesh
		z=np.asarray(mesh.z)[self.meshIdx]	
		
		return min(z) , max(z)
	
	def getMeshIdxYExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in y-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum y-coordinate.
				* (float): Maximum y-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh
		y=np.asarray(mesh.getCellCenters[1])[self.meshIdx]	
		
		return min(y) , max(y)
	
	def getMeshIdxXExtend(self):
		
		"""Returns extend of ROI's ``meshIdx`` in x-coordinate.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum x-coordinate.
				* (float): Maximum x-coordinate.
				
		"""
		
		mesh=self.embryo.simulation.mesh
		x=np.asarray(mesh.getCellCenters[0])[self.meshIdx]	
		
		return min(x) , max(x)
	
	def getMeshIdxExtend(self):
		
		"""Returns extend of ROI's ``meshIdx``.
		
		Returns:
			tuple: Tuple containing:
				
				* (float): Minimum x-coordinate.
				* (float): Maximum x-coordinate.
				* (float): Minimum y-coordinate.
				* (float): Maximum y-coordinate.
				* (float): Minimum z-coordinate.
				* (float): Maximum z-coordinate.
		"""
		
		xmin,xmax=self.getMeshIdxXExtend()
		ymin,ymax=self.getMeshIdxYExtend()
		zmin,zmax=self.getMeshIdxZExtend()
		
		return xmin,xmax,ymin,ymax,zmin,zmax
		
	def getType(self):
		
		"""Returns type of ROI, splitting off all module names etc. .
		
		Returns:
			str: Type of ROI.
		
		"""
		
		typ=str(type(self))
		before,typ,after=typ.split("'")
		typ=typ.replace('pyfrp_ROI.','')
		typ=typ.replace('ROI','')
		
		typ=typ.replace('pyfrp.subclasses.','')
		
		return typ
	
	def setColor(self,color):
		
		"""Sets color of ROI.
		
		Color can be either ``str``, ``float`` or ``tuple``. See also: 
		http://matplotlib.org/api/colors_api.html
		
		Args:
			color (str): New color.
		
		Returns:
			str: New color.
		
		"""
		
		self.color=color
		return self.color
	
	def setUseForRim(self,b):
		
		"""Marks the ROI to be used for rim calculation.
		
		Args:
			b (bool): True if ROI should be used, False else.
		
		Returns:
			bool: Current flag value.
		
		"""
		
		if self.numExt>0 and b==True:
			printWarning('Be careful, region '+self.name+' is set for rim calculation but has indices outside of image.')
		self.useForRim=b
		return self.useForRim
	
	def getUseForRim(self):
		
		"""Returns if the ROI is used for rim calculation.
		
		Returns:
			bool: Current flag value.
		
		"""
		
		return self.useForRim
	
	def getColor(self):
		
		"""Returns color of ROI.
		"""
		
		return self.color
	
	def emptyIdxs(self):
		
		"""Flushes all indices, inserting empty lists for 
		all of them.
		
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				* meshIdx (list): Mesh indices.
		
			
		"""
		
		self.imgIdxX=[]
		self.imgIdxY=[]
		self.meshIdx=[]
		
		return self.getAllIdxs()
	
	def copyIdxs(self,r):
		
		"""Copies indices of other ROI and inserts them into ROI.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to take indices from.
		
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				* meshIdx (list): Mesh indices.
			
		"""
		
		
		self.imgIdxX=r.imgIdxX
		self.imgIdxY=r.imgIdxY
		
		self.extImgIdxX=r.extImgIdxX
		self.extImgIdxY=r.extImgIdxY
		
		self.meshIdx=r.meshIdx
		
		return self.getAllIdxs()
	
	def getAllIdxs(self):
		
		"""Returns all index arrays of ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				* meshIdx (list): Mesh indices.
		
		"""
		
		return self.imgIdxX,self.imgIdxY,self.meshIdx
	
	def getImgMask(self):
		
		"""Returns image mask of ROI.
		
		Returns:
			numpy.ndarray: Image mask.
		
		"""
		
		return self.imgMask
	
	def getExtMask(self):
		
		"""Returns extended mask of ROI.
		
		Returns:
			numpy.ndarray: Extended mask.
		
		"""
		
		return self.extMask
	
	def computeNumExt(self):
		
		"""Computes number of extended pixels of ROI.
		
		Returns:
			int: Number of extended pixels.
			
		"""
		
		self.numExt=len(self.extImgIdxX)
		return self.numExt
	
	def getNumExt(self):
		
		"""Returns number of extended pixels of ROI.
		
		Returns:
			int: Number of extended pixels.
			
		"""
		
		return self.numExt
	
	def setDataVec(self,vec):
		
		"""Sets data vector of ROI.
		
		Args:
			vec (numpy.ndarray): New data vector.
			
		Returns:
			numpy.ndarray: New data vector.
		
		"""
		
		self.dataVec=vec 
		return self.dataVec

	def getDataVec(self):
		
		"""Returns current data vector of ROI.
			
		Returns:
			numpy.ndarray: Current data vector.
		
		"""
		
		return self.dataVec
	
	def setSimVec(self,vec):
		
		"""Sets simulation vector of ROI.
		
		Args:
			vec (numpy.ndarray): New simulation vector.
			
		Returns:
			numpy.ndarray: New simulation vector.
		
		"""
		
		self.simVec=vec 
		return self.simVec

	def getSimVec(self):
		
		"""Returns current simulation vector of ROI.
			
		Returns:
			numpy.ndarray: Current simulation vector.
		
		"""
		
		return self.simVec
	
	def computeImgMask(self):
		
		"""Computes image mask of ROI.
		
		Image mask is a ``dataResPx * dataResPx`` array with the value
		``1`` for pixels inside ROI and ``0`` elsewhere.
		
		Returns:
			numpy.ndarray: Image mask.
		
		"""
		
		vals=np.zeros((self.embryo.dataResPx,self.embryo.dataResPx))
		self.imgMask=pyfrp_idx_module.ind2mask(vals,self.imgIdxX,self.imgIdxY,1)
		return self.imgMask
		
	def computeExtMask(self):
		
		"""Computes mask of extended pixels of ROI.
		
		Mask is a 2D array with the value
		``1`` for pixels inside ROI and ``0`` elsewhere.
		
		.. note:: Returns ``None,None,None`` if there are no extended pixels.
		
		Also returns coordinate arrays, since offset of extended mask is not ``[0,0]``.
		See also http://docs.scipy.org/doc/numpy/reference/generated/numpy.meshgrid.html .
		
		Returns:
			tuple: Tuple containing:
				
				* mX (numpy.ndarray): Coordinate array corresponding to pixels of extended mask.
				* mY (numpy.ndarray): Coordinate array corresponding to pixels of extended mask.
				* extMask (numpy.ndarray): Extended mask.
		
		"""
		
		if len(self.extImgIdxX)==0:
			return None,None,None
		
		minX=min(self.extImgIdxX)
		maxX=max(self.extImgIdxX)
		
		minY=min(self.extImgIdxY)
		maxY=max(self.extImgIdxY)
		
		X=np.arange(minX,maxX+1)
		Y=np.arange(minY,maxY+1)
		
		mX,mY=np.meshgrid(X,Y)
		
		vals=np.zeros((len(X),len(Y)))
		
		idXtemp=np.asarray(self.extImgIdxX)+abs(minX)
		idYtemp=np.asarray(self.extImgIdxY)+abs(minY)
		
		idXtemp=idXtemp.astype('int')
		idYtemp=idYtemp.astype('int')
		
		self.extMask=pyfrp_idx_module.ind2mask(vals,idXtemp,idYtemp,1)
		
		return mX,mY,self.extMask
		
	def showImgIdx(self,ax=None):
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["imgIdx"],sup=self.name+" imgIdx")
			ax=axes[0]
		
		self.computeImgMask()
	
		ax.imshow(self.imgMask)
		plt.draw()
		
		
		return ax
	
	def showExtImgIdx(self,ax=None):
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["extImgIdx"],sup=self.name+" imgIdx")
			ax=axes[0]
		
		mX,mY,self.extMask=self.computeExtMask()
		if mX!=None:
			ax.contourf(mX,mY,self.extMask)
			plt.draw()
			
		return ax
	
	def showMeshIdx(self,ax=None):
		
		x,y,z=self.embryo.simulation.mesh.getCellCenters()
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["MeshIdx"],sup=self.name+" MeshIdx",proj=['3d'])
			ax=axes[0]
		
		#Somehow need to convert to np array since slicing does not work for fipy variables
		x=np.asarray(x)[self.meshIdx]
		y=np.asarray(y)[self.meshIdx]
		z=np.asarray(z)[self.meshIdx]
		
		ax.scatter(x,y,z,c=self.color)
		plt.draw()
		
		return ax
	
	def showMeshIdx2D(self,ax=None):
		
		
		
		x,y,z=self.embryo.simulation.mesh.getCellCenters()
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["MeshIdx"],sup=self.name+" MeshIdx")
			ax=axes[0]
		
		#Somehow need to convert to np array since slicing does not work for fipy variables
		x=np.asarray(x)[self.meshIdx]
		y=np.asarray(y)[self.meshIdx]
	
		ax.scatter(x,y,c=self.color)
		plt.draw()
		
		return ax
	
	
	def computeIdxs(self,matchMesh=False,debug=False):
		
		"""Computes image and mesh indices of ROI. 
	
		Will do this by:
		
			* Compute image indices.
			* Match image indices with master ROI.
			* Compute external indices.
			* Compute mesh indices.
			* Match mesh indices with the ones of master ROI.
			
		.. note:: If no master ROI is defined, will not do anything.
		
		.. note:: If master ROI has not been indexed yet, will first index it, then continue. 
		
		.. note:: Will skip mesh index computation if there is no mesh generated yet.
		
		Keyword Args:
			matchMesh (bool): Match mesh indices with master ROI.
			debug (bool): Print out debugging messages.
			
		Return:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x direction.
				* imgIdxY (list): Image indices in y direction.
				* meshIdx (list): Mesh indices.
		
		"""
		
		if self.embryo.getMasterROIIdx()==None:
			printWarning("No Master ROI has been defined yet. Will not continue compute ROI indices.")
			return self.getAllIdxs()
		else:
			masterROI=self.embryo.getMasterROI()
		
		if self!=masterROI:
			if len(masterROI.imgIdxX)==0:
				printWarning("Idxs of Master ROI have not been computed. Will compute them first.")
				masterROI.computeIdxs(debug=debug)
		
		startIdx=time.clock()
		
		if  type(self) is not customROI:
			
			self.computeImgIdx(debug=debug)
			self.matchImgIdx(masterROI)
			self.computeExtIdx(debug=debug)
			
			if self.embryo.simulation!=None:
				if self.embryo.simulation.mesh.mesh==None:
					printWarning("Mesh has not been generated, will not compute meshIdxs")
				else:
				
					self.computeMeshIdx(self.embryo.simulation.mesh)
				
					if matchMesh:
						if self!=masterROI:
							self.matchMeshIdx(masterROI)
			else:
				printWarning("Simulation object does not exist yet, hence won't index for mesh.")
		else:
			
			self.updateIdxs()
			self.matchImgIdx(masterROI)
			if self.embryo.simulation!=None:
				if self.embryo.simulation.mesh.mesh!=None:
					self.matchMeshIdx(masterROI)
		
		if debug:
			print('Compute Idxs: ', startIdx-time.clock())
		
		return self.getAllIdxs()
	
	def computeExtIdx(self,debug=False):
		
		"""Computes indices of external pixels.
		
		Does this by comparing extended pixels of ``self`` with the one of the master ROI.
		
		Keyword Args:
			debug (bool): Print out debugging messages.
		
		Return:
			tuple: Tuple containing:
			
				* extImgIdxX (list): External image indices in x direction.
				* extImgIdxY (list): External image indices in y direction.
		
		"""
		
		m=self.embryo.getMasterROI()
		rois=[self,m]
		[self.extImgIdxX,self.extImgIdxY]=pyfrp_idx_module.getCommonExtendedPixels(rois,self.embryo.dataResPx,debug=debug)
		self.computeNumExt()
		return self.extImgIdxX,self.extImgIdxY
	
	def matchImgIdx(self,r):
		
		"""Matches image indices of ``self`` with the ones of ROI ``r``.
		
		Does this by generating masks of both ROIs and multiplicating them.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to match with.
		
		Return:
			tuple: Tuple containing:
			
				* imgIdxX (list): Matched image indices in x direction.
				* imgIdxY (list): Matched image indices in y direction.
		
		"""
		
		self.computeImgMask()
		self.imgMask=self.imgMask*r.computeImgMask()
		self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)	
		return self.imgIdxX,self.imgIdxY
	
	def matchMeshIdx(self,r,matchZ=False):
		
		"""Matches mesh indices of ROI with the ones of a different ROI.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): ROI to match with.
			
		Keyword Args:
			matchZ (bool): Also match with respect to z.
			
		Returns:
			numpy.ndarray: Matched list of mesh indices.
		
		"""
		
		x,y,z=self.embryo.simulation.mesh.getCellCenters()
		
		x=np.asarray(x)[self.meshIdx]
		y=np.asarray(y)[self.meshIdx]
		
		ins=r.checkXYInside(x,y)
		
		self.meshIdx=np.asarray(self.meshIdx)
		self.meshIdx=self.meshIdx[np.where(ins)[0]]
		
		return self.meshIdx
	
	def checkZInside(self,z):
		
		"""Checks if z coordinate is within ROIs z-range.
		
		Arg:
			z (float): z-coordinate.
			
		Returns:
			bool: True if inside.
		"""
		
		if self.zmin<= z and z<=self.zmax:
			return True
		else:
			return False
	
	def showIdxs(self,axes=None):
		
		"""Shows all three types of ROI indices.
		
		Shows:
			
			* Mesh indices by calling :py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.showMeshIdx`
			* Image indices by calling :py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.showImgIdx`
			* Extended indices by calling :py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.showExtImgIdx`
		
		Keyword Args:
			axes (list): List of axes to plot in (at least ``len(axes)=3``). 
		
		Returns:
			list: List of ``matplotlib.axes``.	
		
		"""
		
		wereNone=False
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["MeshIdx","imgIdx","extImgIdx"],sup=self.name+" Idx",proj=['3d',None,None])
			wereNone=True
			
		self.showMeshIdx(axes[0])
		self.showImgIdx(axes[1])
		self.showExtImgIdx(axes[2])
		
		if wereNone:
			axes[0].set_title(self.name+" Idx")
			
		return axes
	
	def checkSymmetry(self,debug=False):
			
		img=np.zeros((self.embryo.dataResPx,self.embryo.dataResPx))
		img[self.imgIdxX,self.imgIdxY]=1
		return pyfrp_img_module.symmetryTest(img,debug=debug)
		
	def idxs2Quad(self,debug=False):
		
		if not self.checkSymmetry():
			printWarning('Cannot reduce region '+self.name+' to quadrant. Indices are not symmetric.')
			return self.getAllIdxs()
		
		self.imgIdxX,self.imgIdxY=imgIdx2Quad(self.imgIdxX,self.imgIdxY,self.embryo.dataResPx,debug=debug)
		
		if 'Quad' not in self.embryo.geometry.typ:
			printWarning('Will not update mesh indices, geometry is not set to quad.')
		else:	
			self.computeMeshIdx()
		
		return self.getAllIdxs()
	
	def idxs2Full(self):
		return self.computeIdxs()
	
	def resetDataVec(self):
		
		"""Resets data vector to an empty list"""
		
		self.setDataVec([])
		return self
	
	def resetSimVec(self):
		
		"""Resets simulation vector to an empty list"""
		
		self.setSimVec([])
		return self
	
	def plotData(self,ax=None,color=None,linewidth=1,legend=True,linestyle='-',label=None,legLoc=-1):
		
		"""Plot data vector of ROI.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linestyle (str): Linestyle of plot.
			linewidth (float): Linewidth of plot.
			legend (bool): Show legend.
			legLoc (int): Location of legend.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	
		
		"""
		
		if color==None:
			color=self.color
		
		if label==None:
			label=self.name + ' simulated'
		
		ax = pyfrp_plot_module.plotTS(self.embryo.tvecData,self.dataVec,ax=ax,linewidth=linewidth,color=color,label=self.name + ' data',
		title="Data",sup=self.name+" data",linestyle=linestyle,legend=legend,legLoc=legLoc)
		
		return ax
	
	def plotDataPinned(self,ax=None,color=None,linewidth=1,legend=True,linestyle='-',label=None,legLoc=-1):
		
		"""Plot pinned data vector of ROI.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linestyle (str): Linestyle of plot.
			linewidth (float): Linewidth of plot.
			legend (bool): Show legend.
			legLoc (int): Location of legend.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	
		
		"""
		
		if color==None:
			color=self.color
		
		if label==None:
			label=self.name + ' simulated'
		
		ax = pyfrp_plot_module.plotTS(self.embryo.tvecData,self.dataVecPinned,ax=ax,linewidth=linewidth,color=color,label=self.name + ' data',
		title="Data Pinned",sup=self.name+" data",linestyle=linestyle,legend=legend,legLoc=legLoc)
		
		return ax
	
	def plotSim(self,ax=None,color=None,linewidth=1,legend=True,linestyle='--',label=None,legLoc=-1):
		
		"""Plot simulation vector of ROI.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linestyle (str): Linestyle of plot.
			linewidth (float): Linewidth of plot.
			legend (bool): Show legend.
			legLoc (int): Location of legend.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	
		
		"""
		
		if color==None:
			color=self.color
		
		if label==None:
			label=self.name + ' simulated'
		
		ax = pyfrp_plot_module.plotTS(self.embryo.simulation.tvecSim,self.simVec,ax=ax,linewidth=linewidth,color=color,
		label=label,title="Simulation",sup=self.name+" simulation",linestyle=linestyle,legend=legend,legLoc=legLoc)
			
		return ax
	
	def plotSimPinned(self,ax=None,color=None,linewidth=1,legend=True,linestyle='--',label=None,legLoc=-1):
		
		"""Plot pinned simulation vector of ROI.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linestyle (str): Linestyle of plot.
			linewidth (float): Linewidth of plot.
			legend (bool): Show legend.
			legLoc (int): Location of legend.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	
		
		"""
		
		if color==None:
			color=self.color
		
		if label==None:
			label=self.name + ' simulated'
		
		ax = pyfrp_plot_module.plotTS(self.embryo.simulation.tvecSim,self.simVecPinned,ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + ' simulated',title="Simulation Pinned",sup=self.name+" simulation",linestyle=linestyle,legend=legend,legLoc=legLoc)
			
		return ax

	def findIncluded(self):
		
		"""Returns list of :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI` objects
		in which ROI is included.
		
		Returns:
			list: List of customROIs.
		
		"""
		
		incl=[]
		for r in self.embryo.ROIs:
			if type(r) is customROI:
				if r.roiIncluded(self):
					incl.append(r)
		return incl			
					
	def isMaster(self):
		
		"""Returns if ROI is masterROI.
		
		Returns:
			bool: True if masterROI.
			
		"""	
		
		return self==self.embryo.getMasterROI()
	
	def getMaxExtendPlane(self):
		
		"""Returns in which plane ("xy","xz","yz") the ROI has the biggest extend.
		
		Returns:
			str: Plane with largest extend.
		"""
		
		#Compute extend for all three coordinat4e
		xmin,xmax,ymin,ymax,zmin, zmax = self.getExtend()
		ext=[abs(xmax-xmin),abs(ymax-ymin),abs(zmax-zmin)]
		
		#Get 2 largest values
		mExts,indExts=pyfrp_misc_module.getIdxOfNLargest(ext,2)
		
		plane=""
		if 0 in indExts:
			plane=plane+"x"
		if 1 in indExts:
			plane=plane+"y"
		if 2 in indExts:
			plane=plane+"z"
			
		if len(plane)!=2:
			printWarning("Something went wrong finding plane. Plane is " + plane)
		
		return plane
	
	def getPlaneMidCoordinate(self):
		
		"""Returns midpoint of extend orthogonal to plane of maximum extension.
		
		Returns:
			float: Midpoint.
		
		"""
		
		plane=self.getMaxExtendPlane()
		xmin,xmax,ymin,ymax,zmin, zmax = self.getExtend()
		
		if plane=='xy':
			return (zmax+zmin)/2.
			
		elif plane=='xz':
			return (ymax+ymin)/2.
			
		elif plane=='yz':
			return (xmax+xmin)/2.
	
	def getOrthogonal2Plane(self):
		
		"""Returns orthogonal direction to plane of maximum extension.
		
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.getPlaneMidCoordinate` and
		:py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.getMaxExtendPlane` .
		
		Returns:
			str: Direction.
		
		"""
		
		plane=self.getMaxExtendPlane()
		
		if plane=='xy':
			return 'z'
			
		elif plane=='xz':
			return 'y'
			
		elif plane=='yz':
			return 'x'
		
	def getExtend(self):
		
		"""Returns x-/y-/z-extend of ROI.
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
				* zmin (float): Minimum z-coordinate.
				* zmax (float): Maximum z-coordinate.
					
		"""
		
		[xmin,xmax],[ymin,ymax] = self.computeXYExtend()
		zmin, zmax = self.getZExtend()
		
		return xmin,xmax,ymin,ymax,zmin, zmax 
		
	def plotSolutionVariable(self,phi,ax=None,vmin=None,vmax=None,nlevels=25,colorbar=True,plane='xy',zs=None,zdir=None,mask=True,nPts=1000,mode='normal',title="Solution Variable",
			  typ='contour'):
		
		"""Plots simulation solution variable over all indices of ROI as 2D contour plot.
		
		.. note:: If no ``ax`` is given, will create new one.
		
		``plane`` variable controls in which plane the solution variable is supposed to be plotted. 
		Acceptable input variables are ``"xy","xz","yz"``. See also
		:py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.getMaxExtendPlane`.
		
		See also http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.tricontourf .
		
		
		.. warning:: ``matplotlib.pyplot.tricontourf`` has problems when ``phi`` only is in a single level of contour plot.
		   To avoid this, we currently add some noise in this case just to make it plottable. This is not the most elegant
		   solution.
		
		You can find a more detailed explanation in the documentation of :py:func:`pyfrp.modules.pyfrp_plot_module.plotSolutionVariable`.
		
		Args:
			phi (fipy.CellVariable): Solution variable.
			
		Keyword Args:
			ax (matplotlib.axes): Axes used for plotting.
			vmin (float): Minimum value displayed in contour plot.
			vmax (float): Maximum value displayed in contour plot.
			nlevels (int): Number of contour levels.
			colorbar (bool): Display color bar.
			plane (str): Plane in which solution variable is supposed to be plotted.
			zs (float): In case of a 3D plot, height in direction zdir where to put contour.
			zdir (str): Orthogonal direction to plane.
			nPts (int): Number of points used for interpolating (only if ``mode=normal``).
			mode (str): Which contour function to use.
			title (str): Title of plot.
			typ (str): Type of plot.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		#Get values of phi if necessary
		if hasattr(phi,'value'):
			val=np.asarray(phi.value)
		else:
			val=np.asarray(phi)
		
		#Get x/y/z/values coordinates in ROI
		x,y,z=self.embryo.simulation.mesh.getCellCenters()
		x=x[self.meshIdx]
		y=y[self.meshIdx]
		z=z[self.meshIdx]
		val=val[self.meshIdx]
		
		#Finding distance threshold later used for masking nodes.
		dmaxX,dmaxY,dmaxZ=self.getMaxNodeDistance()
		dFact=1.5
		
		if plane=='xy':
			X=x
			Y=y
			D=dFact*max([dmaxX,dmaxY])
		elif plane=='xz':	
			X=x
			Y=z
			D=dFact*max([dmaxX,dmaxZ])
		elif plane=='yz':
			X=y 
			Y=z
			D=dFact*max([dmaxY,dmaxZ])
		else:
			printError("Don't understand plane="+plane+". Will not plot.")
			return None
		
		if not mask:
			D=None
		
		ax=pyfrp_plot_module.plotSolutionVariable(x,y,val,ax=ax,vmin=vmin,vmax=vmax,nlevels=nlevels,colorbar=colorbar,
					    plane=plane,zs=zs,zdir=zdir,sup=self.name,dThresh=D,nPts=nPts,mode=mode,title=title,typ=typ)
		
		return ax
	
	def getSimConc(self,phi,append=True):
		
		"""Computes the simulation concentration over ROI.
		
		Args:
			phi (fipy.CellVariable): Solution variable.
			
		Keyword Args:
			append (bool): Append result to simulation vector.
		
		Returns:
			float: Simulation concentration over ROI.
		
		
		"""
		
		cvs=self.embryo.simulation.mesh.getCellVolumes()
		
		c=pyfrp_integration_module.getAvgConc(phi,cvs,self.meshIdx)
		
		if append:
			self.simVec.append(c)
			
		return c	
	
	def pinAllTS(self,bkgdVal=None,normVal=None,bkgdValSim=None,normValSim=None,debug=False):
		
		"""Pins both data and simulation timeseries of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_fit_module.pinConc`.
		
		Keyword Args:
			bkgdVal (float): Use this background value instead of newly computing it.
			normVal (float): Use this norming value instead of newly computing it.
			bkgdValSim (float): Use this background value for simulation timeseries instead of newly computing it.
			normValSim (float): Use this norming value for simulation timeseries instead of newly computing it.
			debug (bool): Print debugging messages.
		
		Returns:
			tuple: Tuple containing:
			
				* dataVecPinned (numpy.ndarray): Pinned data vector.
				* simVecPinned (numpy.ndarray): Pinned simulation vector.
		
		"""
		
		bkgdValSim=pyfrp_misc_module.assignIfVal(bkgdValSim,bkgdVal,None)
		normValSim=pyfrp_misc_module.assignIfVal(normValSim,normVal,None)
		
		self.dataVecPinned = self.pinDataTS(bkgdVal=bkgdVal,normVal=normVal,debug=debug)
		self.simVecPinned = self.pinSimTS(bkgdVal=bkgdValSim,normVal=normValSim,debug=debug)
		
		return self.dataVecPinned,self.simVecPinned
	
	def pinDataTS(self,bkgdVal=None,normVal=None,debug=False):
		
		"""Pins data timeseries of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_fit_module.pinConc`.
		
		Keyword Args:
			bkgdVal (float): Use this background value instead of newly computing it.
			normVal (float): Use this norming value instead of newly computing it.
			debug (bool): Print debugging messages.
		
		Returns:	
			numpy.ndarray: Pinned data vector.		
		
		"""
		
		if bkgdVal==None or normVal==None:
			bkgdValTemp,normValTemp = self.embryo.computePinVals(debug=debug)
			bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
			normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
			
		self.dataVecPinned=pyfrp_fit_module.pinConc(self.dataVec,bkgdVal,normVal,axes=None,debug=debug,tvec=self.embryo.tvecData,color=self.color)
		
		return self.dataVecPinned
	
	def pinSimTS(self,bkgdVal=None,normVal=None,debug=False):
		
		"""Pins simulation timeseries of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_fit_module.pinConc`.
		
		Keyword Args:
			bkgdVal (float): Use this background value instead of newly computing it.
			normVal (float): Use this norming value instead of newly computing it.
			debug (bool): Print debugging messages.
		
		Returns:	
			numpy.ndarray: Pinned simulation vector.		
		
		"""
	
		if bkgdVal==None or normVal==None:
			bkgdValTemp,normValTemp = self.embryo.computePinVals(debug=debug)
			bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
			normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
		
		self.simVecPinned=pyfrp_fit_module.pinConc(self.simVec,bkgdVal,normVal,axes=None,debug=debug,tvec=self.embryo.simulation.tvecSim,color=self.color)
		
		return self.simVecPinned
	
	def getFittedVec(self,fit):
		
		"""Returns fitted simulation vector of ROI of given fit.
		
		.. note:: To avoid crashes, function returns empty list
		   if ROI is in ``ROIsFItted`` but has not been fitted yet.
		   Also inserts an empty list at the respective index.
		
		Args:
			fit (pyfrp.subclasses.pyfrp_fit): Fit object.
			
		Returns:
			numpy.ndarray: Fitted simulation vector.
		
		"""
		
		try:
			return fit.fittedVecs[fit.ROIsFitted.index(self)]
		except IndexError:
			fit.fittedVecs.insert(fit.ROIsFitted.index(self),[])
			#return fit.fittedVecs[fit.ROIsFitted.index(self)]	
			#This solves a problem only temporarily.
			return []
		
	def getdataVecFitted(self,fit):
		
		"""Returns fitted data vector of ROI of given fit.
		
		.. note:: To avoid crashes, function returns empty list
		   if ROI is in ``ROIsFItted`` but has not been fitted yet.
		   Also inserts an empty list at the respective index.
		
		Args:
			fit (pyfrp.subclasses.pyfrp_fit): Fit object.
			
		Returns:
			numpy.ndarray: Fitted data vector.
		
		"""
		
		try:
			return fit.dataVecsFitted[fit.ROIsFitted.index(self)]
		except IndexError:
			fit.dataVecsFitted.insert(fit.ROIsFitted.index(self),[])
			return []
		
	def plotFit(self,fit,ax=None,color=None,linewidth=1,legend=True,title=None,linestyles=['-','-.'],show=True):
		
		"""Plot fit for ROI.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linestyles (list): Linestyles of data and simulation.
			linewidth (float): Linewidth of plot.
			legend (bool): Show legend.
			show (bool): Show figure.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	
		
		"""
		
		if color==None:
			color=self.color
		
		if title==None:
			title="Fit "+fit.name
			
		ax = pyfrp_plot_module.plotTS(fit.tvecFit,self.getFittedVec(fit),ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + fit.name,title=title,sup=self.name+" fitted",linestyle=linestyles[1],legend=legend,show=show)
		
		ax = pyfrp_plot_module.plotTS(fit.tvecFit,self.getdataVecFitted(fit),ax=ax,linewidth=linewidth,color=color,
		label=self.name + ' ' + fit.name,title=title,sup=self.name+" fitted",linestyle=linestyles[0],legend=legend,show=show)
		
		return ax
	
	def isAnalyzed(self):
		
		"""Checks if ROI has been analyzed.
		
		Returns:
			bool: True if ROI has been analyzed.
		
		"""
		
		return len(self.embryo.tvecData)==len(self.dataVec)
	
	def isSimulated(self):
		
		"""Checks if ROI has been simulated.
		
		Returns:
			bool: True if ROI has been simulated.
		
		"""
		
		if self.embryo.simulation!=None:
			return len(self.embryo.simulation.tvecSim)==len(self.simVec)
		return False
	
	def isFitted(self):
		
		"""Checks if ROI has been fitted in ALL fits of embryo.
		
		Returns:
			bool: True if ROI has been fitted.
		
		"""
		
		fitted=False
		for fit in self.embryo.fits:
			if self in fit.ROIsFitted and len(self.getFittedVec(fit))>0:
				fitted=True
		return fitted
	
	def getInterpolationError(self):
		
		"""Prints out interpolation error for the volume of this ROI.
		
		Interpolation error is defined as:
		
		``dataVec[0]/simVec[0]``,
		
		That is, by how much does the first simulation value defer from first data
		value.
		
		Returns:
			float: Interpolation error.
		
		"""
		
		if self.isSimulated() and self.isAnalyzed():
			try:
				return self.dataVec[0]/self.simVec[0]
			except ZeroDivisionError:
				printWarning("Dividing by zero. Going to return infinity.")
				return np.inf
		else:
			printWarning("ROI is either not simulated or analyzed. Cannot return interpolation error.")
			return 0.
	
	def getEncapsulatingBox(self):
		
		"""Returns encapsulating box of ROI.
		
		That is, a box defined by ``[xmin,xmax],[ymin,ymax],[zmin,zmax]``
		in which ROI lies fully within.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List describing extend in x-direction (``[xmin,xmax]``).
				* yExtend (list): List describing extend in y-direction (``[ymin,ymax]``).
				* zExtend (list): List describing extend in z-direction (``[zmin,zmax]``).
				
		"""
		
		xExtend,yExtend=self.computeXYExtend()
		zExtend=self.getZExtend()
		return xExtend,yExtend,zExtend
		
	def refineInMeshByField(self,factor=3.,addZ=15.,findIdxs=True,debug=False,run=True,fnOut=None):
		
		"""Refines mesh inside ROI by adding box field to mesh file.
		
		The mesh size inside the box is computed by ``mesh.volSizePx/factor``. To ensure
		that there are enough original nodes inside ROI that then allow refinement from,
		``addZ`` pixels is added in z-direction both below and above the ROI.
		
		See also :py:func:`pyfrp.subclasses.pyfrp_mesh.mesh.addBoxField`.
		
		Keyword Args:
			factor (float): Refinement factor.
			addZ (float): Number of pixels added above and below ROI for box field.
			findIdxs (bool): Find mesh indices of ROI after refinement.
			run (bool): Run Gmsh to generate new mesh after refinement.
			debug (bool): Print debugging messages.
			fnOut (str): Path to output geo file.
			
		Returns:
			str: Path to new .geo file.
			
		"""
		
		xExtend,yExtend,zExtend=self.getEncapsulatingBox()
		zExtend=[zExtend[0]-addZ,zExtend[1]+addZ]
		
		if debug:
			print("Adding Box Field for ROI " + self.name)
			print("Mesh Nodes in ROI before: ", len(self.meshIdx))
			
		fnOut=self.embryo.simulation.mesh.addBoxField(self.embryo.simulation.mesh.volSizePx/factor,xExtend,yExtend,zExtend,comment=self.name+" field",run=run,fnOut=fnOut)
	
		if findIdxs:
			self.computeMeshIdx(self.embryo.simulation.mesh)
		
		if debug and findIdxs:
			print("Mesh Nodes in ROI after: ", len(self.meshIdx))
			
		return fnOut
	
	def adaptRefineInMeshByField(self,nNodesReq,factor=3.,addZ=15.,zIncrement=1.,fIncrement=1.,nNodesMax='inf',debug=False,ROIReq=None,fnOut=None):
		
		"""Refines mesh inside ROI adaptively until a given number of nodes inside ROI 
		is reached.
		
		Does this by:
			
			* Refining through :py:func:`refineInMeshByField`.
			* Computing mesh indices via :py:func:`computeMeshIdx`.
			* If number of nodes did not change, increase ``addZ``, else increase ``factor``.
			* Check if desired number of nodes is reached or not, if not, repeat.
			
		.. note:: If the new number of nodes in the ROI exceeds ``nNodesMax``, will revert the last step
		   and perform the other operation, e.g. increasing ``addZ`` instead of ``factor`` and vice versa.
		
		.. note:: If ``ROIReq`` is given, will try to refine in ``self`` such that ``ROIReq`` has at least ``nNodesReq``
		   mesh nodes. If it is not given, ``nNodesReq`` refers to the nodes in ``self``.
		
		Args:
			nNodesReq (int): Desired number of nodes inside ROI.
		
		Keyword Args:
			factor (float): Refinement factor.
			addZ (float): Number of pixels added above and below ROI for box field.
			zIncrement (float): Number of pixels addZ is increased per adaptive step.
			fIncrement (float): Stepsize of refinement factor.
			nNodesMax (float): Maximum number of nodes allowed in ROI.
			debug (bool): Print debugging messages.
			ROIReq (pyfrp.subclasses.pyfrp_ROI.ROI): The ROI object that is referred to with nNodesReq.
			fnOut (str): Path to output geo file.
			
		Returns:
			int: Final number of nodes in ROI.
			
		"""
		
		#Convert nNodesMax if necessary
		nNodesMax=pyfrp_misc_module.translateNPFloat(nNodesMax)
		
		#Get current node numbers
		if ROIReq==None:
			self.computeMeshIdx(self.embryo.simulation.mesh)
			nNodes=len(self.meshIdx)
			nNodesAll=self.embryo.simulation.mesh.getNNodes()
		else:
			ROIReq.computeMeshIdx(ROIReq.embryo.simulation.mesh)
			nNodes=len(ROIReq.meshIdx)
			nNodesROIReq=len(ROIReq.meshIdx)
			nNodesAll=ROIReq.embryo.simulation.mesh.getNNodes()
		nNodesROI=len(self.meshIdx)
		
		#Init flags
		mode=0
		i=0
		
		#As long as requirement isn't met, refine
		while nNodes<nNodesReq:
			
			self.refineInMeshByField(factor=factor,addZ=addZ,findIdxs=True,debug=False,run=True,fnOut=fnOut)
			
			#Compute updated idxs
			nNodesAllNew=self.embryo.simulation.mesh.getNNodes()
			if ROIReq==None:
				nNodesNew=len(self.meshIdx)
			else:
				ROIReq.computeMeshIdx(ROIReq.embryo.simulation.mesh)
				nNodesNew=len(ROIReq.meshIdx)
				nNodesROIReqNew=len(ROIReq.meshIdx)
			nNodesROINew=len(self.meshIdx)
			
			#Print out current status
			if debug:
				print("Iteration ", i, ". ")
				print("Current parameters: addZ = ", addZ, " factor = ", factor) 
				print("Total mesh nodes: ", nNodesAllNew)
				print("Mesh Nodes in ROI before refinement: " , nNodesROI, " and after ", nNodesROINew, ".")
				if ROIReq!=None:
					print("Mesh Nodes in ROIReq before refinement: " , nNodesROIReq, " and after ", nNodesROIReqNew, ".")
			
			#Check if nNodes requirement is met
			if nNodesNew<nNodesReq:
				if nNodesAllNew==nNodesAll:
					if debug: 
						print("nNodesAll did not change, will increase addZ by ",zIncrement,". \n")
					addZ=addZ+zIncrement
					mode=0
				else:
					if debug:
						print("nNodes not large enough yet, will increase factor by ",fIncrement,". \n")
					factor=factor+fIncrement
					mode=1
			
			#Check if maximum number of nodes was exceeded
			elif nNodesNew>nNodesMax:
				if debug:
					print("Number of nodes exceeded maximum allowed number", nNodesMax, ".")
				
				if mode==0:
					if debug:
						print("Previously tried to increase addZ. Will try old addZ, but increase factor by " ,fIncrement,". \n")
					addZ=addZ-zIncrement
					factor=factor+fIncrement
					mode=1
				elif mode==1:
					if debug:
						print("Previously tried to increase factor. Will try old factor, but increase addZ by " ,zIncrement,". \n")
					addZ=addZ+zIncrement
					factor=factor-fIncrement
					mode=0
						
			i=i+1
			
			#Update old counter
			nNodes=nNodesNew
			nNodesAll=nNodesAllNew
			nNodesROI=nNodesROINew
			if ROIReq!=None:
				nNodesROIReq=nNodesROIReqNew
			
		return nNodes
				
	def printDetails(self):
		
		"""Prints out all attributes of ROI object."""
		
		print("ROI ", self.name, " details:")
		printAllObjAttr(self)
		print() 
	
	def plotSimConcProfile(self,phi,ax=None,direction='x',mode='normal',nbins=20,color=None,label=None,legend=False):
		
		"""Plots concentration profile of solution variable in 
		single direction.
		
		``mode`` can be either ``"normal"`` or ``"hist"``. If ``mode="hist"``, will plot a histogram with ``nbins`` bins using
		:py:func:`pyfrp.modules.pyfrp_misc_module.simpleHist`.
		
		.. note:: ``direction`` sets in which direction the profile should be plotted. if ``direction="r"``, then function
		   will plot a radial profile and uses ``self.embryo.geometry.center`` as center if ROI does not have a center,
		   else uses center of ROI.
		
		.. note:: Will create axes if not given via ``ax``.
		
		Example:
		
		Grab ROI:
		
		>>> sl=emb.getROIByName("Slice")
		
		Make some plot:
		
		>>> fig,axes=pyfrp_plot_module.makeSubplot([2,2])
		
		Plot some concentration profiles:
		
		>>> ax=sl.plotSimConcProfile(emb.simulation.IC,mode='hist',color='g',label='direction = x',nbins=100,ax=axes[0],legend=False)
		>>> ax=sl.plotSimConcProfile(emb.simulation.IC,mode='hist',direction='y',color='r',label='direction = y',nbins=100,ax=axes[1],legend=False)
		>>> ax=sl.plotSimConcProfile(emb.simulation.IC,mode='hist',direction='r',color='b',nbins=100,label='direction = r',ax=axes[2],legend=False)
		>>> ax=sl.plotSimConcProfile(emb.simulation.IC,mode='normal',direction='r',color='b',label='direction = r',ax=axes[3],legend=False)
		
		.. image:: ../imgs/pyfrp_ROI/plotSimConcProfile.png
		
		Args:
			phi (fipy.CellVariable): Solution variable
			
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			direction (str): Direction to be plotted (x/y/z/r).
			color (str): Color of plot.
			legend (bool): Show legend.
			label (str): Label of plot.
			nbins (int): Number of bins of histogram.
			mode (str): Either ``normal`` or ``hist``.
			
		Returns:
			matplotlib.axes: Matplotlib axes used for plotting.
		
		"""
		
		if color==None:
			color=self.color
		if label==None:
			label=self.name
		
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Concentration profile"],sup=self.name+" phi")
			ax=axes[0]
		
		if direction=='x':
			x=self.embryo.simulation.mesh.getCellCenters()[0]
		elif direction=='y':
			x=self.embryo.simulation.mesh.getCellCenters()[1]
		elif direction=='z':
			x=self.embryo.simulation.mesh.getCellCenters()[2]
		elif direction=='r':
			if hasattr(self,'center'):
				center=self.center
			else:
				center=self.embryo.geometry.center
			x=np.sqrt((self.embryo.simulation.mesh.getCellCenters()[0]-center[0])**2+(self.embryo.simulation.mesh.getCellCenters()[1]-center[1])**2)
		else:
			printError('Direction '+ direction+ 'unknown. Will not plot.')
			return ax
		
		x=np.asarray(x)[self.meshIdx]
		if hasattr(phi,'value'):
			v=np.asarray(phi.value)[self.meshIdx]
		else:
			v=np.asarray(phi)[self.meshIdx]
			
		vSorted,xSorted=pyfrp_misc_module.sortListsWithKey(v,x)
		
		if mode=='hist':
			xSorted,vSorted=pyfrp_misc_module.simpleHist(xSorted,vSorted,bins=nbins)
		
		pyfrp_plot_module.plotTS(xSorted,vSorted,color=color,label=label,legend=legend,ax=ax)
		ax.set_xlabel(direction)
		ax.set_ylabel("Concentration")
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def getCopy(self):
		
		"""Returns deepcopy of ROI object.
		
		Uses ``copy.copy`` to generate copy of object, see also https://docs.python.org/2/library/copy.html .
		
		``copy.deepcopy`` also generates copies of other objects, including ``ROI.embryo``.
		
		"""
		
		return copy.copy(self)
	
	def getNMeshNodes(self):
		
		"""Returns number of mesh indices inside ROI.
		
		Returns:
			int: Number of nodes.
			
		"""
		
		return len(self.meshIdx)
	
	def getNImgPxs(self):
		
		"""Returns number of image pixels inside ROI.
		
		Returns:
			int: Number of indices.
			
		"""
		
		return len(self.imgIdxX)
	
		
	def getMaxNodeDistance(self):
		
		"""Returns maximum node distance in x/y/z direction 
		for all nodes in ROI.
		
		Returns:
			tuple: Tuple containing:
			
				* dmaxX (float): Maximum distance in x-direction
				* dmaxY (float): Maximum distance in y-direction
				* dmaxZ (float): Maximum distance in z-direction
				
		"""
		
		distances=self.embryo.simulation.mesh.mesh.cellDistanceVectors
		
		return max(distances[0][self.meshIdx]),max(distances[1][self.meshIdx]),max(distances[2][self.meshIdx])
		
	def genGmshDomain(self,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Translates ROI into gmsh domain object.
		
		This object can then be used to write ROIs to ``.geo`` files.
				
		.. note:: If ``minID==None``, will grab maximum ID via :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getMaxGeoID` and add 1.
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Domain object.
		
		"""
		
		if minID==None:
			minID=self.embryo.geometry.getMaxGeoID()+1
		
		d=pyfrp_gmsh_geometry.domain()
		printWarning("This ROI type does not have genGmshDomain right now. This might change in further versions.")
		return d
	
	def writeToGeoFile(self,fn=None,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .geo file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.geo`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.genGmshDomain`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
			
		Returns:
			str: Path to geo file.
		
		"""
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".geo"
			fn=fn.replace(" ","_")
			
		printWarning("ROI of type "+str(self.getType)+" does not have writeToGeoFile right now. This might change in further versions.")
		
		return fn
	
	def genMeshFile(self,fn=None,volSizePx=20.,debug=False,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .msh file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.msh`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.writeToGeoFile`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
			
		Returns:
			str: Path to mesh file.
		
		"""
		
		#Check filename
		if fn!=None:
			if not fn.endswith(".geo"):
				fn,ext=os.path.splitext(fn)
				fn=fn+".geo"
		
		#Make geo file
		fn=self.writeToGeoFile(fn=fn,volSizePx=volSizePx,minID=minID)
		fnMsh=fn.replace(".geo",".msh")
		
		#Run gmsh
		pyfrp_gmsh_module.runGmsh(fn,fnOut=fnMsh,debug=debug,volSizeMax=volSizePx)
	
		return fnMsh
	
	def genAsOpenscadInGeometry(self):
		
		"""Generates intersection between ROI and geometry as solid python object.
		
		See also :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.genAsOpenscad` and 
		:py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.genAsOpenscad`.
		
		Returns:
			solid.solidpython.cylinder: Solid python object. 
		
		"""
		
		openScadGeo=self.embryo.geometry.genAsOpenscad()
		openScadROI=self.genAsOpenscad()
		
		return openScadGeo*openScadROI
	
	def render2Openscad(self,fn=None,segments=48):
		
		"""Generates .scad file for the ROI.
		
		.. note:: If ``fn`` is not given, will save .scad file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.scad``.
		
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		Returns:
			str: Output filename.
		
		"""
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".scad"
			fn=fn.replace(" ","_")
			
		solid.scad_render_to_file(self.genAsOpenscad(), filepath=fn,file_header='$fn = %s;' % segments, include_orig_code=False)
		
		return fn
		
	def render2OpenscadInGeometry(self,fn=None,segments=48):
		
		"""Generates .scad file for the intersection between ROI and geometry.
		
		.. note:: If ``fn`` is not given, will save .scad file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.scad``.
		
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		Returns:
			str: Output filename.
		
		"""
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".scad"
			fn=fn.replace(" ","_")
			
		solid.scad_render_to_file(self.genAsOpenscadInGeometry(), filepath=fn,file_header='$fn = %s;' % segments, include_orig_code=False)
		
		return fn
	
	def render2Stl(self,fn=None,segments=48):
		
		"""Generates .stl file for the ROI.
		
		Will do this by:
		
			* Generating openscad object via :py:func:`genAsOpenscad`.
			* Rendering this to scad file via :py:func:`render2Openscad`.
			* Calling :py:func:`pyfrp.modules.pyfrp_openscad_module.runOpenscad`.
		
		.. note:: If ``fn`` is not given, will save .stl file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.stl``.
				
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		Returns:
			str: Output filename.
		
		"""
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".scad"
			fn=fn.replace(" ","_")
		
		fnStl=fn
		fnScad=fnStl.replace(".stl",".scad")
		
		self.render2Openscad(fn=fnScad,segments=segments)
		pyfrp_openscad_module.runOpenscad(fnScad,fnOut=fnStl)
		
		return fnStl
	
	def render2StlInGeometry(self,fn=None,segments=48):
		
		"""Generates .stl file for the intersection between ROI and geometry.
		
		Will do this by:
		
			* Generating openscad object via :py:func:`genAsOpenscadInGeometry`.
			* Rendering this to scad file via :py:func:`render2OpenscadInGeometry`.
			* Calling :py:func:`pyfrp.modules.pyfrp_openscad_module.runOpenscad`.
		
		.. note:: If ``fn`` is not given, will save .stl file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.stl``.
				
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		Returns:
			str: Output filename.
		
		"""
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".scad"
			fn=fn.replace(" ","_")
		
		fnStl=fn
		fnScad=fnStl.replace(".stl",".scad")
		
		self.render2OpenscadInGeometry(fn=fnScad,segments=segments)
		pyfrp_openscad_module.runOpenscad(fnScad,fnOut=fnStl)
		
		return fnStl
	
	def addBoundaryLayerAtSurfaces(self,fn=None,segments=48):
		
		"""Adds boundary layer around ROI to the mesh.
		
		Does this by:
		
			* Generating a stl file describing ROI, see also :py:func:`pyfrp.subclasses.pyfrp_ROI.ROI.render2StlInGeometry`.
			* Read in stl file as new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.domain` via
			  :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.readStlFile`.
			* Simplify new geometry via :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.simplifySurfaces`.
			* Extracting selected surfaces via :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.gmshElement.extract`.
			* If selected, surface boundaries are approximated into splines via 
			  :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.gmshElement.extract`.
			* Reading in geometry's .geo file via :py:func:`pyfrp.sublcasses.pyfrp_geometry.geometry.readGeoFile`.
			* Merging ROI geometry into main geometry via :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.merge`.
			* Adding a boundary layer mesh via :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addBoundaryLayerField`.
			* Adding all surfaces of ROI's domain to boundary layer, see 
			  :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField.addFaceListByID`.
			* Writing new .geo file.
			* Setting new .geo file as ``fnGeo``.
			* Running :py:func:`genMesh`.
			* Clean up .stl and .scad files that are not needed anymore.
		
		.. note:: ``volSizeLayer`` only allows a single definition of mesh size in layer. Note that the 
		   :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField` class allows different mesh sizes
		   normal and along surfaces. For more information, see its documentation.
		
		.. note:: If no ``fnOut`` is given, will save a new .geo file in same folder as original ``fnGeo`` with subfix:
		   ``fnGeo_roiName_BL.geo``.
		
		.. note:: :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.simplifySurfaces` is not a simple procedure, 
		   we recommend reading its documentation.
		
		If ``volSizePx`` is given, will overwrite mesh's ``volSizePx`` and set it globally at all nodes.
		
		Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): An ROI.
		
		Keyword Args:
			fnOut (str): Path to new .geo file.
			segments (int): Number of segments used for convex hull of surface.
			simplify (bool): Simplify surfaces of stl file.
			iterations (int): Number of iterations used for simplification.
			triangIterations (int): Number of iterations used for subdivision of surfaces.
			addPoints (bool): Allow adding points inside surface triangles.
			fixSurfaces (bool): Allow fixing of surfaces, making sure they are coherent with Gmsh requirements.
			debug (bool): Print debugging messages.
			volSizePx (float): Global mesh density.
			volSizeLayer (float): Boundary layer mesh size.
			thickness (float): Thickness of boundary layer.
			cleanUp (bool): Clean up temporary files when finished.
			approxBySpline (bool): Approximate curvatures by spline.
			angleThresh (float): Threshold angle under which loops are summarized.
			faces (list): List of faces.
			onlyAbs (bool): Take absolute value of faces into account.
		
		Returns:
			str: Path to new .geo file.
		
		"""
		
		if hasattr(self.embryo,'simulation'):
			printError("addBoundaryLayerAtSurfaces: Embryo does not have simulation yet.")
			return ""
		
		return self.embryo.simulation.mesh.addBoundaryLayerAroundROI(self,fnOut=fnOut,segments=segments,simplify=simplify,iterations=iterations,triangIterations=triangIterations,
							fixSurfaces=fixSurfaces,debug=debug,volSizePx=volSizePx,volSizeLayer=volSizeLayer,thickness=thickness,cleanUp=cleanUp,
							approxBySpline=approxBySpline,angleThresh=angleThresh,faces=faces,onlyAbs=onlyAbs)	
			
		
	
class radialROI(ROI):
	
	"""Radial ROI class.
	
	Inherits from :py:class:`ROI`. 
	
	Main attributes are:
	
		* ``radius``: Radius of ROI.
		* ``center``: Center of ROI.
	
	"""

	def __init__(self,embryo,name,Id,center,radius,color='b'):
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.radius=radius
		self.center=center
		
	def setRadius(self,r):
		
		"""Sets radius of ROI.
		
		Args:
			r (float): New radius
			
		Returns:
			float: New radius.
		
		"""
		
		self.radius=r
		return self.radius	
		
	def getRadius(self):
		
		"""Returns current radius of ROI.
		
		Returns:
			float: Current radius.
		
		"""
		
		return self.radius
	
	def setCenter(self,c):
		
		"""Sets radius of ROI.
		
		Args:
			c (list): New center.
			
		Returns:
			list: New center.
		
		"""
		
		self.center=c
		return self.center
	
	def getCenter(self):
		
		"""Returns current center of ROI.
		
		Returns:
			list: Current center.
		
		"""
		
		return self.center	
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getCircleIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getCircleIdxImg(self.center,self.radius,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getCircleIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		
		self.meshIdx=pyfrp_idx_module.getCircleIdxMesh(self.center,self.radius,radius,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		"""Shows ROI in a 2D plot.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linewidth (float): Linewidth of plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	

		"""
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		patch = ptc.Circle(self.center,self.radius,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def center2Mid(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setCenter([self.embryo.dataResPx/2+0.5,self.embryo.dataResPx/2+0.5])
		else:
			return self.setCenter([self.embryo.dataResPx/2,self.embryo.dataResPx/2])
	
	def makeReducable(self,auto=False,debug=False):
		
		oldCenter=self.getCenter()
		self.center2Mid()
		
		if not auto:
			a=input("Change center of region "+ self.name + " from " + str(oldCenter) + ' to ' + str(self.getCenter()) + ' ? [Y/N]')
		
			if a=='N':
				self.setCenter(oldCenter)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setCenter(oldCenter)
				return False
			return True
		return False
		
		
	def checkCentered(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return bool((self.center[0]==self.embryo.dataResPx/2+0.5) and (self.center[1]==self.embryo.dataResPx/2+0.5))
		else:
			return bool((self.center[0]==self.embryo.dataResPx/2.) and (self.center[1]==self.embryo.dataResPx/2.))
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideCircle`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideCircle(x,y,self.center,self.radius)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.center[0]-self.radius,self.center[0]+self.radius]
		self.yExtend=[self.center[1]-self.radius,self.center[1]+self.radius]
		return self.xExtend, self.yExtend
	
	def getCenterOfMass(self):
		
		"""Returns center of mass of ROI.
		
		For a radial ROI, this is equivalent to the ``center``.
		"""
		
		return self.center
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		.. note:: Will grab extent of geometry to find bounds in z-direction. 
		
		Returns:
			solid.solidpython.cylinder: Solid python object. 
		
		"""
		
		try:
			ext=self.embryo.geometry.getZExtend()
		except AttributeError:
			printError("genAsOpenscad: Cannot grab extend from geometry of type " + self.embryo.geometry.typ)
			
		openScadROI=solid.translate([self.center[0],self.center[1],min(ext)])(solid.cylinder(r=self.radius,h=abs(max(ext)-min(ext))))
		return openScadROI
		
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#slice ROI class
	
class sliceROI(ROI):
	
	def __init__(self,embryo,name,Id,height,width,sliceBottom,color='b'):
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.height=height
		self.width=width
		self.sliceBottom=sliceBottom
		
		self.computeZExtend()
		
	def computeZExtend(self):
		
		if self.sliceBottom:
			self.setZExtend(self.height,self.height+self.width)
		else:
			self.setZExtend(self.height-0.5*self.width,self.height+0.5*self.width)
		
		return self.zmin,self.zmax
		
	def setHeight(self,h):
		self.height=h
		self.computeZExtend()
		return self.height
	
	def getHeight(self):
		return self.height
	
	def setSliceBottom(self,s):
		self.sliceBottom=s
		self.computeZExtend()
		return self.sliceBottom
	
	def getSliceBottom(self):
		return self.sliceBottom
	
	def setWidth(self,w):
		self.width=w
		self.computeZExtend()
		return self.width
	
	def getWidth(self):
		return self.width
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getAllIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getAllIdxImg(self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getSliceIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		x,y,z=mesh.getCellCenters()
		self.meshIdx=pyfrp_idx_module.getSliceIdxMesh(z,self.zmin,self.zmax)
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		Only returns ``True``, since ``sliceROI`` is not limited in
		x/y-direction.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y], all ``True``.
		
		"""
		
		return np.ones(np.asarray(x).shape).astype(bool)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		.. note:: Since sliceROI theoretically is not having any limits in x/y-direction,
		   function returns limits given by input image, that is, ``[0,embryo.dataResPx]``.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[1,self.embryo.dataResPx]
		self.yExtend=[1,self.embryo.dataResPx]
		return self.xExtend, self.yExtend
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		.. note:: Will grab extent of geometry to find bounds in x/y-direction. 
		
		Returns:
			solid.solidpython.cube: Solid python object. 
		
		"""
		
		try:
			ext=self.embryo.geometry.getXYExtend()
		except AttributeError:
			printError("genAsOpenscad: Cannot grab extend from geometry of type " + self.embryo.geometry.typ)
		
		z=self.getOpenscadZExtend()
		
		zmin,zmax=min(z),max(z)
		
		openScadROI=solid.translate([ext[0],ext[2],zmin])(solid.cube([ext[1]-ext[0],ext[3]-ext[2],abs(zmax-zmin)]))
		return openScadROI
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Radial and slice ROI class
	
class radialSliceROI(sliceROI,radialROI):
	
	def __init__(self,embryo,name,Id,center,radius,height,width,sliceBottom,color='b'):
		
		radialROI.__init__(self,embryo,name,Id,center,radius,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
		
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getCircleIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getCircleIdxImg(self.center,self.radius,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getCircleIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
	
		self.meshIdx=pyfrp_idx_module.getCircleIdxMesh(self.center,self.radius,mesh,zmin=self.zmin,zmax=self.zmax)
		
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideCircle`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideCircle(x,y,self.center,self.radius)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.center[0]-self.radius,self.center[0]+self.radius]
		self.yExtend=[self.center[1]-self.radius,self.center[1]+self.radius]
		return self.xExtend, self.yExtend
	
	def genGmshDomain(self,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Translates ROI into gmsh domain object.
		
		This object can then be used to write ROIs to ``.geo`` files.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addCuboidByParameters`.
		
		.. note:: If ``minID==None``, will grab maximum ID via :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getMaxGeoID` and add 1.
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Domain object.
		
		"""
		
		d=pyfrp_gmsh_geometry.domain()
		d.addCylinderByParameters([self.center[0],self.center[1]],self.radius,self.height,self.width,volSizePx,
			   plane="z",genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol)
		
		if minID==None:
			minID=self.embryo.geometry.getMaxGeoID()+1
		
		d.incrementAllIDs(minID)
		
		return d	
	
	def writeToGeoFile(self,fn=None,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .geo file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.geo`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.genGmshDomain`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
				
		Returns:
			str: Path to geo file.
		
		"""
		
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".geo"
			fn=fn.replace(" ","_")
			
			
		d=self.genGmshDomain(volSizePx=volSizePx,genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol,minID=minID)
		d.writeToFile(fn)
		
		return fn	
	
	def genAsOpenscad(self,allowInf=False):
		
		"""Generates ROI  as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Keyword Args:
			allowInf (bool): Allow infinity in bounds of z-direction.
		
		Returns:
			solid.solidpython.cylinder: Solid python object. 
		
		"""
		
			
		z=self.getOpenscadZExtend()
		
		zmin,zmax=min(z),max(z)
		openScadROI=solid.translate([self.center[0],self.center[1],zmin])(solid.cylinder(r=self.radius,h=abs(zmax-zmin)))
		
		return openScadROI
		
	#def plotIn3D(self,domain=None,ax=None):
	
	###NOTE need this function here!!!
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Square ROI class
				
class squareROI(ROI):
	
	def __init__(self,embryo,name,Id,offset,sidelength,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.sidelength=sidelength
		self.offset=offset
		
	def setSideLength(self,s):
		self.sidelength=s
		return self.sidelength

	def getSideLength(self):
		return self.sidelength
	
	def setOffset(self,c):
		self.offset=c
		return self.offset
	
	def getOffset(self):
		return self.offset
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getSquareIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getSquareIdxImg(self.offset,self.sidelength,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getSquareIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getSquareIdxMesh(self.sidelength,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		"""Shows ROI in a 2D plot.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linewidth (float): Linewidth of plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	

		"""
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
		
		patch = ptc.Rectangle(self.offset,self.sidelength,self.sidelength,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def centerOffset(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setOffset([self.embryo.dataResPx/2+0.5-self.sidelength/2,self.embryo.dataResPx/2+0.5-self.sidelength/2])
		else:
			return self.setOffset([self.embryo.dataResPx/2-self.sidelength/2,self.embryo.dataResPx/2-self.sidelength/2])
		
	def makeReducable(self,auto=False,debug=False):
		
		oldOffset=self.getOffset()
		self.centerOffset()
		
		if not auto:
			a=input("Change offset of region "+ self.name + " from " + str(oldOffset) + ' to ' + str(self.getOffset()) + ' ? [Y/N]')
			
			if a=='N':
				self.setOffset(oldOffset)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setOffset(oldCenter)
				return False
			return True
		return False
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideSquare`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideSquare(x,y,self.offset,self.sidelength)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelength]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelength]
		return self.xExtend, self.yExtend
	
	def getCorners(self):
		
		"""Returns corners of square in counter-clockwise order, starting with offset.
		
		Returns:
			list: List of 2D coordinates of corners.
		
		"""
		
		corner1=np.asarray(self.offset)
		corner2=np.asarray([self.offset[0]+self.sidelength,self.offset[1]])
		corner3=np.asarray([self.offset[0]+self.sidelength,self.offset[1]+self.sidelength])
		corner4=np.asarray([self.offset[0],self.offset[1]+self.sidelength])
		
		return [corner1,corner2,corner3,corner4]
	
	def getCenterOfMass(self):
		
		r"""Computes center of mass of ROI.
		
		The center of mass is computed by
		
		.. math:: c = \frac{1}{N} \sum\limits_{i=1}{N} x_i ,
		
		where :math:`c` is the center of mass, :math:`N` the number of corners and :math:`x_i` is the 
		coordinate of corner :math:`i` .
		
		Returns:
			numpy.ndarray: Center of mass.
		
		"""
		
		corners=self.getCorners()
		
		CoM=cornes[0]
		
		for i in range(1,len(corners)):
			CoM=CoM+corners[i]
		
		CoM=CoM/len(corners)
		
		return CoM
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		.. note:: Will grab extent of geometry to find bounds in z-direction. 
		
		Returns:
			solid.solidpython.cube: Solid python object. 
		
		"""
		
		try:
			ext=self.embryo.geometry.getZExtend()
		except AttributeError:
			printError("genAsOpenscad: Cannot grab extend from geometry of type " + self.embryo.geometry.typ)
			
		openScadROI=solid.translate([self.offset[0],self.offset[1],min(ext)])(solid.cube([self.sidelength,self.sidelength,abs(max(ext)-min(ext))]))
		
		return openScadROI
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Square and slice ROI class

class squareSliceROI(squareROI,sliceROI):
	
	def __init__(self,embryo,name,Id,offset,sidelength,height,width,sliceBottom,color='b'):	
		
		squareROI.__init__(self,embryo,name,Id,offset,sidelength,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getSquareIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getSquareIdxImg(self.offset,self.sidelength,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getSquareIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getSquareIdxMesh(self.sidelength,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideSquare`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideSquare(x,y,self.offset,self.sidelength)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelength]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelength]
		return self.xExtend, self.yExtend
	
	def genGmshDomain(self,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Translates ROI into gmsh domain object.
		
		This object can then be used to write ROIs to ``.geo`` files.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addCuboidByParameters`.
		
		.. note:: If ``minID==None``, will grab maximum ID via :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getMaxGeoID` and add 1.
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Domain object.
		
		"""
		
		d=pyfrp_gmsh_geometry.domain()
		d.addCuboidByParameters([self.offset[0],self.offset[1],self.height],self.sidelength,self.sidelength,self.width,volSizePx,
			   plane="z",genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol)
		
		if minID==None:
			minID=self.embryo.geometry.getMaxGeoID()+1
		
		d.incrementAllIDs(minID)
		
		return d	
	
	def writeToGeoFile(self,fn=None,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .geo file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.geo`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.genGmshDomain`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
			
		Returns:
			str: Path to geo file.
		
		"""
		
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".geo"
			fn=fn.replace(" ","_")
			
		
		d=self.genGmshDomain(volSizePx=volSizePx,genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol,minID=minID)
		d.writeToFile(fn)
			
		return fn
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Returns:
			solid.solidpython.cube: Solid python object. 
		
		"""
		
		z=self.getOpenscadZExtend()
		zmin,zmax=min(z),max(z)
		openScadROI=solid.translate([self.offset[0],self.offset[1],zmin])(solid.cube([self.sidelength,self.sidelength,abs(zmax-zmin)]))
		return openScadROI
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Rectangle ROI class

class rectangleROI(ROI):
	
	def __init__(self,embryo,name,Id,offset,sidelengthX,sidelengthY,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.sidelengthX=sidelengthX
		self.sidelengthY=sidelengthY
		self.offset=offset
		
	def setSideLengthX(self,s):
		self.sidelengthX=s
		return self.sidelengthX

	def getSideLengthX(self):
		return self.sidelengthX
	
	def setSideLengthY(self,s):
		self.sidelengthY=s
		return self.sidelengthY

	def getSideLengthY(self):
		return self.sidelengthY
	
	def setOffset(self,c):
		self.offset=c
		return self.offset
	
	def getOffset(self):
		return self.offset
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getRectangleIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getRectangleIdxImg(self.offset,self.sidelengthX,self.sidelengthY,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getRectangleIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getRectangleIdxMesh(self.sidelengthX,self.sidelengthY,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		"""Shows ROI in a 2D plot.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linewidth (float): Linewidth of plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	

		"""
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		patch = ptc.Rectangle(self.offset,self.sidelengthX,self.sidelengthY,fill=False,linewidth=linewidth,color=color)
		ax.add_patch(patch)
		return ax
	
	def centerOffset(self):
		if np.mod(self.embryo.dataResPx,2)==0:
			return self.setOffset([self.embryo.dataResPx/2+0.5-self.sidelengthX/2,self.embryo.dataResPx/2+0.5-self.sidelengthY/2])
		else:
			return self.setOffset([self.embryo.dataResPx/2-self.sidelengthX/2,self.embryo.dataResPx/2-self.sidelengthY/2])
		
	def makeReducable(self,atuo=False,debug=False):
		
		oldOffset=self.getOffset()
		self.centerOffset()
		
		if not auto:
			a=input("Change offset of region "+ self.name + " from " + str(oldOffset) + ' to ' + str(self.getOffset()) + ' ? [Y/N]')
			
			if a=='N':
				self.setOffset(oldOffset)
				return False
			elif a=='Y':
				pass
			
			if not self.checkSymmetry():
				printWarning('Cannot make region '+self.name+' reducable.')
				self.setOffset(oldCenter)
				return False
			return True
		return False
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideRectangle`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideRectangle(x,y,self.offset,self.sidelengthX,self.sidelengthY)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelengthX]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelengthY]
		return self.xExtend, self.yExtend
	
	def getCorners(self):
		
		"""Returns corners of rectangle in counter-clockwise order, starting with offset.
		
		Returns:
			list: List of 2D coordinates of corners.
		
		"""
		
		corner1=np.asarray(self.offset)
		corner2=np.asarray([self.offset[0]+self.sidelengthX,self.offset[1]])
		corner3=np.asarray([self.offset[0]+self.sidelengthX,self.offset[1]+self.sidelengthY])
		corner4=np.asarray([self.offset[0],self.offset[1]+self.sidelengthY])
		
		return [corner1,corner2,corner3,corner4]
	
	def getCenterOfMass(self):
		
		r"""Computes center of mass of ROI.
		
		The center of mass is computed by
		
		.. math:: c = \frac{1}{N} \sum\limits_{i=1}{N} x_i ,
		
		where :math:`c` is the center of mass, :math:`N` the number of corners and :math:`x_i` is the 
		coordinate of corner :math:`i` .
		
		Returns:
			numpy.ndarray: Center of mass.
		
		"""
		
		corners=self.getCorners()
		
		CoM=cornes[0]
		
		for i in range(1,len(corners)):
			CoM=CoM+corners[i]
		
		CoM=CoM/len(corners)
		
		return CoM
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		.. note:: Will grab extent of geometry to find bounds in z-direction. 
		
		Returns:
			solid.solidpython.cube: Solid python object. 
		
		"""
		
		try:
			ext=self.embryo.geometry.getZExtend()
		except AttributeError:
			printError("genAsOpenscad: Cannot greab extend from geometry of type " + self.embryo.geometry.typ)
			
		openScadROI=solid.translate([self.offset[0],self.offset[1],min(ext)])(solid.cube([self.sidelengthX,self.sidelengthY,abs(max(ext)-min(ext))]))
		return openScadROI
	
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Rectangle and slice ROI class

class rectangleSliceROI(rectangleROI,sliceROI):
	
	def __init__(self,embryo,name,Id,offset,sidelengthX,sidelengthY,height,width,sliceBottom,color='b'):
		
		rectangleROI.__init__(self,embryo,name,Id,offset,sidelengthX,sidelengthY,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getRectangleIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getRectangleIdxImg(self.sidelengthX,self.sidelengthY,self.offset,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getRectangleIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getRectangleIdxMesh(self.sidelengthX,self.sidelengthY,self.offset,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsideRectangle`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsideRectangle(x,y,self.offset,self.sidelengthX,self.sidelengthY)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend=[self.offset[0],self.offset[0]+self.sidelengthX]
		self.yExtend=[self.offset[1],self.offset[1]+self.sidelengthY]
		return self.xExtend, self.yExtend
	
	def genGmshDomain(self,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Translates ROI into gmsh domain object.
		
		This object can then be used to write ROIs to ``.geo`` files.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addCuboidByParameters`.
		
		.. note:: If ``minID==None``, will grab maximum ID via :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getMaxGeoID` and add 1.
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Domain object.
		
		"""
		
		d=pyfrp_gmsh_geometry.domain()
		d.addCuboidByParameters([self.offset[0],self.offset[1],self.height],self.sidelengthX,self.sidelengthY,self.width,volSizePx,
			   plane="z",genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol)
		
		if minID==None:
			minID=self.embryo.geometry.getMaxGeoID()+1
		
		d.incrementAllIDs(minID)
		
		return d	
	
	def writeToGeoFile(self,fn=None,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .geo file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.geo`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.genGmshDomain`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
		
		Returns:
			str: Path to geo file.
		
		"""
		
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".geo"
			fn=fn.replace(" ","_")
			
		d=self.genGmshDomain(volSizePx=volSizePx,genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol,minID=minID)
		d.writeToFile(fn)
		
		return fn	
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Returns:
			solid.solidpython.cube: Solid python object. 
		
		"""
		
		
		z=self.getOpenscadZExtend()
		zmin,zmax=min(z),max(z)
		openScadROI=solid.translate([self.offset[0],self.offset[1],zmin])(solid.cube([self.sidelengthX,self.sidelengthY,abs(zmax-zmin)]))
		return openScadROI
	
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon ROI class
			
class polyROI(ROI):
	
	def __init__(self,embryo,name,Id,corners,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.corners=corners
	
	def setCorners(self,corners):
		self.corners=corners
		return corners
	
	def getCorners(self):
		return self.corners
	
	def addCorner(self,c,pos=-1):
		self.corners.insert(pos,c)
		return self.corners
	
	def appendCorner(self,c):
		self.corners.append(c)
		return self.corners
	
	def removeCorner(self,pos):
		self.corners.pop(pos)
		return self.corners
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getPolyIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getPolyIdxImg(self.corners,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
			
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getPolyIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getPolyIdxMesh(self.corners,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx
	
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		"""Shows ROI in a 2D plot.
		
		If no color is specified, will use color specified in ``ROI.color``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linewidth (float): Linewidth of plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	

		"""
		
		if color==None:
			color=self.color
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
		patch = ptc.Polygon(self.corners,closed=True,fill=False,linewidth=linewidth,color=color)
		
		ax.add_patch(patch)
		return ax
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsidePoly`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsidePoly(x,y,self.corners)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		cornersNP=np.array(self.corners)
	
		xmax=cornersNP[:,0].max()
		xmin=cornersNP[:,0].min()
		ymax=cornersNP[:,1].max()
		ymin=cornersNP[:,1].min()
		
		self.xExtend=[xmin,xmax]
		self.yExtend=[ymin,ymax]
		return self.xExtend, self.yExtend
	
	def getCenterOfMass(self):
		
		r"""Computes center of mass of ROI.
		
		The center of mass is computed by
		
		.. math:: c = \frac{1}{N} \sum\limits_{i=1}{N} x_i ,
		
		where :math:`c` is the center of mass, :math:`N` the number of corners and :math:`x_i` is the 
		coordinate of corner :math:`i` .
		
		Returns:
			numpy.ndarray: Center of mass.
		
		"""
		
		CoM=np.asarray(self.cornes[0])
		
		for i in range(1,len(self.corners)):
			CoM=CoM+np.asarray(self.corners[i])
		
		CoM=CoM/len(corners)
		
		return CoM
	
	def moveCorner(self,idx,x,y):
		
		"""Moves corner to new postion.
		
		Args:
			idx (int): Index of corner to be moved.
			x (float): New x-coordinate.
			y (float): New y-coordinate.
			
		Results:
			list: Updated corners list.
			
		"""
		
		if idx==-1:
			idx=len(self.corners)-1
		
		for i,corner in enumerate(self.corners):
			if i==idx:
				
				self.corners[idx]=[x,y]
				
		return self.corners		
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Returns:
			solid.solidpython.linear_extrude: Solid python object. 
		
		"""
		
		try:
			ext=self.embryo.geometry.getZExtend()
		except AttributeError:
			printError("genAsOpenscad: Cannot greab extend from geometry of type " + self.embryo.geometry.typ)
		
		poly=solid.polygon(self.corners)
		extruded=solid.linear_extrude(height = abs(max(ext)-min(ext)), center = False)(poly)
		openScadROI=solid.translate([0,0,min(ext)])(extruded)
		return openScadROI
	
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon and slice ROI class
			
class polySliceROI(polyROI,sliceROI):
	
	def __init__(self,embryo,name,Id,corners,height,width,sliceBottom,color='b'):	
		
		polyROI.__init__(self,embryo,name,Id,corners,color=color)
		sliceROI.__init__(self,embryo,name,Id,height,width,sliceBottom,color=color)
		
		self.corners=corners
	
	def computeImgIdx(self,debug=False):
		
		"""Computes image indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getPolyIdxImg`.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		Returns:
			tuple: Tuple containing:
			
				* imgIdxX (list): Image indices in x-direction.
				* imgIdxY (list): Image indices in y-direction.
				
		"""
		
		[self.imgIdxX,self.imgIdxY]=pyfrp_idx_module.getPolyIdxImg(self.corners,self.embryo.dataResPx,debug=debug)
		return self.imgIdxX,self.imgIdxY
	
	def computeMeshIdx(self,mesh):
		
		"""Computes mesh indices of ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.getPolyIdxMesh`.
		
		Args:
			mesh (fipy.GmshImporter3D): Fipy mesh object.
			
		Returns:
			list: Newly computed mesh indices.
				
		"""
		
		self.meshIdx=pyfrp_idx_module.getPolyIdxMesh(self.corners,mesh,zmin=self.zmin,zmax=self.zmax)
		return self.meshIdx	
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		See also :py:func:`pyfrp.modules.pyfrp_idx_module.checkInsidePoly`.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		return pyfrp_idx_module.checkInsidePoly(x,y,self.corners)
	
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		cornersNP=np.array(corners)
	
		xmax=cornersNP[:,0].max()
		xmin=cornersNP[:,0].min()
		ymax=cornersNP[:,1].max()
		ymin=cornersNP[:,1].min()
		
		self.xExtend=[xmin,xmax]
		self.yExtend=[ymin,ymax]
		
		return self.xExtend, self.yExtend
	
	def genGmshDomain(self,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Translates ROI into gmsh domain object.
		
		This object can then be used to write ROIs to ``.geo`` files.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addPrismByParameters`.
		
		.. note:: If ``minID==None``, will grab maximum ID via :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getMaxGeoID` and add 1.
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Domain object.
		
		"""
		
		d=pyfrp_gmsh_geometry.domain()
		d.addPrismByParameters(self.corners,volSizePx,z=self.zmin,height=self.zmax,plane="z",genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol)
		
		if minID==None:
			minID=self.embryo.geometry.getMaxGeoID()+1
		
		d.incrementAllIDs(minID)
		
		return d
	
	def writeToGeoFile(self,fn=None,volSizePx=20.,genLoops=True,genSurfaces=True,genVol=True,minID=None):
		
		"""Writes ROI to geo file.
		
		.. note:: If ``fn`` is not given, will save .geo file of ROI in same folder as the geometry file of the embryo with the following path:
		   ``path/to/embryos/geo/file/nameOfEmbryo_nameOfROI.geo`` .
		   
		See also :py:func:`pyfrp.subclasses.pyfrp_ROI.polySliceROI.genGmshDomain`.   
		
		Keyword Args:
			volSizePx (float): Mesh size of vertices.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genVol (bool): Generate surface loop and corresponding volume.
			minID (int): Id at which geo IDs should start.
		
		Returns:
			str: Path to geo file.
		
		"""
		
		
		if fn==None:
			folder=os.path.dirname(self.embryo.geometry.fnGeo)
			fn=pyfrp_misc_module.slashToFn(folder)+self.embryo.name+"_"+self.name+".geo"
			fn=fn.replace(" ","_")
			
		d=self.genGmshDomain(volSizePx=volSizePx,genLoops=genLoops,genSurfaces=genSurfaces,genVol=genVol,minID=minID)
		d.writeToFile(fn)
		
		return fn	
			
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Returns:
			solid.solidpython.linear_extrude: Solid python object. 
		
		"""
		
		z=self.getOpenscadZExtend()
		zmin,zmax=min(z),max(z)
		poly=solid.polygon(self.corners)
		extruded=solid.linear_extrude(height = abs(zmax-zmin), center = False)(poly)
		openScadROI=solid.translate([0,0,zmin])(extruded)
		return openScadROI
		
		
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Polygon and slice ROI class
	
class customROI(ROI):
	
	def __init__(self,embryo,name,Id,color='b'):	
		
		ROI.__init__(self,embryo,name,Id,color=color)
		
		self.ROIsIncluded=[]
		self.procedures=[]
		
	def addROI(self,r,p):
		if r not in self.ROIsIncluded:
			self.ROIsIncluded.append(r)
			self.procedures.append(p)
		return self.ROIsIncluded
	
	def removeROI(self,r):
		if r in self.ROIsIncluded:
			i=self.ROIsIncluded.index(r)
			self.ROIsIncluded.remove(r)
			self.procedures.pop(i)
			
		return self.ROIsIncluded
	
	def mergeROIs(self,r):
		
		if len(self.ROIsIncluded)==0:
			self.copyIdxs(r)
		else:
			self.computeImgMask()
			self.imgMask=self.imgMask*r.computeImgMask()
			
			self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)
			
			self.meshIdx=pyfrp_misc_module.matchVals(self.meshIdx,r.meshIdx)
				
			self.addROI(r,1)
			
			self.extImgIdxX,self.extImgIdxY = pyfrp_idx_module.getCommonExtendedPixels(self.ROIsIncluded,self.embryo.dataResPx,procedures=self.procedures,debug=False)
			
		return self.ROIsIncluded
		
	def substractROIs(self,r):
		
		if len(self.ROIsIncluded)==0:
			self.copyIdxs(r)
		else:
			self.computeImgMask()
			self.imgMask=self.imgMask*(1-r.computeImgMask())
			self.imgIdxX,self.imgIdxY=pyfrp_idx_module.mask2ind(self.imgMask,self.embryo.dataResPx)
			
			self.meshIdx=pyfrp_misc_module.complValsSimple(self.meshIdx,r.meshIdx)
				
			self.addROI(r,-1)
			
			self.extImgIdxX,self.extImgIdxY = pyfrp_idx_module.getCommonExtendedPixels(self.ROIsIncluded,self.embryo.dataResPx,procedures=self.procedures,debug=False)
			
		return self.ROIsIncluded
	
	def getROIsIncluded(self):
		return self.ROIsIncluded
		
	def setROIsIncluded(self,l):
		self.ROIsIncluded=l
		return self.ROIsIncluded
	
	def updateIdxs(self):
			
		self.emptyIdxs()

		for i in range(len(self.ROIsIncluded)):
			if i==0:
				self.copyIdxs(self.ROIsIncluded[i])
			else:
				if self.procedures[i]==1:
					self.mergeROIs(self.ROIsIncluded[i])
				elif self.procedures[i]==-1:
					self.substractROIs(self.ROIsIncluded[i])
				else:
					printWarning("Unknown Procedure" + str(self.procedures[i]) + " in Custom ROI " + self.name +". Not going to do anything.")
		
		self.computeNumExt()
		
		return self.getAllIdxs()
		
	def showBoundary(self,color=None,linewidth=3,ax=None):
		
		"""Shows ROI in a 2D plot by plotting all included ROIs.
		
		If no color is specified, will use color specified in ``ROI.color``. If ``color=="each"``,
		will plot each included ROI in its respective color.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes used for plotting. If not specified, will generate new one.
			color (str): Color of plot.
			linewidth (float): Linewidth of plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.	

		"""
		
		if color==None:
			color=self.color
		
		if color=='each':
			color=None
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["boundary"],sup=self.name+" boundary")
			ax = axes[0]
			
			img=np.nan*np.ones((self.embryo.dataResPx,self.embryo.dataResPx))
			ax.imshow(img)
			
		for r in self.ROIsIncluded:
			if hasattr(r,'showBoundary'):
				r.showBoundary(color=color,ax=ax,linewidth=linewidth)
		
		return ax
	
	def checkXYInside(self,x,y):
		
		"""Checks if coordinates are inside ROI.
		
		Does this by looping through all ROIs specified in ``ROIsIncluded``
		and checking if x/y is supposed to lie inside or outside of 
		the respective ROI.
		
		Args:
			x (np.ndarray): Array of x-coordinates.
			y (np.ndarray): Array of y-coordinates.
			
		Returns:
			np.ndarray: Array of booleans with corresponding to [x,y].
		
		"""
		
		b=True
		for i,r in enumerate(self.ROIsIncluded):
			if self.procedures[i]==1:
				b=b and r.checkXYInside(x,y)
			elif self.procedures[i]==-1:
				b=b and not r.checkXYInside(x,y)
		return b
			
	def computeXYExtend(self):
		
		"""Computes extend of ROI in x/y direction.
		
		Returns:
			tuple: Tuple containing:
				
				* xExtend (list): List containing minimum/maximum x-coordinate (``[xmin,xmax]``).
				* yExtend (list): List containing minimum/maximum y-coordinate (``[ymin,ymax]``).
		
		"""
		
		self.xExtend,self.yExtend=pyfrp_idx_module.getCommonXYExtend(self.ROIsIncluded)
		return self.xExtend,self.yExtend
	
	def roiIncluded(self,r):
		
		"""Returns if a ROI is included in customROI.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): A ROI.
			
		Returns:
			bool: ``True`` if included, ``False`` else.
		
		"""
		
		return r in self.ROIsIncluded
	
	def genAsOpenscad(self):
		
		"""Generates ROI as solid python object.
		
		Useful if ROI is used to be passed to openscad.
		
		Returns:
			solid.solidpython.openscad_object: Solid python object. 
		
		"""
		
		for i,r in enumerate(self.ROIsIncluded):
			
			if i==0:
				openScadROI=r.genAsOpenscad()
			else:
				if self.procedures[i]==-1:
					openScadROI=openScadROI-r.genAsOpenscad()
				elif self.procedures[i]==1:
					openScadROI=openScadROI+r.genAsOpenscad()
				
		return openScadROI
		
		
