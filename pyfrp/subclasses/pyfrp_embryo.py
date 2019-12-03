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

"""Essential PyFRAP module containing :py:class:`embryo` class. 
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module 
from pyfrp.modules import pyfrp_img_module 
from pyfrp.modules import pyfrp_fit_module 
from pyfrp.modules import pyfrp_stats_module 
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules.pyfrp_term_module import *

#PyFRAP Objects
from . import pyfrp_geometry
from . import pyfrp_simulation
from . import pyfrp_analysis
from . import pyfrp_ROI
from . import pyfrp_fit

#matplotlib
import matplotlib.pyplot as plt

#Time 
import time

#Misc
import os


#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================


class embryo:
	
	"""Main PyFRAP class, gathering all the data and parameters of FRAP experiment.
	
	The ``embryo`` class basically stores:
	
		* A minimum set of basic FRAP parameters.
		* A list of :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` classes, describing all ROIs used for 
		  for evaluating simulation and analysis results.
		* A :py:class:`pyfrp.subclasses.pyfrp_geometry.geometry` class describing the 3-dimensional 
		  geometry of the experiment.
		* A :py:class:`pyfrp.subclasses.pyfrp_analysis.analysis` class describing how the dataset is 
		  going to be analyzed.
		* A :py:class:`pyfrp.subclasses.pyfrp_simulation.simulation` class describing how the dataset is 
		  going to be simulated.
		* A list of :py:class:`pyfrp.subclasses.pyfrp_fit.fit` classes, storing different fitting options
		  and results.
	
	The ``embryo`` class comes with a comprehensive set of methods aimed at making it as powerful as possible,
	while still keeping it simple. The hierarchical structure should make it easy to navigate through a FRAP
	dataset.
	
	"""
	
	#Creates new embryo object
	def __init__(self,name):
		
		#Dataset name
		self.name = name
		
		#Data Image specifics
		self.dataEnc='uint16'
		self.dataFT='tif'
		self.dataResPx=512
		self.dataResMu=322.34
		
		#DataSet files
		self.fileList=[]
		self.fnDatafolder=""
		
		#Data time specifics
		self.frameInterval=10
		self.nFrames=300
		self.tStart=0
		self.tEnd=self.tStart+self.frameInterval*(self.nFrames-1)
		self.tvecData=np.linspace(self.tStart,self.tEnd,self.nFrames)
		
		#Imaging specifics (In mu)
		self.sliceDepthMu=30.
		self.sliceWidthMu=2.5
		
		#Zoom Imaging (need to call this different)
		self.convFact=self.computeConvFact(updateDim=False)
		
		#Imaging specifics (In pixels)
		self.sliceWidthPx=self.sliceWidthMu/self.convFact
		self.sliceDepthPx=self.sliceDepthMu/self.convFact
		self.sliceHeightPx=-self.sliceDepthPx
		self.sliceBottom=0
		
		#Bleaching specifics (in Mu)
		self.sideLengthBleachedMu=2*79.54
		
		#Bleaching specifics (in pxs) (can take out)
		self.sideLengthBleachedPx=self.sideLengthBleachedMu/self.convFact
		self.offsetBleachedPx=[self.dataResPx/2-self.sideLengthBleachedPx/2, self.dataResPx/2-self.sideLengthBleachedPx/2]
		
		#List of ROIs
		self.ROIs=[]
		
		#Master ROI
		self.masterROIIdx=None
		
		#Geometry
		self.geometry=None
		
		#Simulation
		self.simulation=None
		
		#Analysis
		self.analysis=None
		
		#Fitting 
		self.fits=[]

	def addFit(self,fit):
		
		"""Appends fit object to list of fits
		
		Args:
			fit (pyfrp.subclasses.pyfrp_fit.fit): fit object.
			
		Returns:
			list: Updated ``fits`` list.
		
		"""
		
		self.fits.append(fit)
		return self.fits
	
	def newFit(self,name):
		
		"""Creates new fit object and appends it to list of fits.
		
		Args:
			name (str): Name of fit object.
			
		Returns:
			pyfrp.subclasses.pyfrp_fit.fit: Newly created fit.
		
		"""
		
		fit=pyfrp_fit.fit(self,name)
		self.addFit(fit)
		return fit
		
	def deleteFit(self,i):
		
		"""Deletes fit with index ``i`` from ``fits`` list.
		
		Args:
			i (int): Index of fit to be deleted.
			
		Returns:
			list: Updated ``fits`` list.
		
		"""
		
		self.fits.pop(i)
		return self.fits
	
		
	def save(self,fn=None,copyMeshFiles=True,debug=False):
		
		"""Saves embryo object to pickle file.
		
		If ``fn=None`` will save to ``self.name.emb``.
		
		Keyword Args:
			fn (str): Output filename.
			copyMeshFiles (bool): Copy meshfiles to embryo file destination.
			debug (bool): Print out debugging messages.
			
		Returns:
			str: Output filename.
		
		"""
		
		if fn==None:
			fn=self.name+".emb"
		
		if copyMeshFiles:
			
			print("before copyMeshFiles:", fn)
			print("fnGeo", self.geometry.fnGeo)
			
			fnGeo,fnMesh=pyfrp_IO_module.copyMeshFiles(fn,self.geometry.fnGeo,self.simulation.mesh.fnMesh,debug=debug)
			self.geometry.setFnGeo(fnGeo)
			self.simulation.mesh.setFnMesh(fnMesh)
			
			fnNew=os.path.splitext((os.path.split(fn)[-1]))[0]
			self.renameMeshFiles(fn=fnNew,debug=debug)
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)
		
		print("Saved "+  self.name+ " to " + fn)
		return fn
	
	def renameMeshFiles(self,fn=None,debug=False):
		
		"""Renames meshfiles associated with embryo fo ``fn``.
		
		If ``fn=None`` will rename to ``self.name.ext``, where ``ext`` is 
		either ``.geo`` or ``.msh``.
		
		Example: 
		
		>>> emb.geometry.getFnGeo()
		>>> path/to/meshfiles/dome.geo
		>>> emb.getName()
		>>> myEmbryo
		>>> emb.renameMeshFiles()
		>>> emb.geometry.getFnGeo()
		>>> path/to/meshfiles/myEmbryo.geo
		
		.. note:: Will automatically update ``geometry.fnGeo`` and 
		   ``simulation.mesh.fnMesh`` properties.
		
		Keyword Args:
			fn (str): Desired filename.
			debug (bool): Print out debugging messages.
		
		Returns:
			tuple: Tuple containing:
				
				* fnGeoNew (str): New path to geometry file.
				* fnMeshNew (str): New path to mesh file.
				
		"""
		
		if fn==None:
			fn=self.name+".emb"
		
		if self.geometry!=None:
			fnGeoNew=pyfrp_IO_module.copyAndRenameFile(self.geometry.fnGeo,fn,debug=debug)
			self.geometry.setFnGeo(fnGeoNew)
				
		if self.simulation.mesh.mesh!=None:
			fnMeshNew=pyfrp_IO_module.copyAndRenameFile(self.simulation.mesh.fnMesh,fn,debug=debug)
			self.simulation.mesh.setFnMesh(fnMeshNew)
			
		return self.geometry.fnGeo, self.simulation.mesh.fnMesh	
		
	
	def copy(self):
		
		"""Copies embryo, preserving all attributes and methods.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo.embryo: Embryo copy.
		
		"""
		
		copiedEmbryo=cpy.deepcopy(self)
		return copiedEmbryo
	
	def updateVersion(self):
		
		"""Updates embryo object to current version, making sure that it possesses
		all attributes.
		
		Creates a new embryo object and compares ``self`` with the new embryo object.
		If the new embryo object has a attribute that ``self`` does not have, will
		add attribute with default value from the new embryo object.
		
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo.embryo: ``self``
			
		"""
		
		embtemp=embryo("temp")
		pyfrp_misc_module.updateObj(embtemp,self)
		return self
	
	def computeConvFact(self,updateDim=True):
		
		"""Computes conversion factor between um to px (unit=um/px).
		
		If ``updateDimensions`` is selected, will update all dimensions of
		embryo object data rely on ``convFact``.
		
		Keyword Args:
			updateDim (bool): Automatically updated all convFact relevant attributes.
		
		Returns:
			float: New conversion factor.
			
		"""
		
		self.convFact = self.dataResMu/self.dataResPx
		if updateDim:
			self.updatePxDimensions()
		return self.convFact
	
	def updatePxDimensions(self):
		
		"""Updates all all convFact relevant attributes.
		
		Returns:
			tuple: Tuple containing:
			
				* sliceWidthPx (float): Updated slice width in px.
				* sliceDepthPx (float): Updated slice depth in px.
				* sliceHeightPx (float): Updated slice height in px.
				* sideLengthBleachedPx (float): Updated side length of bleached square in px.
				* offsetBleachedPx (float): Updated offset if bleached square in px.
				
		"""
		
		
		self.sliceWidthPx=self.sliceWidthMu/self.convFact
		self.sliceDepthPx=self.sliceDepthMu/self.convFact
		self.sliceHeightPx=-self.sliceDepthPx
			
		self.updateBleachedRegion()

		return self.sliceWidthPx,self.sliceDepthPx, self.sliceHeightPx, self.sideLengthBleachedPx, self.offsetBleachedPx
	
	def setName(self,n):
		
		"""Sets embryo name."""
		
		self.name=n
		return self.name
	
	def getName(self):
		
		"""Returns embryo name."""
		
		return self.name
	
	def setDataResMu(self,res):
		
		"""Sets resolution of data in :math:`\mu m`."""
		
		self.dataResMu=res 
		self.computeConvFact()
		return self.dataResMu
	
	def setDataResPx(self,res):
		
		"""Sets resolution of data in px. """
		
		self.dataResPx=int(res) 
		self.computeConvFact()
		return self.dataResPx
	
	def setDataFT(self,f):
		
		"""Sets data filetype of datasets, for example *.tif* ."""
		
		self.dataFT=f 
		return self.dataFT
	
	def setDataEnc(self,e):
		
		"""Sets data encoding of datasets, for example ``uint16`` ."""
		
		self.dataEnc=e 
		return self.dataEnc
	
	def getDataFT(self):
		
		"""Returns current data filetype."""
		
		return self.dataFT
	
	def getDataEnc(self):
		
		"""Returns current data encoding."""
		
		return self.dataEnc
	
	def setFrameInterval(self,dt):
		
		"""Sets interval between imaging frames, then updates all time vectors.
		
		Args: 
			dt (float): New frame interval in seconds.
		
		Returns: 
			float: Set frame interval.
		
		"""
		
		self.frameInterval=dt
		self.updateTimeDimensions()
		return self.frameInterval
	
	def setTStart(self,t):
		
		"""Sets start time of experiment, then updates all time vectors.
		
		Args: 
			t (float): Start time in seconds.
		
		Returns: 
			float: Set start time.
		
		"""
		
		self.tStart=t
		self.updateTimeDimensions()
		return self.tStart
	
	def setTEnd(self,t):
		
		"""Sets end time of experiment, then updates all time vectors.
		
		Args: 
			t (float): End time in seconds.
		
		Returns: 
			float: Set end time.
		
		"""
		
		self.tEnd=t
		self.updateTimeDimensions()
		return self.tEnd
	
	def updateNFrames(self):
		
		"""Updates number of frames, then updates all time vectors.
		
		.. note:: Gets number of frames by reading number of files of type ``dataFT`` in ``fnDatafolder``. 
		   So you should make sure that there are no extra files in ``fnDatafolder``.
		
		
		Returns: 
			int: Updated number of frames.
		
		"""
		
		self.nFrames=len(self.fileList)
		self.updateTimeDimensions()
		return self.nFrames
	
	def getNFrames(self):
		
		"""Returns current number of frames."""
		
		return self.nFrames
	
	def setNFrames(self,n):
		
		"""Sets number of frames, then updates all time vectors. """
		
		self.nFrames=n
		self.updateTimeDimensions()
		return self.nFrames
	
	def getDataResPx(self):
		
		"""Returns resolution of data in px."""
		
		return self.dataResPx
	
	def getDataResMu(self):
		
		"""Returns resolution of data in :math:`\\um m` ."""
		
		return self.dataResMu
	
	def getFrameInterval(self):
		
		"""Returns current frame interval."""
		
		return self.frameInterval
	
	def getTStart(self):
		
		"""Returns current exerpiment start time."""
		
		return self.tStart
	
	def getTEnd(self):
		
		"""Returns current exerpiment end time."""
		
		return self.tEnd
	
	def getTvecData(self):
		
		"""Returns current data time vector."""
		
		return self.tvecData
	
	def updateTimeDimensions(self,simUpdate=True):
		
		"""Updates time dimensions using information in ``frameInterval``, ``nFrames`` 
		and ``tStart``.
		
		.. note:: If embryo object already possesses simulation object, will update simulation
		   time dimensions using :py:func:`pyfrp.subclasses.pyfrp_simulation.simulation.toDefaultTvec` .
		   If you have different settings in simulation vector that you want to preserve, select
		   ``simUpdate=False``.
		   
		Keyword Args:
			simUpdate (bool): Update simulation time dimensions.
		
		Returns:
			numpy.ndarray: Updated data time vector.
		
		"""
		
		self.tEnd=self.tStart+self.frameInterval*(self.nFrames-1)
		self.tvecData=np.linspace(self.tStart,self.tEnd,self.nFrames)
		if self.simulation!=None:
			if simUpdate:
				self.simulation.toDefaultTvec()
		return self.tvecData
	
	def setSliceDepthMu(self,d,updateGeometry=True):
		
		"""Sets resolution of data in :math:`\mu m` and then updates all dimensions of
		embryo object data rely on ``sliceDepthMu``.
		
		If ``updateGeometry`` is selected, will automatically update geometry.
		
		.. note:: Not every geometry depends on slice depth. If the geometry object has
		   a method ``restoreDefault``, this will be called. Otherwise geometry is not updated.
		
		Args:
			d (float): New slice depth in :math:`\mu m`.
			
		Keyword Args:
			updateGeometry (bool): Update geometry.
			
		Returns:
			float: New slice depth.
			
		"""
		
		self.sliceDepthMu=d
		self.updatePxDimensions()
		
		if updateGeometry:
			if hasattr(self.geometry,"restoreDefault"):
				self.geometry.restoreDefault()
			
		
		return self.sliceDepthMu
	
	def setMasterROIIdx(self,idx):
		
		"""Define ROI with index ``idx`` as master ROI.
		"""
		
		self.masterROIIdx=idx
		return self.masterROIIdx
	
	def setSideLengthBleachedMu(self,s):
		
		"""Sets sidelength of bleached square in :math:`\mu m`, will then 
		update bleached region parameters by calling :py:func:`updateBleachedRegion`.
		
		.. warning:: Attributes ``sideLengthBleachedMu`` and ``offsetBleachedPx`` will 
		   deprecated in further versions. Bleached region definition will then solely rely 
		   on ROI definitions. 
		
		"""
		
		self.sideLengthBleachedMu=s
		self.updateBleachedRegion()
		return self.sideLengthBleachedMu
		
	def updateBleachedRegion(self):
		
		r"""Updates sidelength and offset of bleached square in px.
		
		Bleached region is defined around center of image. That is, the offset is set to
		   
		.. math:: x_{\mathrm{offset}} = \frac{1}{2}(r_{\mathrm{px}},r_{\mathrm{px}})^T - (s_{\mathrm{bleached,px}},s_{\mathrm{bleached,px}})^T  
		   
		where :math:`r_{\mathrm{px}}` is the resolution of the data in pixel, and :math:`s_{\mathrm{bleached,px}` is
		the sidelength of the bleached square in pixel.
		
		.. warning:: Attributes ``sideLengthBleachedMu`` and ``offsetBleachedPx`` will 
		   deprecated in further versions. Bleached region definition will then solely rely 
		   on ROI definitions. 
		
		"""
		
		self.sideLengthBleachedPx=self.sideLengthBleachedMu/self.convFact
		self.offsetBleachedPx=[self.dataResPx/2-self.sideLengthBleachedPx/2, self.dataResPx/2-self.sideLengthBleachedPx/2]
		return self.sideLengthBleachedPx, self.offsetBleachedPx
		
	def getMasterROIIdx(self):
		
		"""Returns index of master ROI"""
		
		return self.masterROIIdx
	
	def getROIIdx(self,r):
		
		"""Returns index of ROI in ``ROIs`` list.
		
		Returns ``-1`` if ROI is not in list.
		
		Args:
			r (pyfrp.subclasses.pyfrp_ROI.ROI): Some ROI.
			
		Returns:
			int: Index of ROI.
			
		"""
		
		try:
			return self.ROIs.index(r)
		except ValueError:
			printError("ROI "+ r.name+ " is not in ROIs list of embryo "+self.name+".")
			return -1
	
	def addROI(self,roi):
		
		"""Adds ROI to ``ROIs`` list.
		
		Args:
			roi (pyfrp.subclasses.pyfrp_ROI.ROI): A PyFRAP ROI.
			
		Returns:
			list: Updated ROIs list.
		
		"""
	
		self.ROIs.append(roi)
		
		return self.ROIs
	
	def newROI(self,name,Id,zmin='-inf',zmax='inf',color='b',asMaster=False):
		
		"""Creates new simple :py:class:`pyfrp.subclasses.pyfrp_ROI.ROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			
		Keyword Args:
			zmin (float): Lower boundary in z-direction.
			zmax (float): upper boundary in z-direction.
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.ROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.ROI(self,name,Id,zmin=zmin.inf,zmax=zmax,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
		
	def newRadialROI(self,name,Id,center,radius,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.radialROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			center (list): Center of radial ROI.
			radius (float): Radius of radial ROI.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.radialROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.radialROI(self,name,Id,center,radius,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
		
	def newSliceROI(self,name,Id,height,width,sliceBottom,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.sliceROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		.. note:: If ``sliceBottom==True``, slice ranges from ``height<=z<=height+width``, otherwise
		   from ``height-width/2<=z<=height+width/2``.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			height (float): z-coordinate of slice.
			width (float): width in z-direction of slice.
			sliceBottom (bool): Put origin of slice at botton of slice.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.sliceROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.sliceROI(self,name,Id,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRadialSliceROI(self,name,Id,center,radius,height,width,sliceBottom,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		.. note:: If ``sliceBottom==True``, slice ranges from ``height<=z<=height+width``, otherwise
		   from ``height-width/2<=z<=height+width/2``.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			center (list): Center of radial ROI.
			radius (float): Radius of radial ROI.
			height (float): z-coordinate of slice.
			width (float): width in z-direction of slice.
			sliceBottom (bool): Put origin of slice at botton of slice.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.radialSliceROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.radialSliceROI(self,name,Id,center,radius,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newSquareROI(self,name,Id,offset,sidelength,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.squareROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: Offset is set to be left-bottom corner of square.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			offset (list): Offset of of square.
			sidelength (float): Sidelength of square.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.squareROI: Newly created ROI object.
		
		"""
		
		
		roi=pyfrp_ROI.squareROI(self,name,Id,offset,sidelength)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newSquareSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.squareSliceROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: Offset is set to be left-bottom corner of square.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		.. note:: If ``sliceBottom==True``, slice ranges from ``height<=z<=height+width``, otherwise
		   from ``height-width/2<=z<=height+width/2``.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			center (list): Center of radial ROI.
			offset (list): Offset of of square.
			sidelength (float): Sidelength of square.
			width (float): width in z-direction of slice.
			sliceBottom (bool): Put origin of slice at botton of slice.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.squareSliceROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.squareSliceROI(self,name,Id,offset,sidelength,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRectangleROI(self,name,Id,offset,sidelengthX,sidelengthY,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.rectangleROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: Offset is set to be left-bottom corner of rectangle.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			offset (list): Offset of of square.
			sidelengthX (float): Sidelength of rectangle in x-direction.
			sidelengthY (float): Sidelength of rectangle in y-direction.
				
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.rectangleROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.rectangleROI(self,name,Id,offset,sidelengthX,sidelengthY,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newRectangleSliceROI(self,name,Id,offset,sidelengthX,sidelengthY,height,width,sliceBottom,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.rectangleSliceROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: Offset is set to be left-bottom corner of rectangle.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		.. note:: If ``sliceBottom==True``, slice ranges from ``height<=z<=height+width``, otherwise
		   from ``height-width/2<=z<=height+width/2``.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			center (list): Center of radial ROI.
			offset (list): Offset of of rectangle.
			sidelengthX (float): Sidelength of rectangle in x-direction.
			sidelengthY (float): Sidelength of rectangle in y-direction.
			width (float): width in z-direction of slice.
			sliceBottom (bool): Put origin of slice at botton of slice.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.rectangleSliceROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.rectangleSliceROI(self,name,Id,offset,sidelengthX,sidelengthY,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newPolyROI(self,name,Id,corners,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.polyROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		Each corner given in ``corners`` list must be a list of form ``[x,y]``.
		
		.. note:: Polygon is automatically closed.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			corners (list): List of ``[x,y]`` coordinates describing corners.
				
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.polyROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.polyROI(self,name,Id,corners,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newPolySliceROI(self,name,Id,corners,height,width,sliceBottom,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.polySliceROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		Each corner given in ``corners`` list must be a list of form ``[x,y]``.
		
		.. note:: Polygon is automatically closed.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		.. note:: If ``sliceBottom==True``, slice ranges from ``height<=z<=height+width``, otherwise
		   from ``height-width/2<=z<=height+width/2``.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			center (list): Center of radial ROI.
			corners (list): List of ``[x,y]`` coordinates describing corners.
			width (float): width in z-direction of slice.
			sliceBottom (bool): Put origin of slice at botton of slice.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.polySliceROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.polySliceROI(self,name,Id,corners,height,width,sliceBottom,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def newCustomROI(self,name,Id,color='b',asMaster=False):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI` object
		and adds it to the ``ROIs`` list of embryo object.
		
		.. note:: You can use :py:meth:`getFreeROIId` to find a unused ID for the newROI.
		
		Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			
		Keyword Args:
			color (str): Color of ROI.
			asMaster (bool): Set ROI as master ROI?
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.customROI: Newly created ROI object.
		
		"""
		
		roi=pyfrp_ROI.customROI(self,name,Id,color=color)
		self.ROIs.append(roi)
		if asMaster:
			self.masterROIIdx=self.ROIs.index(roi)
		return roi
	
	def getFreeROIId(self):
		
		"""Returns first free ID of ROIs.
		"""
		
		IdList=pyfrp_misc_module.objAttrToList(self.ROIs,'Id')
		notFound=True
		i=0
		while notFound:
			if i not in IdList:
				return i
			i=i+1
			
	def removeROI(self,i):
		
		"""Removes ROI of index ``i``.
		"""
		
		self.ROIs.pop(i)
		if i==self.masterROIIdx:
			printWarning('You are deleting the masterROI. You need to select a new one before you can analyze or simulate.')
		return self.ROIs
	
	def listROIs(self):
		
		"""Prints out all ROIs and their respective type."""
		
		for r in self.ROIs:
			print(r.name , type(r))
		return True
	
	def genDefaultROIs(self,center,radius,rimFactor=0.66,masterROI=None,bleachedROI=None,rimROI=None,sliceHeightPx=None,clean=True):
		
		"""Creates a standard set of ROI objects and adds them to ``ROIs`` list.
		
		The set of ROIs covers the generally most useful ROIs used for FRAP experiments, providing already all the settings
		that PyFRAP uses for rim computation etc.
		
		ROIs contain:
		
			* All: ROI covering all pixels and mesh nodes. :py:class:`pyfrp.subclasses.pyfrp_ROI.sliceROI`
			* All Square: ROI covering all pixels and mesh nodes inside bleached region. :py:class:`pyfrp.subclasses.pyfrp_ROI.squareROI`
			* All Out: ROI covering all pixels and mesh nodes outside bleached region. :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI`
			* Slice: ROI covering all pixels and mesh nodes inside recorded field of view. :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
			* Slice rim: ROI covering all pixels that are not used for rim concentration computation. :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
			* Rim: ROI covering all pixels that are used for rim concentration computation. :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI`
			* Bleached Square: ROI covering all pixels and mesh nodes inside bleached region and imaging slice. :py:class:`pyfrp.subclasses.pyfrp_ROI.squareSliceROI`
			* Out: ROI covering all pixels and mesh nodes outside bleached region but inside imaging slice. :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI`
			
		.. note:: Will automatically set ``Slice`` as ``masterROI``.
		
		.. note:: If ``masterROI`` is given, will use this ROI instead of *Slice* as master ROI via :py:func:`setMasterROIIdx`. Will not create *Slice* at all. 
		   If ``masterROI`` is not in ``ROIs`` list yet, it will be automatically added to the list. 
		
		.. note:: If ``bleachedROI`` is given, will use this ROI instead of *Bleached Square*. Will generate copy of ``bleachedROI`` to generate *All Square* type 
		   ROI.
		
		.. note:: If ``rimROI`` is given, *Slice rim* is not created, and ``rimROI`` is used instead. 
		
		.. note:: If ``sliceHeightPx`` is given, will create ROIs *Slice*, *Out*, *Bleached Square* and *Rim* at this height, otherwise will use value stored in ``embryo.sliceHeightPx``.
		
		Args:	
			center (list): Center of circle defining imaging slice.
			radius (float): Radius of circle defining imaging slice.
			
		Keyword Args:	
			rimFactor (float): Factor describing percentage of imaging slice excluded from rim computation.
			sliceHeightPx (float): Height of slice ROI in px.
			clean (bool): Will remove all ROIs from embryo object before creating new ones.
			masterROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is supposed to be used as a masterROI.
			bleachedROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is supposed to be used to indicate the bleached region.
			rimROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is substracted from *Slice*. Should lie withing *Slice*.
			
		Returns:
			list: Updated list of ROIs.
		
		
		"""
		
		#Clean up if necessary
		if clean:
			self.ROIs=[]
			
		#Add masterROI and bleachedROI to ROIs list if they are given (necessary such that getFreeROIId works)
		if masterROI!=None and masterROI not in self.ROIs:
			self.ROIs.append(masterROI)
		if bleachedROI!=None and bleachedROI not in self.ROIs:
			self.ROIs.append(bleachedROI)
		if rimROI!=None and rimROI not in self.ROIs:
			self.ROIs.append(rimROI)
		
		#Check if sliceHeight is given
		if sliceHeightPx==None:
			sliceHeightPx=self.sliceHeightPx
		
		#Create All ROI
		allsl=self.newSliceROI("All",self.getFreeROIId(),0,np.inf,False,color=(0.5,0.3,0.4))
		
		#Create All Square
		if bleachedROI!=None:
			zmin,zmax=allsl.getZExtend()
			allsqu=bleachedROI.getCopy() 
			allsqu.setName("All "+bleachedROI.getName())
			allsqu.setZExtend(zmin,zmax)
			allsqu.setId(self.getFreeROIId())
			allsqu.setColor((1.,140/255.,0.))
			self.ROIs.append(allsqu)
		else:	
			allsqu=self.newSquareROI('All Square',self.getFreeROIId(),self.offsetBleachedPx,self.sideLengthBleachedPx,color=(1.,140/255.,0.))
		
		#Create All Out
		allout=self.newCustomROI("All Out",self.getFreeROIId(),color=(0.3,0.4,0.5))
		allout.addROI(allsl,1)
		allout.addROI(allsqu,-1)
		
		#Create Slice
		if masterROI!=None:
			sl=masterROI
			self.setMasterROIIdx(self.getROIIdx(sl))
			if sl.getName()!="Slice":
				printWarning("masterROI is not named Slice, this might lead to problems.")
				
		else:
			sl=self.newRadialSliceROI("Slice",self.getFreeROIId(),center,radius,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='g',asMaster=True)
		
		#Create Rim
		if rimROI!=None:
			sl2=rimROI
		else:
			sl2=self.newRadialSliceROI("Slice rim",self.getFreeROIId(),center,rimFactor*radius,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='y')
		 
		rim=self.newCustomROI("Rim",self.getFreeROIId(),color='m')
		rim.addROI(sl,1)
		rim.addROI(sl2,-1)
		rim.setUseForRim(True)
		
		#Create Bleached Square
		if bleachedROI!=None:
			squ=bleachedROI
			if squ.getName()!="Bleached Square":
				printWarning("Bleached ROI is not named Bleached Square, this might lead to problems.")
		else:	
			squ=self.newSquareSliceROI("Bleached Square",self.getFreeROIId(),self.offsetBleachedPx,self.sideLengthBleachedPx,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='b')
		
		#Create Out
		out=self.newCustomROI("Out",self.getFreeROIId(),color='r')
		out.addROI(sl,1)
		out.addROI(squ,-1)
		
		return self.ROIs
	
	def genMinimalROIs(self,center,radius,rimFactor=0.66,masterROI=None,bleachedROI=None,rimROI=None,sliceHeightPx=None,clean=True):
		
		"""Creates a minimal standard set of ROI objects and adds them to ``ROIs`` list.
		
		The set of ROIs covers the generally most useful ROIs used for FRAP experiments, providing already all the settings
		that PyFRAP uses for rim computation etc.
		
		ROIs contain:
		
			* Slice: ROI covering all pixels and mesh nodes inside recorded field of view. :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
			* Slice rim: ROI covering all pixels that are not used for rim concentration computation. :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
			* Rim: ROI covering all pixels that are used for rim concentration computation. :py:class:`pyfrp.subclasses.pyfrp_ROI.customROI`
			* Bleached Square: ROI covering all pixels and mesh nodes inside bleached region and imaging slice. :py:class:`pyfrp.subclasses.pyfrp_ROI.squareSliceROI`
			
		.. note:: Will automatically set ``Slice`` as ``masterROI``.
		
		.. note:: If ``masterROI`` is given, will use this ROI instead of *Slice* as master ROI via :py:func:`setMasterROIIdx`. Will not create *Slice* at all. 
		   If ``masterROI`` is not in ``ROIs`` list yet, it will be automatically added to the list. 
		
		.. note:: If ``bleachedROI`` is given, will use this ROI instead of *Bleached Square*. Will generate copy of ``bleachedROI`` to generate *All Square* type 
		   ROI.
		
		.. note:: If ``rimROI`` is given, *Slice rim* is not created, and ``rimROI`` is used instead. 
		
		.. note:: If ``sliceHeightPx`` is given, will create ROIs *Slice*, *Out*, *Bleached Square* and *Rim* at this height, otherwise will use value stored in ``embryo.sliceHeightPx``.
		
		Args:	
			center (list): Center of circle defining imaging slice.
			radius (float): Radius of circle defining imaging slice.
			
		Keyword Args:	
			rimFactor (float): Factor describing percentage of imaging slice excluded from rim computation.
			sliceHeightPx (float): Height of slice ROI in px.
			clean (bool): Will remove all ROIs from embryo object before creating new ones.
			masterROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is supposed to be used as a masterROI.
			bleachedROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is supposed to be used to indicate the bleached region.
			rimROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI that is substracted from *Slice*. Should lie withing *Slice*.
			
		Returns:
			list: Updated list of ROIs.
		

		"""
		
		#Clean up if necessary
		if clean:
			self.ROIs=[]
			
		#Add masterROI and bleachedROI to ROIs list if they are given (necessary such that getFreeROIId works)
		if masterROI!=None and masterROI not in self.ROIs:
			self.ROIs.append(masterROI)
		if bleachedROI!=None and bleachedROI not in self.ROIs:
			self.ROIs.append(bleachedROI)
		if rimROI!=None and rimROI not in self.ROIs:
			self.ROIs.append(rimROI)
		
		#Check if sliceHeight is given
		if sliceHeightPx==None:
			sliceHeightPx=self.sliceHeightPx
					
		#Create Slice
		if masterROI!=None:
			sl=masterROI
			self.setMasterROIIdx(self.getROIIdx(sl))
			if sl.getName()!="Slice":
				printWarning("masterROI is not named Slice, this might lead to problems.")
				
		else:
			sl=self.newRadialSliceROI("Slice",self.getFreeROIId(),center,radius,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='g',asMaster=True)
		
		#Create Rim
		if rimROI!=None:
			sl2=rimROI
		else:
			sl2=self.newRadialSliceROI("Slice rim",self.getFreeROIId(),center,rimFactor*radius,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='y')
		 
		rim=self.newCustomROI("Rim",self.getFreeROIId(),color='m')
		rim.addROI(sl,1)
		rim.addROI(sl2,-1)
		rim.setUseForRim(True)
		
		#Create Bleached Square
		if bleachedROI!=None:
			squ=bleachedROI
			if squ.getName()!="Bleached Square":
				printWarning("Bleached ROI is not named Bleached Square, this might lead to problems.")
		else:	
			squ=self.newSquareSliceROI("Bleached Square",self.getFreeROIId(),self.offsetBleachedPx,self.sideLengthBleachedPx,sliceHeightPx,self.sliceWidthPx,self.sliceBottom,color='b')
		
		return self.ROIs
	
	def getROIs(self):
		
		"""Returns ``ROIs`` list."""
		
		return self.ROIs
	
	def getROIByName(self,name):
		
		"""Returns ROI in ``ROIs`` list with specified name.
		
		If ROI with ``name`` does not exist, returns ``None``.
		"""
		
		for r in self.ROIs:
			if r.name==name:
				return r
		
		printWarning('Cannot find ROI with name '+ name + '. Will return None. This can lead to further problems.')
		
		return None

	def getROIById(self,Id):
		
		"""Returns ROI in ``ROIs`` list with specified ID.
		
		If ROI with ``Id`` does not exist, returns ``None``.
		"""
		
		for r in self.ROIs:
			if r.Id==Id:
				return r
			
		
		printWarning('Cannot find ROI with Id '+ Id + '. Will return None. This can lead to further problems.')
		
		return None
	
	def updateFileList(self):
		
		"""Updates file list containing all names of recovery images.
		
		If new ``fileList`` is not empty, will update number of frames
		to match number of files in ``fileList`` by calling :py:meth:`updateNFrames`.
		
		Returns:
			list: Updated file list.
		
		"""
		
		self.fileList=pyfrp_misc_module.getSortedFileList(self.fnDatafolder,self.dataFT)
		if len(self.fileList)==0:
			printWarning("There are no files of type " + self.dataFT + " in " + self.fnDatafolder +" . This can lead to problems.")
		else:
			self.updateNFrames()
		return self.fileList
	
	def setFileList(self,l):
		
		"""Sets file list to ``l``."""
		
		self.fileList=l
		return self.fileList
	
	def getFileList(self):
		
		"""Returns list of recovery data files."""
		
		return self.fileList
	
	def setDataFolder(self,fn):
		
		"""Set folder containing recovery data files.
		
		Will automatically try to update ``fileList`` by calling
		:py:meth:`updateFileList`.
		
		Args:
			fn (str): Path to folder containing data files.
		
		Returns:
			str: New data folder path.
		
		"""
		
		fn=fn.replace('//','/')
		self.fnDatafolder=fn
		self.updateFileList()
		return self.fnDatafolder
	
	def getDataFolder(self):
		
		"""Returns folder containing recovery data files."""
		
		return self.fnDatafolder
	
	def getGeometry(self):
		
		"""Returns embryo's geometry object.
		
		Returns:
			pyfrp.subclasses.pyfrp_geometry.geometry: PyFRAP geometry object.
		"""
		
		return self.geometry
	
	def setGeometry2ZebraFishDomeStage(self,center,imagingRadius,radiusScale=1.1):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.zebrafishDomeStage`.
		
		Args:
			center (list): Center of geometry.
			imagingRadius (float): Radius of embryo at imaging slice.
			
		Keyword Args:
			radiusScale (float): Scaling factor defining how much bigger outer radius is to inner radius.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.zebrafishDomeStage: New zebrafish geometry.
		"""
		
		self.geometry=pyfrp_geometry.zebrafishDomeStage(self,center,imagingRadius,radiusScale=radiusScale)
		return self.geometry
	
	def setGeometry2Cylinder(self,center,radius,height):
			
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.cylinder`.
		
		Args:
			center (list): Center of geometry.
			radius (float): Radius of cylinder.
			height (float): Height of cylinder.
				
		Returns:
			pyfrp.subclasses.pyfrp_geometry.cylinder: New cylinder geometry.
		"""
		
		
		self.geometry=pyfrp_geometry.cylinder(self,center,radius,height)
		return self.geometry
	
	def setGeometry2Cone(self,center,upperRadius,lowerRadius,height):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.cone`.
		
		Args:
			center (list): Center of geometry.
			upperRadius (float): Radius at upper end of cone.
			lowerRadius (float): Radius at lower end of cone.
			height (float): Height of cylinder.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.cone: New cone geometry.
		"""
		
		self.geometry=pyfrp_geometry.cone(self,center,upperRadius,lowerRadius,height)
		return self.geometry
	
	def setGeometry2Ball(self,center,imagingRadius):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.ball`.
		
		Args:
			center (list): Center of geometry.
			imagingRadius (float): Radius of embryo in imaging slice.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.ball: New ball geometry.
		"""
		
		self.geometry=pyfrp_geometry.xenopusBall(self,center,imagingRadius)
		return self.geometry
	
	def setGeometry2ZebraFishDomeStageQuad(self,center,imagingRadius,radiusScale=1.1):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.zebrafishDomeStageQuad`.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Args:
			center (list): Center of geometry.
			imagingRadius (float): Radius of embryo at imaging slice.
			
		Keyword Args:
			radiusScale (float): Scaling factor defining how much bigger outer radius is to inner radius.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.zebrafishDomeStageQuad: New zebrafish quadrant geometry.
		"""
		
		
		self.geometry=pyfrp_geometry.zebrafishDomeStageQuad(self,center,imagingRadius,radiusScale=radiusScale)
		return self.geometry
	
	def setGeometry2CylinderQuad(self,center,radius,height):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.cylinderQuad`.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Args:
			center (list): Center of geometry.
			radius (float): Radius of cylinder.
			height (float): Height of cylinder.
				
		Returns:
			pyfrp.subclasses.pyfrp_geometry.cylinder: New cylinder quadrant geometry.
		"""
		
		self.geometry=pyfrp_geometry.cylinderQuad(self,center,radius,height)
		return self.geometry
	
	def setGeometry2BallQuad(self,center,imagingRadius):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.ballQuad`.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Args:
			center (list): Center of geometry.
			imagingRadius (float): Radius of embryo in imaging slice.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.ball: New ball quadrant geometry.
		"""
		
		self.geometry=pyfrp_geometry.xenopusBallQuad(self,center,imagingRadius)
		return self.geometry
	
	def setGeometry2Custom(self,center,fnGeo="",dim=3):
		
		"""Sets embryo's geometry to :py:class:`pyfrp.subclasses.pyfrp_geometry.custom`.
		
		Args:
			center (list): Center of geometry.
			fnGeo (str): Path to geometry file.
			
		Returns:
			pyfrp.subclasses.pyfrp_geometry.custom: New custom geometry.
		"""
		
		self.geometry=pyfrp_geometry.custom(self,center,fnGeo,dim=dim)
		return self.geometry
	
	def newSimulation(self):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_simulation.simulation` object and sets it as
		embryos's ``simulation``.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation.simulation: New simulation object.
		"""
		
		self.simulation=pyfrp_simulation.simulation(self)
		return self.simulation
	
	def newAnalysis(self):
		
		"""Creates new :py:class:`pyfrp.subclasses.pyfrp_analysis.analysis` object and sets it as
		embryos's ``analysis``.
		
		Returns:
			pyfrp.subclasses.pyfrp_analysis.analysis: New analysis object.
		"""
		
		self.analysis=pyfrp_analysis.analysis(self)
		return self.analysis
	
	def computeROIIdxs(self,signal=None,debug=True):
		
		"""Computes image, extended and mesh indices of all ROIs in embryo's ``ROIs`` list.
		
		Keyword Args:
			signal (PyQt5.QtCore.pyqtSignal): PyQT signal to send progress to GUI.
			debug (bool): Print final debugging messages and show debugging plots.
		
		Returns:
			list: Updated list of ROIs.
		"""
		
		for i,r in enumerate(self.ROIs):
			startInit=time.clock()
			r.computeIdxs()
			
			if debug:
				print(r.name, time.clock()-startInit)
			if signal:
				signal.emit(int(100.*i)/float(len(self.ROIs)))
			
		return self.ROIs
	
	def geometry2Quad(self):
		
		"""Converts current geometry to quadrant reduced version if available,
		keeping essential parameters the same.
		
		If geometry is already reduced or there is no quadrant version of the geometry,
		will do nothing and return unchanged geometry.
		
		.. warning:: Quadrant reduction is still experimental.
		
		.. note:: Will set ``fnGeo`` back to default value in ``meshfiles`` folder.
	
		Returns:
			pyfrp.subclasses.pyfrp_geometry.geometry: Updated geometry.
		"""
		
		if self.geometry==None:
			printError("No geometry selected yet.")
			return self.geometry	
		
		if 'Quad' in self.geometry.typ:
			printError("Geometry already set to quad.")
			return self.geometry	
		
		if self.geometry.typ=='zebrafishDomeStage':
			self.setGeometry2BallQuad(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		elif self.geometry.typ=='cylinder':
			self.setGeometry2CylinderQuad(self,self.geometry.getCenter(),self.geomtry.getRadius())
		elif self.geometry.typ=='xenopusBall':
			self.setGeometry2BallQuad(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		else:
			printError("Geometry of type " + self.geometry.typ + " has no quadrant reduced version.") 
		
		
		return self.geometry	
			
	
	def geometry2Full(self):
		
		"""If current geometry was in quadrant reduced version, converts it to 
		full version, keeping essential parameters the same.
		
		.. warning:: Quadrant reduction is still experimental.
		
		.. note:: Will set ``fnGeo`` back to default value in ``meshfiles`` folder.
	
		Returns:
			pyfrp.subclasses.pyfrp_geometry.geometry: Updated geometry.
		"""
		
		if self.geometry==None:
			printError("No geometry selected yet.")
			return
		
		if 'Quad' not in self.geometry.typ:
			printError("Geometry already set to full.")
			return
		
		if self.geometry.typ=='zebrafishDomeStageQuad':
			self.setGeometry2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		elif self.geometry.typ=='cylinderQuad':
			self.setGeometry2Cylinder(self,self.geometry.getCenter(),self.geomtry.getRadius())
		elif self.geometry.typ=='xenopusBallQuad':
			self.setGeometry2Ball(self,self.geometry.getCenter(),self.geomtry.getImagingRadius())
		
		return self.geometry	
	
	def setEmbryo2Quad(self):
		
		"""Reduces both geometry and ROIs of embryo to quadrant, such that 
		embryo object is fully quadrant reduced.
		
		Returns:
			tuple: Tuple Containing:
			
				* ``self.geometry`` (pyfrp.subclasses.pyfrp_geometry.geometry): Updated geometry.
				* ``self.ROIs`` (list): Updated list of ROIs.
		"""
		
		self.geometry2Quad()
		self.ROIs2Quad()
		return self.geometry,self.ROIs
	
	def setEmbryo2Full(self):	
		
		"""Sets both geometry and ROIs to full mode.
		
		Returns:
			tuple: Tuple Containing:
			
				* ``self.geometry`` (pyfrp.subclasses.pyfrp_geometry.geometry): Updated geometry.
				* ``self.ROIs`` (list): Updated list of ROIs.
		"""
		
		self.geometry2Full()
		self.ROIs2Full()
		return self.geometry,self.ROIs
	
	def ROIs2Quad(self):
		
		"""Reduces ROIs to quadrant reduced mode, mapping all their indices in first quadrant.
		
		.. note:: Use :py:func:setEmbryo2Quad to make sure that both geometry and ROIs are reduced.
		
		.. warning:: Quadrant reduction is still experimental.
		
		Returns:
			list: Updated list of ROIs.
		"""
		
		for r in self.ROIs:
			r.idxs2Quad()
		return self.ROIs
	
	def ROIs2Full(self):
		
		"""Sets ROIs to full mode.

		Returns:
			list: Updated list of ROIs.
		"""
		
		for r in self.ROIs:
			r.idxs2Full()
		return self.ROIs
	
	def showAllROIBoundaries(self,ax=None,withImg=False,idx=0):
		
		"""Shows boundaries of all ROIs in ``ROIs`` list.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			withImg (bool): Shows data image inside same axes.
			idx (int): Index of data image to be shown.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIs:
			if hasattr(r,'showBoundary'):
				ax=r.showBoundary(ax=ax)
		
		if withImg:
			ax=self.showDataImg(ax=ax,idx=idx)
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def loadDataImg(self,idx):
		
		"""Loads data image in ``fnDatafolder`` of index ``idx``.
	
		Args:
			idx (int): Index of data image to be loaded.
			
		Returns:
			numpy.ndarray: Loaded image.
		"""
		
		return pyfrp_img_module.loadImg(self.fnDatafolder+self.fileList[idx],self.dataEnc)
	
	def showDataImg(self,ax=None,idx=0):
		
		"""Shows data image in ``fnDatafolder`` of index ``idx``.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			idx (int): Index of data image to be shown.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],sup=["Embryo" + self.name],titles=["DataImg "+str(idx)],tight=False)
			ax=axes[0]
			
		img=self.loadDataImg(idx)	
		ax.imshow(img)
		
		return ax
		
	def showAllROIIdxs(self,axes=None):
		
		"""Shows image, extended and mesh indices of all ROIs in ``ROIs`` list.
		
		If no axes are given via ``axes``, will create new list of matplotlib axes. 
		If axes are given, they should have ``len(axes)=3*len(ROIs)``.
		
		.. note:: If the mesh has a large number of nodes, this can take up a lot
		   of memory due to every node being plotted in a matplotlib scatter plot.
		
		Keyword Args:
			axes (list): List of matplotlib.axes to be plotted in.
			
		Returns:
			matplotlib.axes: List of axes used for plotting.
		
		"""
		
		
		if axes==None:
			proj=len(self.ROIs)*['3d']+len(self.ROIs)*[None]+len(self.ROIs)*[None]
			fig,axes = pyfrp_plot_module.makeSubplot([3,len(self.ROIs)],sup=["Embryo" + self.name + " ROI Indices"],tight=True,proj=proj)
			
		for i,r in enumerate(self.ROIs):
			print(r.name)
			currAxes=[axes[0+i],axes[len(self.ROIs)+i],axes[2*len(self.ROIs)+i]]
			r.showIdxs(axes=currAxes)
		
		return axes
			
	
	def plotAllData(self,ax=None,legend=True):
		
		"""Plots all data timeseries for all ROIs in ``ROIs`` list.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			legend (bool): Show legend in plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIs:
			ax=r.plotData(ax=ax,legend=legend)
	
		return ax	
	
	def plotAllSim(self,ax=None,legend=True):
		
		"""Plots all simulation timeseries for all ROIs in ``ROIs`` list.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			legend (bool): Show legend in plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIs:
			ax=r.plotSim(ax=ax)
		
		if legend:
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	def plotAllDataPinned(self,ax=None,legend=True):
	
		"""Plots all pinned data timeseries for all ROIs in ``ROIs`` list.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			legend (bool): Show legend in plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIs:
			ax=r.plotDataPinned(ax=ax,legend=legend)
		
		return ax	
	
	def plotAllSimPinned(self,ax=None,legend=True):
		
		"""Plots all pinned simulation timeseries for all ROIs in ``ROIs`` list.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to be plotted in.
			legend (bool): Show legend in plot.
			
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		for r in self.ROIs:
			ax=r.plotSimPinned(ax=ax,legend=legend)
		
		if legend:
			ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
		
		return ax	
	
	
	def checkQuadReducable(self,tryFix=False,auto=False,debug=False):
		
		"""Checks if embryo is reducable to quadrant by checking if all
		ROIs are either point or axis symmetric around geometry center.
		
		.. note:: You want to call :py:meth:showAllROIBoundaries and 
		   :py:meth:computeROIIdxs afterwards to make sure everything went
		   properly.
		   
		.. warning:: Quadrant reduction is still experimental.   
		   
		Keyword Args:
			tryFix (bool): Tries to readjust ROIs into reducable form.
			auto (bool): Readjust ROIs automatically.
			debug (bool): Print debugging messages.
			
		Returns:
			bool: True if embryo is reducable.
		
		"""
		
		reducable=True
		for r in self.ROIs:
			if not r.checkSymmetry(debug=debug):
				if tryFix:
					if hasattr(r,'makeReducable'):
						reducable=r.makeReducable(auto=auto,debug=debug)
				else:		
					printWarning('Cannot reduce region '+ self.name+' to quadrant. Indices are not symmetric.')
					reducable=False
				
		return reducable		
	
	def makeQuadReducable(self,auto=False,debug=False):
		
		"""Makes embryo quadrant reducable by:
		
			* Checks if embryo is reducable to quadrant by checking if all
			  ROIs are either point or axis symmetric around geometry center.
			* Tries to readjust non-symmetric ROIs to make them quad-reducable.
			* If ROIs are reducable, will center geometry at center of data image.
		
		.. note:: You want to call :py:meth:showAllROIBoundaries and 
		   :py:meth:computeROIIdxs afterwards to make sure everything went
		   properly.
		   
		.. warning:: Quadrant reduction is still experimental.   
		   
		Keyword Args:
			tryFix (bool): Tries to readjust ROIs into reducable form.
			auto (bool): Readjust ROIs automatically.
			debug (bool): Print debugging messages.
			
		Returns:
			bool: True if embryo is reducable.
		
		"""
		
		reducable=self.checkQuadReducable(tryFix=True,auto=auto,debug=debug)

		if reducable:
			self.geometry.centerInImg()
		return reducable
	
	def getMasterROI(self):
		
		"""Returns master ROI."""
		
		return self.ROIs[self.masterROIIdx]
	
	def getOptimalAllROI(self,name='All',makeNew=False):
		
		"""Readjusts ROI with name *All*  (if existent) to cover whole geometry."""
		
		self.geometry.setAllROI(name=name,makeNew=makeNew)
	
	def computePinVals(self,useMin=True,useMax=True,bkgdVal=None,debug=False):
		
		"""Compute overall pinning values over all ROIs.
		
		Keyword Args:
			useMin (bool): Use minimum value for background computation.
			useMax (bool): Use maximum value for norm value computation.
			bkgdVal (float): Use this background value instead of newly computing it.
			debug (bool): Print debugging messages.
	
		Returns:
			tuple: Tuple containing:
				
				* bkgdVal (float): Background value.
				* normVal (float): Norming value.
			
		"""
		
		if bkgdVal==None:
			bkgdVal=self.computeBkgd(useMin=useMin,debug=debug)
		normVal=self.computeNorm(bkgdVal,useMax=useMax,debug=debug)
		
		return bkgdVal,normVal
	
	def computeBkgd(self,useMin=False,fromTS='both',debug=False):
		
		"""Computes background value over all ROIs.
		
		If ``useMin==False``, will use value at first index of data/simulation
		vectors as norming value.
		
		.. note:: Use ``fromTS='both'`` to use values from both simulation and 
		   data for background computation. If ``fromTS='data'``, will only use 
		   data, of ``fromTS='sim'`` will only use simulation vectors.
		
		Keyword Args:
			useMin (bool): Use minimum value for background computation.
			fromTS (bool): Which time series to use for background computation.
			debug (bool): Print debugging messages.
		
		Returns:
			float: Background value. 
			
		"""
		
		bkgds=[]
		for r in self.ROIs:
			if fromTS in ['data','both']:
				if r.isAnalyzed():
					bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.dataVec,useMin=useMin,useMax=False,debug=debug)
					bkgds.append(bkgdTemp)
				else:
					if debug:
						printWarning("ROI " + r.name + " is not analyzed yet, cannot use it for background computation.")
						
			if fromTS in ['sim','both']:
				if r.isSimulated():
					bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.simVec,useMin=useMin,useMax=False,debug=debug)
					bkgds.append(bkgdTemp)
				else:
					if debug:
						printWarning("ROI " + r.name + " is not simulated yet, cannot use it for background computation.")
		
		if len(bkgds)==0:
			if debug:
				printWarning("No bkgd value could be found. Will return 0.")
			return 0.
		
		return min(bkgds)	
	
	def computeNorm(self,bkgdVal,useMax=True,fromTS='both',debug=False):
		
		"""Computes background value over all ROIs.
		
		If ``useMax==False``, will use value at last index of data/simulation
		vectors as norming value.
	
		.. note:: Use ``fromTS='both'`` to use values from both simulation and 
		   data for background computation. If ``fromTS='data'``, will only use 
		   data, of ``fromTS='sim'`` will only use simulation vectors.
		
		Args:
			bkgdVal (float): Use this background value instead of newly computing it.
		
		Keyword Args:
			useMax (bool): Use maximum value for norm value computation.
			fromTS (bool): Which time series to use for background computation.
			debug (bool): Print debugging messages.
		
		Returns:
			float: Norming value. 
			
		"""
		
		norms=[]
		for r in self.ROIs:
			if fromTS in ['data','both']:
				if r.isAnalyzed():
					bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.dataVec,bkgdVal=bkgdVal,useMin=False,useMax=useMax,debug=debug)
					norms.append(norm)
				else:
					if debug:
						printWarning("ROI " + r.name + " is not analyzed yet, cannot use it for norming computation.")
						
			if fromTS in ['sim','both']:
				if r.isSimulated():
					bkgdTemp,norm=pyfrp_fit_module.computePinVals(r.simVec,bkgdVal=bkgdVal,useMin=False,useMax=useMax,debug=debug)
					norms.append(norm)
				else:
					if debug:
						printWarning("ROI " + r.name + " is not simulated yet, cannot use it for norming computation.")
		if len(norms)==0:
			if debug:
				printWarning("No norm value could be found. Will return 1.")
			return 1.
			
		return max(norms)	
			
	def computeIdealFRAPPinVals(self,bkgdName='Bleached Square',normName='Slice',debug=False,useMin=False,useMax=False,sepSim=True,switchThresh=0.95):
		
		"""Computes background and norming value using optimized settings.
		
		Idea: Instead using values from all ROIs to compute pinning values, select 
		two ROIs that should lead to optimal pinning values. In the default case, this
		is the ROI describing the bleached region (here we expect the lowest intensities),
		and the slice (here we expect the overall end concentration the experiment is converging to).
		
		If ``useMin==False``, will use value at first index of data/simulation
		vectors as norming value.
			
		If ``useMax==False``, will use value at last index of data/simulation
		vectors as norming value.
	
		.. note:: If ``bkgdName`` and ``normName`` are not set differently, will look for ROIs
		   with these names for background and norming computation, respectively. If they don't exist,
		   will return ``None``. 
		   :py:meth:`genDefaultROIs` will make sure that both those ROIs exist.
		   
		.. warning:: Not all ROIs are suitable for pinning value computation. Generally ROIs that have the
		   least extended volume proof most suitable.
		
		.. note:: ``switchThresh`` checks if recovery is complete. This is important if recovery curves are 
		   very slow and bleached region does not reach full recovery. If this happens, we divide by a number 
		   that is smaller than 1 and hence boost all timeseries way above one instead of limiting it below one
		
		.. warning:: ``sepSim==True`` makes sure that we never get negative intensities in the pinned simulation 
		   vectors. Since interpolation is never perfect, this can happen if :math:`bkgdValue(data)>bkgdValue(sim)`.
		
		Keyword Args:
			bkgdName (str): Name of ROI used for background computation.
			normName (str): Name of ROI used for norming computation.
			useMin (bool): Use minimum value for background computation.
			useMax (bool): Use maximum value for norm value computation.
			fromTS (bool): Which time series to use for background computation.
			debug (bool): Print debugging messages.
			sepSim (bool): Use seperate pinning values for simulation vectors.
		
		Returns:
			tuple: Tuple containing:
				
				* bkgdVal (float): Background value for data vectors.
				* normVal (float): Norming value for data vectors.
				* bkgdValSim (float): Background value for simulation vectors.
				* normValSim (float): Norming value for simulation vectors.
		
		"""
		
		bkgdROI=self.getROIByName(bkgdName)
		normROI=self.getROIByName(normName)
		
		if bkgdROI==None or normROI==None:
			return None,None,None,None
		
		bkgdValSim=None
		normValSim=None
		
		bkgds=[]
		norms=[]
		
		bkgdBdata,normBdata=pyfrp_fit_module.computePinVals(bkgdROI.dataVec,useMin=useMin,useMax=useMax,debug=debug)
		bkgdBsim,normBsim=pyfrp_fit_module.computePinVals(bkgdROI.simVec,useMin=useMin,useMax=useMax,debug=debug)
		
		bkgdVal=min([bkgdBdata,bkgdBsim])
		
		if sepSim:
			bkgdValSim=bkgdBsim
			bkgdVal=bkgdBdata
		
		"""Need to check if we need to switch ROI. This is important if recovery curves are very slow and bleached region
		does not reach full recovery. If this happens, we divide by a number that is smaller than 1 and hence boost all 
		timeseries way above one instead of limiting it below one."""
		
		if normBdata<switchThresh or normBsim<switchThresh:
			bkgdNdata,normNdata=pyfrp_fit_module.computePinVals(normROI.dataVec,useMin=useMin,useMax=useMax,bkgdVal=bkgdVal,debug=debug)
			bkgdNsim,normNsim=pyfrp_fit_module.computePinVals(normROI.simVec,useMin=useMin,useMax=useMax,bkgdVal=bkgdVal,debug=debug)
			
			normVal=max([normNsim,normNdata])
			if sepSim:
				normValSim=normNsim
				normVal=normNdata
			
			print("in switch")
			
			if debug:
				printNote('Switched to region ' +  normName + ' for computation of normalization value.' )
		else:	
			
			print(normBdata+bkgdBdata, normBsim+bkgdBsim)
			
			print("not in switch")
			normVal=max([normBdata,normBsim])
			if sepSim:
				normValSim=normBsim
				normVal=normBdata
		
		return bkgdVal, normVal, bkgdValSim, normValSim
		
	def pinAllROIs(self,bkgdVal=None,normVal=None,bkgdValSim=None,normValSim=None,useMin=False,useMax=False,debug=False):
		
		"""Pins both simulation and data vectors of all ROIs.
		
		.. note:: If no bkgdVal or normVal is giving, will try to compute it using
		   :py:meth:`computePinVals`. Only then input of ``useMax`` and ``useMin``
		   are relevant.
		
		Keyword Args:
			bkgdVal (float): Background value used for data pinning.
			normVal (float): Norming value used for data pinning.
			bkgdValSim (float): Background value used for simulation pinning.
			normValSim (float): Norming value used for simulation pinning.
			useMin (bool): Use minimum value for background computation.
			useMax (bool): Use maximum value for norm value computation.
			debug (bool): Print debugging messages.
			
		Returns:
			list: Updated list of ROIs.
		
		"""
		
		bkgdValTemp,normValTemp = self.computePinVals(useMin=useMin,useMax=useMax,debug=debug)
		
		bkgdVal = pyfrp_misc_module.assignIfVal(bkgdVal,bkgdValTemp,None)
		normVal = pyfrp_misc_module.assignIfVal(normVal,normValTemp,None)
		
		for r in self.ROIs:
			r.pinAllTS(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=debug)
		return self.ROIs	
			
	def printAllAttr(self,full=False):
		
		"""Prints out all attributes of embryo object.""" 
		
		printAllObjAttr(self)
		
		if full:
			self.geometry.printDetails()
			self.analysis.printAllAttr()
			self.simulation.printAllAttr()
			self.simulation.mesh.printAllAttr()
				
		
	def isAnalyzed(self):
		
		"""Returns ``True`` if all ROIs have been analyzed."""
		
		b=True
		for r in self.ROIs:
			b=b and r.isAnalyzed()
		return b

	def isSimulated(self):
		
		"""Returns ``True`` if all ROIs have been simulated."""
		
		b=True
		for r in self.ROIs:
			b=b and r.isSimulated()
		return b
	
	def isFitted(self):
		
		"""Returns ``True`` if all fits are fitted."""
		
		b=True
		for fit in self.fits:
			b=b and fit.isFitted()
		return b	
			
	def getFitByName(self,name):
		
		"""Returns fit in ``fits`` list with name ``name``. 
		
		If it doesn't exists, returns ``None``.
		"""
		
		for fit in self.fits:
			if fit.name==name:
				return fit
			
		printWarning('Cannot find fit with name '+ name + '. Will return None. This can lead to further problems.')
		
		return None
	
	def getInterpolationError(self):
		
		"""Prints out interpolation error by ROI.
		"""
		
		for r in self.ROIs:
			print(r.name, r.getInterpolationError())
	
	def checkROIIdxs(self,debug=False):
		
		"""Checks if all ROIs have their mesh and image indices computed.
		
		Keyword Args:
			debug (bool): Print debugging messages.
		
		Returns:
			tuple: Tuple containing:
			
				* img (bool): True if all ROIs have up-to-date image indices.
				* mesh (bool): True if all ROIs have up-to-date mesh indices.
		"""
		
		img=True
		mesh=True
		for r in self.ROIs:
			if len(r.imgIdxX)==0:
				if debug:
					printWarning("ROI " + r.name+ " needs its img indices computed. Will not be able to perform image analysis otherwise.")
				img=False
			if len(r.meshIdx)==0:
				if debug:
					printWarning("ROI " + r.name+ " needs its mesh indices computed. Will not be able to perform simulation otherwise.")
				mesh=False
		return img,mesh
	
	def quickAnalysis(self,maxDExpPx=None,timeScale='log',runIdxs=True,runAnalysis=True,runSimulation=True,runPin=True,runFit=True,bkgdName='Bleached Square',normName='Slice',useMin=False,useMax=False,sepSim=True,switchThresh=0.95):
			
		"""Performs complete FRAP analysis of embryo object including:
		
			* Finding ROI indices by calling :py:meth:`computeROIIdxs` .
			* Running image analysis through :py:meth:`pyfrp.subclasses.pyfrp_analysis.analysis.run` .
			* Generating optimal time simulation vector through  :py:meth:`pyfrp.subclasses.pyfrp_simulation.simulation.getOptTvecSim`
			* Running simulation through :py:meth:`pyfrp.subclasses.pyfrp_simulation.simulation.run` .
			* Pin simulation and data timeseries using :py:meth:`computeIdealFRAPPinVals` and :py:meth:`pinAllROIs`.
			* Run all fits in ``fits`` list by calling :py:meth:`pyfrp.subclasses.pyfrp_fit.fit.run`
			
		Keyword Args:
			maxDExpPx (float): Maximum expected diffusion coefficient.
			timeScale (str): Linear (``'lin'``) or logarithmic (``'log'``) time scaling.
			runIdxs (bool): Run indexing.
			runAnalysis (bool): Run image analysis.
			runSimulation (bool): Run simulation.
			runPin (bool): Run pinning.
			runFit (bool): Run fitting.
			bkgdName (str): Name of ROI used for background computation.
			normName (str): Name of ROI used for norming computation.
			useMin (bool): Use minimum value for background computation.
			useMax (bool): Use maximum value for norm value computation.
			sepSim (bool): Use seperate pinning values for simulation vectors.
			
			
		"""
		
		if runIdxs:
			self.computeROIIdxs()
		
		#Image analysis
		if runAnalysis:
			self.analysis.run(showProgress=True)

		#Simulation
		if maxDExpPx!=None:
			self.simulation.getOptTvecSim(maxDExpPx)
		
		if timeScale=='log':
			self.simulation.toLogTimeScale()
		
		if runSimulation:
			self.simulation.run(showProgress=True)

		#Pin Concentrations
		if runPin:
			bkgdVal,normVal,bkgdValSim,normValSim=self.computeIdealFRAPPinVals(debug=True,useMin=False,useMax=False,switchThresh=0.95,bkgdName=bkgdName,normName=normName)
			self.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=False)

		#Run all fits
		if runFit:
			for fit in self.fits:
				fit.run(debug=False)
			
		
	def clearAllAttributes(self):
		
		"""Replaces all attribute values of embryo object with ``None``, except ``name``.
		
		Useful if embryos are seperated and molecule file needs to be compressed.
		
		Returns:
			bool: True if success, False else.
		
		"""
		
		try:
			for item in vars(self):
				if item!="name":
					vars(self)[str(item)]=None
			return True				
		except:
			printError("Failed to clearAllAttributes in embryo" + self.name +" .")
			return False
				
	def sliceEmbryo(self,nSlices,direction='z'):
		
		"""Slices embryo in z-direction in ``nSlices``.
		
		Creates ``nSlices`` new :py:class:`pyfrp.subclasses.pyfrp_ROI.sliceROI` ROIs
		and appends them to ``ROIs`` list.
		
		.. note:: Slice ROIs will be created with ``sliceBottom=False``.
		
		Args:
			nSlices (int): Number of slices embryo is supposed to get cut.
		
		Returns:
			list: Updates ``ROIs`` list.
			
		"""
		
		if direction=='z':
			
			#Get zExtend and height of geometry
			zmin,zmax = self.geometry.getZExtend()
			height=abs(zmax-zmin)
			
			sliceWidth = height/nSlices
			
			for i in range(nSlices):
				self.newSliceROI("Slice"+str(i),self.getFreeROIId(),zmin+(i+0.5)*sliceWidth,sliceWidth,False,color=[1,1/(i+2.),0.])
		
		elif direction=='x':
			xmin,xmax,ymin,ymax = self.geometry.getXYExtend()
			xwidth=abs(xmax-xmin)
			ywidth=abs(ymax-ymin)
			sliceWidth = xwidth/nSlices
			
			for i in range(nSlices):
				self.newRectangleROI("Slice"+str(i),self.getFreeROIId(),[xmin+i*sliceWidth,ymin],sliceWidth,ywidth,color=[1,1/(i+2.),0.])
				
		elif direction=='y':
			xmin,xmax,ymin,ymax = self.geometry.getXYExtend()
			xwidth=abs(xmax-xmin)
			ywidth=abs(ymax-ymin)
			sliceWidth = ywidth/nSlices
			
			for i in range(nSlices):
				self.newRectangleROI("Slice"+str(i),self.getFreeROIId(),[xmin,ymin+i*sliceWidth],xwidth,sliceWidth,color=[1,1/(i+2.),0.])
				
				
		return self.ROIs
	
	
	def compareFitsByAIC(self,ROIs=None,sigma=1,fromSSD=True,thresh=None,printOut=True):
		
		"""Compares all fits of embryo using the Akaike information criterion (AIC).
		
		For a detailed explanation of the model selection procedure, please refer to
		:py:func:`pyfrp.modules.pyfrp_stats_module.compareFitsByAIC`.
		
		If ``printOut`` is selected, will print a list of selected fits and the 
		underlying Akaike values in a readable table.
		
		If the AIC or the AICc should be used can be determined using 
		:py:func:`pyfrp.modules.pyfrp_stats_module.useAIC`.
		
		Keyword Args:	
			ROIs (list): List of ROIs to be considered for computation.
			sigma (float): Standard deviation of normal distribution if known.
			fromSSD (bool): Simply use SSD as maximum likelihood.
			thresh (float): Probability range for model selection.
			
		Returns:
			tuple: Tuple containing:
			
				* AICs (list): List of AIC values of the respective fits.
				* deltaAICs (numpy.ndarray): List of Akaike difference values of the respective fits.
				* weights (numpy.ndarray): List of Akaike difference weights of the respective fits.
				* acc (list): List of acceptable fits by model selection.
				* ks (list): List of number of parameters fitted of the respective fits.
				* ns (list): List of number of datapoints fitted of the respective fits.
		"""
		
		AICs, deltaAICs,weights,acc,ks,ns = pyfrp_stats_module.compareFitsByAIC(self.fits,ROIs=ROIs,sigma=sigma,fromSSD=fromSSD,thresh=thresh)
		
		
		if printOut:
			fitNames = pyfrp_misc_module.objAttrToList(self.fits,'name')
			printNote("Acceptable fits under threshold "+ str(thresh) +":")
			for fit in acc:
				print(fit.name)
			
			print("Details of AIC analysis can be found in the following table:")
			printTable([fitNames,AICs,deltaAICs,weights,ks,ns],["fitNames","AICs","deltaAICs","weights","k","n"],col=True)
			
		return AICs, deltaAICs,weights, acc,ks,ns
			
	def compareFitsByCorrAIC(self,ROIs=None,sigma=1,fromSSD=True,thresh=None,printOut=True):
		
		"""Compares all fits of embryo using the corrected Akaike information criterion (AICc).
		
		For a detailed explanation of the model selection procedure, please refer to
		:py:func:`pyfrp.modules.pyfrp_stats_module.compareFitsByCorrAIC`.
		
		If ``printOut`` is selected, will print a list of selected fits and the 
		underlying Akaike values in a readable table.
		
		If the AIC or the AICc should be used can be determined using 
		:py:func:`pyfrp.modules.pyfrp_stats_module.useAIC`.
		
		Keyword Args:	
			ROIs (list): List of ROIs to be considered for computation.
			sigma (float): Standard deviation of normal distribution if known.
			fromSSD (bool): Simply use SSD as maximum likelihood.
			thresh (float): Probability range for model selection.
			
		Returns:
			tuple: Tuple containing:
			
				* AICs (list): List of AICc values of the respective fits.
				* deltaAICs (numpy.ndarray): List of Akaike difference values of the respective fits.
				* weights (numpy.ndarray): List of Akaike difference weights of the respective fits.
				* acc (list): List of acceptable fits by model selection.
				* ks (list): List of number of parameters fitted of the respective fits.
				* ns (list): List of number of datapoints fitted of the respective fits.
		
		"""
		
		AICs, deltaAICs,weights,acc,ks,ns = pyfrp_stats_module.compareFitsByCorrAIC(self.fits,ROIs=ROIs,sigma=sigma,fromSSD=fromSSD,thresh=thresh)
		
		if printOut:
			fitNames = pyfrp_misc_module.objAttrToList(self.fits,'name')
			printNote("Acceptable fits under threshold "+ str(thresh) +":")
			for fit in acc:
				print(fit.name)
			
			print("Details of AIC analysis can be found in the following table:")
			printTable([fitNames,AICs,deltaAICs,weights,ks,ns],["fitNames","AICcs","deltaAICs","weights","k","n"],col=True)
		
		return AICs, deltaAICs,weights, acc,ks,ns	
	
	def getRimFactorByPx(self,radius,px):
		
		r"""Returns the correct rimFactor if a rim of with ``px`` in a 
		ROI of radius ``radius`` is desired.
		
		Rim factor :math:`f` is then given by:
		
		.. math:: f=1-\frac{p}{r},
		
		where :math:`p` is ``px`` and :math:`r` is ``radius``.
		
		Args:
			radius (float): Radius of ROI.
			px (float): Number of pixels that rim should be wide.
			
		Returns:
			float: Calculated rim factor.
		
		"""
		
		return 1-px/radius
	
	def fixFilePaths(self):
		
		"""Fixes paths to geometry/meshfiles.
		
		Returns:
			bool: True if success for all paths.
			
		
		"""
		
		b=True
		try:
			self.geomtry.fnGeo=pyfrp_misc_module.fixPath(self.geomtry.fnGeo)	
		except:
			b=False
		try:
			self.simulation.mesh.fnMesh=pyfrp_misc_module.fixPath(self.geomtry.fnMesh)	
		except:
			b=False
			
		return b
		
		
	
	
	#def grabDataDetails(self):
	###NOTE make a function that automatically grabs filetype and what not			
