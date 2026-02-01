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

"""PyFRAP module containing geometry classes. The :py:class:`geometry` class is 
a simple geometry class providing basic parameters and methods, parenting different
more specific geometries such as:

	* :py:class:`zebrafishDomeStage`: Describing a zebrafish embryo in dome stage.
	* :py:class:`cylinder`: A simple cylinder.
	* :py:class:`xenopusBall`: A simple ball geometry (looking like a xenopus embryo).
	* :py:class:`cone`: A (truncated) cone.

For most of the geometries, this module also provides quadrant reduced versions of the geometry,
reducing the geometry to the first quadrant around the center. 

.. note:: Rules for adding new geometries:
	
	   * Always subclass from  :py:class:`geometry`
	   * Always center geometry around ``geometry.center``. That includes having defining the center
	     in the .geo file by ``center_x`` and ``center_y``.
	   * Unit is pixels.
	   * Be careful with method overwrites. Use them wisely.
	   * Include the geometries into the GUI.
	   * Make them accessable by sharing them.

"""


#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_gmsh_module
from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules.pyfrp_term_module import *

#OS
import os,os.path
import shutil

#Solid/Opescad
import solid
import solid.utils

#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class geometry(object):
	
	"""Basic PyFRAP geometry class. 
	
	Stores all the necessary information to describe a geometry. Comes with helpful methods
	for 
	
		* Centering (see :py:func:`centerInImg`)
		* Defining optimal ROIs (see :py:func:`optimalAllROI`)
		* Geometry file parsing (see :py:func:`readGeoFile`)
		* Geometry plotting (see :py:func:`plotGeometry`)

	Args:
		embryo (pyfrp.subclasses.pyfrp_emrbyo.embryo): Embryo class that geometry belongs to.
		typ (str): Type of geometry.
		fnGeo (str): Path to gmsh .geo file describing the geometry.
		center (numpy.ndarray): Center of geometry.
		
	
	
	"""
	
	def __init__(self,embryo,typ,fnGeo,center,dim=3):
		
		self.typ=typ
		self.embryo=embryo
		self.fnGeo=pyfrp_misc_module.fixPath(fnGeo)
		self.center=center
		self.dim=dim
		
		self.geoFileParameters={}
	
	def getDim(self):
		
		"""Returns dimension of geometry.
		
		Returns:
			int: Dimension.
		"""
		
		return self.dim
			
	def getEmbryo(self):
		
		"""Returns :py:class:`pyfrp.subclasses.pyfrp_embryo.embryo` instance 
		that geometry belongs to."""
		
		return self.embryo
	
	def getTyp(self):
		
		"""Returns type of geometry.
		
		Returns:
			str: Type of geometry.
		"""
		
		
		return self.typ
	
	def setFnGeo(self,fn):
		
		"""Sets path to .geo file.
		
		Args:
			fn (str): Path to file.
		
		Returns:
			str: Path to file.
		"""
		
		self.fnGeo=pyfrp_misc_module.fixPath(fn)
		return self.fnGeo
	
	def getFnGeo(self):
		
		"""Returns path to .geo file.
		
		Returns:
			str: Path to file.
		"""
		
		return self.fnGeo
	
	def setCenter(self,c,updateInFile=True):
		
		"""Sets geometry center.
		
		.. note:: The geometry ``center`` attribute is then also set in
		   .geo file if ``updateInFile`` is selected.
		
		Args:
			c (numpy.ndarray): New center ``[x,y]``.
			
		Keyword Args:
			updateInFile (bool): Update center in .geo file.
			
		Returns:
			numpy.ndarray: New center.
		
		"""
		
		self.center=c
		
		if updateInFile:
			if hasattr(self,'updateGeoFile'):
				self.updateGeoFile()
		return self.center
	
	def getCenter(self):
		
		"""Returns geometry center.
		
		Returns:
			numpy.ndarray: New center.
		
		"""
		
		return self.center
	
	def center2Mid(self,updateInFile=True):
		
		"""Sets geometry center to center of image.
		
		Uses ``embryo.dataResPx`` to calculate image center.
		
		.. note:: The geometry ``center`` attribute is then also set in
		   .geo file if ``updateInFile`` is selected.
			
		Keyword Args:
			updateInFile (bool): Update center in .geo file.
			
		Returns:
			numpy.ndarray: New center.
		
		"""
		
		if mod(self.embryo.dataResPx,2)==0:
			return self.setCenter([self.embryo.dataResPx/2+0.5,self.embryo.dataResPx/2+0.5],updateInFile=updateInFile)
		else:
			return self.setCenter([self.embryo.dataResPx/2,self.embryo.dataResPx/2],updateInFile=updateInFile)
	
	def centerInImg(self):
		
		"""Sets geometry center to center of image, updates .geo file and 
		if avaialable remeshes.
		
		Uses ``embryo.dataResPx`` to calculate image center.
		
		.. note:: The geometry ``center`` attribute is then also set in
		   .geo file.
				
		Returns:
			numpy.ndarray: New center.
		
		"""
		
		oldCenter=self.getCenter()
		self.centerMid()
		a=input("Change center of geometry from " + oldCenter + ' to ' + self.getCenter() + ' ? [Y/N]')
		if a=='Y':
			self.updateGeoFile()
			if None!=self.embryo.simulation:
				self.embryo.simulation.mesh.genMesh()
		
		else:
			self.setCenter(oldCenter)
		
		return 	self.getCenter()
	
	def setAllROI(self,name='All',makeNew=False,updateIdxs=False):
		
		"""Tries to set the optimal *All* ROI for a specific geometry.
		
		Keyword Args:
			name (str): Name of All ROI to look for.
			makeNew (bool): Generate a new All ROI. 
			updateIdxs (bool): Update indices of new ROI.
			
		Returns:
			pyfrp.subclasses.pyfrp_ROI.ROI: New ROI instance.
		
		"""
		
		if not hasattr(self,'optimalAllROI'):
			printWarning('Cannot set AllROI. Geometry of type ' + self.typ + ' has no method for this!')
			return None
		
		if makeNew:
			rnew=self.optimalAllROI(name)
		else:	
			r=self.embryo.getROIByName(name)
			if r==None:
				printWarning('Cannot set AllROI. ROI with name ' + name + ' does not exist!')
				return None
			
			
			incl=r.findIncluded()
				
			name=r.getName()
			Id=r.getId()
			color=r.getColor()
			asMaster=r.isMaster()
			
			rnew=self.optimalAllROI(name=name,Id=Id,color=color,asMaster=asMaster)
			
			for rincl in incl:
				ind=rincl.ROIsIncluded.index(r)
				rincl.ROIsIncluded.pop(ind)
				rincl.ROIsIncluded.insert(ind,rnew)
				
				if updateIdxs:
					rincl.updateIdxs()
			
			ind=self.embryo.ROIs.index(r)
			self.embryo.removeROI(ind)
			self.embryo.ROIs.insert(ind,self.embryo.ROIs.pop(-1))
				
		return rnew
	
	def readGeoFile(self):
		
		"""Reads the .geo file and parses it into a :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.domain` 
		instance. 
		
		Returns:
			:py:class:`pyfrp.modules.pyfrp_gmsh_geometry.domain`: Domain containing geometry. 
		
		"""
		
		try:
			domain,self.geoFileParameters=pyfrp_gmsh_IO_module.readGeoFile(self.fnGeo)
		except:
			printWarning("Cannot read geo file " + self.fnGeo)
			return None
		
		return domain
	
	def plotGeometry(self,ax=None,color='k',ann=False,backend='mpl',asSphere=False,linewidth=1):
		
		"""Plots geometry in 3D.
		
		Reads the .geo file and parses it into a :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.domain` 
		instance. Then draws the domain.
		
		If no axes are given via ``ax``, will create new matplotlib axes.
		
		Keyword Args:
			ax (matplotlib.axes): Axes to draw in.
			color (str): Color of plot.
			ann (bool): Show annotations.
			backend (str): Backend used for drawing.
			asSphere (bool): Draw vertices as spheres.
			linewidth (float): Line width.
		
		Returns:
			matplotlib.axes: Axes used for plotting.
		
		"""
		
		domain=self.readGeoFile()
		
		if domain==None:
			return ax
		
		#if ax==None:
			#fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Geometry " +self.typ],proj=['3d'])
			#ax=axes[0]
		
		ax=domain.draw(ax=ax,color=color,ann=ann,backend=backend,asSphere=asSphere,linewidth=linewidth)
		
		return ax
	
	def getZExtend(self):
		
		"""Returns extend in z-direction by reading out vertices from 
		.geo file and returning maximum and minimum z-coordinates.
		
		Returns:
			tuple: Tuple containing:
				
				* zmin (float): Minimum z-coordinate.
				* zmax (float): Maximum z-coordinate.
				
		"""
		
		domain=self.readGeoFile()
		
		z=[]
		for v in domain.vertices:
			z.append(v.x[2])	
		
		return min(z), max(z) 
	
	def getXYExtend(self):
		
		"""Returns extend in x/y-direction by reading out vertices from 
		.geo file and returning maximum and minimum x/y-coordinates.
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
				
		"""
		
		domain=self.readGeoFile()
		
		x=[]
		y=[]
		for v in domain.vertices:
			x.append(v.x[0])	
			y.append(v.x[1])	
		
		return min(x), max(x) , min(y), max(y)
	
	def getExtend(self):
		
		"""Returns extend in x/y/z-direction. 
		
		Will call :py:func:`getXYExtend` and :py:func:`getZExtend` for it.
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
				* zmin (float): Minimum z-coordinate.
				* zmax (float): Maximum z-coordinate.
		
		"""
		
		return list(self.getXYExtend())+list(self.getZExtend())
		
	def printDetails(self):
		
		"""Prints out all details of geometry object.
		"""
		
		print("Geometry of embryo ", self.embryo.name, " Details.")
		printAllObjAttr(self)
		print() 
	
	def moveGeoFile(self,fn):
		
		"""Moves geometry file to different directory.
		
		.. note:: This function actually copies the file so that files in ``pyfrp/meshfiles/`` 
		   will not be removed.
		
		Will update ``geometry.fnGeo`` to the new file location.
		
		.. note:: If existent, will also copy the corresponding mesh file.
		
		Args:
			fn (str): Path of folder where geo file is supposed to go.
		
		Returns:
			str: New file location.
		
		"""
		
		shutil.copy(self.fnGeo,fn)
		
		if os.path.isfile(self.fnGeo.replace('.geo','.msh')):
			shutil.copy(self.fnGeo.replace('.geo','.msh'),fn)
			
			if hasattr(self.embryo,'simulation'):
				self.embryo.simulation.mesh.fnMesh=os.path.join(fn,os.path.basename(self.fnGeo).replace('.geo','.msh'))
				
		self.fnGeo=os.path.join(fn,os.path.basename(self.fnGeo))
	
		return self.fnGeo
	
	def getMaxGeoID(self):
		
		"""Returns maximum ID over all elements in .geo file.
		
		Sell also :py:func:`readGeoFile` and :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.getAllMaxID`.
		
		Returns:
			int: Maximum ID.
		
		"""
		
		domain=self.readGeoFile()
		
		return domain.getAllMaxID()
	
	def render2Openscad(self,fn=None,segments=48):
		
		"""Generates .scad file for the geometry.
		
		.. note:: If ``fn=None``, then will use the same filename and path as .geo file.
		
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		"""
		
		if fn==None:
			fn=self.fnGeo.replace(".geo",".scad")
		
		solid.scad_render_to_file(self.genAsOpenscad(), filepath=fn,file_header='$fn = %s;' % segments, include_orig_code=False)
		
	def render2Stl(self,fn=None,segments=48):
		
		"""Generates .stl file for the geometry.
		
		.. note:: If ``fn=None``, then will use the same filename and path as .geo file.
		
		Keyword Args:
			fn (str): Output filename.
			segments (int): Number of segments used for convex hull of surface.
		
		"""
		
		if fn==None:
			fn=self.fnGeo.replace(".geo",".stl")
		
		fnStl=fn
		fnScad=fnStl.replace(".stl",".scad")
		
		self.render2Openscad(fn=fnScad,segments=segments)
		pyfrp_openscad_module.runOpenscad(fnScad,fnOut=fnStl)
		
		return fnStl
	
class zebrafishDomeStage(geometry):
	
	"""Geometry describing a zebrafish embryo in dome stage.
	
	For information about zebrafish stages, see 
	http://onlinelibrary.wiley.com/doi/10.1002/aja.1002030302/abstract;jsessionid=1EAD19FE5563DAA94E3C22C5D5BEEC85.f01t03.
	
	.. image:: ../imgs/pyfrp_geometry/dome.png
	
	The zebrafish geometry is basically described by two half-balls with different radii piled on top of each other.
	The crucial geometrical parameters are:
		
		* ``outerRadius``: Radius of the outer ball.
		* ``innerRadius``: Radius of the inner ball.
		* ``centerDist``: Distance between the centers of the two balls.
	
	PyFRAP automatically computes these parameters given 
	
		* ``imagingRadius``: Radius of embryo at imaging depth.
		* ``imagingHeight``: Imaging depth.
		* ``radiusScale``: Scaling factor between radii.
		
	For details of this computations, see :py:func:`computeDome`.
		

	"""
	
	def __init__(self,embryo,center,imagingRadius,radiusScale=1.1):
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		
		super(zebrafishDomeStage, self).__init__(embryo,"zebrafishDomeStage",mdir+"dome.geo",center)
		
		#How much bigger is the inner radius than the outerRadius
		self.radiusScale=radiusScale
		self.imagingRadius=imagingRadius
	
		#Compute geometry properties
		self.restoreDefault()
		
		self.updateGeoFile()
	
	def setOuterRadius(self,r):
		
		"""Sets outer radius and updates zebrafish geometry.
		
		Args:
			r (float): New radius.
		
		Returns:
			float: New radius.
		
		"""
		
		self.outerRadius=r
		self.computeDome()
		return self.outerRadius
	
	def setRadiusScale(self,s):
		
		"""Sets scaling factor between outer and inner radius 
		and updates zebrafish geometry.
		
		Args:
			s (float): New scaling factor.
		
		Returns:
			float: New scaling factor.
		
		"""
		
		self.radiusScale=s
		self.computeDome()
		return self.radiusScale
	
	def setImagingRadius(self,r):
		
		"""Sets imaging radius and updates zebrafish geometry.
		
		Args:
			r (float): New radius.
		
		Returns:
			float: New radius.
		
		"""
		
		self.imagingRadius=r
		self.computeDome()
		return self.imagingRadius
	
	def setImagingHeight(self,h):
		
		"""Sets imaging height and updates zebrafish geometry.
		
		Args:
			h (float): New height.
		
		Returns:
			float: New height.
		
		"""
		
		self.imagingHeight=h
		if self.embryo.sliceHeightPx!=self.imagingHeight:
			printWarning("imagingHeight of geometry is not identical with imaging depth of embryo. This can lead to problems.")
		
		self.computeDome()
		return self.imagingHeight
	
	def getImagingRadius(self):
		
		"""Returns imaging radius.
		"""
		
		return self.imagingRadius
	
	def getImagingHeight(self):
		
		"""Returns imaging height.
		"""
		
		return self.imagingHeight
	
	def getOuterRadius(self):
		
		"""Returns outer radius.
		"""
		
		return self.outerRadius
	
	def getInnerRadius(self):
		
		"""Returns inner radius.
		"""
		
		return self.innerRadius
	
	def getRadiusScale(self):
		
		"""Returns radius scaling factor.
		"""
		
		return self.radiusScale
	
	def computeDome(self):
		
		r"""Updates zebrafish geometry.
		
		Computes zebrafish geometry as follows:
		
		.. math:: r_{\mathrm{outer}}=\frac{r_{\mathrm{imaging}}^2+h_{\mathrm{imaging}}^2}{-2 h_{\mathrm{imaging}}},\\
		   r_{\mathrm{inner}}=s_{\mathrm{radius}}r_{\mathrm{outer}},\\
		   d_{\mathrm{center}}=\sqrt{r_{\mathrm{inner}}^2-r_{\mathrm{outer}}^2},
		   
		where :math:`r_{\mathrm{outer}}` is the ``outerRadius``, :math:`r_{\mathrm{inner}}` is the ``innerRadius``, :math:`h_{\mathrm{imaging}}` 
		is the ``imagingHeight`` and :math:`r_{\mathrm{imaging}}` is the ``imagingRadius``.
		
		Returns:
			tuple: Tuple containing:
			
				* innerRadius (float): New inner radius.
				* centerDist (float): New distance between centers.
				
		"""
		
		self.outerRadius=(self.imagingRadius**2+self.imagingHeight**2)/(2*(-self.imagingHeight))
		self.innerRadius=self.outerRadius*self.radiusScale
		self.centerDist=np.sqrt(self.innerRadius**2-self.outerRadius**2) 
		return self.innerRadius,self.centerDist
	
	def restoreDefault(self):
		
		"""Restores default values.
		
		Only default value of zebrafish geometry is that ``imagingHeight`` is set to
		be equal to ``sliceHeightPx`` of ``embryo``.
		
		"""
		
		self.imagingHeight=self.embryo.sliceHeightPx
		self.computeDome()
	
	def updateGeoFile(self,debug=False):
		
		"""Updates .geo file of geometry.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		"""
		
		pyfrp_gmsh_module.updateDomeGeo(self.fnGeo,self.imagingRadius,self.imagingHeight,self.center,run=False,debug=debug)
		
	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False):	
		
		"""Sets optimal ROI to a :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
		with radius ``outerRadius``, center ``center``, covering the whole z-range of 
		geometry.
		
		Keyword Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			color (str): Color that ROI is going to be associated with.
			asMaster (bool): Make new ROI masterROI.
			
		
		"""
		
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getOuterRadius(),0,np.inf,False,color=color,asMaster=asMaster)
	
	def getZExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getZExtend`.
		
		By default, zebrafish geometry is set to range from ``-outerRadius`` to ``0``. 
		"""
		
		z=[-self.outerRadius,0]
		
		return min(z),max(z)
		
	def getXYExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getXYExtend`.
		
		By default, zebrafish geometry is set to range from ``center[i]-outerRadius`` to ``center[i]+outerRadius``. 
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
		
		"""
		
		return self.center[0]-self.outerRadius,self.center[0]+self.outerRadius,self.center[1]-self.outerRadius,self.center[1]+self.outerRadius
	
	def getVolume(self):
		
		r"""Returns volume of geometry.
		
		Volume is computed by:
		
		.. math:: V_{dome}=\frac{\pi}{6}(4 r_{\mathrm{outer}}^3 -  (d_{\mathrm{center}}-r_{\mathrm{inner}})(3r_{\mathrm{outer}}^2+(d_{\mathrm{center}}-r_{\mathrm{inner}})^2))	
		
		Returns:
			float: Volume of zebrafish dome.
		
		"""
		
		Vouter=4./6.*np.pi * self.outerRadius**2
		Vinner=1/6.*np.pi * (self.centerDist-self.innerRadius)*(3*self.outerRadius+(self.centerDist -  self.innerRadius)**2)
		
		return Vouter-Vinner
	
	def genAsOpenscad(self):
		
		"""Generates zebrafish geometry as solid python object.
		
		Useful if geometry is used to be passed to openscad.
		
		Returns:
			solid.solidpython.openscad_object: Solid python object. 
		
		"""
		
		outerBall=solid.translate([self.center[0],self.center[1],-self.outerRadius])(solid.sphere(r=self.outerRadius))
		innerBall=solid.translate([self.center[0],self.center[1],-self.outerRadius-self.centerDist])(solid.sphere(r=self.innerRadius))
	
		return outerBall-innerBall
		
class zebrafishDomeStageQuad(zebrafishDomeStage):
	
	"""Geometry describing a zebrafish embryo in dome stage, reduced to first quadrant.
	
	Inherits from :py:class:`zebrafishDomeStage`. Please refer to its documentation for
	further details.
	
	.. image:: ../imgs/pyfrp_geometry/domequad.png
	
	"""
	
	def __init__(self,embryo,center,imagingRadius,radiusScale=1.1):
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		zebrafishDomeStage.__init__(self,embryo,center,imagingRadius,radiusScale=1.1)
		geometry.__init__(self,embryo,"zebrafishDomeStageQuad",mdir+"quad_dome.geo",center)
	
	
class cylinder(geometry):
	
	"""Geometry describing a cylinder.
	
	.. image:: ../imgs/pyfrp_geometry/cylinder.png
	
	The crucial geometrical parameters are:
		
		* ``radius``: Radius of the cylinder.
		* ``height``: Height of the cylinder.
		
	"""
	
	def __init__(self,embryo,center,radius,height):
		mdir=pyfrp_misc_module.getMeshfilesDir()
		super(cylinder, self).__init__(embryo,"cylinder",mdir+"cylinder.geo",center)

		self.radius=radius
		self.height=height
		self.updateGeoFile()
		
	def setHeight(self,h):
		
		"""Sets cylinder height.
		
		Args:
			h (float): New height.
		
		Returns:
			float: New height.
		
		"""
		
		self.height=h
		return self.height
	
	def setRadius(self,r):
		
		"""Sets cylinder radius.
		
		Args:
			h (float): New radius.
		
		Returns:
			float: New radius.
		
		"""
		
		self.radius=r
		return self.radius
	
	def getHeight(self):
		
		"""Returns cylinder height.
		
		Returns:
			float: Height.
		
		"""
		
		return self.height
	
	def getRadius(self):
		
		"""Returns cylinder radius.
		
		Returns:
			float: Radius.
		
		"""
		
		return self.radius
	
	def updateGeoFile(self,debug=False):
		
		"""Updates .geo file of geometry.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		"""
		
		pyfrp_gmsh_module.updateCylinderGeo(self.fnGeo,self.radius,self.height,self.center,run=False,debug=debug)

	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		
		"""Sets optimal ROI to a :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
		with radius ``radius``, center ``center``, covering the whole z-range of 
		geometry.
		
		Keyword Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			color (str): Color that ROI is going to be associated with.
			asMaster (bool): Make new ROI masterROI.
			
		
		"""
		
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
	
	def getZExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getZExtend`.
		
		By default, cylinder geometry is set to range from ``-height`` to ``0``. 
		"""
		z=[-self.height,0]
		return min(z),max(z)
	
	def getXYExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getXYExtend`.
		
		By default, cylinder geometry is set to range from ``center[i]-radius`` to ``center[i]+radius``. 
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
		
		"""
		
		return self.center[0]-self.radius,self.center[0]+self.radius,self.center[1]-self.radius,self.center[1]+self.radius
	
	def genAsOpenscad(self):
		
		"""Generates cylinder geometry as solid python object.
		
		Useful if geometry is used to be passed to openscad.
		
		Returns:
			solid.solidpython.openscad_object: Solid python object. 
		
		"""
		
		z=self.getZExtend()
		cylinder=solid.translate([self.center[0],self.center[1],min(z)])(solid.cylinder(r=self.radius,h=abs(self.height)))
		
		return cylinder
	
class cylinderQuad(cylinder):
	
	"""Geometry describing a cylinder, reduced to first quadrant.
	
	Inherits from :py:class:`cylinder`. Please refer to its documentation for
	further details.
	
	.. image:: ../imgs/pyfrp_geometry/cylinderquad.png
	
	"""
	
	def __init__(self,embryo,center,radius,height):
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		cylinder.__init__(self,embryo,center,radius,height)
		geometry.__init__(self,embryo,"cylinderQuad",mdir+"quad_cylinder.geo",center)
		
	
class xenopusBall(geometry):
	
	"""Geometry describing a ball.
	
	This geometry is similar to a xenopus in stage 7-10, see http://www.xenbase.org/anatomy/alldev.do.
	
	.. image:: ../imgs/pyfrp_geometry/ball.png
	
	The crucial geometrical parameters are:
		
		* ``radius``: Radius of the ball.
		
	PyFRAP automatically computes these parameters given 
	
		* ``imagingRadius``: Radius of embryo at imaging depth.
		* ``imagingHeight``: Imaging depth.
		
	For details of this computations, see :py:func:`computeBall`.
		
	"""
	
	def __init__(self,embryo,center,imagingRadius):		
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		super(xenopusBall, self).__init__(embryo,"xenopusBall",mdir+"ball.geo",center)
		
		self.imagingRadius=imagingRadius
		
		self.restoreDefault()
		self.updateGeoFile()
	
	def setImagingRadius(self,r):
		
		"""Sets imaging radius and updates ball geometry.
		
		Args:
			h (float): New radius.
		
		Returns:
			float: New radius.
		
		"""
		
		self.imagingRadius=r
		self.computeBall()
		return self.imagingRadius
	
	def setImagingHeight(self,h):
		
		"""Sets imaging height and updates ball geometry.
		
		Args:
			h (float): New height.
		
		Returns:
			float: New height.
		
		"""
		
		
		self.imagingHeight=h
		self.computeBall()
		return self.imagingHeight
	
	def getImagingRadius(self,r):
		
		"""Returns imaging radius."""
		
		return self.imagingRadius
	
	def getImagingHeight(self,h):
		
		"""Returns imaging height."""
		
		return self.imagingHeight
	
	def computeBall(self):
		
		r"""Computes ball geometry from ``imagingRadius`` and ``imagingHeight``.
		
		Computes ball geometry as follows:
		
		.. math:: r=\frac{r_{\mathrm{imaging}}^2+h_{\mathrm{imaging}}^2}{-2 h_{\mathrm{imaging}}},
		   
		where :math:`r` is the ``radius``, :math:`h_{\mathrm{imaging}}` 
		is the ``imagingHeight`` and :math:`r_{\mathrm{imaging}}` is the ``imagingRadius``.
		
		The center of the ball is set to ``[center[0],center[1],-radius]`` (Only in .geo file).
		
		"""
		
		self.radius=(self.imagingRadius**2+self.imagingHeight**2)/(2*(-self.imagingHeight))
	
	def restoreDefault(self):
		
		"""Restores default values.
		
		Only default value of ball geometry is that ``imagingHeight`` is set to
		be equal to ``sliceHeightPx`` of ``embryo``.
		
		"""
		
		self.imagingHeight=self.embryo.sliceHeightPx
		self.computeBall()
		
	def updateGeoFile(self,debug=False):
		
		"""Updates .geo file of geometry.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		"""
		
		pyfrp_gmsh_module.updateBallGeo(self.fnGeo,self.radius,self.center,run=False,debug=debug)
	
	def getRadius(self):
		
		"""Returns ball radius."""
		
		return self.radius
	
	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		
		"""Sets optimal ROI to a :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
		with radius ``radius``, center ``center``, covering the whole z-range of 
		geometry.
		
		Keyword Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			color (str): Color that ROI is going to be associated with.
			asMaster (bool): Make new ROI masterROI.
			
		
		"""
		
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
	
	def getZExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getZExtend`.
		
		By default, ball geometry is set to range from ``-2*imagingRadius`` to ``0``. 
		"""
		
		z=[-2*self.imagingRadius,0]
		
		return min(z),max(z)
	
	def getXYExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getXYExtend`.
		
		By default, ball geometry is set to range from ``center[i]-radius`` to ``center[i]+radius``. 
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
		
		"""
		
		return self.center[0]-self.radius,self.center[0]+self.radius,self.center[1]-self.radius,self.center[1]+self.radius
	
class xenopusBallQuad(xenopusBall):
	
	"""Geometry describing a ball, reduced to first quadrant.
	
	Inherits from :py:class:`xenopusBall`. Please refer to its documentation for
	further details.
	
	.. warning:: ``meshfiles/quad_ball.geo`` does not exist yet. 
	
	"""
	
	def __init__(self,embryo,center,imagingRadius):
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		xenopusBall.__init__(self,embryo,center,imagingRadius)	
		geometry.__init__(self,embryo,"xenopusBallQuad",mdir+"quad_ball.geo",center)	

class cone(geometry):
	
	"""Geometry describing a cut-off cone.
	
	.. image:: ../imgs/pyfrp_geometry/cone.png
	
	The crucial geometrical parameters are:
		
		* ``upperRadius``: Upper radius of cone.
		* ``lowerRadius``: Lower radius of cone.
		* ``height``: Height of cone.
		
	.. note:: Can also be extended to a real cone by setting 
	   ``lowerRadius=0``.
	
	"""
	
	def __init__(self,embryo,center,upperRadius,lowerRadius,height):
		
		mdir=pyfrp_misc_module.getMeshfilesDir()
		super(cone, self).__init__(embryo,"cone",mdir+"cone.geo",center)

		self.upperRadius=upperRadius
		self.lowerRadius=lowerRadius
		
		self.height=height
		
		self.updateGeoFile()
		
	def setHeight(self,h):
		
		"""Sets cone height.
		
		Args:
			h (float): New height.
		
		Returns:
			float: New height.
		
		"""
		
		self.height=h
		return self.height
	
	def setLowerRadius(self,r):
		
		"""Sets cone lower radius.
		
		Args:
			h (float): New lower radius.
		
		Returns:
			float: New lower radius.
		
		"""
		
		self.lowerRadius=r
		return self.lowerRadius
	
	def getLowerRadius(self):
		
		"""Returns lower radius of cone."""
		
		return self.lowerRadius
	
	def setUpperRadius(self,r):
		
		"""Sets cone upper radius.
		
		Args:
			h (float): New upper radius.
		
		Returns:
			float: New upper radius.
		
		"""
		
		self.upperRadius=r
		return self.upperRadius
	
	def getUpperRadius(self):
		
		"""Returns upper radius of cone."""
		
		return self.upperRadius
	
	def getHeight(self):
		
		"""Returns cone radius."""
		
		return self.height
	
	def computeSliceHeightFromRadius(self,radius):
		
		r"""Returns the slice height given a slice radius.
		
		Slice height is computed by
		
		.. math:: s(r) = \frac{h}{l-u} (r-u)
		
		where :math:`l,u` are lower and upper radius respectively and :math:`h` is cone height.
		
		Args:
			radius (float): Slice radius.
			
		Returns:
			float: Slice height.
		
		"""
		
		sliceHeight=self.height/(self.lowerRadius-self.upperRadius)*(radius-self.upperRadius) 
		return sliceHeight
	
	def computeRadiusFromSliceHeight(self,height):
		
		r"""Returns the slice radius given a slice height.
		
		Slice radius is computed by
		
		.. math:: r(s) = \frac{l-u}{h} s +u
		
		where :math:`l,u` are lower and upper radius respectively and :math:`h` is cone height.
		
		Args:
			height (float): Slice height.
			
		Returns:
			float: Slice radius.
		
		"""
		
		radius=(self.upperRadius-self.lowerRadius)/self.height*height+self.upperRadius
		return radius
	
	def updateGeoFile(self,debug=False):
		
		"""Updates .geo file of geometry.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			
		"""
		
		pyfrp_gmsh_module.updateConeGeo(self.fnGeo,self.upperRadius,self.lowerRadius,abs(self.height),self.center,run=False,debug=debug)

	def optimalAllROI(self,name='',Id=0,color='b',asMaster=False,roi=None):	
		
		"""Sets optimal ROI to a :py:class:`pyfrp.subclasses.pyfrp_ROI.radialSliceROI`
		with radius ``upperRadius``, center ``center``, covering the whole z-range of 
		geometry.
		
		Keyword Args:
			name (str): Name of new ROI.
			Id (int): ID of new ROI.
			color (str): Color that ROI is going to be associated with.
			asMaster (bool): Make new ROI masterROI.
			
		
		"""
		
		return self.embryo.newRadialSliceROI(name,Id,self.getCenter(),self.getUpperRadius(),0.,np.inf,False,color=color,asMaster=asMaster)
	
	def getZExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getZExtend`.
		
		By default, cone geometry is set to range from ``-height`` to ``0``. 
		"""
		
		z=[-self.height,0]
		
		return min(z),max(z)
		
	def getXYExtend(self):
		
		"""Overwrites :py:func:`pyfrp.subclasses.pyfrp_geometry.geometry.getXYExtend`.
		
		By default, cone geometry is set to range from ``center[i]-max([self.upperRadius,self.lowerRadius])`` 
		to ``center[i]+max([self.upperRadius,self.lowerRadius])``. 
		
		Returns:
			tuple: Tuple containing:
				
				* xmin (float): Minimum x-coordinate.
				* xmax (float): Maximum x-coordinate.
				* ymin (float): Minimum y-coordinate.
				* ymax (float): Maximum y-coordinate.
		
		"""
		
		maxRadius=max([self.upperRadius,self.lowerRadius])
		
		return self.center[0]-maxRadius,self.center[0]+maxRadius,self.center[1]-maxRadius,self.center[1]+maxRadius
	
	def genAsOpenscad(self):
		
		"""Generates cone geometry as solid python object.
		
		Useful if geometry is used to be passed to openscad.
		
		Returns:
			solid.solidpython.openscad_object: Solid python object. 
		
		"""
		
		z=self.getZExtend()
		cone=solid.translate([self.center[0],self.center[1],min(z)])(solid.cylinder(r1=self.lowerRadius,h=self.height,r2=self.upperRadius))
		
		return cone
	
class custom(geometry):
	
	"""Custom geometry class for custom geometry configurations.

	"""
	
	def __init__(self,embryo,center,fnGeo,dim=3):
		
		geometry.__init__(self,embryo,"custom",fnGeo,center,dim=dim)
	
	