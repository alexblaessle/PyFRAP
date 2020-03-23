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

"""Essential PyFRAP module containing :py:class:`mesh` class. 
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================


#PyFRAP
from pyfrp.modules import pyfrp_gmsh_module
from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules import pyfrp_integration_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_gmsh_geometry
from pyfrp.modules.pyfrp_term_module import *

#FiPy
import fipy

#Numpy/Scipy
import numpy as np

#Misc
import os
import os.path


#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================

class mesh(object):
	
	"""Mesh class for PyFRAP.
	
	The mesh class stores all information about location and creation of the mesh used for a simulation. 
	It is directly associated with the :py:class:`pyfrp.subclasses.pyfrp_simulation.simulation` object
	that uses it.
	
	Meshes can either be created via running Gmsh onto the *.geo* file of the :py:class:`pyfrp.subclasses.pyfrp_geometry.geometry`,
	or by running Gmsh internally from FiPy using some predefined functions (limited geometry support). See also 
	:py:func:`genMesh`.
	
	The most important attributes are:
	
		* ``mesh``: The actual mesh as a ``fipy.GmshImporter3D`` object.
		* ``fromFile``: Flag that controls if mesh should be created from *.geo* file or not.
		* ``volSizePx``: Mesh element size in px.
	
	Besides mesh storage and creation, the mesh class contains useful functions such as:
		
		* Mesh refinement, see :py:func:`refine`, :py:func:`addBoxField` and :py:func:`forceMinMeshDensityInROI`.
		* Information output, see :py:func:`printStats`.
		* Plotting, see :py:func:`plotMesh` and :py:func:`plotDensity`.
	
	Args:
		simulation (pyfrp.subclasses.pyfrp_simulation.simulation): Simulation object.
	
	"""

	def __init__(self,simulation):
		
		#Naming/ID
		self.simulation=simulation
		self.mesh=None
		self.restoreDefaults()
	
	def setVolSizePx(self,v,remesh=True,fnOut=None):
		
		"""Sets volSize of mesh in px.
		
		.. note:: If ``fnOut=None``, then either ``fnMesh`` will be used, or,
		   if ``fnMesh`` is not set yet, will use ``geometry.fnGeo``.
		
		Args:
			v (float): New volSize.
			
		Keyword Args:
			remesh (bool): Generate mesh with new volSize.
			fnOut (str): Output filepath for meshfile.
			dim (int): Dimension of mesh.
		
		Returns:
			float: New volSize.
		
		"""
		
		
		self.volSizePx=v
		self.updateGeoFile()
		if remesh:
			self.fnMesh=self.genMesh(fnOut=fnOut)
		return self.volSizePx
	
	def getVolSizePx(self):
		
		"""Returns mesh volSize in px.
		
		Returns:
			float: VolSize.
		
		"""
		
		return self.volSizePx
		
	def getSimulation(self):
		
		"""Returns :py:class:`pyfrp.subclasses.pyfrp_simulation` that mesh belongs to.
		
		Returns:
			pyfrp.subclasses.pyfrp_simulation: Simulation object.
		"""
		
		return self.simulation
	
	def restoreDefaults(self):	
		
		"""Restores default parameters of mesh.
		
		Default parameters are:
		
			* ``mesh.geometry=mesh.simulation.embryo.geometry``
			* ``mesh.fromFile=True``
			* ``mesh.volSizePx=20``
			* ``mesh.fnMesh=""``
	
		"""
		
		#if self.simulation.embryo.geometry==None:
			#printWarning("Geometry of embryo not specified, mesh generation will not work!")
		
		
		self.fromFile=True
		self.volSizePx=20.
		self.fnMesh=""
		
	def genMesh(self,fnOut=None,debug=False,gmshBin=None):
		
		"""Main mesh generation function.
		
		.. note:: If ``fnOut=None``, will use ``geometry.fnGeo``.
		
		.. note:: If ``fromFile=True``, will generate from .geo file running
		   gmsh directly on the file. If not, will try to run hard coded FiPy
		   version for mesh generation via :py:func:`runFiPyMeshGenerator` .
		
		Keyword Args:
			fnOut (str): Output filepath for meshfile.
			debug (bool): Print debugging messages.
			dim (int): Dimension of mesh.
			
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""
		
		if fnOut==None:
			fnOut=self.simulation.embryo.geometry.fnGeo.replace(".geo",".msh")
		
		dim=self.simulation.embryo.geometry.getDim()
		
		if self.fromFile:
			self.fnMesh=pyfrp_misc_module.fixPath(fnOut)
		
			pyfrp_gmsh_module.runGmsh(self.simulation.embryo.geometry.fnGeo,fnOut=fnOut,debug=debug,volSizeMax=self.volSizePx,dim=dim,gmshBin=gmshBin)
			
			self.importMeshFromFile(self.fnMesh)
		else:
			self.runFiPyMeshGenerator(self.simulation.embryo.geometry.typ)

		return self.mesh
	
	def runFiPyMeshGenerator(self,typ):
		
		"""Runs gmsh on the via FiPy internally defined meshes.
		
		Available meshes:
			
			* cylinder
			* zebrafishDomeStage
			* xenopusBall
			
		.. note:: Any refinement method will not work if mesh is created this way.
		
		Args:
			typ (str): Type of mesh to be created (see list above).
				
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
			
		"""
		
		if typ=="cylinder":
			self.mesh=pyfrp_gmsh_module.genFiPyCylinderMesh(self.volSizePx,self.simulation.embryo.geometry.radius,self.simulation.embryo.geometry.height,self.simulation.embryo.geometry.center)
		elif typ=="zebrafishDomeStage":
			self.mesh=pyfrp_gmsh_module.genFiPyDomeMesh(self.volSizePx,self.simulation.embryo.geometry.innerRadius,self.simulation.embryo.geometry.outerRadius,self.simulation.embryo.geometry.centerDist,self.simulation.embryo.geometry.center)
		elif typ=="xenopusBall":
			self.mesh=pyfrp_gmsh_module.genFiPyBallMesh(self.volSizePx,self.simulation.embryo.geometry.radius,self.simulation.embryo.geometry.center)
		else:
			printWarning("Geometry type " + typ + "unknown for runFiPyMeshGenerator. Check pyfrp_gmsh_module for available geometries.")
		
		return self.mesh
		
	def importMeshFromFile(self,fn):
		
		"""Imports mesh from a Gmsh .msh file.
		
		See also http://www.ctcms.nist.gov/fipy/fipy/generated/fipy.meshes.html. 
		
		Args:
			fn (str): Filepath to meshfile.
		
		Keyword Args:
			dim (int): Dimension of mesh.
		
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
			
		"""
		
		try:
			dim=self.simulation.embryo.geometry.getDim()
		except AttributeError:
			printWarning("Was not able to receive geometry's dimension. Will assume dim=3.")
			dim=3 
			
		if dim==3:
			try:
				self.mesh=fipy.GmshImporter3D(fn)
			except AttributeError:
				self.mesh=fipy.Gmsh3D(fn)
		elif dim==2:
			try:
				self.mesh=fipy.GmshImporter2D(fn)
			except AttributeError:
				self.mesh=fipy.Gmsh2D(fn)
		
		else:
			printError("Unknown dimensionality dim = "+str(dim))
		self.fnMesh=fn
		return self.mesh
		
	def setMesh(self,m):
		
		"""Sets ``mesh`` attribute to a new mesh.
		
		Args:
			m (fipy.GmshImporter3D): New mesh.
				
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""
		
		self.mesh=m
		return self.mesh
	
	def setFromFile(self,v):
		
		"""Sets flag if mesh is supposed to be created from file (recommended)
		or from internally defined mesh creation method.
		
		Args:
			v (bool): New flag value.
				
		Returns:
			bool: New flag value.
		
		"""
		
		self.fromFile=v
		return self.fromFile
	
	def getMesh(self):
		
		"""Returns mesh that is used for simulation.
		
		Returns:
			fipy.GmshImporter3D: Gmsh mesh object.
		
		"""

		return self.mesh
	
	def getFnMesh(self):
		
		"""Returns the filepath of meshfile.
		"""
		
		return self.fnMesh
	
	def setFnMesh(self,fn):
		
		
		"""Sets the filepath of meshfile.
		
		Imports the new mesh right away using :py:func:`importMeshFromFile`.
		
		"""
		
		
		self.fnMesh=pyfrp_misc_module.fixPath(fn)
		self.importMeshFromFile(self.fnMesh)
		return self.fnMesh
	
	def updateGeoFile(self,debug=False):
		
		"""Updates geometry file by writing new ``volSizePx``
		into ``embryo.geometry.fnGeo``.
		
		Keyword Args:
			debug (bool): Print debugging messages.
		
		Returns:
			str: Path to ouput meshfile.
		
		"""
		
		if os.path.isfile(self.simulation.embryo.geometry.fnGeo):
			self.fnMesh=pyfrp_gmsh_module.updateVolSizeGeo(self.simulation.embryo.geometry.fnGeo,self.volSizePx,run=False,debug=debug)
		else:
			printWarning("Cannot update meshfile, since geometry.fnGeo does not exist or is invalid.")
		return self.fnMesh
	
	def refine(self,debug=False):
		
		"""Refines mesh by splitting.
		
		See also http://gmsh.info/doc/texinfo/gmsh.html .
		
		Keyword Args:
			debug (bool): Print debugging messages.
		
		"""
		
		pyfrp_gmsh_module.refineMsh(self.fnMesh,debug=debug)
		self.importMeshFromFile(self.fnMesh)
		return self.fnMesh
	
	def forceMinMeshDensityInROI(self,ROI,density,stepPercentage=0.1,debug=False,findIdxs=True,method='refine',maxCells=100000):
		
		"""Forces global mensh density such that a certain density is reached in a
		given ROI.
		
		Tries to achive a mesh density ``density`` in ``ROI`` by globally refining mesh either 
		through decreasing ``volSizePx`` by ``stepPercentage`` percent (``method=volSize``), or 
		by using Gmsh's ``-refine`` option (``method=refine``). If maximum number of cells is 
		exceeded, will use the last mesh that did not exceed ``maxCells``.
		
		Args:
			ROI (pyfrp.subclasses.pyfrp_ROI.ROI): ROI object.
			density (float): Desired density.
			
		Keyword Args:
			stepPercentage (float): If method is ``volSize``, percentage of ``volSize`` decrease.
			method (str): Refinement method (``refine``/``volSize``).
			maxCells (int): Total maximum number of mesh cells allowed.
			findIdxs (bool): Find ROI indices after refinement.
			debug (bool): Print debugging messages.
		
		Returns:
			float: New ``volSizePx``
		
		"""
		
		#Set counter to 0
		j=0
		
		#Loop until mesh density has been found
		while ROI.getMeshDensity()<density:
			
			#Backup old stuff
			volSizeBackup=self.getVolSizePx()
			meshBackup=self.getMesh()
			
			#Refine
			if method=='volSize':
				self.setVolSizePx((1-stepPercentage)*self.getVolSizePx())
			elif method=='refine':
				self.refine()
			else:
				printError("Unknown method: ", method)
				return
			
			#Recompute Idxs
			ROI.computeMeshIdx(self.mesh)
			
			#Debugging output
			if debug:
				if method=='volSize':
					print("Tried volSizePx ", self.getVolSizePx())
				elif method=='refine':
					print("Refinement step, ", j)
				
				print("Current density ", ROI.getMeshDensity(), " desired density ", density)
				print("Current number of nodes in ROI ", len(ROI.meshIdx))
				print(" Mesh now has ", np.shape(self.mesh.x)[0], "cells.") 
			
			#Check if we reached maxCells
			if np.shape(self.mesh.x)[0]>maxCells:
				if debug:
					printWarning("Reached maxCells, will keep previous iteration.")
				
				self.mesh=meshBackup
				self.setVolSizePx(volSizeBackup,remesh=False)
				
			#Increment counter
			j=j+1
			
		if debug:
			print("volSizePx = ", self.getVolSizePx(), "is sufficient.")
		
		#Recompute idxs for all ROIs
		if findIdxs:
			self.simulation.embryo.computeROIIdxs()
					
		return self.getVolSizePx()	
	
	def getNNodes(self):
		
		"""Returns number of nodes in mesh. 
		
		If no mesh has been generated yet, will return 0.
		
		Returns:
			int: Number of nodes.
		
		"""
		
		if self.mesh==None:
			return 0
		else:
			return len(self.getCellCenters()[0])
	

		
	def writeVTKFile(self,fn="",sub=False):
		
		"""Writes mesh into vtk file.
		
		Uses *meshIO* (https://github.com/nschloe/meshio), to convert the mesh saved
		in ``fnMesh`` to a .vtk file. 
		
		If ``sub==True``, will start a seperate subprocess and submit 
		``pyfrp_meshIO_script.py`` to it. This can be sometimes useful, since PyFRAP
		sometimes tends to crash otherwise.
		
		If no output path is given via ``fn``, will use same path as ``fnMesh``.
		
		.. note:: *meshIO* only gets imported inside this function, making PyFRAP running
		   even without the package installed. However, this feature will only run with
		   *meshIO*.
		
		Keyword Args:
			fn (str): Optional output path.
			sub (bool): Subprocess flag.
		
		Returns:
			str: Used output path.
		
		"""
		
		
		if not os.path.isfile(self.fnMesh):
			printWarning("Filepath to meshfile has not been specified yet. Cannot write VTK file.")
		
		if fn=="":
			fn=self.fnMesh.replace('.msh','.vtk')
		
		if sub:
			
			cmd = "python pyfrp_meshIO_script.py "+ self.fnMesh
			
			import shlex
			import subprocess
			args = shlex.split(cmd)
			
			p = subprocess.Popen(args)
			p.wait()
		else:
			#MeshIO
			import meshio
			
			points, cells, point_data, cell_data, field_data = meshio.read(self.fnMesh)
			meshio.write(fn,points,cells,point_data=point_data,cell_data=cell_data,field_data=field_data)
			
		return fn
	
	def importVTKFile(self,fnVTK="",sub=False):
		
		"""Imports a .vtk file into a vtk renderer.
		
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		If ``sub==True``, will start a seperate subprocess and submit 
		``pyfrp_meshIO_script.py`` to it. This can be sometimes useful, since PyFRAP
		sometimes tends to crash otherwise.
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
	
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			sub (bool): Subprocess flag.
			
		Returns:
			vtk.vtkRenderer: Renderer object.
		
		"""
		
		#vtk
		from pyfrp.modules import pyfrp_vtk_module
		
		if not os.path.isfile(self.fnMesh):
			printWarning("Filepath to meshfile has not been specified yet. Cannot plot.")
		
		if fnVTK=="":
			fnVTK=self.writeVTKFile(sub=sub)
		
		renderer = pyfrp_vtk_module.importVTKMeshFile(fnVTK,bkgdColor=[1,1,1],color=[0,0,0])
	
		return renderer
	
	def plotMesh(self,fromFile=False,fnVTK="",bkgd=[255,255,255],color=[0,0,0],renderer=None,renderWindow=None):
		
		"""Plots the mesh using VTK.
		
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			bkgd (list): Background color of renderer in RGB values.
			color (list): Color of mesh in RGB values.
			renderer (vtk.vtkRenderer): Some renderer.
			renderWindow (vtk.RenderWindow): Some renderWindow or dummy value.
			
		Returns:
			vtk.vtkRenderer: Renderer object.
		
		"""
		
		from pyfrp.modules import pyfrp_vtk_module
	
		# Import vtk file
		if fromFile:
			renderer=self.importVTKFile(fnVTK=fnVTK)
		
		# Generate everything directly from mesh
		else:
			grid=pyfrp_vtk_module.meshToUnstructeredGrid(self.mesh)
			renderer = pyfrp_vtk_module.unstructedGridToWireframe(grid,bkgdColor=bkgd,color=color,renderer=renderer)
			renderer, renderWindow, renderWindowInteractor = pyfrp_vtk_module.makeVTKCanvas(offScreen=False,bkgd=bkgd,renderer=renderer)
			
			#rendererWindow=renderer.GetRenderWindow()
			#print rendererWindow
			
		# Create renderWindow
		if renderWindow==None:
			renderer, renderWindow, renderWindowInteractor = pyfrp_vtk_module.makeVTKCanvas(offScreen=False,bkgd=bkgd,renderer=renderer)
		
			# This seems to be the key to have a renderer without displaying. However it returns floating 
			# Errors.
			##renderWindow.OffScreenRenderingOn() 
				
		#Start
		renderWindow.GetInteractor().Initialize()
		renderWindow.GetInteractor().Start()
			
			#print rendererWindow
			
		return renderWindow
	
	def saveMeshToPS(self,fnOut,fnVTK="",renderer=None):
		
		"""Saves mesh to postscript file.
		
		Supported extensions are:
		
			* '.ps'  (PostScript)
			* '.eps' (Encapsualted PostScript)
			* '.pdf' (Portable Document Format)
			* '.tex' (LaTeX)
			* '.svg' (Scalable Vector Graphics)
			
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		If no ``renderer`` is given, will create one using :py:func:`plotMesh`. 
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
		Some code taken from http://www.programcreek.com/python/example/23102/vtk.vtkGL2PSExporter .
		
		Args:
			fnOut (str): Path to output file.
		
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			renderer (vtk.vtkOpenGLRenderer): Renderer.
			magnification (int): Degree of magnification.
			
		Returns:
			vtk.vtkGL2PSExporter: Exporter object.
		
		"""
		
		from pyfrp.modules import pyfrp_vtk_module
		
		#Plot again if necessary
		if renderer==None:
			renderer=self.plotMesh(fnVTK=fnVTK)
		
		return pyfrp_vtk_module.saveRendererToPS(renderer,fnOut)
	
	def saveMeshToImg(self,fnOut,fnVTK="",renderer=None,magnification=10,show=True):
		
		"""Saves mesh to image file.
		
		Supported extensions are:
			
			* '.ps'  (PostScript)
			* '.eps' (Encapsualted PostScript)
			* '.pdf' (Portable Document Format)
			* '.jpg' (Joint Photographic Experts Group)
			* '.png' (Portable Network Graphics)
			* '.pnm' (Portable Any Map)
			* '.tif' (Tagged Image File Format)
			* '.bmp' (Bitmap Image)
		
		If ``fnVTK`` is not given, will generate .vtk file from meshfile stored in ``fnMesh``
		using :py:func:`writeVTKFile`.
		
		If no ``renderer`` is given, will create one using :py:func:`plotMesh`. 
		
		.. note:: This function imports *vtk*. *vtk* is only necessary in a few functions,
		   hence only imported when needed. This should make PyFRAP more portable.
		
		Some code taken from http://www.programcreek.com/python/example/23102/vtk.vtkGL2PSExporter .
		
		Args:
			fnOut (str): Path to output file.
		
		Keyword Args:
			fnVTK (str): Path to input vtk file.
			renderer (vtk.vtkOpenGLRenderer): Renderer.
			magnification (int): Degree of magnification.
			show (bool): Show vtk render window.
			
		Returns:
			vtk.vtkExporter: Exporter object.
		
		"""
		
		from pyfrp.modules import pyfrp_vtk_module
		
		#Plot again if necessary
		if renderer==None:
			if show:
				rendererWindow=self.plotMesh(fnVTK=fnVTK)
				renderer=rendererWindow.GetRenderers().GetFirstRenderer()
			else:
				renderer=self.importVTKFile(fnVTK=fnVTK)
				
		return pyfrp_vtk_module.saveRendererToImg(renderer,fnOut,magnification=magnification)
	
	def printStats(self,tetLenghts=False):
		
		"""Prints out statistics of mesh.
		
		Also calculates all tetraheder lengths if ``tetLenghts`` is selected. 
		This might take some time depending on mesh size.
		
		Keyword Args:
			tetLenghts (bool): Also calculate and print out tetrahedra sidelengths.
		
		"""
	
		
		x,y,z = self.getCellCenters()
		
		print("-------------------------------------------")
		print("Mesh Statistics:")
		print("-------------------------------------------")
		print("Mesh has ", np.shape(x)[0] , " cells")
		print("Mesh has ", self.mesh._numberOfVertices, " vertices")
		print("Mesh has ", self.mesh.numberOfFaces, " faces")
		print("min x=", min(x), "max x=", max(x))
		print("min y=", min(y), "max y=", max(y))
		print("min z=", min(z), "max z=", max(z))
		print("Maximum cell volume= ", max(self.getCellVolumes()))
		print("Minimum cell volume= ", min(self.getCellVolumes()))
			
		print("Maximum cell volume is", max(self.getCellVolumes()), "in cell number=", np.argmax(self.getCellVolumes()))
		print("Minimum cell volume is", min(self.getCellVolumes()), "in cell number=", np.argmin(self.getCellVolumes()))
		print("Average cell volume is", np.mean(self.getCellVolumes()))

		if tetLenghts:
			slsVec=self.calcAllTetSidelenghts()
			print("Average sidelength of tetrahedron in self.mesh:", np.mean(slsVec))
			print("Maximum sidelength of tetrahedron in self.mesh:", max(slsVec))
			print("Minimum sidelength of tetrahedron in self.mesh:", min(slsVec))	
		
		print()
		
		return
	
	def calcAllTetSidelenghts(self):
		
		"""Calculates sidelengths of all tetrahedra.
		
		See also :py:func:`pyfrp.modules.pyfrp_integration_module.calcTetSidelengths`.
		
		Returns:
			list: List of all sidelengths.
		"""
		
		#Calculating sidelengths of tetrahedra
		slsVec=[]
		for i in range(np.shape(self.mesh._getOrderedCellVertexIDs())[1]):
			currVert=self.mesh._getOrderedCellVertexIDs()[:,i]
			
			pt1=self.mesh.vertexCoords[:,currVert[0]]
			pt2=self.mesh.vertexCoords[:,currVert[1]]
			pt3=self.mesh.vertexCoords[:,currVert[2]]
			pt4=self.mesh.vertexCoords[:,currVert[3]]
			
			sl1,sl2,sl3=pyfrp_integration_module.calcTetSidelengths(pt1,pt2,pt3,pt4)
			
			slsVec.append(sl1)
			slsVec.append(sl2)
			slsVec.append(sl3)
			
		return slsVec
	
	def plotDensity(self,axes=None,hist=True,bins=100,color='b'):
		
		"""Plots the mesh density in x/y/z-direction.
		
		``hist=True`` is recommended, since otherwise plots generally appear fairly 
		noisy. 
		
		.. note:: If no ``axes`` are given or they do not have the necessary size, 
		   will create new ones.
		
		.. image:: ../imgs/pyfrp_mesh/density_plot.png
		
		Keyword Args:
			axes (list): List of ``matplotlib.axes``.
			hist (bool): Summarize densities in bins.
			bins (int): Number of bins used for hist.
			color (str): Color of plot.
			
		Returns:
			list: List of ``matplotlib.axes``.
			
		"""
		
		x,y,z=self.getCellCenters()
		
		volSortedByX,xSorted=pyfrp_misc_module.sortListsWithKey(self.getCellVolumes(),x)
		volSortedByY,ySorted=pyfrp_misc_module.sortListsWithKey(self.getCellVolumes(),y)
		volSortedByZ,zSorted=pyfrp_misc_module.sortListsWithKey(self.getCellVolumes(),z)
		
		if axes==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["Density(x)","Density(y)","Density(z)"])
		else:
			if len(axes)<3:
				printWarning("axes do not have right have, will create new ones.")
				fig,axes = pyfrp_plot_module.makeSubplot([1,3],titles=["Density(x)","Density(y)","Density(z)"])
		
		if hist:
			xSorted,volSortedByX=pyfrp_misc_module.simpleHist(xSorted,volSortedByX,bins)
			ySorted,volSortedByY=pyfrp_misc_module.simpleHist(ySorted,volSortedByY,bins)
			zSorted,volSortedByZ=pyfrp_misc_module.simpleHist(zSorted,volSortedByZ,bins)
			
		axes[0].plot(xSorted,volSortedByX,color=color)
		axes[1].plot(ySorted,volSortedByY,color=color)
		axes[2].plot(zSorted,volSortedByZ,color=color)
		
		for ax in axes:
			pyfrp_plot_module.redraw(ax)
		
		return axes
	
	def plotCellCenters(self,ax=None,proj=None,color='k',indicateHeight=False,s=5.,roi=None):
		
		"""Plots location of cell centers of mesh.
		
		.. note:: If no ``ax`` are given will create new ones.
		
		If ``proj=[3d]``, will create 3D scatter plot, otherwise project cell centers in 
		2D.
		
		Example:
		
		Create figure
		
		>>> fig,axes = pyfrp_plot_module.makeSubplot([2,2],titles=['2D','2D indicate','3D','3D indicate'],proj=[None,None,'3d','3d'])

		Plot in 4 different ways
		
		>>> mesh.plotCellCenters(ax=axes[0],s=1.)
		>>> mesh.plotCellCenters(ax=axes[1],indicateHeight=True,s=5.)
		>>> mesh.plotCellCenters(ax=axes[2],s=3.)
		>>> mesh.plotCellCenters(ax=axes[3],indicateHeight=True,s=3.)
		
		.. image:: ../imgs/pyfrp_mesh/plotCellCenters.png
		
		Keyword Args:
			ax (matplotlib.axes): Axes to plot in.
			proj (list): List of projections.
			color (str): Color of mesh nodes.
			indicateHeight (bool): Indicate height by color.
			s (float): Size of marker.
			roi (pyfrp.subclasses.pyfrp_ROI): ROI.
			
		Returns:
			matplotlib.axes: Matplotlib axes.
		
		
		"""
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeSubplot([1,1],titles=["Cell Centers"],proj=proj)
			ax=axes[0]
		
		
		x,y,z = self.getCellCenters()
		
		if roi!=None:
			x=x[roi.meshIdx]
			y=y[roi.meshIdx]
			z=z[roi.meshIdx]
			
		if pyfrp_plot_module.is3DAxes(ax):
			if indicateHeight:
				#color=cm.jet()
				ax.scatter(x,y,z,c=z,s=s)
			else:	
				ax.scatter(x,y,z,c=color,s=s)
		else:
			if indicateHeight:
				ax.scatter(x,y,c=z,s=s)
			else:
				ax.scatter(x,y,c=color,s=s)
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
		
	def printAllAttr(self):
		
		"""Prints out all attributes of mesh object.""" 
		
		print("Mesh of embryo ", self.simulation.embryo.name, " Details.")
		printAllObjAttr(self)
	
	def addBoxField(self,volSizeIn,rangeX,rangeY,rangeZ,newFile=True,fnAppendix="_box",comment="newField",run=False,fnOut=None):
		
		"""Adds box field to mesh.
		
		Box fields allow to refine certain areas of the mesh, see also 
		http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes .
		
		.. note:: Will keep ``volSizePx`` as volSize outside of the box.
		
		.. note:: If ``fnOut`` is not specified, will do the following: 
		   
		   - If ``newFile=True``, will create new file with path ``fnGeo+/field/custom/fnGeo+fnAppendix.geo``
		   - Else writes into ``fnGeo``.
		
		Args:
			volSizeIn (float): volSize in px inside the box.
			rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
			rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
			rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
			
		Keyword Args:	
			newFile (bool): Write new mesh into a new .geo file.
			fnAppendix (str): Append this to new file name.
			comment (str): Comment in .geo file before definition of box field.
			run (bool): Run Gmsh on new .geo file afterwards.
			fnOut (str): Path to output geo file.
			dim (int): Dimension of mesh.
			
		Returns:
			str: Path to new .geo file.
			
		"""
		
		if fnOut==None:
			if newFile:
				if "_box" not in self.simulation.embryo.geometry.fnGeo:	
					pyfrp_misc_module.mkdir(self.simulation.embryo.geometry.fnGeo+"/field/")
					pyfrp_misc_module.mkdir(self.simulation.embryo.geometry.fnGeo+"/field/custom/")
					fnOut=os.path.dirname(self.simulation.embryo.geometry.fnGeo)+"/field/custom/"+os.path.basename(self.simulation.embryo.geometry.fnGeo).replace(".geo",fnAppendix+"_"+self.simulation.embryo.geometry.embryo.name+".geo")
				else:
					fnOut=fnOut=self.simulation.embryo.geometry.fnGeo
			else:
				fnOut=self.simulation.embryo.geometry.fnGeo
			
		pyfrp_gmsh_IO_module.addBoxField(self.simulation.embryo.geometry.fnGeo,volSizeIn,self.volSizePx,rangeX,rangeY,rangeZ,comment=comment,fnOut=fnOut)
		self.simulation.embryo.geometry.setFnGeo(fnOut)
		
		if run:
			self.genMesh()
		
		return fnOut
	
	def addBoundaryLayerAroundROI(self,roi,fnOut=None,segments=48,simplify=True,iterations=3,triangIterations=2,
			       fixSurfaces=True,debug=False,volSizePx=None,volSizeLayer=10,thickness=15.,cleanUp=True,
			       approxBySpline=True,angleThresh=0.95,faces='all',onlyAbs=True):
		
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
			dim (int): Dimension of mesh.
			
		Returns:
			str: Path to new .geo file.
		
		"""
		
		if fnOut==None:
			fnOut=self.simulation.embryo.geometry.fnGeo.replace(".geo","_"+roi.name+"_BL.geo")
		
		#Build stl file for ROI
		fnStl=roi.render2StlInGeometry(segments=segments,fn=fnOut.replace(".geo",".stl"))
		
		#Read in stl file and simplify
		dROI=pyfrp_gmsh_IO_module.readStlFile(fnStl)
		
		if simplify: 
			dROI.simplifySurfaces(iterations=iterations,triangIterations=triangIterations,fixSurfaces=fixSurfaces,debug=debug,addPoints=bool(triangIterations>0))		
	
		# If we extract faces, we have to put them in new ROI
		if faces!='all' and len(faces)>0:
			
			# Grab surfaces
			sfs=dROI.getRuledSurfacesByNormal(faces,onlyAbs=onlyAbs)
			
			# Insert in new domain
			d2=pyfrp_gmsh_geometry.domain()
			for sf in sfs:
				sf.extract(d=d2,debug=debug)
			
			# Overwrite sfs and dROI
			sfs=d2.ruledSurfaces
			dROI=d2
			
		#Approximate complicated curves by splines
		if approxBySpline:
			for sf in sfs:
				sf.lineLoop.approxBySpline(angleThresh=angleThresh)	
		dROI.fixAllLoops()
		
		#Read in geometry and merge 
		dGeo=self.simulation.embryo.geometry.readGeoFile()
		dGeo.merge(dROI)
		
		#Create boundary layer
		blf=dGeo.addBoundaryLayerField(hfar=volSizeLayer,hwall_n=volSizeLayer,hwall_t=volSizeLayer,thickness=thickness,Quads=0.)
		
		#Add surfaces to boundary layer mesh
		blf.addFaceListByID(pyfrp_misc_module.objAttrToList(sfs,'Id'))
		blf.setAsBkgdField()
		
		#Set global volSize
		if volSizePx!=None:
			dGeo.setGlobalVolSize(volSizePx)
			self.volSizePx=volSizePx
		else:
			dGeo.setGlobalVolSize(self.getVolSizePx())
		
		#Write and set as new geometry
		dGeo.writeToFile(fnOut)
		self.simulation.embryo.geometry.setFnGeo(fnOut)
		
		#Generate mesh
		self.genMesh()
		
		#Clean up temporary files
		if cleanUp:
			os.remove(fnStl)
			os.remove(fnStl.replace(".stl",".scad"))
		
		return fnOut
		
	
	def getMaxNodeDistance(self):
		
		"""Returns maximum node distance in x/y/z direction.
		
		Returns:
			tuple: Tuple containing:
			
				* dmaxX (float): Maximum distance in x-direction
				* dmaxY (float): Maximum distance in y-direction
				* dmaxZ (float): Maximum distance in z-direction
				
		"""
		
		distances=self.mesh.cellDistanceVectors
		
		return max(distances[0]),max(distances[1]),max(distances[2])
		
	def getCellCenters(self):
		
		"""Returns cell centers of mesh.
		
		If the mesh is 2-dimensional, places all nodes at ``z=self.simulation.embryo.sliceHeight``.
		
		Returns:
			tuple: Tuple containing:
			
				* x (list): x-coordinates of cells.
				* y (list): y-coordinates of cells.
				* z (list): z-coordinates of cells.
				
		"""
		
		if self.mesh==None:
			return [],[],[]
		
		try:
			cellCenters=self.mesh.getCellCenters()
		except AttributeError:
			cellCenters=self.mesh.cellCenters
		
		if len(cellCenters)==3:
			return cellCenters
		else:
			z=self.simulation.embryo.sliceHeightPx*np.ones((len(cellCenters[0]),))
			return cellCenters[0],cellCenters[1],z
		
	
	def getCellVolumes(self):
		
		try:
			cellVolumes=self.mesh.getCellVolumes()
		except AttributeError:
			cellVolumes=self.mesh.cellVolumes
		return 	cellVolumes
			
		
	
	#def getCellCenters():
		
		#try:
			#cellCenters=self.mesh.getCellCenters()
		#except AttributeError:
			#cellCenters=self.mesh.cellCenters
		#return cellCenters	
