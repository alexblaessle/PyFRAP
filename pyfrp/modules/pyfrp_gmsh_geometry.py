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

"""PyFRAP module for creating/extracting gmsh geometries for PyFRAP toolbox. Module mainly has the following classes:

	* A ``domain`` class, acting as a canvas.
	* A ``vertex`` class, substituting gmsh's *Point*.
	* A ``edge`` class, parenting all different kind of edges.
	* A ``line`` class, substituting gmsh's *Line*.
	* A ``arc`` class, substituting gmsh's *Circle*.
	* A ``bSpline`` class, substituting gmsh's *bSpline*.
	* A ``lineLoop`` class, substituting gmsh's *Line Loop*.
	* A ``ruledSurface`` class, substituting gmsh's *Ruled Surface*.
	* A ``surfaceLoop`` class, substituting gmsh's *Surface Loop*.
	* A ``volume`` class, substituting gmsh's *Volume*.
	* A ``field`` class, parenting all different kind of fields.
	* A ``attractorField`` class, substituting gmsh's *Attractor* field.
	* A ``boundaryLayerField`` class, substituting gmsh's *Boundary Layer* field.
	* A ``thresholdField`` class, substituting gmsh's *Threshold* field.
	* A ``minField`` class, substituting gmsh's *Min* field.
	
This module together with pyfrp.pyfrp_gmsh_IO_module and pyfrp.pyfrp_gmsh_module works partially as a python gmsh wrapper, however is incomplete.
If you want to know more about gmsh, go to http://gmsh.info/doc/texinfo/gmsh.html .
"""
#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#String
import string

#PyFRAP modules
from . import pyfrp_plot_module
from .pyfrp_term_module import *
from . import pyfrp_misc_module
from . import pyfrp_gmsh_IO_module
from . import pyfrp_idx_module
from . import pyfrp_geometry_module
from . import pyfrp_IO_module
from . import pyfrp_vtk_module

#Matplotlib
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import copy as cpy

		
#===========================================================================================================================================================================
#Class definitions
#===========================================================================================================================================================================
	
class domain:
	
	"""Domain class storing embryo geometry entities.
	
	Args:
		edges (list): List of edges.
		vertices (list): List of vertices.
		arcs (list): List of arcs.
		lines (list): List of lines.
		bSplines (list): List of bSplines.
		lineLoops (list): List of lineLoops.
		surfaceLoops (list): List of surfaceLoops.
		ruledSurfaces (list): List of ruledSurfaces.
		volumes (list): List of volumes.
		fields (list): List of fields.
		annXOffset (float): Offset of annotations in x-direction.
		annYOffset (float): Offset of annotations in y-direction.
		annZOffset (float): Offset of annotations in z-direction.
			
	"""
	
	def __init__(self):
		
		
		#Lists to keep track of all geometrical entities.
		self.edges=[]
		self.vertices=[]
		self.arcs=[]
		self.lines=[]
		self.bSplines=[]
		self.lineLoops=[]
		self.ruledSurfaces=[]
		self.surfaceLoops=[]
		self.volumes=[]
		self.fields=[]
		self.bkgdField=None
		
		#Some settings for plotting
		self.annXOffset=3.
		self.annYOffset=3.
		self.annZOffset=3.
		
	def addVertex(self,x,Id=None,volSize=None,checkExist=False):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` instance
		at point ``x`` and appends it to ``vertices`` list.
		
		.. note:: ``volSize`` does not have any effect on the geometry itself but is simply 
		   stored in the vertex object for further usage.
		
		If ``checkExist=True``, checks if a vertex at same location already exists. If so, will return 
		that vertex instead.
		
		Args:
			x (numpy.ndarray): Coordinate of vertex.
			
		Keyword Args:
			Id (int): ID of vertex.
			volSize (float): Element size at vertex.
			checkExist (bool): Checks if a vertex at same location already exists.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.vertex: New vertex instance.
		
		"""
		
		if checkExist:
			for v in self.vertices:
				if (v.x==x).all():
					return v
		
		newId=self.getNewId(self.vertices,Id)
		
		v=vertex(self,x,newId,volSize=volSize)
		self.vertices.append(v)	
		
		return v
	
	def addLine(self,v1,v2,Id=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.line` instance
		at point ``x`` and appends it to ``edges`` and ``lines`` list.
		
		Args:
			v1 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
			v2 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
			
		Keyword Args:
			Id (int): ID of line.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.line: New line instance.
		
		"""
		
		newId=self.getNewId(self.edges,Id)
		
		e=line(self,v1,v2,newId)
		self.lines.append(e)
		self.edges.append(e)
		
		return e
	
	def addArc(self,vstart,vcenter,vend,Id=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.arc` instance
		at point ``x`` and appends it to ``edges`` and ``arcs`` list.
		
		Args:
			vstart (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
			vcenter (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Center vertex.
			vend (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
			
		Keyword Args:
			Id (int): ID of arc.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.arc: New line instance.
		
		"""
		
		
		newId=self.getNewId(self.edges,Id)
			
		a=arc(self,vstart,vcenter,vend,newId)
		self.arcs.append(a)
		self.edges.append(a)
		
		return a
	
	def addBSpline(self,vertices,Id=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.line` instance
		at point ``x`` and appends it to ``edges`` and ``lines`` list.
		
		Args:
			vertices (list): List of vertex objects.
			
		Keyword Args:
			Id (int): ID of spline.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.bSpline: New spline instance.
		
		"""
		
		newId=self.getNewId(self.edges,Id)
		
		e=bSpline(self,vertices,newId)
		self.bSplines.append(e)
		self.edges.append(e)
		
		return e
	
	def addCircleByParameters(self,center,radius,z,volSize,plane="z",genLoop=False,genSurface=False,checkExist=True):
		
		"""Adds circle to domain by given center and radius.
		
		Will create 5 new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` objects 
		``[vcenter,v1,v2,v3,v4]`` and four new `pyfrp.modules.pyfrp_gmsh_geometry.arc` objects
		[a1,a2,a3,a4] and builds circle.
		
		Circle  will be at ``z=z`` and vertices will have mesh size ``volSize``.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addCircleByParameters([256,256],100,50,30.)
		>>> d.addCircleByParameters([256,256],100,50,30.,plane="x")
		>>> d.addCircleByParameters([256,256],100,50,30.,plane="y")
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addCircleByParameters.png
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		Args:
			center (numpy.ndarray): Center of circle.
			radius (float): Radius of the circle.
			z (float): Height at which circle is placed.
			volSize (float): Mesh size of vertices.
		
		Keyword Args:
			plane (str): Plane in which circle is placed.
			genLoop (bool): Create lineLoop.
			genSurface (bool): Create ruledSurface.
			checkExist (bool): Checks if a vertex at same location already exists.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* arcs (list): List of arcs.
				* loop (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): Line loop.
				* surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Ruled Surface.
				
		
		"""
		
		# Define coordinates
		xcenter=pyfrp_geometry_module.flipCoordinate([center[0],center[1],z],plane,origAxis="z")
		x1=pyfrp_geometry_module.flipCoordinate([center[0]+radius,center[1],z],plane,origAxis="z")
		x2=pyfrp_geometry_module.flipCoordinate([center[0],center[1]+radius,z],plane,origAxis="z")
		x3=pyfrp_geometry_module.flipCoordinate([center[0]-radius,center[1],z],plane,origAxis="z")
		x4=pyfrp_geometry_module.flipCoordinate([center[0],center[1]-radius,z],plane,origAxis="z")
		
		# Add vertices
		vcenter=self.addVertex(xcenter,volSize=volSize)
		v1=self.addVertex(x1,volSize=volSize,checkExist=checkExist)
		v2=self.addVertex(x2,volSize=volSize,checkExist=checkExist)
		v3=self.addVertex(x3,volSize=volSize,checkExist=checkExist)
		v4=self.addVertex(x4,volSize=volSize,checkExist=checkExist)
		
		# Add Arcs
		a1=self.addArc(v1,vcenter,v2)
		a2=self.addArc(v2,vcenter,v3)
		a3=self.addArc(v3,vcenter,v4)
		a4=self.addArc(v4,vcenter,v1)
		
		if genLoop or genSurface:
			loop=self.addLineLoop(edgeIDs=[a1.Id,a2.Id,a3.Id,a4.Id])
		else:
			loop=None
			
		if genSurface:
			surface=self.addRuledSurface(lineLoopID=loop.Id)
		else:
			surface=None
		
		return [vcenter,v1,v2,v3,v4],[a1,a2,a3,a4],loop,surface
	
	def addPolygonByParameters(self,coords,volSize,z=0.,plane="z",genLoop=False,genSurface=False):
		
		"""Adds polygon to domain by given vertex coordinates.
		
		Will create a list of new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` objects 
		and a list of new `pyfrp.modules.pyfrp_gmsh_geometry.line` objects
		connecting the vertices.
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		.. note:: Vertices can be given either as a 
		
			* list of coordinate triples ``[[x1,y1,z1],[x2,y2,z2],...]``.
			* list of x-y-coordinates and a given z-coordinate ``[[x1,y1,z],[x2,y2,z],...]``.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addPolygonByParameters([[100,100,100],[200,200,100],[200,100,100]],30.)
		>>> d.addPolygonByParameters([[100,100,100],[200,200,100],[200,100,100]],30.,plane="x")
		>>> d.addPolygonByParameters([[100,100,100],[200,200,100],[200,100,100]],30.,plane="y")
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addPolygonByParameters.png
		
		
		.. note:: Vertices are created in the order of the coordinates and connected in the same order.
		

		Args:
			coords (list): List of coordinates.
			volSize (float): Mesh size of vertices.
			
		Keyword Args:
			plane (str): Plane in which polygon is placed.
			z (float): Height at which polygon is placed.
			genLoop (bool): Create lineLoop.
			genSurface (bool): Create ruledSurface.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* lines (list): List of connecting lines.
				* loop (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): Line loop.
				* surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Ruled Surface.
				
		"""
		
		# Define coordinates
		xs=[]
		for c in coords:
			if len(c)==3:
				xs.append(pyfrp_geometry_module.flipCoordinate([c[0],c[1],c[2]],plane,origAxis="z"))
			else:
				xs.append(pyfrp_geometry_module.flipCoordinate([c[0],c[1],z],plane,origAxis="z"))
		
		# Add vertices
		vertices=[]
		for x in xs:
			vertices.append(self.addVertex(x,volSize=volSize))
		
		# Add Lines
		lines=[]
		for i in range(len(vertices)):
			lines.append(self.addLine(vertices[i],vertices[pyfrp_misc_module.modIdx(i+1,vertices)]))
		
		# Add LineLoop
		if genLoop or genSurface:
			loop=self.addLineLoop(edgeIDs=pyfrp_misc_module.objAttrToList(lines,'Id'))
		else:
			loop=None
		
		# Add surface
		if genSurface:
			surface=self.addRuledSurface(lineLoopID=loop.Id)
		else:
			surface=None
		
		return vertices,lines,loop,surface
	
	def addRectangleByParameters(self,offset,sidelengthX,sidelengthY,z,volSize,plane="z"):
		
		"""Adds rectangle to domain by given offset and sidelengths.
		
		Will create a list of four :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` objects 
		and a list of four `pyfrp.modules.pyfrp_gmsh_geometry.line` objects
		connecting the vertices.
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		.. note:: The ``offset`` is defined as the bottom left corner.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addRectangleByParameters([256,256],100,200,50,30.)
		>>> d.addRectangleByParameters([256,256],100,200,50,30.,plane="x")
		>>> d.addRectangleByParameters([256,256],100,200,50,30.,plane="y")
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addRectangleByParameters.png
		
		Args:
			offset (numpy.ndarray): Offset of rectangle.
			sidelengthX (float): Sidelength in x-direction.
			sidelengthY (float): Sidelength in y-direction.
			z (float): Height at which rectangle is placed.
			volSize (float): Mesh size of vertices.
			
		Keyword Args:
			plane (str): Plane in which rectangle is placed.
			genLoop (bool): Create lineLoop.
			genSurface (bool): Create ruledSurface.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* lines (list): List of connecting lines.
				* loop (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): Line loop.
				* surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Ruled Surface.
				
		"""
		
		coords=[[offset[0],offset[1],z],[offset[0]+sidelengthX,offset[1],z],
		[offset[0]+sidelengthX,offset[1]+sidelengthY,z],[offset[0],offset[1]+sidelengthY,z]]
		
		return self.addPolygonByParameters(coords,volSize,plane=plane)
	
	def addSquareByParameters(self,offset,sidelength,z,volSize,plane="z"):
		
		"""Adds square to domain by given offset and sidelength.
		
		Will create a list of four :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.vertex` objects 
		and a list of four `pyfrp.modules.pyfrp_gmsh_geometry.line` objects
		connecting the vertices.
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		.. note:: The ``offset`` is defined as the bottom left corner.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addSquareByParameters([256,256],100,50,30.)
		>>> d.addSquareByParameters([256,256],100,50,30.,plane="x")
		>>> d.addSquareByParameters([256,256],100,50,30.,plane="y")
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addSquareByParameters.png
		
		Args:
			offset (numpy.ndarray): Offset of square.
			sidelength (float): Sidelength of square.
			
			z (float): Height at which square is placed.
			volSize (float): Mesh size of vertices.
			
		Keyword Args:
			plane (str): Plane in which square is placed.
			genLoop (bool): Create lineLoop.
			genSurface (bool): Create ruledSurface.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* lines (list): List of connecting lines.
				* loop (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): Line loop.
				* surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Ruled Surface.
				
		"""
		
		return self.addRectangleByParameters(offset,sidelength,sidelength,z,volSize,plane=plane)
	
	def addPrismByParameters(self,coords,volSize,height=1.,z=0.,plane="z",genLoops=True,genSurfaces=True,genSurfaceLoop=True,genVol=True):
		
		r"""Adds prism to domain by given vertex coordinates.
		
		Will create:
		
			* 2 new polygons, see :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addPolygonByParameters`.
			* n :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.line` objects connecting the two polyogns.
		
		If selected, will create:
			
			* n+2 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` objects around the 6 surfaces.
			* n+2 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` objects.
			* 1 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop`.
			* 1 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.volume`.
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		.. note:: Vertices can be given either as a 
		
			* list of coordinate triples ``[[x1,y1,z1],[x2,y2,z2],...]``. Then the list of vertices needs to be of length :math:`2n`, where
			  where :math:`n` is the number of corners of the top and lower polygon. Otherwise :py:func:`addPrismByParameters` will crash.
			* list of x-y-coordinates, a given z-coordinate and height. This will place the vertices at ``[[x1,y1,z],[x2,y2,z],...]`` and 
			  ``[[x1,y1,z+height],[x2,y2,z+height],...]``.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addPrismByParameters([[256,256],[200,220],[200,200],[210,210],[220,200]],30.,z=50.,height=40.,plane="z",genLoops=True,genSurfaces=True,genVol=True)
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addPrismByParameters.png
		
		.. note:: Vertices are created in the order of the coordinates and connected in the same order.
		
		Args:
			coords (list): List of coordinates.
			volSize (float): Mesh size of vertices.
			
		Keyword Args:
			plane (str): Plane in which prism is placed.
			z (float): Height at which first polygon is placed.
			height (float): Height of prism.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genSurfaceLoop (bool): Generate surface loop.
			genVol (bool): Generate corresponding volume.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* lines (list): List of lines.
				* loops (list): List of loops.
				* surfaces (list): List of surfaces.
				* surfaceLoop (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): Generated surface loop.
				* vol (pyfrp.modules.pyfrp_gmsh_geometry.volume): Generated volume.
		
		"""
		
		# Create upper and lower polygons
		if len(coords[0])==3:
			if np.mod(len(coords),2)!=0:
				printError("addPrismByParameters: You gave a list of 3-dimensional vertex coordinates. However,the number of coordinates is odd, will not be able to continue.")
				return
				
			vertices,lines,ltemp,stemp = self.addPolygonByParameters(coords,volSize,z=0.,plane="z",genSurface=genSurface,genLoop=genLoop)
			vertices1=vertices[:len(vertices)/2]
			vertices2=vertices[len(vertices)/2:]
			lines1=lines[:len(lines)/2]
			lines2=lines[len(lines)/2:]
				
		else:	
			vertices1,lines1,ltemp,stemp = self.addPolygonByParameters(coords,volSize,z=z,plane="z")
			vertices2,lines2,ltemp,stemp = self.addPolygonByParameters(coords,volSize,z=z+height,plane="z")
			
			
		# Connect them with lines
		lines3=[]
		for i in range(len(vertices1)):
			lines3.append(self.addLine(vertices1[i],vertices2[i]))
		
		# Add loops
		loops=[]
		if genLoops:
			
			# Loops of upper and lower polygon
			loops.append(self.addLineLoop(edgeIDs=pyfrp_misc_module.objAttrToList(lines1,"Id")))
			loops.append(self.addLineLoop(edgeIDs=pyfrp_misc_module.objAttrToList(lines2,"Id")))
			
			# Loops of side faces
			for i in range(len(lines1)):
				loops.append(self.addLineLoop(edgeIDs=[-lines1[i].Id,lines3[i].Id,lines2[i].Id,-lines3[pyfrp_misc_module.modIdx(i+1,lines1)].Id]))
				
		# Add surfaces
		surfaces=[]
		if genSurfaces:
			for loop in loops:
				surfaces.append(self.addRuledSurface(lineLoopID=loop.Id))
				
		# Make surface loop
		if genSurfaceLoop:
			surfaceLoop=self.addSurfaceLoop(surfaceIDs=pyfrp_misc_module.objAttrToList(surfaces,'Id'))
		else:
			surfaceLoop=None
		
		# Make volume
		if genVol:	
			vol=self.addVolume(surfaceLoopID=surfaceLoop.Id)
		else:
			vol=None	
		
		return [vertices1,vertices2],[lines1,lines2,lines3],loops,surfaces,surfaceLoop,vol
	
	def addCuboidByParameters(self,offset,sidelengthX,sidelengthY,height,volSize,plane="z",genLoops=True,genSurfaces=True,genSurfaceLoop=True,genVol=True):
		
		"""Adds Cuboid to domain by given offset, sidelengths in x- and y-direction and height.
		
		Will define vertices and then call :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addPrismByParameters`.
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		
		>>> d.draw()
		
		will generate:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addCuboidByParameters.png
		
		Args:
			offset (numpy.ndarray): Offset of cuboid.
			sidelengthX (float): Sidelength in x-direction.
			sidelengthY (float): Sidelength in y-direction.
			height (float): Height of cuboid.
			volSize (float): Mesh size of vertices.
			
		Keyword Args:
			plane (str): Plane in which prism is placed.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genSurfaceLoop (bool): Generate surface loop.
			genVol (bool): Generate corresponding volume.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* lines (list): List of lines.
				* loops (list): List of loops.
				* surfaces (list): List of surfaces.
				* surfaceLoop (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): Generated surface loop.
				* vol (pyfrp.modules.pyfrp_gmsh_geometry.volume): Generated volume.
		
		"""
		
		# Define coordinates
		coords=[[offset[0],offset[1]],[offset[0]+sidelengthX,offset[1]],
		[offset[0]+sidelengthX,offset[1]+sidelengthY],[offset[0],offset[1]+sidelengthY]]
		
		return self.addPrismByParameters(coords,volSize,height=height,z=offset[2],plane="z",genLoops=genLoops,genSurfaces=genSurfaces,genSurfaceLoop=genSurfaceLoop,genVol=genVol)
	
	def addBallByParameters(self,center,radius,z,volSize,genLoops=True,genSurfaces=True,genSurfaceLoop=True,genVol=True,checkExist=True):
		
		"""Adds ball to domain by given center and radius.
		
		Will create.
		
			* 3 new circles, see :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addCircleByParameters`.
			
		If selected, will create:
			
			* 8 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` objects around the 8 surfaces.
			* 8 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` objects.
			* 1 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop`.
			* 1 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.volume`.
		
		For example:
		
		>>> center=[256,50]
		>>> radius=100
		>>> Z=0
		>>> volSize=20
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addBallByParameters(center,radius,Z,volSize,genLoops=True,genSurfaces=True,genVol=True,checkExist=True)
		>>> d.draw()
		
		would return:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addBallByParameters.png
		
		Args:
			center (numpy.ndarray): Center of cylinder.
			radius (float): Radius of the cylinder.
			z (float): Height at which cylinder is placed.
			volSize (float): Mesh size of vertices.
		
		Keyword Args:
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genSurfaceLoop (bool): Generate surface loop.
			genVol (bool): Generate volume.
			checkExist (bool): Checks if a vertex at same location already exists.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* arcs (list): List of arcs.
				* loops (list): List of loops.
				* surfaces (list): List of surfaces.
				* surfaceLoop (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): Generated surface loop.
				* vol (pyfrp.modules.pyfrp_gmsh_geometry.volume): Generated volume.
		
		"""
	
		# Add 3 circles
		v1,a1,ll1,s1=self.addCircleByParameters(center,radius,z,volSize,genLoop=False)
		v2,a2,ll2,s2=self.addCircleByParameters([z,center[0]],radius,center[1],volSize,genLoop=False,plane='x')
		v3,a3,ll3,s3=self.addCircleByParameters([center[1],z],radius,center[0],volSize,genLoop=False,plane='y')
		vertices=v1+v2+v3
		arcs=a1+a2+a3
		
		# Define line loops
		ll=[]
		ll.append([a1[0],a2[0],a3[0]])
		ll.append([a1[1],a3[0],a2[3]])
		ll.append([a1[2],a2[3],a3[1]])
		ll.append([a1[3],a3[1],a2[0]])

		ll.append([a1[0],a2[1],a3[3]])
		ll.append([a1[1],a3[3],a2[2]])
		ll.append([a1[2],a2[2],a3[2]])
		ll.append([a1[3],a3[2],a2[1]])

		# Generate line loops
		lineLoops=[]
		if genLoops:
			for l in ll:
				lnew=self.addLineLoop(edgeIDs=pyfrp_misc_module.objAttrToList(l,'Id'))
				lnew.fix()
				lineLoops.append(lnew)
			
		# Generate surfaces
		surfaces=[]
		if genSurfaces:
			for l in lineLoops:
				surfaces.append(self.addRuledSurface(lineLoopID=l.Id))
			
		# Make surface loop
		if genSurfaceLoop:
			surfaceLoop=self.addSurfaceLoop(surfaceIDs=pyfrp_misc_module.objAttrToList(surfaces,'Id'))
		else:
			surfaceLoop=None
		
		# Make volume
		if genVol:	
			vol=self.addVolume(surfaceLoopID=surfaceLoop.Id)
		else:
			vol=None
		
		return vertices,arcs,lineLoops,surfaces,surfaceLoop,vol
			
	def addCylinderByParameters(self,center,radius,z,height,volSize,plane="z",genLoops=True,genSurfaces=True,genSurfaceLoop=True,genVol=True,checkExist=True):
		
		"""Adds cylinder to domain by given center and radius and height.
		
		Will create.
		
			* 2 new circles at ``z=z`` and ``z=z+height``, see :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addCircleByParameters`.
			* 4 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.line` objects connecting the two circles.
		
		If selected, will create:
			
			* 6 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` objects around the 6 surfaces.
			* 6 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` objects.
			* 1 :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop`.
			* 1 corresponding :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.volume`.
		
		For example:
		
		>>> d=pyfrp_gmsh_geometry.domain()
		>>> d.addCylinderByParameters([256,256],100,50,100,30.,plane="z",genLoops=True,genSurfaces=True,genVol=True)
		>>> d.draw()
		
		would return:
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/addCylinderByParameters.png
		
		.. note:: Plane can be given as ``"x","y","z"``. See also :py:func:`pyfrp.modules.pyfrp_geometry_module.flipCoordinate`.
		
		Args:
			center (numpy.ndarray): Center of cylinder.
			radius (float): Radius of the cylinder.
			z (float): Height at which cylinder is placed.
			height (float): Height of cylinder.
			volSize (float): Mesh size of vertices.
		
		Keyword Args:
			plane (str): Plane in which cylinder is placed.
			genLoops (bool): Generate line loops.
			genSurfaces (bool): Generate surfaces.
			genSurfaceLoop (bool): Generate surface loop.
			genVol (bool): Generate volume.
			checkExist (bool): Checks if a vertex at same location already exists.
			
		Returns:
			tuple: Tuple containing:
			
				* vertices (list): List of vertices.
				* arcs (list): List of arcs.
				* lines (list): List of lines.
				* loops (list): List of loops.
				* surfaces (list): List of surfaces.
				* surfaceLoop (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): Generated surface loop.
				* vol (pyfrp.modules.pyfrp_gmsh_geometry.volume): Generated volume.
		
		"""
		
		# Check input
		if genVol and not genSurfaces:
			printError("Cannot create volume when there are no surfaces.")
		if genSurfaces and not genLoops:
			printError("Cannot create surfaces when there are no loops.")
			
		# Create circles
		vertices1,arcs1,ltemp,stemp=self.addCircleByParameters(center,radius,z,volSize,plane=plane)
		vertices2,arcs2,ltemp,stemp=self.addCircleByParameters(center,radius,z+height,volSize,plane=plane)
		
		# Create connecting lines
		lines=[]
		lines.append(self.addLine(vertices1[1],vertices2[1]))
		lines.append(self.addLine(vertices1[2],vertices2[2]))
		lines.append(self.addLine(vertices1[3],vertices2[3]))
		lines.append(self.addLine(vertices1[4],vertices2[4]))
		
		# Generate loops
		loops=[]
		if genLoops:
			loops.append(self.addLineLoop(edgeIDs=[arcs1[0].Id,arcs1[1].Id,arcs1[2].Id,arcs1[3].Id]))
			loops.append(self.addLineLoop(edgeIDs=[arcs2[0].Id,arcs2[1].Id,arcs2[2].Id,arcs2[3].Id]))
			loops.append(self.addLineLoop(edgeIDs=[-lines[0].Id,arcs1[0].Id,lines[1].Id,-arcs2[0].Id]))
			loops.append(self.addLineLoop(edgeIDs=[-lines[1].Id,arcs1[1].Id,lines[2].Id,-arcs2[1].Id]))
			loops.append(self.addLineLoop(edgeIDs=[-lines[2].Id,arcs1[2].Id,lines[3].Id,-arcs2[2].Id]))
			loops.append(self.addLineLoop(edgeIDs=[-lines[3].Id,arcs1[3].Id,lines[0].Id,-arcs2[3].Id]))
				  
		# Generate surfaces
		surfaces=[]
		surfaceIds=[]
		if genSurfaces:
			for loop in loops:
				surfaces.append(self.addRuledSurface(lineLoopID=loop.Id))
				surfaceIds.append(surfaces[-1].Id)
				
		# Generate surface loop and volume
		if genSurfaceLoop:
			surfaceLoop=self.addSurfaceLoop(surfaceIDs=surfaceIds)
		else:	
			surfaceLoop=None
		
		if genVol:
			vol=self.addVolume(surfaceLoopID=surfaceLoop.Id)
		else:
			vol=None
		
		return [vertices1,vertices2],[arcs1,arcs2],lines,loops,surfaces,surfaceLoop,vol
	
	def insertVertex(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts vertex into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.vertex): A vertex.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated edges list.
		
		"""
		
		if self.getVertexByX(obj.x)[0]!=False:
			if debug:
				printWarning("Vertex with x=" +str(obj.x) + " already exists.")
		
		return self.insertElement("vertices",obj,copy=copy,strict=strict,debug=debug)
		
	def insertEdge(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts edge into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.edge): A edge.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated edges list.
		
		"""
		
		LOld=len(self.edges)
		l=self.insertElement("edges",obj,copy=copy,strict=strict,debug=debug)
		b=(LOld<len(self.edges))
		
		if b:
			if obj.typ==0:
				self.lines.append(obj)
			if obj.typ==1:
				self.arcs.append(obj)
			if obj.typ==2:
				self.bSpline.append(obj)
				
		return l
	
	def insertLineLoop(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts line loop into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): A line loop.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated lineLoops list.
		
		"""
		
		return self.insertElement("lineLoops",obj,copy=copy,strict=strict,debug=debug)
	
	def insertRuledSurface(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts ruled surface into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): A ruled surface.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated ruledSurfaces list.
		
		"""
		
		return self.insertElement("ruledSurfaces",obj,copy=copy,strict=strict,debug=debug)
	
	def insertSurfaceLoop(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts surface loop into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): A surface loop.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated surfaceLoops list.
		
		"""
		
		return self.insertElement("ruledSurfaces",obj,copy=copy,strict=strict,debug=debug)
	
	def insertVolume(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts volume into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.volume): A volume.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updates volumes list.
		
		"""
		
		return self.insertElement("volumes",obj,copy=copy,strict=strict,debug=debug)
	
	def insertField(self,obj,copy=False,strict=True,debug=False):
		
		"""Inserts field into domain.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.insertElement`.
		
		Args:
			obj (pyfrp.modules.pyfrp_gmsh_geometry.field): A field.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updates fields list.
		
		"""
		
		return self.insertElement("fields",obj,copy=copy,strict=strict,debug=debug)
	
	def insertElement(self,element,obj,copy=False,strict=True,debug=False):
		
		"""Inserts gmshElement into domain.
		
		Checks if there is already a element with ID.
		
		.. note:: If ``copy=True``, will generate copy of element. This might mess
		   with some connection between elements. Thus ``copy=False`` as default.
		
		Possible values for ``element`` are:
			
			* vertices
			* lines
			* arcs
			* lineLoops
			* bSplines
			* ruledSurfaces
			* surfaceLoops
			* volumes
			* fields
			* auto
			
		.. note:: ``element='auto'`` will automatically detect the type of element and insert it at the right 
		   point.
		
		Will automatically set ``self`` as element's domain.
		
		.. note:: If ``strict=True``, will not allow double IDs.
		
		Args:
			element (str): Name of element list where object belongs.
			obj (pyfrp.modules.pyfrp_gmsh_geometry.gmshElement): Element to insert.
		
		Keyword Args:
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: Updated respective element list.
			
		"""
		
		if element=='auto':
			element = obj.getTypeListName()
			if element in ['edges','bSplines','arcs','lines']:
				return self.insertEdge(obj,strict=strict,copy=copy,debug=debug)
			if element=='vertices':
				return self.insertVertex(obj,strict=strict,copy=copy,debug=debug)
				
		if self.checkIdExists(obj.Id,getattr(self,element),debug=debug):
			if debug:
				printWarning(obj.getType() + " with Id=" +str(obj.getID()) + " already exits.")
			if strict:
				return getattr(self,element)
				
		if copy:
			e=obj.getCopy()
		else:
			e=obj
			
		e.domain=self
		getattr(self,element).append(e)
		
		return getattr(self,element)
	
	def getRuledSurfacesByNormal(self,direction,onlyAbs=True):
		
		"""Returns all surfaces in domain that have given normal vector.
		
		The direction can be given in multiple ways:
		
			* A ``numpy.ndarray``: The method will look for all surfaces with same normal vector than array.
			* A ``str``: The method will first check if ``direction='all'`` is given. If so, return all surfaces.
			  Otherwise the method will decode the string ("x"/"y","z") using 
			  :py:func:`pyfrp.modules.pyfrp_geometry_module.decodeEuclideanBase`, then proceed the same way as 
			  with the ``numpy.ndarray``.
			* A ``list`` of the previous options: Will find all surface matching each of them.
			
		.. note:: If ``onlyAbs=True``, will only look for matches in terms of absolute value. If a list of directions is 
		   given, then one can also specifiy a list of ``onlyAbs`` values.
			
		Args:
			direction (numpy.ndarray): Direction to be matched.
			
		Keyword Args:
			onlyAbs (bool): Only try to match in terms of absolute value.
			
		Returns:
			list: List of matching surfaces.
		
		"""
		
		# Check if we need to return all of them
		if direction=='all':
			return list(self.ruledSurfaces)
		
		# Result list
		sfs=[]
		
		# If list, call recursively
		if isinstance(direction,list):
			
			if not isinstance(onlyAbs,list):
				onlyAbs=len(direction)*[onlyAbs]
		
			for i,d in enumerate(direction):
				sfs=sfs+self.getRuledSurfacesByNormal(d,onlyAbs=onlyAbs[i])
			return sfs	
		
		# If given as string, decode
		if isinstance(direction,str):
			d=pyfrp_geometry_module.decodeEuclideanBase(direction)
		else:
			d=direction
		
		# Look for matching surfaces
		for sf in self.ruledSurfaces:
			if onlyAbs:
				if pyfrp_misc_module.compareVectors(abs(sf.getNormal().astype('float')),d.astype('float')):
					sfs.append(sf)
			else:
				if pyfrp_misc_module.compareVectors(sf.getNormal().astype('float'),d.astype('float')):	
					sfs.append(sf)
					
		return sfs			
					
	def checkIdExists(self,Id,objList,debug=False):
		
		"""Checks if any object in ``objList`` already has ID ``Id``.
		
		Args:
			Id (int): ID to be checked.
			objList (list): List of objects, for example ``edges``.
		
		Keyword Args:
			debug (bool): Print debugging output.
		
		Returns:
			bool: True if any object has ID ``Id``.
		
		"""
		
		IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
		if Id in IdList:
			if debug:
				printWarning("Object with Id " + str(Id) + " already exists.")
			return True
		return False
	
	def getNewId(self,objList,Id=None):
		
		"""Returns free ID for object type.
		
		Args:
			objList (list): List of objects, for example ``edges``.
			
		Keyword Args:
			Id (int): ID to be checked.
			
		Returns:
			int: New free ID.
		
		"""
		
		if Id==None:
			newId=self.incrementID(objList)
		else:
			if self.checkIdExists(Id,objList):
				newId=self.incrementID(objList)
			else:
				newId=Id
		
		return newId
		
	def incrementID(self,objList):
		
		"""Returns ID that is by one larger for a specific 
		object type.
		
		Args:
			objList (list): List of objects, for example ``edges``.
				
		Returns:
			int: Incremented ID.
		
		"""
		
		if len(objList)==0:
			newId=1
		else:
			IdList=pyfrp_misc_module.objAttrToList(objList,'Id')
			newId=max(IdList)+1		
		return newId
		
	def getEdgeById(self,ID):
		
		"""Returns edge with ID ``ID``.
		
		Returns ``(False,False)`` if edge cannot be found.
		
		Args:
			ID (int): ID of edge.
				
		Returns:
			tuple: Tuple containing:
				
				* e (pyfrp.modules.pyfrp_gmsh_geometry.edge): Edge.
				* i (int): Position in ``edges`` list.
		
		"""
		
		for i,e in enumerate(self.edges):
			if e.Id==ID:
				return e,i
		return False,False
	
	def getEdgeByVertices(self,v1,v2):
		
		"""Returns edge between vertex ``v1`` and ``v2``.
		
		Returns ``(False,False)`` if edge cannot be found.
		
		Args:
			v1 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Vertex 1.
			v2 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Vertex 2.
			
		Returns:
			tuple: Tuple containing:
				
				* e (pyfrp.modules.pyfrp_gmsh_geometry.edge): Edge.
				* i (int): Position in ``edges`` list.
		
		"""
		
		for i,e in enumerate(self.edges):
			vertices=[e.getFirstVertex(1),e.getLastVertex(1)]
		
			if v1 in vertices and v2 in vertices:
				return e,i
		return False,False	
		
	
	def getLineLoopById(self,ID):
		
		"""Returns lineLoop with ID ``ID``.
		
		Returns ``(False,False)`` if lineLoop cannot be found.
		
		Args:
			ID (int): ID of lineLoop.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): lineLoop.
				* i (int): Position in ``lineLoops`` list.
		
		"""
		
		for i,l in enumerate(self.lineLoops):
			if l.Id==ID:
				return l,i
	
		return False,False
	
	def getRuledSurfaceById(self,ID):
		
		"""Returns ruledSurface with ID ``ID``.
		
		Returns ``(False,False)`` if ruledSurface cannot be found.
		
		Args:
			ID (int): ID of ruledSurface.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): ruledSurface.
				* i (int): Position in ``ruledSurfaces`` list.
		
		"""
		
		for i,l in enumerate(self.ruledSurfaces):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getSurfaceLoopById(self,ID):
		
		"""Returns surfaceLoop with ID ``ID``.
		
		Returns ``(False,False)`` if surfaceLoop cannot be found.
		
		Args:
			ID (int): ID of surfaceLoop.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop): surfaceLoop.
				* i (int): Position in ``surfaceLoops`` list.
		
		"""
		
		for i,l in enumerate(self.surfaceLoops):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getVolumeById(self,ID):
		
		"""Returns volume with ID ``ID``.
		
		Returns ``(False,False)`` if volume cannot be found.
		
		Args:
			ID (int): ID of volume.
				
		Returns:
			tuple: Tuple containing:
				
				* l (pyfrp.modules.pyfrp_gmsh_geometry.volume): volume.
				* i (int): Position in ``volumes`` list.
		
		"""
		
		for i,l in enumerate(self.volumes):
			if l.Id==ID:
				return l,i
		return False,False
	
	def getFieldById(self,ID):
		
		"""Returns field with ID ``ID``.
		
		Returns ``(False,False)`` if field cannot be found.
		
		Args:
			ID (int): ID of field.
				
		Returns:
			tuple: Tuple containing:
				
				* f (pyfrp.modules.pyfrp_gmsh_geometry.field): Field.
				* i (int): Position in ``fields`` list.
		
		"""
		
		for i,f in enumerate(self.fields):
			if f.Id==ID:
				return f,i
		return False,False
	
	def getVertexById(self,ID):
		
		"""Returns vertex with ID ``ID``.
		
		Returns ``(False,False)`` if vertex cannot be found.
		
		Args:
			ID (int): ID of vertex.
				
		Returns:
			tuple: Tuple containing:
				
				* v (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Vertex.
				* i (int): Position in ``vertices`` list.
		
		"""
		
		for i,v in enumerate(self.vertices):
			if v.Id==ID:
				return v,i
		return False,False
	
	def getVertexByX(self,x):
		
		"""Returns vertex at coordinate ``x``.
		
		Returns ``(False,False)`` if vertex cannot be found.
		
		Args:
			x (numpy.ndarry): Coordinate of vertex.
				
		Returns:
			tuple: Tuple containing:
				
				* v (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Vertex.
				* i (int): Position in ``vertices`` list.
		
		"""
		
		for i,v in enumerate(self.vertices):
			if (np.array(x)==v.x).sum()==len(v.x):
				return v,i
		return False,False
		
	def draw(self,ax=None,color='k',ann=None,drawSurfaces=False,surfaceColor='b',alpha=0.2,backend='mpl',asSphere=True,size=5,annElements=[True,True,True],linewidth=1):
		
		"""Draws complete domain.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		With ``annElements`` the user has the possibility to only annotate given elements. For example
		``annElements=[False,True,False]`` only annotates edges. 
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
			asSphere (bool): Draws vertex as sphere (only in vtk mode).
			size (float): Size of vertex (only in vtk mode).
			annElements (list): Only annotate some element types.
			linewidth (float): Line width.
			
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		if ann==None:
			ann=False
		
		for v in self.vertices:
			ax=v.draw(ax=ax,color=color,ann=ann*annElements[0],backend=backend,size=size,asSphere=asSphere,render=False)
		for e in self.edges:
			ax=e.draw(ax=ax,color=color,ann=ann*annElements[1],backend=backend,render=False,linewidth=linewidth)
		if drawSurfaces:
			for s in self.ruledSurfaces:
				ax=s.draw(ax=ax,color=surfaceColor,alpha=alpha,backend=backend,ann=ann*annElements[2])
		
		if backend=="vtk":
			ax=pyfrp_vtk_module.renderVTK(ax,start=False)
		
		return ax
		
	def getExtend(self):
		
		"""Returns extend of domain in all 3 dimensions.
		
		Returns: 
			tuple: Tuple containing:
				
				* minx (float): Minimal x-coordinate.
				* maxx (float): Maximal x-coordinate.
				* miny (float): Minimal y-coordinate.
				* maxy (float): Maximal y-coordinate.
				* minz (float): Minimal z-coordinate.
				* maxz (float): Maximal z-coordinate.
	
		"""
			
		x=[]
		y=[]
		z=[]
		for v in self.vertices:
			x.append(v.x[0])
			y.append(v.x[1])
			z.append(v.x[2])
		return min(x), max(x), min(y),max(y), min(z),max(z)
	
	def verticesCoordsToList(self):
		
		"""Returns list of coordinates from all vertrices.
		
		Returns:
			list: List of (x,y,z) coordinates.
		
		"""
		
		l=[]
		for v in self.vertices:
			l.append(v.x)
		return l
	
	def setGlobalVolSize(self,volSize):
		
		"""Sets volSize for all nodes in geometry.
		
		"""
		
		for v in self.vertices:
			v.volSize=volSize
			
	
	def addLineLoop(self,Id=None,edgeIDs=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` instance
		with given edgeIDs. 
			
		Keyword Args:
			edgeIDs (list): List of edge IDs included in line loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.lineLoop: New lineLoop instance.
		
		"""
		
		newId=self.getNewId(self.lineLoops,Id)
		
		l=lineLoop(self,edgeIDs,newId)
		self.lineLoops.append(l)
		
		return l
	
	def addAllSurfacesToLoop(self):
		
		"""Adds all surfaces in domain to a single surfaceLoop.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop: New surfaceLoop instance.
		
		"""
		
		surfaceIDs=pyfrp_misc_module.objAttrToList(self.ruledSurfaces,'Id')
		
		return self.addSurfaceLoop(surfaceIDs=surfaceIDs)
	
	def addEnclosingVolume(self):
		
		"""Adds volume enclosing all surfaces. 
		
		See also :py:func:`addAllSurfacesToLoop`.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.volume: New volume instance.
		
		"""
		
		s=self.addAllSurfacesToLoop()
		return self.addVolume(surfaceLoopID=s.Id)
		
	def addSurfaceLoop(self,Id=None,surfaceIDs=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop` instance
		with given surfaceIDs. 
			
		Keyword Args:
			surfaceIDs (list): List of surface IDs included in surface loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoop: New surfaceLoop instance.
		
		"""
		
		newId=self.getNewId(self.surfaceLoops,Id)
		
		l=surfaceLoop(self,surfaceIDs,newId)
		self.surfaceLoops.append(l)
		
		return l
	
	def addRuledSurface(self,Id=None,lineLoopID=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` instance
		with given lineLoop. 
			
		Keyword Args:
			lineLoopID (ID): ID of line loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface: New ruledSurface instance.
		
		"""
		
		newId=self.getNewId(self.ruledSurfaces,Id)
		
		l=ruledSurface(self,lineLoopID,newId)
		self.ruledSurfaces.append(l)
		
		return l
	
	def addVolume(self,Id=None,surfaceLoopID=None):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.volume` instance
		with given surfaceLoop. 
			
		Keyword Args:
			surfaceLoopID (ID): ID of surface loop.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.volume: New volume instance.
		
		"""
		
		newId=self.getNewId(self.volumes,Id)
		
		l=volume(self,surfaceLoopID,newId)
		self.volumes.append(l)
		
		return l
	
	def addBoxField(self,Id=None,volSizeIn=10.,volSizeOut=20.,xRange=[],yRange=[],zRange=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.boxField` instance. 
			
		Keyword Args:
			Id (int): ID of field.
			volSizeIn (float): Mesh element volume inside box.
			volSizeOut (float): Mesh element volume outside box.
			xRange (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
			yRange (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
			zRange (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
			
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.boxField: New boxField instance.
		
		"""
		
		newId=self.getNewId(self.fields,Id)
		l=boxField(self,newId,volSizeIn=volSizeIn,volSizeOut=volSizeOut,xRange=xRange,yRange=yRange,zRange=zRange)
		self.fields.append(l)
		
		return l
	
	def addThresholdField(self,Id=None,IField=None,LcMin=5.,LcMax=20.,DistMin=30.,DistMax=60.):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.thresholdField` instance. 
			
		.. image:: ../imgs/pyfrp_gmsh_geometry/thresholdField.png
	
		Keyword Args:
			Id (int): ID of field.
			IField (int): ID of vertex that is center to threshold field.
			LcMin (float): Minimum volSize of threshold field.
			LcMax (float): Maximum volSize of threshold field.
			DistMin (float): Minimun density of field.
			DistMax (float): Maximum density of field.
				
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.thresholdField: New thresholdField instance.
			
		"""
			
		newId=self.getNewId(self.fields,Id)
		l=thresholdField(self,newId,IField=IField,LcMin=LcMin,LcMax=LcMax,DistMin=DistMin,DistMax=DistMax)
		self.fields.append(l)
		
		return l
	
	def addAttractorField(self,Id=None,NodesList=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.attractorField` instance. 
			
		Keyword Args:
			Id (int): ID of field.
			NodesList (list): List of IDs of the Nodes that attractor field centers around.
				
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.attractorField: New attractorField instance.
			
		"""
			
		newId=self.getNewId(self.fields,Id)
		l=attractorField(self,newId,NodesList=NodesList)
		self.fields.append(l)
		
		return l
	
	def addMinField(self,Id=None,FieldsList=[]):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.minField` instance. 
			
		Keyword Args:
			Id (int): ID of field.
			NodesList (list): List of IDs of the Nodes that attractor field centers around.
				
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.minField: New attractorField instance.
			
		"""
			
		newId=self.getNewId(self.fields,Id)
		l=minField(self,newId,FieldsList=FieldsList)
		self.fields.append(l)
		
		return l
	
	def addBoundaryLayerField(self,Id=None,AnisoMax=10000000000,hwall_n=1.,hwall_t=1,ratio=1.1,thickness=10.,hfar=1.,IntersectMetrics=1,Quads=0.):
		
		"""Adds new :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField` instance. 
		
		Keyword Args:
			Id (int): ID of field.
			AnisoMax (float): Threshold angle for creating a mesh fan in the boundary layer.
			IntersectMetrics (int): Intersect metrics of all faces.
			Quad (int): Generate recombined elements in the boundary layer.
			har (float): Element size far from the wall.
			hwall_n (float): Mesh Size Normal to the The Wall.
			hwall_t (float): Mesh Size Tangent to the Wall.
			ratio (float): Size Ratio Between Two Successive Layers.
			thickness (float): Maximal thickness of the boundary layer.
			List (list): List of field IDs.
				
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField: New boundaryLayerField instance.
			
		"""
		
		newId=self.getNewId(self.fields,Id)
		l=boundaryLayerField(self,newId,AnisoMax=AnisoMax,hwall_n=hwall_n,hwall_t=hwall_t,ratio=ratio,thickness=thickness,hfar=hfar,IntersectMetrics=IntersectMetrics,Quads=Quads)
		self.fields.append(l)
		
		return l
	
	def setAnnOffset(self,offset):
		
		"""Sets annotation offset for plotting.
		
		Args:
			offset (numpy.ndarray): New offset.
		
		"""
		
		self.annXOffset=offset[0]
		self.annYOffset=offset[1]
		self.annZOffset=offset[2]
		
	
	def writeToFile(self,fn):
		
		"""Writes domain to file.
		
		Args:
			fn (str): File path to write to.
			
		"""
		
		with open(fn,'w') as f:
			
			self.writeElements("vertices",f)
			self.writeElements("lines",f)
			self.writeElements("arcs",f)
			self.writeElements("bSplines",f)
			self.writeElements("lineLoops",f)
			self.writeElements("ruledSurfaces",f)
			self.writeElements("surfaceLoops",f)
			self.writeElements("volumes",f)
			self.writeElements("fields",f)

	def writeElements(self,element,f):
			
		"""Writes all entities of a specific element type to file.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* bSplines
			* ruledSurfaces
			* surfaceLoops
			* volumes
			* fields
		
		Args:
			element (str): Element type to write.
			f (file): File to write to.
			
		"""
	
		pyfrp_gmsh_IO_module.writeComment(f,element)
		for v in getattr(self,element):
			f=v.writeToFile(f)
		f.write("\n")
	
	def incrementAllIDs(self,offset):
		
		"""Adds offset to all entity IDs.
		
		Args:
			offset (int): Offset to be added.
			
		"""
		
		self.incrementIDs(offset,"vertices")
		self.incrementIDs(offset,"lines")
		self.incrementIDs(offset,"arcs")
		self.incrementIDs(offset,"bSplines")
		self.incrementIDs(offset,"lineLoops")
		self.incrementIDs(offset,"ruledSurfaces")
		self.incrementIDs(offset,"surfaceLoops")
		self.incrementIDs(offset,"volumes")
		
	def incrementIDs(self,offset,element):
		
		"""Adds offset to all entity IDs.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
			* fields
		
		Args:
			offset (int): Offset to be added.
			element (str): Element type to increment.
		
		"""
		
		for e in getattr(self,element):		
			e.Id=e.Id+offset
	
	def setDomainGlobally(self):
		
		"""Makes sure that ``self`` is domain for all 
		elements.
			
		"""
		
		self.setDomainForElementType("vertices")
		self.setDomainForElementType("lines")
		self.setDomainForElementType("arcs")
		self.setDomainForElementType("bSplines")
		self.setDomainForElementType("lineLoops")
		self.setDomainForElementType("ruledSurfaces")
		self.setDomainForElementType("surfaceLoops")
		self.setDomainForElementType("volumes")
		self.setDomainForElementType("fields")
		
		
	def setDomainForElementType(self,element):
		
		"""Makes sure that ``self`` is domain for all 
		elements of given type.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
			* fields
			
		Args:
			offset (int): Offset to be added.
			element (str): Element type to increment.
		
		"""
		
		for e in getattr(self,element):
			e.domain=self
		
	def getMaxID(self,element):
		
		"""Returns maximum ID for a specific element.
		
		Possible elements are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
		
		Args:
			element (str): Element type.
			
		Returns:
			int: Maximum ID.
		"""
		
		IDs=[]
		
		for e in getattr(self,element):
			IDs.append(e.Id)
		
		try:
			return max(IDs)
		except ValueError:
			0.
			
	def getAllMaxID(self):
		
		"""Returns maximum ID over all elements.
			
		Returns:
			int: Maximum ID.
		"""
		
		IDs=[]
		
		IDs.append(self.getMaxID("vertices"))
		IDs.append(self.getMaxID("lines"))
		IDs.append(self.getMaxID("arcs"))
		IDs.append(self.getMaxID("lineLoops"))
		IDs.append(self.getMaxID("ruledSurfaces"))
		IDs.append(self.getMaxID("surfaceLoops"))
		IDs.append(self.getMaxID("volumes"))
		
		return max(IDs)
	
	def getAllFieldsOfType(self,typ):
		
		"""Returns all fields of domain with specific typ.
		
		Returns:
			list: List of :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.field` objects.
		
		"""
		
		fs=[]
		for f in self.fields:
			if f.typ==typ:
				fs.append(f)
		return fs		
	
	def getBkgdField(self):
		
		"""Returns background field of domain.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.field: Background field.
			
		"""
		
		return self.bkgdField
	
	def hasBkgdField(self):
		
		"""Checks if domain already has a background field.
		
		Returns:
			bool: True if background field already exists.
		"""
	
		return self.bkgdField!=None
	
	def genMinBkgd(self,FieldsList=[]):
		
		"""Generates minimum field as background field.
		
		If domain already has minimum field, will take it and set it 
		as background field. If domain has multiple minimum fields, will take 
		the first one that appears in ``fields`` list.
		
		Keyword Args:
			FieldsList (list): List of field IDs included in minField.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.minField: Minimum field. 
		
		"""
		
		#Generate minField if not existent
		minFields=self.getAllFieldsOfType("min")
		if len(minFields)==0:
			if len(FieldsList)==0:
				FieldsList=pyfrp_misc_module.objAttrToList(self.fields,'Id')
			minField=self.addMinField(FieldsList=FieldsList)
		else:
			if self.hasBkgdField():
				if self.getBkgdField() in minFields:
					minField=self.getBkgdField()
				else:
					minField=minFields[0]
					minField.setAsBkgdField()
		
		if not self.hasBkgdField():
			minField.setAsBkgdField()
		
		return self.getBkgdField()
	
	def getAllObjectsWithProp(self,objName,attr,val):
		
		"""Filters all objects of type objName given attribute value.
		
		Possible objects names are:
		
			* vertices
			* lines
			* arcs
			* lineLoops
			* ruledSurfaces
			* surfaceLoops
			* volumes
			* fields
		
		.. note:: ``val`` can have any datatype.
		
		Args:
			objName (str): Name of object list.
			attr (str): Name of attribute.
			val (str): Value of attribute.
		
		Returns:
			list: List of objects that fulfill requirement.
		"""
		
		objects=getattr(self,objName)
		filteredObjects=pyfrp_misc_module.getAllObjWithAttrVal(objects,attr,val)
		
		return filteredObjects
	
	def simplifySurfaces(self,iterations=3,triangIterations=2,addPoints=False,fixSurfaces=True,debug=False):
		
		"""Tries to simplify surfaces inside the domain.
		
		Does this by:
			
			* For ``iterations`` iterations, do: 
				* Find all surfaces with the same normal vector.
				* Try to fuse this surfaces, see also :py:func:`pyfrp.modules.pyfrp_geometry_module.ruledSurface.fuse`.
				* Clean up edges via :py:func:`pyfrp.modules.pyfrp_geometry_module.domain.cleanUpUnusedEdges`.
				
			* Fixing loops via :py:func:`pyfrp.modules.pyfrp_geometry_module.domain.fixAllLoops`.
			* Fixing surfaces via :py:func:`pyfrp.modules.pyfrp_geometry_module.domain.fixAllSurfaces`.
		
		Keyword Args:
			iterations (int): Number of iterations used for simplification.
			triangIterations (int): Number of iterations used for subdivision of surfaces.
			addPoints (bool): Allow adding points inside surface triangles.
			fixSurfaces (bool): Allow fixing of surfaces, making sure they are coherent with Gmsh requirements.
			debug (bool): Print debugging messages.
		
		
		"""
		
		#Remember the stats from the start
		x=len(self.ruledSurfaces)
		y=len(self.lineLoops)
		z=len(self.edges)
		
		#Compute the normal of all surfaces
		for surface in self.ruledSurfaces:		
			surface.getNormal()
				
		#Loop through iterations
		for k in range(iterations):
				
			#Loop through surfaces
			for i,surface in enumerate(self.ruledSurfaces):
				
				#Get all surfaces with same normal vector
				sameNormal=self.getAllObjectsWithProp("ruledSurfaces","normal",surface.normal)
				sameNormal=sameNormal+self.getAllObjectsWithProp("ruledSurfaces","normal",-surface.normal)
				
				#Loop through all with same normal
				for j,sN in enumerate(sameNormal):
					if sN==surface:
						continue
					
					#Fuse
					if surface.fuse(sN,debug=debug):
						if debug:
							print("Successfully fused ", surface.Id, sN.Id)
							
			#Clean up edges
			self.cleanUpUnusedEdges(debug=debug)
			
		#Print some final statistics
		if debug:
			print("Surfaces: Before =" , x , " After:" , len(self.ruledSurfaces)) 
			print("lineLoops: Before =" , y , " After:" , len(self.lineLoops)) 
			print("Edges: Before =" , z , " After:" , len(self.edges)) 
		
		#raw_input()
		
		#Fix loops and surfaces
		self.fixAllLoops(debug=debug)
		if fixSurfaces:
			self.fixAllSurfaces(iterations=triangIterations,addPoints=addPoints)
			
	def cleanUpUnusedEdges(self,debug=False):
		
		"""Cleans up all unused edges in domain.
		
		See also: :py:func:`pyfrp.pyfrp_modules.pyfrp_gmsh_geometry.edge.delete`.
		
		Keyword Args:
			debug (bool): Print debugging output.
		
		"""
		
		for edge in self.edges:
			edge.delete(debug=debug)
	
	def fixAllLoops(self,debug=False):
		
		"""Tries to fix all loops in domain.
		
		See also: :py:func:`pyfrp.pyfrp_modules.pyfrp_gmsh_geometry.lineLoop.fix`.
		
		Keyword Args:
			debug (bool): Print debugging output.
		
		"""
		
		for loop in self.lineLoops:
			loop.fix()
	
	def fixAllSurfaces(self,debug=False,iterations=2,addPoints=False):
		
		"""Tries to fix all surfaces in domain.

		Does this by reiniating all ``lineLoop``. 
		
		See also: :py:func:`pyfrp.pyfrp_modules.pyfrp_gmsh_geometry.ruledSurface.initLineLoop`.
		
		Keyword Args:
			iterations (int): Number of iterations used for subdivision of surfaces.
			addPoints (bool): Allow adding points inside surface triangles.
			debug (bool): Print debugging messages.
				
		"""
		
		for surface in self.ruledSurfaces:
			surface.initLineLoop(surface.lineLoop.Id,debug=debug,iterations=iterations,addPoints=addPoints)
	
	def save(self,fn):
		
		"""Saves domain to pickle file.
		
		Args:
			fn (str): Output filename.
		"""
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)
	
	def merge(self,d):
		
		"""Merges domain d into this domain.
		
		Does this by:
		
			* Incrementing all IDs in ``d`` such that there is no overlap with ``self``.
			* Merging all element lists.
			* Making sure that all elements refer to ``self`` as domain.
			
		See also :py:func:`incrementAllIDs` and :py:func:`setDomainGlobally`.
		
		Args:
			d (pyfrp.modules.pyfrp_geometry_module.domain): Domain to merge.
		
		"""
		
		d.incrementAllIDs(self.getAllMaxID()+1)
		
		self.edges=self.edges+d.edges
		self.vertices=self.vertices+d.vertices
		self.arcs=self.arcs+d.arcs
		self.bSplines=self.bSplines+d.bSplines
		self.lines=self.lines+d.lines
		self.lineLoops=self.lineLoops+d.lineLoops
		self.ruledSurfaces=self.ruledSurfaces+d.ruledSurfaces
		self.surfaceLoops=self.surfaceLoops+d.surfaceLoops
		self.volumes=self.volumes+d.volumes
		self.fields=self.fields+d.fields
		
		self.setDomainGlobally()
	
	def removeDuplicates(self,debug=False):
		
		self.removeDuplicateEdgeIDs(debug=debug)
		self.removeDuplicateVerticesIDs()
		
	def removeDuplicateVerticesIDs(self):
		
		"""Checks if multiple vertices have the same ID and tries to remove one of them.
		
		Checks if vertices with same ID have the same coordinate. If so, remove all but one. Otherwise fixes
		index.
		
		Returns:
			list: Updated vertices list.
			
		"""
		
		# Loop through edges
		for i,v in enumerate(self.vertices):
			for j in range(i+1,len(self.vertices)):
				
				# Check if same ID
				if self.vertices[j].Id==self.vertices[i].Id:
					
					# Check if same start/end vertex
					if pyfrp_misc_module.compareVectors(self.vertices[j].x,self.vertices[i].x):			
						self.vertices.remove(self.vertices[j])
					else:
						newId=self.getNewId(self.vertices,None)
						self.vertices[j].setID(newId)
		
		return self.vertices
		
		
	def removeDuplicateEdgeIDs(self,debug=False):
		
		"""Checks if multiple edges have the same ID and tries to remove one of them.
		
		Checkss if edges with same ID have the same start and end vertex. If so, removes it all but one. Otherwise fixes
		index.
		
		Returns:
			list: Updated edges list.
			
		"""
		
		edgeIDs=pyfrp_misc_module.objAttrToList(self.edges,'Id')
		print(edgeIDs)
		# Loop through edges
		for i,e in enumerate(self.edges):
			for j in range(i+1,len(self.edges)):
				print(i,j,self.edges[i].Id,self.edges[j].Id)
				# Check if same ID
				if self.edges[j].Id==self.edges[i].Id:
					
					print("same id",self.edges[j].Id) 
					
					# Check if same start/end vertex
					if (self.edges[j].getFirstVertex()==self.edges[i].getFirstVertex()) and (self.edges[j].getLastVertex()==self.edges[i].getLastVertex()):			
						self.edges[j].delete(debug=debug)
					else:
						newId=self.getNewId(self.edges,None)
						self.edges[j].setID(newId)
		
		return self.edges
	

class gmshElement(object):
	
	def __init__(self,domain,Id):
		
		self.domain=domain
		self.Id=Id
		
	def getID(self):
		
		"""Returns ID of element.
		
		Returns:
			int: ID of element.
		
		"""
		
		return self.Id
		
	def setID(self,Id):
		
		"""Sets ID of element.
		
		Args:
			Id (int): New ID of element.
			
		Returns:
			int: New ID of element.
		
		"""
		
		self.Id=Id
		return self.Id
		
	def getCopy(self):
		
		"""Returns copy of element.
		
		Uses `copy.copy` to generate copy.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.gmshElement: Copy of element.
		
		"""
		
		return cpy.copy(self)
		
	def getType(self):
		
		"""Returns type of element.
		
		Returns:
			str: Type of element.
		
		"""
		
		t=str(type(self))
		t=t.split("'")[1]
		t=t.replace("pyfrp.modules.pyfrp_gmsh_geometry.","")
		
		return t
	
	def getTypeListName(self):
		
		"""Returns the element lists name.
		
		Returns:
			str: Name of element list.
		
		"""
		
		if self.getType()=="vertex":
			return "vertices"
		if self.getType()=="line":
			return "lines"
		if self.getType()=="arc":
			return "arcs"
		if self.getType()=="edge":
			return "edges"
		if self.getType()=="lineLoop":
			return "lineLoops"
		if self.getType()=="ruledSurface":
			return "ruledSurfaces"
		if self.getType()=="surfaceLoop":
			return "surfaceLoops"
		if self.getType()=="volume":
			return "volumes"
		if self.getType()=="field":
			return "fields"
		
	def getTypeList(self):
		
		"""Returns the element list of domain for this element.
		
		Returns:
			list: Element list.
		"""
		
		return self.getDomain().getattr(self.getTypeListName)
		
	
	def getDomain(self):
		
		"""Returns element's domain.
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: Element's domain.
		
		"""
		
		return self.domain
		
	def setDomain(self,d):	
		
		"""Sets element's domain.
		
		Args:
			d (pyfrp.modules.pyfrp_gmsh_geometry.domain): New domain
		
		Returns:
			pyfrp.modules.pyfrp_gmsh_geometry.domain: New domain.
		
		"""
		
		self.domain=d
		
		return self.domain
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return []
	
	def getAllSubElements(self,elements=[]):
		
		"""Finds all elements that are necessary to define this element recursively.
		
		Returns:
			list: List of elements.
		
		"""
		
		elements=list(elements)
		
		if len(self.getSubElements())==0:
			return elements
		else:	
			
			for el in self.getSubElements():
				elements.append(el)
				elements=el.getAllSubElements(elements=elements)
				
			return elements	
		
	def extract(self,d=None,strict=True,copy=False,debug=False):
		
		"""Extracts element and all elements necessary to define it.
		
		.. note:: If ``d`` is specified, then all extracted elements are inserted into ``d`` using
		   :py:func:`insertElement`.
		
		Keyword Args:
			d (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain to insert element
			copy (bool): Inserts copy of object.
			strict (bool): Don't allow IDs to be assigned to multiple elements.
			debug (bool): Print debugging output.
			
		Returns:
			list: List of elements.
		"""
		
		elmts=[self]+self.getAllSubElements()
		
		if d!=None:
			for el in elmts:
				d.insertElement('auto',el,strict=strict,copy=copy,debug=debug)
		
		return elmts
		
class vertex(gmshElement):
	
	"""Vertex class storing information from gmsh .geo Points.
	
	.. note:: ``volSize`` does not have any effect on the geometry itself but is simply 
		stored in the vertex object for further usage.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain vertex belongs to.
		x (numpy.ndarray): Coordinate of vertex.
		Id (int): ID of vertex.
			
	Keyword Args:
		volSize (float): Element size at vertex.
		
	"""
	
	def __init__(self,domain,x,Id,volSize=None):
		
		gmshElement.__init__(self,domain,Id)
		
		self.x=np.array(x)
		self.volSize=volSize
			
	def draw(self,ax=None,color=None,ann=None,backend="mpl",asSphere=True,size=10,render=False):
		
		"""Draws vertex.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of vertex.
			ann (bool): Show annotations.
			asSphere (bool): Draws vertex as sphere (only in vtk mode).
			size (float): Size of vertex (only in vtk mode).
			render (bool): Render in the end (only in vtk mode).
		
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		if backend=="mpl":
			ax=self.drawMPL(ax=ax,color=color,ann=ann)
		if backend=="vtk":
			ax=self.drawVTK(color=color,ann=ann,ax=ax,asSphere=asSphere,size=size,render=render)
		
		return ax
	
	def drawMPL(self,ax=None,color=None,ann=None):
		
		"""Draws vertrex into matplotlib axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot`.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of domain.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeGeometryPlot()
			ax=axes[0]
			
		ax.scatter(self.x[0],self.x[1],self.x[2],c=color)
		if ann:
			ax.text(self.x[0]+self.domain.annXOffset, self.x[1]+self.domain.annYOffset, self.x[2]+self.domain.annZOffset, "p"+str(self.Id), None)
			
		pyfrp_plot_module.redraw(ax)
		
		return ax
			
	def drawVTK(self,size=10,asSphere=True,ax=None,ann=None,color=[0,0,0],render=False):
		
		"""Draws vertrex into VTK renderer.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new ``vtkRenderer``, 
		   see also :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		Keyword Args:
			ax (vtk.vtkRenderer): Renderer to draw in.
			color (str): Color of vertex.
			ann (bool): Show annotations.
			asSphere (bool): Draws vertex as sphere.
			size (float): Size of vertex.
			render (bool): Render in the end.
				
		Returns:
			vtk.vtkRenderer: Updated renderer.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			ax,renderWindow,renderWindowInteractor=pyfrp_vtk_module.makeVTKCanvas()

		pyfrp_vtk_module.drawVTKPoint(self.x,asSphere=asSphere,color=color,size=size,renderer=ax)
		
		if ann:
			printWarning("Annotations don't properly work with backend=vtk .")
			pyfrp_vtk_module.drawVTKText("p"+str(self.Id),[self.x[0]+self.domain.annXOffset, self.x[1]+self.domain.annYOffset, self.x[2]+self.domain.annZOffset],renderer=ax)
		
		if render:
			ax=pyfrp_vtk_module.renderVTK(ax,start=False)
		
		return ax
		
	def setX(self,x):
		
		"""Sets coordinate if vertex to ``x``.
		
		Returns:
			numpy.ndarray: New vertex coordinate.
		"""
		
		self.x=x
		return self.x
	
	def writeToFile(self,f):
		
		"""Writes vertex to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Point("+str(self.Id)+")= {" + str(self.x[0]) + ","+ str(self.x[1])+ "," + str(self.x[2]) + ',' + str(self.volSize) + "};\n" )
		
		return f
	
	def addToAttractor(self,attrField=None,LcMin=5.,LcMax=20.,DistMin=30.,DistMax=60.):
		
		"""Adds vertex to a attractor field. 
		
		If no field is given, will create new one with given parameters. Will also create 
		a new threshhold field around attractor and add fields to minField. If no minField exists,
		will create a new one too and set it as background field.
		
		See also :py:func:`addAttractorField`, :py:func:`addThresholdField`, :py:func:`addMinField` and :py:func:`genMinBkgd`.
		
		Keyword Args:
			attrField (pyfrp.modules.pyfrp_gmsh_geometry.attractorField): Attractor field object.
			LcMin (float): Minimum volSize of threshold field.
			LcMax (float): Maximum volSize of threshold field.
			DistMin (float): Minimun density of field.
			DistMax (float): Maximum density of field.
			
		Returns:	
			pyfrp.modules.pyfrp_gmsh_geometry.attractorField: Attractor field around vertex.
		
		"""
		
		#Generate attractor field if not given
		if attrField==None:
			attrField=self.domain.addAttractorField(NodesList=[self.Id])
		else:
			attrField.addNodeByID(self.Id)
		
		#Generate threshhold field if not already existent
		threshFields=attrField.includedInThresholdField()
		if len(threshFields)==0:
			threshField=self.domain.addThresholdField(IField=attrField.Id,LcMin=LcMin,LcMax=LcMax,DistMin=DistMin,DistMax=DistMax)
		else:
			threshField=threshFields[0]
		
		self.domain.genMinBkgd(FieldsList=[threshField.Id])
		
		return attrField	
	
	def addToBoundaryLayer(self,boundField=None,**fieldOpts):
		
		"""Adds vertex to a boundary layer field. 
		
		If no field is given, will create new one with given parameters and add it to a minField. If no minField exists,
		will create a new one too and set it as background field.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addBoundaryLayerField`
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addMinField` and 
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.genMinBkgd`.
		
		Keyword Args:
			boundField (pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField): Boundary layer field object.
			fieldOpts (dict): See documentation of boundary layer field of all available options.
			
		Returns:	
			pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField: Boundary layer field around vertex.
		
		"""
		
		#Generate attractor field if not given
		if boundField==None:
			boundField=self.domain.addBoundaryLayerField()
		
		#Add Vertex
		boundField.addNodeByID(self.Id)
		
		#Set options
		boundField.setFieldAttributes(**fieldOpts)
		
		#Generate background field
		self.domain.genMinBkgd(FieldsList=[boundField.Id])
			
		return boundField	
	
class edge(gmshElement):
	
	"""Edge class storing information from gmsh .geo circles and lines.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain vertex belongs to.
		Id (int): ID of edge.
		typ (int): Type of edge (1=arc/0=line).
			
	"""
	
	def __init__(self,domain,Id,typ):
		
		gmshElement.__init__(self,domain,Id)
		self.typ=typ
		
	def getDomain(self):
		
		"""Returns domain edge belongs to."""
		
		return self.domain
	
	def getTyp(self):
		
		"""Returns Type of edge."""
		
		return self.typ
	
	def decodeTyp(self):
		
		"""Decodes type of edge into string."""
		
		if typ==1:
			return "arc"
		elif typ==0:
			return "line"
		elif typ==2:
			return "bSpline"
	
	def addToBoundaryLayer(self,boundField=None,**fieldOpts):
		
		"""Adds edge to a boundary layer field. 
		
		If no field is given, will create new one with given parameters and add it to a minField. If no minField exists,
		will create a new one too and set it as background field.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addBoundaryLayerField`
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addMinField` and 
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.genMinBkgd`.
		
		Keyword Args:
			boundField (pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField): Boundary layer field object.
			fieldOpts (dict): See documentation of boundary layer field of all available options.
			
		Returns:	
			pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField: Boundary layer field around edge.
		
		"""
		
		#Generate attractor field if not given
		if boundField==None:
			boundField=self.domain.addBoundaryLayerField()
		
		#Add Vertex
		boundField.addEdgeByID(self.Id)
		
		#Set options
		boundField.setFieldAttributes(**fieldOpts)
		
		#Generate background field
		self.domain.genMinBkgd(FieldsList=[boundField.Id])
			
		return boundField	
	
	def includedInLoop(self):
		
		"""Checks if edge is included in a loop.
		
		Returns:
			tuple: Tuple containing:
			
				* included (bool): True if included.
				* loops (list): List of :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.lineLoop` objects that include edge. 
			
		"""
		
		loops=[]
		
		for i,loop in enumerate(self.domain.lineLoops):
			if self in loop.edges:
				loops.append(loop)
		return len(loops)>0,loops
	
	def includedInField(self):
		
		"""Checks if edge is included in a field.
		
		.. note:: Only checks for boundary layer fields, since they are the only ones who can evolve around edge.
		
		Returns:
			tuple: Tuple containing:
			
				* included (bool): True if included.
				* fields (list): List of :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.fields` objects that include edge. 
			
		"""
		
		fields=[]
		
		for field in self.domain.fields:
			if field.typ=="boundaryLayer":
				if self in field.EdgesList:
					fields.append(fields)
		return len(fields)>0,fields
	
	def delete(self,debug=False):
		
		"""Deletes edge if it is not used in any loop or field.
		
		Returns:
			bool: True if deletion was successful.
		"""
		
		incl,loops=self.includedInLoop()
		if incl:
			if debug:
				printWarning("Was not able to delete edge with ID " + str(self.Id) +". Still part of loops" + str(pyfrp_misc_module.objAttrToList(loops,'Id')) + " .")
			return False
		
		incl,fields=self.includedInField()
		if incl:
			if debug:
				printWarning("Was not able to delete edge with ID " + str(self.Id) +". Still part of field with ID " + str(pyfrp_misc_module.objAttrToList(fields,'Id')) )
			return False
		
		try:
			if self.typ==0:
				self.domain.lines.remove(self)
			if self.typ==1:
				self.domain.arcs.remove(self)
			if self.typ==2:
				self.domain.bSplines.remove(self)
				
			self.domain.edges.remove(self)
		except ValueError:
			if debug:
				printWarning("Could not remove edge " + str(self.Id)+" from elements list. Already seems to be removed.")
			return False
		
		return True
	
class line(edge):
	
		
	"""Line class storing information from gmsh .geo lines.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain line belongs to.
		v1 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
		v2 (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
		Id (int): ID of line.
		
	"""
	
	
	def __init__(self,domain,v1,v2,Id):
		
		edge.__init__(self,domain,Id,0)

		self.v1=v1
		self.v2=v2
	
	def getMiddle(self):
		
		r"""Returns midpoint of line.
		
		.. math:: m = \frac{x(v_1) + x(v_2)}{2}
		     
		Returns:
			numpy.ndarray: Midpoint.
			
		"""
		
		return (self.v1.x+self.v2.x)/2.
	
	def draw(self,ax=None,color=None,ann=None,backend="mpl",render=False,drawVertices=False,linewidth=1):
			
		"""Draws line.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end (only in vtk mode).
			drawVertices (bool): Also draw vertices.
			
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		if backend=="mpl":
			ax=self.drawMPL(ax=ax,color=color,ann=ann,linewidth=linewidth)
		if backend=="vtk":
			ax=self.drawVTK(color=color,ann=ann,ax=ax,render=render,linewidth=linewidth)
		
		if drawVertices:
			ax=self.v1.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
			ax=self.v2.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
			
		return ax
			
	def drawMPL(self,ax=None,color=None,ann=None,linewidth=1):
		
		"""Draws line into matplotlib axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot`.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeGeometryPlot()
			ax=axes[0]
		ax.plot([self.v1.x[0],self.v2.x[0]],[self.v1.x[1],self.v2.x[1]],zs=[self.v1.x[2],self.v2.x[2]],color=color,linestyle='-',linewidth=linewidth)
		if ann:
			m=self.getMiddle()
			ax.text(m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset, "l"+str(self.Id), None)
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def drawVTK(self,ax=None,color=None,ann=None,render=False,linewidth=1):
		
		"""Draws line into VTK renderer.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new ``vtkRenderer``, 
		   see also :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		See also :py:func:`pyfrp.modules.pyfrp_vtk_module.drawVTKLine`.
		
		Keyword Args:
			ax (vtk.vtkRenderer): Renderer to draw in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end.
				
		Returns:
			vtk.vtkRenderer: Updated renderer.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			ax,renderWindow,renderWindowInteractor=pyfrp_vtk_module.makeVTKCanvas()
			
		pyfrp_vtk_module.drawVTKLine(self.v1.x,self.v2.x,color=color,renderer=ax,linewidth=linewidth)
		
		if ann:
			printWarning("Annotations don't properly work with backend=vtk .")
			m=self.getMiddle()
			pyfrp_vtk_module.drawVTKText("p"+str(self.Id),[m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset],renderer=ax)
		
		if render:
			ax=pyfrp_vtk_module.renderVTK(ax,start=False)
		
		return ax
		
	def getLastVertex(self,orientation):
		
		"""Returns last vertex of line given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of line.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.v2
		elif orientation==-1:
			return self.v1
		else:
			printError("Cannot return last vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getFirstVertex(self,orientation):
		
		"""Returns first vertex of line given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of line.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.v1
		elif orientation==-1:
			return self.v2
		else:
			printError("Cannot return first vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def writeToFile(self,f):
		
		"""Writes line to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Line("+str(self.Id)+")= {" + str(self.v1.Id) + "," + str(self.v2.Id) + "};\n" )
		
		return f
	
	def getDirection(self,orientation):
		
		"""Returns direction of line.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of line.
			
		Returns:
			numpy.ndarray: Direction of line.
		
		"""
			
		return self.getLastVertex(orientation).x-self.getFirstVertex(orientation).x
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.getFirstVertex(1),self.getLastVertex(1)]
	
class arc(edge):
	
	"""Arc class storing information from gmsh .geo cicle.
	
	Will compute ``angleOffset``, ``angle`` and ``pOffset`` on creation.
	
	.. image:: ../imgs/pyfrp_gmsh_geometry/arc.png
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain arc belongs to.
		vstart (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Start vertex.
		vcenter (pyfrp.modules.pyfrp_gmsh_geometry.vertex): Center vertex.
		vend (pyfrp.modules.pyfrp_gmsh_geometry.vertex): End vertex.
		Id (int): ID of arc.
		
	"""		
	
	def __init__(self,domain,vstart,vcenter,vend,Id):
		
		edge.__init__(self,domain,Id,1)
		
		self.vcenter=vcenter
		self.vstart=vstart
		self.vend=vend
		
		self.radius=self.computeRadius()
		
		self.pOffset=self.computePOffset()
		
		self.angleOffset=self.computeAngleOffset()
		self.angle=self.computeAngle()
		
	

	def computeAngleOffset(self):
		
		"""Computes and returns offset angle of arc.
		"""
		
		self.angleOffset=pyfrp_geometry_module.getAngle(self.pOffset,self.vstart.x-self.vcenter.x)
		
		return self.angleOffset
	
	def computeAngle(self):
		
		"""Computes and returns angle of arc.
		"""
		
		self.angle=pyfrp_geometry_module.getAngle(self.vstart.x-self.vcenter.x,self.vend.x-self.vcenter.x)
		return self.angle
	
	def computePOffset(self):
		
		"""Computes and returns offset point of arc.
		"""
		
		v1n,v2nb = self.getNormVec()
		
		self.pOffset=self.radius*v2nb
		self.pOffset=self.pOffset/np.linalg.norm(self.pOffset)
		
		return self.pOffset
	
	def getNormVec(self):
		
		"""Computes and returns vectors normal to arc.
		
		Returns:
			tuple: Tuple containing:
				
				* v1n (numpy.ndarray): Normal vector to ``vstart-vcenter``.
				* v2n (numpy.ndarray): Normal vector to ``vend-vcenter``.
					
		"""
		
		v1=self.vstart.x-self.vcenter.x
		v2=self.vend.x-self.vcenter.x
		
		self.v1n = v1/np.linalg.norm(v1)
		v2n = v2/np.linalg.norm(v2)
		v2nb = v2n-np.dot(v2n,self.v1n)*self.v1n
		
		self.v2nb = v2nb/np.linalg.norm(v2nb)
		
		return self.v1n,self.v2nb
	
	def getPlotVec(self):
		
		"""Returns vectors for plotting arc.
		
		Returns:
			tuple: Tuple containing:
				
				* x (numpy.ndarray): x-array.
				* y (numpy.ndarray): y-array.
				* z (numpy.ndarray): z-array.
						
		"""
		
		self.getNormVec()
			
		if np.mod(self.angle,np.pi/2.)<0.01:
			a = np.linspace(0,self.angle,1000)
		else:
			a = np.linspace(self.angleOffset-self.angle,self.angleOffset,1000)
				
		x,y,z=self.getPointOnArc(a)
		
		return x,y,z
	
	def getPointOnArc(self,a):
		
		"""Returns point on arc at angle ``a``.
		
		Returns:
			tuple: Tuple containing:
				
				* x (float): x-coordinate.
				* y (float): y-coordinate.
				* z (float): z-coordinate.
						
		"""
			
		x = self.vcenter.x[0]+np.sin(a)*self.radius*self.v1n[0]+np.cos(a)*self.radius*self.v2nb[0]
		y = self.vcenter.x[1]+np.sin(a)*self.radius*self.v1n[1]+np.cos(a)*self.radius*self.v2nb[1]
		z = self.vcenter.x[2]+np.sin(a)*self.radius*self.v1n[2]+np.cos(a)*self.radius*self.v2nb[2]	
		
		return x,y,z

	def computeRadius(self):
		
		"""Computes and returns radius of arc.
		
		Returns:
			float: Radius of arc.
		"""
		
		self.radius=np.linalg.norm(self.vstart.x-self.vcenter.x)
		return self.radius
		
	def inArc(self,x,debug=False):
		
		"""Tells if coordinate ``x`` is on arc or not.
		
		Returns:
			bool: ``True`` if on arc, ``False`` otherwise.
		"""
		
		a=self.computeAngle(array([self.radius,0])-self.vcenter.x,x-self.vcenter.x)
				
		if np.mod(a,2*np.pi)<self.angle+self.angleOffset and self.angleOffset<=np.mod(a,2*np.pi):
			return True
		else:
			return False
	
	def getRadius(self):
		
		"""Returns radius of arc."""
		
		return self.radius
	
	def getAngle(self):
		
		"""Returns angle of arc."""
		
		return self.angle
	
	def getAngleOffset(self):
		
		"""Returns offset angle of arc."""
		
		return self.angleOffset
	
	def getVstart(self):
		
		"""Returns start vertex of arc."""
		
		return self.vstart
	
	def getVend(self):
		
		"""Returns end vertex of arc."""
		
		return self.vend
	
	def getXstart(self):
		
		"""Returns start coordinate of arc."""
		
		return self.vstart.x
	
	def getXend(self):
		
		"""Returns end coordinate of arc."""
		
		return self.vend.x
	
	def getVcenter(self):
		
		"""Returns center vertex of arc."""
		
		return self.vcenter
	
	def getXcenter(self):
		
		"""Returns center coordinate of arc."""
		
		return self.vcenter.x
	
	def draw(self,ax=None,color=None,ann=None,backend="mpl",render=False,drawVertices=True,linewidth=1):
		
		"""Draws arc.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end (only in vtk mode).
			drawVertices (bool): Also draw vertices.
			
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		if backend=="mpl":
			ax=self.drawMPL(ax=ax,color=color,ann=ann,linewidth=linewidth)
		if backend=="vtk":
			ax=self.drawVTK(color=color,ann=ann,ax=ax,render=render,linewidth=linewidth)
		
		if drawVertices:
			ax=self.vcenter.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
			ax=self.vstart.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
			ax=self.vend.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
			
		return ax
		
	def drawMPL(self,ax=None,color=None,ann=None,render=False,linewidth=1):
		
		"""Draws arc into matplotlib axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot`.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeGeometryPlot()
			ax=axes[0]
			
		x,y,z=self.getPlotVec()
		
		ax.plot(x,y,zs=z,color=color,linestyle='-',linewidth=linewidth)
		
		if ann:
			x,y,z=self.getPointOnArc(self.angle/2.)
			ax.text(x+self.domain.annXOffset, y+self.domain.annYOffset, z+self.domain.annZOffset, "c"+str(self.Id), None)
				
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def drawVTK(self,ax=None,color=None,ann=None,render=False,linewidth=1):
		
		"""Draws arc into VTK renderer.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new ``vtkRenderer``, 
		   see also :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		See also :py:func:`pyfrp.modules.pyfrp_vtk_module.drawVTKArc`.
		
		Keyword Args:
			ax (vtk.vtkRenderer): Renderer to draw in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end.
			
		Returns:
			vtk.vtkRenderer: Updated renderer.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			ax,renderWindow,renderWindowInteractor=pyfrp_vtk_module.makeVTKCanvas()	
		
		pyfrp_vtk_module.drawVTKArc(self.vstart.x,self.vcenter.x,self.vend.x,color=color,renderer=ax,linewidth=linewidth)
		
		if ann:
			printWarning("Annotations don't properly work with backend=vtk .")
			m=self.getMiddle()
			pyfrp_vtk_module.drawVTKText("p"+str(self.Id),[m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset],renderer=ax)
		
		if render:
			ax=pyfrp_vtk_module.renderVTK(ax,start=False)
		
		return ax
	
	def getLastVertex(self,orientation):
		
		"""Returns last vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.getVend()
		elif orientation==-1:
			return self.getVstart()
		else:
			printError("Cannot return last vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getFirstVertex(self,orientation):
		
		"""Returns first vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==-1:
			return self.getVend()
		elif orientation==1:
			return self.getVstart()
		else:
			printError("Cannot return first vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def writeToFile(self,f):
		
		"""Writes arc to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Circle("+str(self.Id)+")= {" + str(self.vstart.Id) + ","+ str(self.vcenter.Id)+ "," + str(self.vend.Id) + "};\n" )
		
		return f
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.vcenter,self.vstart,self.vend]
	
class bSpline(edge):
	
	"""Bspline class storing information from gmsh .geo BSpline.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain arc belongs to.
		vertices (list): List of vertex objects.
		Id (int): ID of spline.
		
	"""		
	
	def __init__(self,domain,vertices,Id):
		
		edge.__init__(self,domain,Id,2)
		
		self.initVertices(vertices)
	
	def initVertices(self,vertices):
		
		"""Initiates list of vertices.
		
		If vertex is given by Id, will use :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.getVertexById`
		to identify vertex. 
		
		Args:
			vertices (list): List of vertex objects.
			
		Returns:
			list: List of vertex objects.
		
		"""
		
		self.vertices=[]
		for v in vertices:
			if isinstance(v,int):
				self.vertices.append(self.domain.getVertexById(v)[0])
			else:
				self.vertices.append(v)
				
		return self.vertices		
	
	def writeToFile(self,f):
		
		"""Writes bSpline to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("BSpline("+str(self.Id)+")= {" )
			
		for i,v in enumerate(self.vertices):
			f.write(str(v.Id))
			if i!=len(self.vertices)-1:
				f.write(",")
			else:
				f.write("};\n")
	
		return f
	
	def getMiddle(self):
		
		r"""Returns midpoint of bSpline.
		
		Midpoint in this case is defined as the coordinate of the mid vertex.
		
		Returns:
			numpy.ndarray: Midpoint.
			
		"""
	
		return self.vertices[int(np.floor(len(self.vertices)/2.))].x
	
	def draw(self,ax=None,color=None,ann=None,backend="mpl",render=False,drawVertices=False,linewidth=1):
			
		"""Draws spline.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end (only in vtk mode).
			drawVertices (bool): Also draw vertices.
			
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		printWarning("Spline drawing currently just draws lines inbetween interpolation points.")
		
		if backend=="mpl":
			ax=self.drawMPL(ax=ax,color=color,ann=ann,linewidth=linewidth)
		if backend=="vtk":
			ax=self.drawVTK(color=color,ann=ann,ax=ax,render=render,linewidth=linewidth)
		
		if drawVertices:
			for v in self.vertices:
				ax=v.draw(ax=ax,color=color,ann=ann,backend=backend,render=render,asSphere=False)
				
		return ax
			
	def drawMPL(self,ax=None,color=None,ann=None,linewidth=1):
		
		"""Draws spline into matplotlib axes.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot`.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line.
			ann (bool): Show annotations.
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
			fig,axes = pyfrp_plot_module.makeGeometryPlot()
			ax=axes[0]
		
		for i in range(len(self.vertices)-1):
			ax.plot([self.vertices[i].x[0],self.vertices[i+1].x[0]],[self.vertices[i].x[1],self.vertices[i+1].x[1]],zs=[self.vertices[i].x[2],self.vertices[i+1].x[2]],color=color,linestyle='-',linewidth=linewidth)
		
		if ann:
			m=self.getMiddle()
			ax.text(m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset, "l"+str(self.Id), None)
		
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def drawVTK(self,ax=None,color=None,ann=None,render=False,linewidth=1):
		
		"""Draws spline into VTK renderer.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new ``vtkRenderer``, 
		   see also :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		See also :py:func:`pyfrp.modules.pyfrp_vtk_module.drawVTKLine`.
		
		Keyword Args:
			ax (vtk.vtkRenderer): Renderer to draw in.
			color (str): Color of line.
			ann (bool): Show annotations.
			render (bool): Render in the end.
				
		Returns:
			vtk.vtkRenderer: Updated renderer.
		
		"""
		
		if ann==None:
			ann=False
		
		if ax==None:
				
			ax,renderWindow,renderWindowInteractor=pyfrp_vtk_module.makeVTKCanvas()
		
		pyfrp_vtk_module.drawVTKPolyLine(pyfrp_misc_module.objAttrToList(self.vertices,'x'),color=color,renderer=ax,linewidth=linewidth)
		
		if ann:
			printWarning("Annotations don't properly work with backend=vtk .")
			m=self.getMiddle()
			pyfrp_vtk_module.drawVTKText("p"+str(self.Id),[m[0]+self.domain.annXOffset, m[1]+self.domain.annYOffset, m[2]+self.domain.annZOffset],renderer=ax)
		
		if render:
			ax=pyfrp_vtk_module.renderVTK(ax,start=False)
		
		return ax
		
	def getLastVertex(self,orientation):
		
		"""Returns last vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.vertices[-1]
		elif orientation==-1:
			return self.vertices[0]
		else:
			printError("Cannot return last vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getFirstVertex(self,orientation):
		
		"""Returns first vertex of arc given a orientation.
		
		Orientation can be either forward (1), or reverse (-1).
		
		Args:
			orientation (int): Orientation of arc.
			
		Returns:
			pyfrp.pyfrp_gmsh_geometry.vertex: Vertex.
		
		"""
		
		if orientation==1:
			return self.vertices[0]
		elif orientation==-1:
			return self.vertices[-1]
		else:
			printError("Cannot return first vertex. Orientation " + str(orientation) + " unknown.")
			return None
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return self.vertices
	
class lineLoop(gmshElement):
	
	"""Lineloop class storing information from gmsh .geo.

	Object has two major attributes:
	
		* edges (list): List of pyfrp.moduels.pyfrp_gmsh_geometry.edge objects.
		* orientations (list): List of orientations of each element, either ``1`` or ``-1`` 
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain loop belongs to.
		edgeIDs (list): List of edge IDs.
		Id (int): ID of loop.
		
	"""		
	
	def __init__(self,domain,edgeIDs,ID):
		
		gmshElement.__init__(self,domain,ID)
		
		self.edges,self.orientations=self.initEdges(edgeIDs)
		
	def initEdges(self,IDs):
		
		"""Constructs ``edges`` and ``orientations`` list at object initiations 
		from list of IDs.
		
		Args:
			IDs (list): List of IDs
			
		Returns:
			tuple: Tuple containing:
			
				* edges (list): List of pyfrp.moduels.pyfrp_gmsh_geometry.edge objects.
				* orientations (list): List of orientations of each element, either ``1`` or ``-1`` 
		
		"""
		
		self.edges=[]
		self.orientations=[]
		
		for ID in IDs:
			self.addEdgeByID(ID)
		
		return self.edges,self.orientations
		
	def addEdgeByID(self,ID):
		
		"""Adds edge to lineloop.
		
		Args:
			ID (int): ID of edge to be added.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		
		self.edges.append(self.domain.getEdgeById(abs(ID))[0])
		self.orientations.append(np.sign(ID))
		
		return self.edges
	
	def insertEdgeByID(self,ID,pos):
		
		"""Inserts edge to lineloop at position.
		
		Args:
			ID (int): ID of edge to be inserted.
			pos (int): Position at which ID to be inserted.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		self.edges.insert(pos,self.domain.getEdgeById(abs(ID))[0])
		self.orientations.insert(pos,np.sign(ID))
		
		return self.edges
	
	def removeEdgeByID(self,ID):
		
		"""Remove edge from lineloop.
		
		Args:
			ID (int): ID of edge to be removed.
			
		Returns:
			list: Updated edgeIDs list.
			
		"""
		
		idx=self.edges.index(abs(ID))
		self.edges.remove(abs(ID))
		self.orientations.pop(idx)
		
		return self.edges
	
	def reverseEdge(self,ID):
		
		"""Reverses the orientation of an edge in the line loop.
		
		Args:
			ID (int): ID of edge to be reversed.
		
		Returns:
			list: Updated orientations list.
			
		"""
		
		e=self.domain.getEdgeById(abs(ID))[0]
		self.orientations[self.edges.index(e)]=-self.orientations[self.edges.index(e)]
		
		return self.orientations
	
	def draw(self,ax=None,color='k',ann=None,backend='mpl',drawVertices=False,linewidth=1):
		
		"""Draws complete line loop.
		
		There are two different backends for drawing, namely 
			
			* Matplotlib (``backend='mpl'``)
			* VTK (``backend='vtk'``)
		
		Matplotlib is easier to handle, but slower. VTK is faster for complex
		geometries.
		
		.. note:: If ``backend=mpl``, ``ax`` should be a ``matplotlib.axes``, if ``backend='vtk'``,
		   ``ax`` should be a ``vtk.vtkRenderer`` object.
		
		.. note:: If no axes is given, will create new one,
		   see also :py:func:`pyfrp.modules.pyfrp_plot_module.makeGeometryPlot` 
		   or :py:func:`pyfrp.modules.pyfrp_vtk_module.makeVTKCanvas`.
		
		.. warning:: Annotations are not properly working with ``backend='vtk'``.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of line loop.
			ann (bool): Show annotations.
			drawVertices (bool): Also draw vertices.
		
		Returns:
			matplotlib.axes: Updated axes.
			
		"""
		
		for e in self.edges:
			ax=e.draw(ax=ax,color=color,ann=ann,backend=backend,drawVertices=drawVertices,linewidth=linewidth)
		
		return ax
	
	def printLoop(self):
		
		"""Prints loop.
		"""
		
		ids=np.array(pyfrp_misc_module.objAttrToList(self.edges,"Id"))
		orients=np.array(self.orientations)
		
		print("Line Loop with ID = "+ str(self.Id)+": "+str(ids*orients))
	
	def fix(self):
		
		"""Fixes loop.
		"""
		
		edgesNew=[self.edges[0]]
		orientationsNew=[self.orientations[0]]
		
		for i in range(1,len(self.edges)):
			
			lastEdge=edgesNew[i-1]
			vLast=lastEdge.getLastVertex(orientationsNew[i-1])
			
			for j in range(len(self.edges)):
				
				currEdge=self.edges[j]
				currOrient=self.orientations[j]
						
				if currEdge==lastEdge:
					continue
				
				if vLast == currEdge.getFirstVertex(currOrient):
					edgesNew.append(currEdge)
					orientationsNew.append(currOrient)
					break 
				
				elif vLast == currEdge.getLastVertex(currOrient):
					edgesNew.append(currEdge)
					orientationsNew.append(-currOrient)
					break
					
				if j==len(self.edges)-1:
					printWarning("Could not fix loop with ID" + str(self.Id))
					print("Edge with ID " +str(lastEdge.Id) + " is not matching with any other edge.")
					return False
				
		self.edges=edgesNew
		self.orientations=orientationsNew
		
		return True
		
	def checkClosed(self,fix=False,debug=False):
		
		"""Checks if lineLoop is closed.
		
		Keyword Args:
			debug (bool): Print debugging messages.
			fix (bool): Close if necessary.
			
		Returns:
			bool: True if closed.
		
		"""
		
		b=True
		
		for i in range(len(self.edges)):
			
			#Get ID of edge
			edge1Temp=self.edges[i]
			orient1=self.orientations[i]
			
			#Get ID of next edge
			edge2Temp=self.edges[pyfrp_misc_module.modIdx(i+1,self.edges)]
			orient2=self.orientations[pyfrp_misc_module.modIdx(i+1,self.edges)]
			
			#Get ID of first/last vertex
			firstVertexId=edge1Temp.getFirstVertex(orient1).Id
			lastVertexId=edge2Temp.getLastVertex(orient2).Id
			
			#Check if vertices are matching
			if firstVertexId!=lastVertexId:
				b=False
				
				if fix:
					self.reverseEdge(edge2Temp.Id)
					b=True
			
					if debug:
						print("Edge with ID " +str(edge1Temp.Id) + " was not matching edge with ID " + str(edge2Temp.Id) + ". \n Fixed this.")
				
					
				if debug:
					printWarning("lineLoop with ID " + str(self.Id) + " does not close." )
					print("Edge with ID " +str(edge1Temp.Id) + " is not matching edge with ID " + str(edge2Temp.Id))
						
		return b
	
	def writeToFile(self,f):
		
		"""Writes line loop to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Line Loop("+str(self.Id)+")= {" )
			
		for i,s in enumerate(self.edges):
			f.write(str(self.orientations[i]*s.Id))
			if i!=len(self.edges)-1:
				f.write(",")
			else:
				f.write("};\n")
	
		return f
	
	def getVertices(self):
		
		"""Returns all vertices included in loop."""
		
		vertices=[]
		for i,edge in enumerate(self.edges):
			vertices.append(edge.getFirstVertex(self.orientations[i]))
		return vertices
	
	def hasCommonEdge(self,loop):
		
		"""Checks if lineLoop has common edges with other lineLoop.
		
		Args:
			loop (pyfrp.modules.pyfrp_gmsh_geometry.lineLoop): lineLoop object.
		
		Returns:
			tuple: Tuple containing:
			
				* hasCommon (bool): True if loops have common edge.
				* edges (list): List of common edges.
						
		"""
		
		edges=[]
		for e in self.edges:
			if e in loop.edges:
				edges.append(e)
		
		return len(edges)>0,edges
	
	def fuse(self,loop,maxL=1000,debug=False,surface=None):
		
		"""Fuses lineLoop with other loop.
		
		"""
		
		#Find common edge
		b,commonEdges=self.hasCommonEdge(loop)
		
		#print "commonEdges", pyfrp_misc_module.objAttrToList(commonEdges,'Id')
		
		if not b:
			printWarning("Cannot fuse lineLoop with ID " + str(self.Id) + " and lineLoop with ID "+ str(loop.Id) +" . Loops do not have common edge.")
			return False

		#Sort edges of loop and pop edges that are in common
		idx=loop.edges.index(commonEdges[0])
		idxLast=loop.edges.index(commonEdges[-1])+1
		
		edgeTemp1,loop.edges=pyfrp_misc_module.popRange(loop.edges,idx,idxLast)
		orientTemp1,loop.orientations=pyfrp_misc_module.popRange(loop.orientations,idx,idxLast)
		
		#print "popped ",pyfrp_misc_module.objAttrToList(edgeTemp1,'Id')
		#print "remain ",pyfrp_misc_module.objAttrToList(loop.edges,'Id')
		
		edges=list(np.roll(loop.edges,len(loop.edges)-idx))
		orientations=list(np.roll(loop.orientations,len(loop.edges)-idx))
		
		#Pop common edge out of this loop
		idx=self.edges.index(commonEdges[0])
		idxLast=self.edges.index(commonEdges[-1])+1
		edgeTemp2,self.edges=pyfrp_misc_module.popRange(self.edges,idx,idxLast)
		orientTemp2,self.orientations=pyfrp_misc_module.popRange(self.orientations,idx,idxLast)
		
		#print "popped 2 ",pyfrp_misc_module.objAttrToList(edgeTemp2,'Id')
		#print "remain ",pyfrp_misc_module.objAttrToList(self.edges,'Id')
		
		#Figure out if edges of other loop are in right order or need to be reversed
		if orientTemp1[0]==orientTemp2[0]:
			edges.reverse()
			orientations.reverse()
			
		if len(edges)>maxL:
			if debug:
				printWarning("Cannot fuse lineLoop with ID " + str(self.Id) + " and lineLoop with ID "+ str(loop.Id) +" . Resulting loop exceeds maxL.")
			return False
		
		#print "inserting ",pyfrp_misc_module.objAttrToList(edges,'Id')
		
		#Insert edges of second loop into loop
		self.edges[idx:idx]=edges
		self.orientations[idx:idx]=orientations
		
		#Check if closed in the end
		#self.checkClosed(fix=True,debug=False)
		self.fix()
		
		#Delete second lineLoop
		if surface!=None:
			loop.removeFromSurface(surface)
		b=loop.delete(debug=debug)
		
		bs=[]
		#Delete common edge
		for edge in edgeTemp1:
			bs.append(edge.delete(debug=debug))
		
		return True
	
	def approxBySpline(self,angleThresh=0.314,debug=False):
		
		"""Approximates parts of line loop by spline.
		
		Summarizes all consecutive lines in loop that have a small angle inbetween to 
		a spline. 
		
		.. note:: The choice of ``angleThresh`` is crucial for this function to work. It should be chosen 
		   on a by-case basis if necessary.
		
		Example: 
		
		Load test file:
		
		>>> d,dd = pyfrp_gmsh_IO_module.readGeoFile("pyfrp/meshfiles/examples/splineTest.geo")
		
		Draw:
		
		>>> d.setAnnOffset([0.1,0.1,0.00])
		>>> ax=d.draw(backend='mpl',asSphere=False,ann=True,annElements=[False,True,False])
		
		returns the following 
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/approxBySpline1.png
		
		Approximate by spline and draw again
		
		>>> d.lineLoops[0].approxBySpline(angleThresh=0.1*np.pi)
		>>> ax=d.draw(backend='mpl',asSphere=False,ann=True,annElements=[False,True,False])
		
		returns
		
		.. image:: ../imgs/pyfrp_gmsh_geometry/approxBySpline2.png
		
		And write to file
		
		>>> d.writeToFile("pyfrp/meshfiles/examples/approximated.geo")
		
		Keyword Args:
			angleThresh (float): Angular threshold in radians.
			debug (bool): Print debugging messages.
			
		Returns:
			bool: True if approximated.
		
		"""
		
		if len(self.edges)<=4:
			return False
		
		nEdges=len(self.edges)
		
		# Bookkeeping lists
		subst=[]
		edgesSubst=[]
		vertices=[]
		appending=False
		for i in range(len(self.edges)):
			
			# Get indices of edges
			idx1=pyfrp_misc_module.modIdx(i,self.edges)
			idx2=pyfrp_misc_module.modIdx(i+1,self.edges)
			
			# Get edges
			e1=self.edges[idx1]
			e2=self.edges[idx2]

			# If either e1 or e2 is not a line, skip right away.
			if e1.typ>0 or e2.typ>0:
				continue
		
			# Compute angle
			angle=pyfrp_geometry_module.getAngle(e1.getDirection(self.orientations[idx1]),e2.getDirection(self.orientations[idx2]))
	
			# If angle satisfies threshold criteria, append edge
			if angle<=angleThresh:
				
				vertices=vertices+[e1.getFirstVertex(self.orientations[idx1]),e1.getLastVertex(self.orientations[idx1])]	
				edgesSubst.append(e1)
				appending=True
			
			# If angle is too large or its the end of the loop, close spline
			if angle>angleThresh or i==len(self.edges)-1:
				if appending==True:
					
					if e1 not in edgesSubst:
						edgesSubst.append(e1)
					
					"""IDEA: Check that includedInLoop gives same return for all edges in edgesSubst. Then finally add 
					spline.
					"""
					noSpline=False
					for j,e in enumerate(edgesSubst):
						inLoop,loops=e.includedInLoop()
							
						if j>0:
							if loops.sort()!=oldLoops.sort():
								noSpline=True
								print("cannot turn into spline because ", edgesSubst[j-1].Id , " and ",edgesSubst[j].Id) 
							
						oldLoops=list(loops)
						
					if not noSpline:
						
						
						""" Check if e2 is already about to be substituted by a spline. Then we can simply add edges together to one spline.
						Otherwise, just create new spline.
						"""
						if len(subst)>0:
							if e2 in subst[0][0]:
								subst[0][0]=edgesSubst+subst[0][0]
								subst[0][2]=list(set(loops+subst[0][2]))
								subst[0][1].vertices=vertices+subst[0][1].vertices		
							else:	
								spline=self.domain.addBSpline(vertices+[e1.getLastVertex(self.orientations[i])])
								subst.append([edgesSubst,spline,loops])
						else:
							spline=self.domain.addBSpline(vertices+[e1.getLastVertex(self.orientations[i])])
							subst.append([edgesSubst,spline,loops])
				
				# Set back Bookkeeping variables
				appending=False
				edgesSubst=[]
				vertices=[]
				
		# Replace edges from loops with spline.
		for sub in subst:
			for loop in sub[2]:
				
				try:
					idx1=loop.edges.index(sub[0][0])
					idx2=loop.edges.index(sub[0][-1])
				except IndexError:
					printWarning("approxBySpline: Cannot find index of either first or last edge.")
				
				if idx1>idx2:
					idxInsert=idx2
				else:
					idxInsert=idx1
				
				for e in sub[0]:
					try:
						loop.edges.remove(e)
					except ValueError:
						printWarning("approxBySpline: Cannot remove edge "+str(e.Id)+" from loop"+loop.Id +".")
				
				loop.edges.insert(idxInsert,sub[1])
				
				if debug:
					print("Substituted edges ", sub[0][0].Id , "-", sub[0][-1].Id, " with spline ", sub[1].Id)
				
				
		# Remove edges from domain.
		for sub in subst:
			for e in sub[0]:
				e.delete(debug=debug)
		
		return nEdges>len(self.edges)
			
	def includedInSurface(self):
		
		"""Checks if loop is included in a surface.
		
		Returns:
			tuple: Tuple containing:
			
				* included (bool): True if included.
				* loops (list): List of :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface` objects that include loop. 
		
		"""
		
		surfaces=[]
		
		for surface in self.domain.ruledSurfaces:
			if self==surface.lineLoop:
				surfaces.append(surface)
		
		return len(surfaces)>0,surfaces

	def delete(self,debug=False):
		
		"""Deletes loop if it is not used in any surface.
		
		Returns:
			bool: True if deletion was successful.
		"""
		
		incl,surfaces=self.includedInSurface()
		if incl:
			printWarning("Was not able to delete loop with ID " + str(self.Id) +". Still part of " + str(len(surfaces)) + " surfaces.")
			return False
		
		self.domain.lineLoops.remove(self)
		
		return True
	
	def removeFromSurface(self,surface):
		
		"""Removes lineLoop from surface.
		"""
		
		if self==surface.lineLoop:
			surface.lineLoop=None
	
	
	def removeFromAllSurfaces(self):
		
		"""Removes lineLoop from all surfaces.
		"""
		
		for surface in self.domain.ruledSurfaces:
			self.removeFromSurface(surface)	
	
	def isCoplanar(self):
		
		"""Returns if all edges lie in single plane.
		
		Does this by 
		
			* picking the first two vertices as first vector ``vec1 = v1 - v0``
			* looping through vertices and computung the normal vector  
			  between ``vec1`` and ``vec2=v[i]-v0``.
			* Checking if all normal vectors are colinear.  
		
		Returns:
			bool: True if coplanar.
		
		"""
		
		#Get vertex coordinates
		coords=pyfrp_misc_module.objAttrToList(self.getVertices(),'x')
		
		#Compute normals
		normals=[]
		for i in range(2,len(coords)):	
			n=pyfrp_geometry_module.computeNormal([coords[0],coords[1],coords[i]])
			normals.append(n)
		
		#Check if normals are all colinear
		b=[]	
		for i in range(1,len(normals)):
			
			# Make sure to skip normal vectors produced from colinear vectors
			if sum(normals[i])==0:
				continue
			
			b.append(pyfrp_geometry_module.checkColinear(normals[0],normals[i]))
			
		return sum(b)==len(b)
	
	def getCenterOfMass(self):
		
		"""Computes center of mass of surface.
		
		Returns:
			numpy.ndarray: Center of mass.
		"""
		
		coords=np.array(pyfrp_misc_module.objAttrToList(self.getVertices(),'x'))
		
		return pyfrp_idx_module.getCenterOfMass(coords)
	
	def getEdges(self):
		
		"""Returns list of edges included in lineLoop."""
		
		return self.edges
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return self.getEdges()

	
class ruledSurface(gmshElement):
	
	"""ruledSurface class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		loopID (int): ID of surrounding loop.
		Id (int): ID of surface.
		
	"""		
	
	def __init__(self,domain,loopID,ID):
		
		gmshElement.__init__(self,domain,ID)
		self.initLineLoop(loopID)
	
	def initLineLoop(self,loopID,debug=False,addPoints=False,iterations=2):
	
		"""Checks length of lineLoop and if length of lineLoop is greater
		than 4, will perform triangulation so Gmsh can handle surface."""
		
		#Get lineLoop
		self.lineLoop=self.domain.getLineLoopById(loopID)[0]
		
		#Check if there is a line loop
		if self.lineLoop==False:
			return False,[]
		
		#Check length
		if len(self.lineLoop.edges)<=4:
			return False,[]
	
		#Compute normal vector
		oldNormal=self.getNormal()
		
		#Create triangulation
		rmat=self.rotateToPlane('xy')
		newNormal=self.getNormal()
		
		#Get vertices
		vertices=self.getVertices()
		coords=pyfrp_misc_module.objAttrToList(vertices,'x')
		
		#Get maximum volSize of vertices
		maxVolSize=max(pyfrp_misc_module.objAttrToList(vertices,'volSize'))
		
		#Get coordinates in plane
		coordsPlane=np.asarray(coords)[:,np.where(abs(self.normal)!=1)[0]]
		
		#Triangulate
		tri,coordsTri = pyfrp_idx_module.triangulatePoly(coordsPlane,addPoints=addPoints,iterations=iterations,debug=True)
		
		#Add 3D dimension to coordsTri	
		coordsTri=np.concatenate((coordsTri,coords[0][2]*np.ones((coordsTri.shape[0],1))),axis=1)
		
		#Loop through each triangle 
		surfacesCreated=[]
		vertices=[]
		for i in range(len(tri)):
			
			edges=[]
			
			#Loop through each vertex 
			for j in range(len(tri[i])):
				
				#Get first vertex, create it if necessary
				v1=self.domain.getVertexByX(coordsTri[tri[i][j]])[0]
				if v1==False:
					v1=self.domain.addVertex(coordsTri[tri[i][j]],volSize=maxVolSize)
				
				#Get second vertex, create it if necessary
				v2=self.domain.getVertexByX(coordsTri[tri[i][pyfrp_misc_module.modIdx(j+1,tri[i])]])[0]
				if v2==False:
					v2=self.domain.addVertex(coordsTri[tri[i][pyfrp_misc_module.modIdx(j+1,tri[i])]],volSize=maxVolSize)
							
				#Check if edge already exists
				if not self.domain.getEdgeByVertices(v1,v2)[0]:
					edges.append(self.domain.addLine(v1,v2))
				else:
					edges.append(self.domain.getEdgeByVertices(v1,v2)[0])
				
				#Remember vertices so we can use them later for turning everything back.
				vertices=vertices+[v1,v2]
				
			#Add line loop
			edgeIDs=pyfrp_misc_module.objAttrToList(edges,"Id")
			loop=self.domain.addLineLoop(edgeIDs=edgeIDs)
			loop.checkClosed(fix=False,debug=False)
			loop.fix()
			
			#Add ruledSurface if necessary
			if i==0:
				self.lineLoop=loop
			else:
				snew=self.domain.addRuledSurface(lineLoopID=loop.Id)	
				surfacesCreated.append(snew)
		
		#Delete original loop
		self.domain.lineLoops.remove(self.domain.getLineLoopById(loopID)[0])
		
		#Remove duplicates in vertices list.
		vertices=pyfrp_misc_module.remRepeatsList(vertices)
		
		#Rotate back
		for v in vertices:
			v.x=np.dot(v.x,rmat.T)
			
		return True,surfacesCreated
		
	def normalToPlane(self):
		
		"""Checks if surface lies within either x-y-/x-z-/y-z-plane.
		
		Does this by checking if ``1.`` is in the normal vector.
		
		Returns:
			bool: True if in plane.
		
		"""
		
		return 1. in self.normal
	
	def isCoplanar(self):
		
		"""Returns if surface lies in single plane.
		
		Returns:
			bool: True if coplanar.
		
		"""
		
		return self.lineLoop.isCoplanar()
	
	def getCenterOfMass(self):
		
		"""Computes center of mass of surface.
		
		Returns:
			numpy.ndarray: Center of mass.
		"""
		
		return self.lineLoop.getCenterOfMass()
		
	def getNormal(self,method='cross'):
		
		"""Computes normal to surface.
		
		First checks if surface is coplanar using :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface.isCoplanar`.
		Then finds two independent vectors that span surface and passes them on to 
		:py:func:`pyfrp.modules.pyfrp_geometry_module.computeNormal`.
		
		Currently there are two methods available:
		
			* ``cross``, see also :py:func:`normalByCross`.
			* ``newells``, see also :py:func:`newells`.
	
		If method is unknown, will fall back to ``cross``.
		
		Keyword Args:
			method (str): Method of normal computation.
		
		Returns:
			numpy.ndarray: Normal vector to surface.
		"""
		
		if not self.isCoplanar():
			printWarning("Surface " + str(self.Id) + " is not coplanar. The resulting normal vector might thus be not correct.")
		
		#Get vertices
		vertices=self.lineLoop.getVertices()
		
		#Find non-colinear vertices
		vec1=vertices[1].x-vertices[0].x
		idx=None
		for i in range(2,len(vertices)):
			tempVec=vertices[i].x-vertices[0].x
			if not pyfrp_geometry_module.checkColinear(vec1,tempVec):
				idx=i
				break
		if idx==None:
			printError("All points in surface "+str(self.Id) + " seem to be colinear. Will not be able to compute normal.")
			print(self.Id)
			self.draw(ann=True)
			
			print(pyfrp_misc_module.objAttrToList(self.lineLoop.getVertices(),'Id'))
			
			input()
			
			return np.zeros((3,))
		
		#Compute normal
		coords=[vertices[0].x,vertices[1].x,vertices[idx].x]
		
		self.normal=pyfrp_geometry_module.computeNormal(coords,method=method)
		
		return self.normal

	def writeToFile(self,f):
		
		"""Writes ruled surface to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Ruled Surface("+str(self.Id)+")= {"+str(self.lineLoop.Id)+ "};\n" )
	
		return f
	
	def addToBoundaryLayer(self,boundField=None,**fieldOpts):
		
		"""Adds surface to a boundary layer field. 
		
		If no field is given, will create new one with given parameters and add it to a minField. If no minField exists,
		will create a new one too and set it as background field.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addBoundaryLayerField`
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.addMinField` and 
		:py:func:`pyfrp.modules.pyfrp_gmsh_geometry.domain.genMinBkgd`.
		
		Keyword Args:
			boundField (pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField): Boundary layer field object.
			fieldOpts (dict): See documentation of boundary layer field of all available options.
			
		Returns:	
			pyfrp.modules.pyfrp_gmsh_geometry.boundaryLayerField: Boundary layer field around edge.
		
		"""
		
		#Generate boundary field if not given
		if boundField==None:
			boundField=self.domain.addBoundaryLayerField()
		
		#Add Vertex
		boundField.addFaceByID(self.Id)
		
		#Set options
		boundField.setFieldAttributes(**fieldOpts)
		
		#Generate background field
		self.domain.genMinBkgd(FieldsList=[boundField.Id])
			
		return boundField	
	
	def hasCommonEdge(self,surface):
		
		"""Checks if surface has common edge with other surface.
		
		Args:
			surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Surface object.
		
		Returns:
			tuple: Tuple containing:
			
				* hasCommon (bool): True if loops have common edge.
				* e (pyfrp.modules.pyfrp_gmsh_geometry.edge): Edge that is in common.
						
		"""
		
		return self.lineLoop.hasCommonEdge(surface.lineLoop)
		
	def fuse(self,surface,maxL=1000,debug=False,sameNormal=False):
		
		"""Fuses surface with another surface.
		
		Will not do anything if surfaces do not have an edge in common. 
		
		"""
	
		if not self.hasCommonEdge(surface)[0]:
			if debug:
				printWarning("Cannot fuse surface with ID " + str(self.Id) + " and surface with ID "+ str(surface.Id) +" . Surfaces do not have common edge.")
			return False
			
		if not self.hasSameNormal(surface):
			if sameNormal:
				if debug:
					printWarning("Cannot fuse surface with ID " + str(self.Id) + " and surface with ID "+ str(surface.Id) +" . Not same normal, but sameNormal="+str(sameNormal))
				return False
			
			if debug:
				printWarning("Fusing surface with ID " + str(self.Id) + " and surface with ID "+ str(surface.Id) +" will alter surface normal.")
		
		b=self.lineLoop.fuse(surface.lineLoop,maxL=maxL,debug=debug,surface=surface)
		if b:
			surface.removeFromAllLoops()
			surface.delete()
		
		return b
	
	def removeFromAllLoops(self):
		
		"""Removes surface from all surface loops.
		"""
		
		for loop in self.domain.surfaceLoops:
			if self in loop.surfaces:
				loop.surfaces.remove(self)
			
	def hasSameNormal(self,surface,sameOrientation=False):
		
		"""Checks if sufrace has the same normal vector as another surface.
		
		Args:
			surface (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): Surface object.
		
		Keyword Args:
			sameOrientation (bool): Forces surfaces to also have same orientation.
			
		Returns:
			bool: True if same normal vector.
			
		"""
		
		# Compute Normals
		self.getNormal()
		surface.getNormal()
		
		if sameOrientation:
			if pyfrp_misc_module.compareVectors(self.normal,surface.normal):
				return True
		else:
			if pyfrp_misc_module.compareVectors(self.normal,surface.normal) or pyfrp_misc_module.compareVectors(self.normal,-surface.normal):
				return True
		return False
			
	def includedInLoop(self):
		
		"""Checks if surface is included in a surfaceLoop.
		
		Returns:
			tuple: Tuple containing:
			
				* included (bool): True if included.
				* loops (list): List of :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.surfaceLoops` objects that include surface. 
		
		"""
		
		loops=[]
		
		for loop in self.domain.surfaceLoops:
			if self in loop.surfaces:
				loops.append(loop)
		
		return len(loops)>0,loops

	def delete(self):
		
		"""Deletes surface if it is not used in any surfaceLoop.
		
		Returns:
			bool: True if deletion was successful.
		"""
		
		incl,loops=self.includedInLoop()
		if incl:
			printWarning("Was not able to delete loop with ID " + str(self.Id) +". Still part of " + str(len(loops)) + " loops.")
			return False
		
		self.domain.ruledSurfaces.remove(self)
		
		return True
	
	def draw(self,ax=None,color='b',edgeColor='k',drawLoop=True,ann=None,alpha=0.2,backend='mpl',linewidth=1):
		
		"""Draws surface and fills it with color.
		
		.. note:: If ``ann=None``, will set ``ann=False``.
		
		.. note:: If no axes is given, will create new one.
		
		.. warning:: Does not work for surfaces surrounded by arcs yet.
		
		Keyword Args:
			ax (matplotlib.axes): Matplotlib axes to be plotted in.
			color (str): Color of surface.
			ann (bool): Show annotations.
			edgeColor (str): Color of lineLoop around.
			alpha (float): Transparency of surface. 
				
		Returns:
			matplotlib.axes: Axes.
		
		"""
		
		if backend!='mpl':
			printError("Cannot draw surface with backend="+backend+". Currently not supported")
			return ax
		
		if ann==None:
			ann=False
			
		if ax==None:
			fig,axes = pyfrp_plot_module.makeGeometryPlot()
			ax=axes[0]
		
		if drawLoop:
			ax=self.lineLoop.draw(ax=ax,color=edgeColor,ann=False,linewidth=linewidth)
		
		for e in self.lineLoop.edges:
			if e in self.domain.arcs:
				printWarning("Cannot draw surface " + str(self.Id) + " yet. Surfaces including arcs are not supported yet.")
				return ax
		
		#Get Vertex coordinates in the form we want (this is probably unnecessarily complicated)
		vertices=self.lineLoop.getVertices()
		coords=pyfrp_misc_module.objAttrToList(vertices,'x')
		coords=np.asarray(coords)
		coords = list(zip(coords[:,0], coords[:,1], coords[:,2]))
		coordsNew=[]
		coordsNew.append(list(coords))	
		
		#Add collection
		coll=Poly3DCollection(coordsNew,alpha=alpha)
		coll.set_facecolor(color)
		ax.add_collection3d(coll)
		
		#annotation
		if ann:
			com=pyfrp_idx_module.getCenterOfMass(np.array(pyfrp_misc_module.objAttrToList(vertices,'x')))
			ax.text(com[0],com[1],com[2], "s"+str(self.Id), None)
		
		#Redraw
		pyfrp_plot_module.redraw(ax)
		
		return ax
	
	def getVertices(self):
		
		"""Returns all vertices included in surface."""
		
		if self.lineLoop==None:
			return []
		else:
			return self.lineLoop.getVertices()
	
	def getEdges(self):
		
		"""Returns all edges included in surface."""
		
		if self.lineLoop==None:
			return []
		else:
			return self.lineLoop.getEdges()
	
	
	def rotateToNormal(self,normal,ownNormal=None):
		
		"""Rotates surface such that it lies in the plane
		with normal vector ``normal``.
		
		See also :py:func:`pyfrp.modules.pyfrp_geometry_module.getRotMatrix`.
		
		Args:
			normal (numpy.ndarray): Normal vector.
		
		Returns: 
			numpy.ndarray: Rotation matrix.
		
		"""
		
		if ownNormal==None:
			ownNormal=self.getNormal()
		
		if pyfrp_geometry_module.checkColinear(normal,ownNormal):
			rmat=np.identity(3)
		else:
			rmat=pyfrp_geometry_module.getRotMatrix(normal,ownNormal)
		
		# Rotate
		for v in self.getVertices():
			v.x=np.dot(v.x,rmat)
			
		return rmat	
			
	def rotateToSurface(self,s):
			
		"""Rotates surface such that it lies in the same plane
		as a given surface.
		
		See also :py:func:`pyfrp.modules.pyfrp_geometry_module.getRotMatrix`.
		
		Args:
			s (pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface): A surface.
		
		Returns: 
			numpy.ndarray: Rotation matrix.
		
		"""
		
		return self.rotateToNormal(s.getNormal())
		
	def rotateToPlane(self,plane):
		
		"""Rotates surface such that it lies in plane.
		
		See also :py:func:`pyfrp.modules.pyfrp_geometry_module.getRotMatrix`.
		
		Possible planes are:
		
			* ``xy``
			* ``xz``
			* ``yz``
			
		Args:
			plane (str): Plane to rotate to.
		
		Returns: 
			numpy.ndarray: Rotation matrix.
		
		"""
		
		if plane=="xz":
			normal=np.array([0,0,1.])
		elif plane=="yz":
			normal=np.array([1.,0,0])
		elif plane=="xy":
			normal=np.array([0,0,1.])
		else:
			printError("Do not know the plane " +plane +". Will not rotate plane")
			return
		
		return self.rotateToNormal(normal)
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.lineLoop]
	
class surfaceLoop(gmshElement):
	
	"""surfaceLoop class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain loop belongs to.
		surfaceIDs (list): List of surfaces.
		Id (int): ID of loop.
		
	"""		
	
	def __init__(self,domain,surfaceIDs,ID):
		
		gmshElement.__init__(self,domain,ID)
		
		self.surfaces=self.initSurfaces(surfaceIDs)	
		
	def initSurfaces(self,IDs):
		
		"""Constructs ``surfaces`` list at object initiations 
		from list of IDs.
		
		Args:
			IDs (list): List of IDs.
			
		Returns:
			list: List of pyfrp.modules.pyfrp_gmsh_geometry.ruledSurface objects.
		
		"""
		
		self.surfaces=[]
		
		for ID in IDs:
			self.addSurfaceByID(ID)
		
		return self.surfaces
	
	def addSurfaceByID(self,ID):
		
		"""Adds surface to surfaceloop.
		
		Args:
			ID (int): ID of surface to be added.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.append(self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def insertSurfaceByID(self,ID,pos):
		
		"""Inserts surface to surfaceloop at position.
		
		Args:
			ID (int): ID of surface to be inserted.
			pos (int): Position at which ID to be inserted.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.insert(pos,self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def removeSurfaceByID(self,ID):
		
		"""Remove surface from surfaceloop.
		
		Args:
			ID (int): ID of surface to be removed.
			
		Returns:
			list: Updated surfaceIDs list.
			
		"""
		
		self.surfaces.remove(self.domain.getRuledSurfaceById(ID)[0])
		
		return self.surfaces
	
	def writeToFile(self,f):
		
		"""Writes surface loop to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Surface Loop("+str(self.Id)+")= {" )
		
		for i,s in enumerate(self.surfaces):
			f.write(str(s.Id))
			if i!=len(self.surfaces)-1:
				f.write(",")
			else:
				f.write("};\n")
	
		return f
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return self.surfaces
	
class volume(gmshElement):

	"""Volume class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		surfaceLoopID (int): ID of surrounding surface loop.
		Id (int): ID of surface loop.
		
	"""		
	
	def __init__(self,domain,surfaceLoopID,ID):
		
		gmshElement.__init__(self,domain,ID)
		self.surfaceLoop=self.domain.getSurfaceLoopById(surfaceLoopID)[0]
	
	def writeToFile(self,f):
		
		"""Writes Volume to file.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f.write("Volume("+str(self.Id)+")= {"+str(self.surfaceLoop.Id)+ "};\n" )
	
		return f
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.surfaceLoop]
	
class field(gmshElement):

	"""Field class storing information from gmsh .geo.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
		typ (str): Type of field.
		
	"""		
	
	def __init__(self,domain,typ,Id):
		
		gmshElement.__init__(self,domain,Id)
		self.typ=typ
	
	def setAsBkgdField(self):
		
		"""Sets this mesh as background field for the whole domain.
		"""
		
		self.domain.bkgdField=self
	
	def isBkgdField(self):
		
		"""Returns true if field is background field.
		
		Returns:
			bool: True if background field.
		"""
		
		return self.domain.bkgdField==self
	
	def setFieldAttr(self,name,val):
		
		"""Sets attribute of field.
		
		.. note:: Value can have any data type.
		
		Args:
			name (str): Name of attribute.
			val (str): Value.
		
		"""
		
		setattr(self,name,val)
	
	def setFieldAttributes(self,**kwargs):
		
		"""Sets multiple field attributes.
		"""
		
		for key, value in kwargs.items():
			self.setFieldAttributes(key,value)
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return []
	
class boxField(field):
	
	"""Box field class storing information from gmsh .geo.
	
	Subclasses from :py:class:`field`.
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
		
	Keyword Args:
		volSizeIn (float): Mesh element volume inside box.
		volSizeOut (float): Mesh element volume outside box.
		xRange (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
		yRange (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
		zRange (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
		
	"""
	
	def __init__(self,domain,Id,volSizeIn=10.,volSizeOut=20.,xRange=[],yRange=[],zRange=[]):
	
		field.__init__(self,domain,"box",Id)
		
		self.VOut=volSizeIn
		self.VIn=volSizeOut
		self.initBox(xRange,yRange,zRange)
		
	def initBox(self,xRange,yRange,zRange):
		
		"""Initializes bounding box.
		"""
		
		self.setRange('X',xRange)
		self.setRange('Y',yRange)
		self.setRange('Z',zRange)
		
	def setRange(self,coord,vec):
		
		"""Sets the bounding box range along a given axis.
		
		Args:
			coord (str): Axis along range is set (``"X","Y","Z"``) 
			vec (list): Range of box ``[minVal,maxVal]``
		
		Returns:
			tuple: Tuple containing:
			
				* coordMin (float): New minimum value.
				* coordMax (float): New maximum value.
				
		"""
		
		try:
			setattr(self,coord+"Min",vec[0])
			setattr(self,coord+"Max",vec[1])
		except IndexError:
			setattr(self,coord+"Min",None)
			setattr(self,coord+"Max",None)
		
		return getattr(self,coord+"Min"),getattr(self,coord+"Max")
		
	def writeToFile(self,f):
			
		"""Writes box field to file.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.writeBoxField`.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
			
		f=pyfrp_gmsh_IO_module.writeBoxField(f,self.Id,self.VIn,self.VOut,[self.XMin,self.XMax],[self.YMin,self.YMax],[self.ZMin,self.ZMax])
		
		if self.isBkgdField():
			f=pyfrp_gmsh_IO_module.writeBackgroundField(f,self.Id)
			
		return f	
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return []
	
	
class attractorField(field):
	
	"""Attractor field class storing information from gmsh .geo.

	Subclasses from :py:class:`field`.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
		
	Keyword Args:
		NodesList (list): List of IDs of the Nodes that attractor field centers around.
			
	"""

	def __init__(self,domain,Id,NodesList=[]):
	
		field.__init__(self,domain,"attractor",Id)
		
		
		self.NodesList=self.initNodesList(NodesList)
	
	def initNodesList(self,NodesList):
		
		"""Adds a list of vertices to NodesList.
		
		See also :py:func:`addNodeByID`.
		
		Args:
			NodesList (list): List of vertex IDs.
			
		Returns:
			list: Updated NodesList.
		
		"""
		
		self.NodesList=[]
		
		for Id in NodesList:
			self.addNodeByID(Id)
			
		return self.NodesList	
			
	def addNodeByID(self,ID):
		
		"""Adds vertex object to NodesList given the ID of the vertex.
		
		Args:
			ID (int): ID of vertex to be added.
			
		Returns:
			list: Updated NodesList.
		
		"""
		
		v,b=self.domain.getVertexById(ID)
		
		if isinstance(b,int):
			self.NodesList.append(v)
		return self.NodesList
	
	def setFieldAttr(self,name,val):
		
		"""Sets field attribute.
		
		.. note:: Value can have any data type.
		
		Args:
			name (str): Name of attribute.
			val (float): Value of attribute.
			
			
		"""
		
		if name=="NodesList":
			self.initNodesList(val)
		else:
			setattr(self,name,val)
	
	def writeToFile(self,f):
			
		"""Writes attractor field to file.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.writeAttractorField`.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
			
		f=pyfrp_gmsh_IO_module.writeAttractorField(f,self.Id,pyfrp_misc_module.objAttrToList(self.NodesList,'Id'))
		
		if self.isBkgdField():
			f=pyfrp_gmsh_IO_module.writeBackgroundField(f,self.Id)
			
		return f	
	
	def includedInThresholdField(self):
		
		"""Returns all the threshholdFields where attractorField is included in.
		
		Returns:
			list: List of threshholdField objects.
		
		"""
		
		threshFields=self.domain.getAllFieldsOfType("threshold")
		
		included=[]
		for tField in threshFields:
			if tField.IField==self.Id:
				included.append(tField)
		
		return included
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return self.NodesList	
	
class thresholdField(field):		
	
	"""Threshold field class storing information from gmsh .geo.
	
	Subclasses from :py:class:`field`.
	
	
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
	
	Keyword Args:
		IField (int): ID of vertex that is center to threshold field.
		LcMin (float): Minimum volSize of threshold field.
		LcMax (float): Maximum volSize of threshold field.
		DistMin (float): Minimun density of field.
		DistMax (float): Maximum density of field.
		
	"""
	
	def __init__(self,domain,Id,IField=None,LcMin=5.,LcMax=20.,DistMin=30.,DistMax=60.):
		
		field.__init__(self,domain,"threshold",Id)
		
		self.IField=IField
		self.LcMin=LcMin
		self.LcMax=LcMax
		self.DistMin=DistMin
		self.DistMax=DistMax
		
	def writeToFile(self,f):
			
		"""Writes threshold field to file.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.writeThresholdField`.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
			
		f=pyfrp_gmsh_IO_module.writeThresholdField(f,self.Id,self.IField,self.LcMin,self.LcMax,self.DistMin,self.DistMax)
		
		if self.isBkgdField():
			f=pyfrp_gmsh_IO_module.writeBackgroundField(f,self.Id)
			
		return f	
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.getDomain().getVertexById(self.IField)]
	
class minField(field):
	
	"""Minimum field class storing information from gmsh .geo.
	
	Subclasses from :py:class:`field`.

	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
		
	Keyword Args:
		FieldsList (list): List of field IDs.
		
	"""
	
	def __init__(self,domain,Id,FieldsList=[]):
		
		field.__init__(self,domain,"min",Id)
		
		self.FieldsList=self.initFieldsList(FieldsList)
	
	def setFieldAttr(self,name,val):
	
		if name=="FieldsList":
			self.initFieldsList(val)
		else:
			setattr(self,name,val)
				
	def initFieldsList(self,FieldsList):
		
		"""Adds a list of vertices to NodesList.
		
		See also :py:func:`addNodeByID`.
		
		Args:
			FieldsList (list): List of field IDs.
			
		Returns:
			list: Updated FieldsList.
		
		"""
		
		self.FieldsList=[]
		
		for Id in FieldsList:
			self.addFieldByID(Id)
		
		return self.FieldsList
		
	def addFieldByID(self,ID):
		
		"""Adds field object to FieldsList given the ID of the field.
		
		Args:
			ID (int): ID of field to be added.
			
		Returns:
			list: Updated FieldsList.
		
		"""
		
		f,b=self.domain.getFieldById(ID)
		
		if isinstance(b,int):
			self.FieldsList.append(f)
		return self.FieldsList
		
	def writeToFile(self,f):
			
		"""Writes minimum field to file.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.writeMinField`.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		f=pyfrp_gmsh_IO_module.writeMinField(f,self.Id,pyfrp_misc_module.objAttrToList(self.FieldsList,'Id'))
		
		if self.isBkgdField():
			f=pyfrp_gmsh_IO_module.writeBackgroundField(f,self.Id)
			
		return f	
	
	def addAllFields(self):
		
		"""Adds all fields in domain to FieldsList if not 
		already in there.
		
		Returns:
			list: Updated FieldsList.
		"""
		
		for f in self.domain.fields:
			if f not in self.FieldsList and f!=self:
				self.FieldsList.append(f)
				
		return self.FieldsList
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return self.FieldsList
		
class boundaryLayerField(field):
	
	r"""Boundary Layer field class storing information from gmsh .geo.
	
	Creates boundary layer mesh around vertices, edges or surfacesin geometry. Boundary layer density 
	is given by
	
	.. math:: h_{wall} * ratio^{(dist/h_{wall})}.
	
	Subclasses from :py:class:`field`.
	
	Example: Adding a box surrounded with a boundary layer to a geometry:
	
	>>> vertices,lines,loops,surfaces,sloops,vols=d.addCuboidByParameters([256-50,256-50,-160],100,100,120,10,genVol=False)

	Adjust volSize:
	
	>>> d.setGlobalVolSize(30.)

	Add boundary layer:
	
	>>> volSizeLayer=10.
	>>> blf=d.addBoundaryLayerField(hfar=volSizeLayer,hwall_n=volSizeLayer,hwall_t=volSizeLayer,thickness=30.,Quads=0.)
	>>> blf.addFaceListByID(pyfrp_misc_module.objAttrToList(surfaces,'Id'))
	>>> blf.setAsBkgdField()
	>>> d.draw()
	
	.. image:: ../imgs/pyfrp_gmsh_geometry/boundaryLayerField_geometry.png
	
	Write to file:
	
	>>> d.writeToFile("dome_boundary.geo")
	
	Generate mesh:
	
	>>> fnMesh=pyfrp_gmsh_module.runGmsh("dome_boundary.geo")
	>>> m=pyfrp_mesh.mesh(None)
	>>> m.setFnMesh(fnMesh)
	>>> m.plotMesh()
	
	.. image:: ../imgs/pyfrp_gmsh_geometry/boundaryLayerField_mesh.png
	
	See also http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes .
	
	Args:
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain surface belongs to.
		Id (int): ID of field.
		
	Keyword Args:
		AnisoMax (float): Threshold angle for creating a mesh fan in the boundary layer.
		IntersectMetrics (int): Intersect metrics of all faces.
		Quad (int): Generate recombined elements in the boundary layer.
		har (float): Element size far from the wall.
		hwall_n (float): Mesh Size Normal to the The Wall.
		hwall_t (float): Mesh Size Tangent to the Wall.
		ratio (float): Size Ratio Between Two Successive Layers.
		thickness (float): Maximal thickness of the boundary layer.
		List (list): List of field IDs.
		
	"""
	
	def __init__(self,domain,Id,AnisoMax=10000000000,hwall_n=1.,hwall_t=1,ratio=1.1,thickness=10.,hfar=1.,IntersectMetrics=1,Quads=0.):
		
		field.__init__(self,domain,"boundaryLayer",Id)
		
		self.AnisoMax=AnisoMax
		self.Quads=Quads
		self.hfar=hfar
		self.hwall_n=hwall_n
		self.hwall_t=hwall_t
		self.IntersectMetrics=IntersectMetrics
		self.ratio=ratio
		self.thickness=thickness
	
		self.EdgesList=[]
		self.FacesList=[]
		self.FanNodesList=[]
		self.FansList=[]
		self.NodesList=[]
	
	
	def addEdgeListByID(self,IDs):
		
		"""Adds a list of edge objects to EdgesList given the ID of the edges.
		
		Args:
			IDs (list): List of IDs of edges to be added.
			
		Returns:
			list: Updated EgesList.
		
		"""
		
		for ID in IDs:
			self.addEdgeByID(ID)
			
		return self.EdgesList	
	
	def addFaceListByID(self,IDs):
		
		"""Adds a list of surfaces objects to FacesList given the ID of the surfaces.
		
		Args:
			IDs (list): List of IDs of surfaces to be added.
			
		Returns:
			list: Updated FacesList.
		
		"""
		
		for ID in IDs:
			self.addFaceByID(ID)
			
		return self.FacesList	
	
	def addNodeListByID(self,IDs):
		
		"""Adds a list of vertex objects to NodesList given the ID of the vertex.
		
		Args:
			IDs (list): List of IDs of vertices to be added.
			
		Returns:
			list: Updated NodesList.
		
		"""
		
		for ID in IDs:
			self.addNodeByID(ID)
			
		return self.NodesList	
	
	def addEdgeByID(self,ID):
		
		"""Adds edge object to EdgesList given the ID of the edge.
		
		Args:
			ID (int): ID of edge to be added.
			
		Returns:
			list: Updated EgesList.
		
		"""
		
		v,b=self.domain.getEdgeById(ID)
		
		if isinstance(b,int):
			self.EdgesList.append(v)
		return self.EdgesList
		
	def addFaceByID(self,ID):
		
		"""Adds surface object to FacesList given the ID of the surface.
		
		Args:
			ID (int): ID of surface to be added.
			
		Returns:
			list: Updated FacesList.
		
		"""
		
		v,b=self.domain.getRuledSurfaceById(ID)
		
		if isinstance(b,int):
			self.FacesList.append(v)
		return self.FacesList
	
	def addNodeByID(self,ID):
		
		"""Adds vertex object to NodesList given the ID of the vertex.
		
		Args:
			ID (int): ID of vertex to be added.
			
		Returns:
			list: Updated NodesList.
		
		"""
		
		v,b=self.domain.getVertexById(ID)
		
		if isinstance(b,int):
			self.NodesList.append(v)
		return self.NodesList
	
	def buildElementDict(self):
		
		"""Builds element dictionary for writing to file.
		"""
		
		elements={}
		for elmnt in ["EdgesList","FacesList","NodesList"]:
			if len(getattr(self,elmnt))>0:
				elements[elmnt]=pyfrp_misc_module.objAttrToList(getattr(self,elmnt),'Id')
		
		return elements
	
	def writeToFile(self,f):
			
		"""Writes boundaryLayerField to file.
		
		See also :py:func:`pyfrp.modules.pyfrp_gmsh_IO_module.writeBoundaryLayerField`.
		
		Args:
			f (file): File to write to.
			
		Returns:
			file: File.
		
		"""
		
		#elements=pyfrp_misc_module.objAttr2Dict(self,attr=["EdgesList","FacesList","NodesList"])
		
		elements=self.buildElementDict()
		fieldOpts=pyfrp_misc_module.objAttr2Dict(self,attr=["AnisoMax","hwall_n","hwall_t","ratio","thickness","hfar","IntersectMetrics","Quads"])
		
		f=pyfrp_gmsh_IO_module.writeBoundaryLayerField(f,self.Id,elements,fieldOpts)
		
		if self.isBkgdField():
			f=pyfrp_gmsh_IO_module.writeBackgroundField(f,self.Id)
			
		return f	
	
	def getSubElements(self):
		
		"""Returns all elements that define this element.
		
		Returns:
			list: List of elements.
		"""
		
		return [self.FacesList,self.EdgesList,self.NodesList]
	
	def setFieldAttr(self,name,val):
		
		"""Sets field attribute.
		
		.. note:: Value can have any data type.
		
		Args:
			name (str): Name of attribute.
			val (float): Value of attribute.
			
			
		"""
		
		
		self.EdgesList=[]
		self.FacesList=[]
		self.FanNodesList=[]
		self.FansList=[]
		self.NodesList=[]
	
		
		if name=="NodesList":
			self.addNodeListByID(val)
		elif name=="FacesList":
			self.addFaceListByID(val)
		elif name=="EdgesList":
			self.addEdgeListByID(val)			
		else:
			setattr(self,name,val)
	
		
