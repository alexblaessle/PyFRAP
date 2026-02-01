"""This module imports all tests/unittests for the
pyfrp_gmsh_geometry."""

from pyfrp.modules import pyfrp_gmsh_geometry
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_gmsh_IO_module

import numpy as np

def test_surfaceFuse():

	"""Test surface fuse function. 

	Reads in .stl, tries to fuse a few surfaces and checks if new surface 
	lineLoop has proper number of vertices."""

	d=pyfrp_gmsh_IO_module.readStlFile(pyfrp_misc_module.fixPath(pyfrp_misc_module.getMeshfilesDir()+"tests/surfaceFuse.stl"))
	
	sfID=1
	toFuseIDs=[2,3,4,5,6,7]
	
	# Grab first surface
	sf1=d.getRuledSurfaceById(sfID)[0]

	# Try to successively fuse
	for ID in toFuseIDs:
		sf1.fuse(d.getRuledSurfaceById(ID)[0],debug=True)
		
	assert pyfrp_misc_module.objAttrToList(sf1.lineLoop.getVertices(),'Id') == [1, 2, 4, 5, 6, 7, 8, 9, 3]
	
def test_domainSimplifySurfaces():

	"""Test domain's simplify surfaces method. 

	Reads in .stl, tries to the geometry and checks if 
	it has been simplified properly."""
	
	d=pyfrp_gmsh_IO_module.readStlFile(pyfrp_misc_module.fixPath(pyfrp_misc_module.getMeshfilesDir()+"tests/surfaceFuse.stl"))
	
	d.simplifySurfaces(triangIterations=0,addPoints=False,fixSurfaces=False,debug=False,iterations=3)
	sameNormal=d.getAllObjectsWithProp("ruledSurfaces","normal",np.array([0,-1.,0]))
	
	assert len(sameNormal) == 1
	
def test_getAllSubElements():	
	
	"""Test gmshElement's getAllSubElements.
	
	Creates cuboid domain and checks that sub-elements of first surface
	have the right number.
	"""
	
	# Create domain
	d=pyfrp_gmsh_geometry.domain()

	# Add Cuboid
	d.addCuboidByParameters([0,0,0],100,150,50,30.,plane="z",genLoops=True,genSurfaces=True,genVol=True)

	assert len(d.ruledSurfaces[0].getAllSubElements())==13