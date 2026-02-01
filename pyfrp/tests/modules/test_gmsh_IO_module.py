"""This module imports all tests/unittests for the
pyfrp_gmsh_IO_module."""

from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules import pyfrp_misc_module


def test_readStlFile():

	"""Test function for readStlFile. 

	Reads in .stl and checks if there are enough
	vertices read in."""

	d=pyfrp_gmsh_IO_module.readStlFile(pyfrp_misc_module.fixPath(pyfrp_misc_module.getMeshfilesDir()+"tests/readStlFile.stl"))
	
	assert len(d.vertices) == 5
	
	