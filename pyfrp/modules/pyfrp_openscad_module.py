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

"""PyFRAP module for handling building geometries via openscad.

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

#Matplotlib
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

#OS
import os
		
#===========================================================================================================================================================================
#Module functions
#===========================================================================================================================================================================

def runOpenscad(fn,fnOut=None,debug=False,fnPath=None):
	
	"""Runs openscad to convert scad file to stl file.
	
	If ``fnOut=None``, then the output file will have the same filename
	as the input file.
	
	Args:
		fn (str): Path to scad file.
		
	Keyword Args:
		fnOut (str): Output filename.
		debug (bool): Show debugging output.
		fnPath (str): Filepath of path file/
	
	Return:
		str: Output filename.
	
	"""
	
	openscadBin=pyfrp_misc_module.getOpenscadBin(fnPath=fnPath)
	addOpenscadToPATHs(fnPath=fnPath)
	
	if fnOut==None:
		fnOut=fn.replace(".scad",".stl")
	
	cmd=openscadBin+" -o " + fnOut+" "+fn
	
	if debug:
		printNote("Running: "+cmd)
	
	os.system(cmd)
	
	return fnOut
	
def addOpenscadToPATHs(fnPath=None):
	
	"""Exports openscad binary defined in path file to current session.
	
	Uses ``getOpenscadBin`` to get path to openscad binary and then sets environment variable.
	
	Keyword Args:
		fnPath (str): Filepath of path file.
	
	"""
	
	# Get openscad bin as defined in path file
	openscadBin=pyfrp_misc_module.getOpenscadBin(fnPath=fnPath)
	fn=os.path.dirname(openscadBin)
	
	return pyfrp_misc_module.addToPATHs(fn)
	
	
	
	
	