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

"""PyFRAP module for running Gmsh on .geo files. Module mainly has the following features:

	* Functions for updating parameters in standard .geo files.
	* Mesh refinement.
	* Running Gmsh

This module together with pyfrp.pyfrp_gmsh_geometry and pyfrp.pyfrp_gmsh_IO_module works partially as a python gmsh wrapper, however is incomplete.
If you want to know more about gmsh, go to http://gmsh.info/doc/texinfo/gmsh.html .
	
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
from numpy import *

#Misc
import os
import shutil
from tempfile import mkstemp
import subprocess
import time
import shlex
import platform

#PyFRAP
from . import pyfrp_gmsh_IO_module
from . import pyfrp_misc_module
from .pyfrp_term_module import *
           
#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def updateCylinderGeo(fn,radius,height,center,run=True,debug=False):
	
	"""Upates parameters in default *cylinder.geo* file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: For function to work, parameters in geo file need to be defined as follows:
	   
		* radius -> radius
		* height -> height
		* center -> center_x , center_y
	
	Args:
		fn (str): Filepath.
		radius (float): New cylinder radius.
		height (float): New cylinder height.
		center (list): New cylinder center.
	
	Keyword Args:
		run (bool): Run gmsh on updated file.
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"height",height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

def updateConeGeo(fn,upperRadius,lowerRadius,height,center,run=True,debug=False):
	
	"""Upates parameters in default *cone.geo* file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: For function to work, parameters in geo file need to be defined as follows:
	   
		* upperRadius -> upper_radius
		* lowerRadius -> lower_radius
		* slice_height -> slice_height
		* center -> center_x , center_y
	
	Args:
		fn (str): Filepath.
		upperRadius (float): New upper cone radius.
		lowerRadius (float): New lower cone radius.
		height (float): New cone height.
		center (list): New cone center.
	
	Keyword Args:
		run (bool): Run gmsh on updated file.
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"upper_radius",upperRadius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"lower_radius",lowerRadius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"height",height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh


def updateBallGeo(fn,radius,center,run=True,debug=False):
	
	"""Upates parameters in default *ball.geo* file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: For function to work, parameters in geo file need to be defined as follows:
	   
		* radius -> radius
		* center -> center_x , center_y
	
	Args:
		fn (str): Filepath.
		radius (float): New ball radius.
		center (list): New ball center.
	
	Keyword Args:
		run (bool): Run gmsh on updated file.
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	
	v=5*int(debug)
		
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])

	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

def updateDomeGeo(fn,radius,slice_height,center,run=False,debug=False):
	
	"""Upates parameters in default *dome.geo* file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: The way that *dome.geo* is written, gmsh will automatically compute 
	   the dome geometry from slice_height and radius.
	
	.. note:: For function to work, parameters in geo file need to be defined as follows:
	   
		* radius -> radius
		* slice_height -> slice_height
		* center -> center_x , center_y
	
	
	Args:
		fn (str): Filepath.
		radius (float): New dome imaging radius.
		slice_height (float): Height of imaging slice.
		center (list): New dome center.
	
	Keyword Args:
		run (bool): Run gmsh on updated file.
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	
	v=5*int(debug)
	
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"radius",radius)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"slice_height",slice_height)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_x",center[0])
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"center_y",center[1])
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh

def updateVolSizeGeo(fn,volSize_px,run=False,debug=False):
	
	"""Upates parameter that defines mesh element volume in *.geo* file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: For function to work, parameters in geo file need to be defined as follows:
	   
		* volSize_px -> volSize_px
	
	Args:
		fn (str): Filepath.
		volSize_px (float): New mesh element size.
	
	Keyword Args:
		run (bool): Run gmsh on updated file.
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	v=5*int(debug)
	pyfrp_gmsh_IO_module.updateParmGeoFile(fn,"volSize_px",volSize_px)
	
	if run:
		gmshBin=getGmshBin()
		os.system(gmshBin + "  -v " + str(v) +" -3 " + fn)
	
	fn_msh=fn.replace(".geo",".msh")
	
	return fn_msh
	


def refineMsh(fn,debug=False):
	
	"""Refines mesh by splitting elements.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	Args:
		fn (str): Filepath.
		
	Keyword Args:
		debug (bool): Print debugging messages.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	v=5*int(debug)
	
	gmshBin=getGmshBin()
	os.system(gmshBin+" -v "+ str(v) + " -refine " + fn)
	return fn

def runGmsh(fn,fnOut=None,debug=False,redirect=False,fnStout=None,fnSterr=None,volSizeMax=None,dim=3):
	
	"""Runs Gmsh generating mesh from .geo file.
	
	.. note:: Debug will also activate full debugging output of gmsh. See also
	   http://gmsh.info/doc/texinfo/gmsh.html#Command_002dline-options .
	
	.. note:: If ``redirect=True``, but ``fnStout`` or  ``fnSterr`` is not specified,
	   will dump stout/sterr into ``meshfiles/gmshLogs/``.
	
	.. note:: Gmsh is run with the following settings (if all flags are activated):
	   ``gmsh -v -3 -optimize -algo3d -clmax volSizeMax -o fnOut fn``
	   This requires that Gmsh was compiled with TetGen algorithm. PyFRAP can be installed with
	   Gmsh + TetGen included by choosing the ``--gmsh`` flag. 
	   See also http://pyfrap.readthedocs.org/en/latest/setup.html#pyfrap-setup-py-api .
	
	Args:
		fn (str): Filepath.
		
	Keyword Args:
		fnOut (str): Output filepath.
		debug (bool): Print debugging messages.
		redirect (bool): Redirect gmsh stout/sterr into seperate files.
		fnStout (str): File for gmsh stout.
		fnSterr (str): File for gmsh sterr.
		volSizeMax (float): Maximum allowed mesh element size.
		dim (int): Dimension of mesh.
		
	Returns:
		str: Path to mesh file.
		
	""" 
	
	#Define where to put log files if necessary
	if fnStout==None:
		fnStout=pyfrp_misc_module.getMeshfilesDir()+'gmshLogs/gmsh.stout'
		
	if fnSterr==None:
		fnSterr=pyfrp_misc_module.getMeshfilesDir()+'gmshLogs/gmsh.sterr'
		
	v=5*int(debug)
	
	#POpen needs to have paths with `/` as seperator, so we need to change in the case of 
	#Windows the path
	if platform.system() in ["Windows"]:
		fn=pyfrp_misc_module.win2linPath(fn)
		if fnOut!=None:
			fnOut=pyfrp_misc_module.win2linPath(fnOut)
	
	# Add Gmsh to PATH and grab gmshBin location
	gmshBin=getGmshBin()
	addGmshToPATHs()
	
	#Define which command to execute
	cmd = gmshBin + " -v " + str(v) +" -"+str(dim)+ " -optimize -algo del3d"
	if volSizeMax!=None:
		cmd = cmd + " -clmax " + str(volSizeMax)
	if fnOut!=None:
		cmd = cmd + " -o " + fnOut + " "
	else:
		fnOut=fn.replace('.geo','.msh')
	cmd = cmd+ " " + fn
	
	#Print out what will be done
	if debug:
		print("Will execute:")
		print(cmd)
	
	#Split command in list for subprocess
	args = shlex.split(cmd)
	
	#redirect stdout and stderr if selected
	if redirect:
		stoutFile = open(fnStout,'wb')
		sterrFile = open(fnSterr,'wb')
	else:	
		stoutFile = None
		sterrFile = None
		
	#Call gmsh via subprocess and wait till its done
	try:
		p = subprocess.Popen(args,stdout=stoutFile,stderr=sterrFile)
		p.wait()
	except:
		printError("Gmsh is not running properly, something is wrong.")
	
	#Fix path to make sure that it exists and is consistent with OS
	fnOut=pyfrp_misc_module.fixPath(fnOut)
	
	return fnOut

def getGmshBin(fnPath=None):
	
	"""Returns path to Gmsh binary defined in *path* file.	
	""" 
	
	gmshBin=pyfrp_misc_module.getPath("gmshBin",fnPath=fnPath)
	
	# Just in case, replace backslashes with slashes
	gmshBin=gmshBin.replace('\\','/')
	
	return gmshBin

def addGmshToPATHs(fnPath=None):
	
	"""Exports gmsh binary defined in path file to current session.
	
	Uses ``getGmshBin`` to get path to gmsh binary and then sets environment variable.
	
	"""
	
	# Get gmsh bin as defined in path file
	fn=os.path.dirname(getGmshBin(fnPath=fnPath))
	
	return pyfrp_misc_module.addToPATHs(fn)
	
	
	
	
	
	
	
	
	
	
	
	