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

"""PyFRAP module for reading/writing gmsh .geo files. Module mainly has the following features:

	* Read .geo files.
	* Translate geometric entities and variables defined in .geo files.
	* Construct :py:class:`pyfrp.pyfrp_gmsh_geometry.domain` object describing complete geometry.
	* Update parameters in .geo files.
	* Add/Remove some geometric entities.
	* Add/update box fields to allow refinement of certain ROIs in mesh.

This module together with pyfrp.pyfrp_gmsh_geometry and pyfrp.pyfrp_gmsh_module works partially as a python gmsh wrapper, however is incomplete.
If you want to know more about gmsh, go to http://gmsh.info/doc/texinfo/gmsh.html .
	
"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy (use indirect import here, so convertMathExpr relates to numpy functions automatically when translating)
from numpy import *

#PyFRAP Modules
from . import pyfrp_gmsh_geometry
from . import pyfrp_misc_module
from . import pyfrp_gmsh_module

from .pyfrp_term_module import *

#Misc
import shutil
from tempfile import mkstemp
import os

from stl import mesh as meshstl

                   
#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================


def splitLine(line,delim="=",closer=";"):
	
	"""Splits line at ``delim``, trimming ``closer``.
	
	Example:
	
	>>> splitLine("Point(3)={1,3,1};")
	>>> ("Point(3)","{1,3,1}")
	
	Args:
		line (str): Line to be splitted.
	
	Keyword Args: 
		delim (str): Delimiter at which to be splitted.
		closer (str): Closing character to be trimmed.
		
	Returns:
		tuple: Tuple containing:
			
			* var (str): Name of variable.
			* val (str): Value of variable
		
	"""
	
	line=line.strip()
	var,val=line.split(delim)
	
	val=val.strip('\n')
	val=val.strip(closer)
	
	var=var.strip()
	val=val.strip()
	
	return var,val

def getId(var,delimOpen="(",delimClose=")"):
	
	"""Returns ID of object that is given between the delimiters ``delimOpen``
	and ``delimClose``.
	
	Example:
	
	>>> getId("Point(3)")
	>>> ("Point",3) 
	
	Args:
		var (str): String describing geoFile variable name.
	
	Keyword Args: 
		delimOpen (str): Openening delimiter of ID.
		delimClose (str): Closing delimiter of ID.
		
	Returns:
		tuple: Tuple containing:
			
			* typ (str): Type of geometric variable.
			* Id (str): Id of geometric variable.
		
	"""
	
	typ,Id=var.split(delimOpen)
	Id=int(Id.split(delimClose)[0])
	return typ,Id

def getVals(val,parmDic,openDelim="{",closeDelim="}",sep=","):
	
	"""Translates value of parameter into list of floats.
	
	Uses parameter dictionary to translate predefined variables into
	floats.
	
	Example:
	
	>>> getVals("{10,3,5}")
	>>> [10.,3.,5.]
		
	Args:
		val (str): Value string of geometric variable.
		parmDic (dict): Parameter dictionary.
	
	Keyword Args:
		openDelim (str): Opening delimiter.
		closeDelim (str): Closing delimiter.
		sep (str): Seperator between values.
	
	Returns:
		rList (list): List of translated values.
	
	"""
	
	if openDelim in val:
		valList,i=pyfrp_misc_module.str2list(val,dtype="str",openDelim=openDelim,closeDelim=closeDelim,sep=sep)
		wasList=True
	else:
		valList=[val]
		wasList=False
	
	rList=[]
	
	for v in valList:
		v=v.strip()
		rList.append(applyParmDic(v,parmDic))
	if wasList:	
		return rList
	else:
		return rList[0]
	
def convertMathExpr(val):
	
	"""Converts math expressions from .geo syntax into python
	syntax.
	
	.. note:: Not all translations have been implemented yet. You can 
	   simply add here expressions by adding a translation to the 
	   translations list (``translations.append([CExpression,PythonExpression])``).
	
	"""
	
	translations=[]
	translations.append(["^","**"])
	translations.append(["Sqrt","sqrt"])
	
	for translation in translations:
		val=val.replace(translation[0],translation[1])
	
	return val

def applyParmDic(val,parmDic):
	
	"""Applies parameter dictionary to variable value.
	
	Example: 
	
	>>> parmDic={'radius',3}
	>>> applyParmDic('radius',parmDic)
	>>> 3
	
	And also applies mathemtical expressions:
	
	>>> parmDic={'radius',3}
	>>> applyParmDic('radius^2-radius',parmDic)
	>>> 6
	
	Args:
		val (str): Value string of geometric variable.
		parmDic (dict): Parameter dictionary.
		
	Returns:
		val (float): Evaluated value.
	
	"""
	
	keys,lengths=sortKeysByLength(parmDic)
	keys.reverse()
	
	val=convertMathExpr(val)

	for key in keys:
		val=val.replace(key,"("+str(parmDic[key])+")")	

	try:
		return eval(val)
	except NameError:
		printWarning("applyParmDic: Could not evaluate value" + str(val))
		return val
	
	
def sortKeysByLength(dic):
	
	"""Sorts dictionary by length of keys. 
	"""
	
	lengths=[]
	for key in list(dic.keys()):
		lengths.append(len(key))
	
	keys,lengths=pyfrp_misc_module.sortListsWithKey(list(dic.keys()),lengths)
	return keys,lengths
		
def readParameter(line,parmDic):
	
	"""Reads in parameter from line and translates values using ``parmDic``.
	
	Args:
		line (str): Line to be splitted.
		parmDic (dict): Parameter dictionary.
		
	Returns:
		tuple: Tuple containing:
		
			* var (str): Name of variable.
			* val (float): Value of variable
	
	"""
	
	var,val = splitLine(line)
	
	val=applyParmDic(val,parmDic)
	return var,val

def readGeoLine(line,parmDic,domain):
	
	"""Reads in line from .geo file. 
	
	Tries to extract type of geometric object and its parameters 
	and uses this to append a geomtric entity to ``domain``. 
	
	If ``line`` describes a parameter, stores parameter name and its value
	in ``parmDic``.
	
	Args:
		line (str): Line to be splitted.
		parmDic (dict): Parameter dictionary.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object, storing all geometric entities.
		
	Returns:
		tuple: Tuple containing:
		
			* parmDic (dict): Updated parameter dictionary.
			* typ (str): Object type type.
			* Id (int): ID of object.
			* vals (list): Values of object.
			* domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Updated domain object.
			
	"""
	
	if line.startswith('//'):
		#This line is a comment, return parmDic
		return parmDic, "comment", -1, [],domain
	
	if "[" in line or "]" in line:
		#This is a field line
		
		domain = readFieldLine(line,domain,parmDic)
		
		return parmDic, "field", -3, [], domain
		
	if "Background Field" in line:
		domain=initBkgdField(line,domain)
		
	if "{" in line or "}" in line:
		#This line is some sort of object
		
		var,val = splitLine(line)
		typ,Id = getId(var)
		vals=getVals(val,parmDic)
		
		if typ=="Point":
			domain.addVertex([vals[0],vals[1],vals[2]],Id=Id,volSize=vals[3])
		elif typ=="Line":
			v1,idx=domain.getVertexById(vals[0])
			v2,idx=domain.getVertexById(vals[1])
			domain.addLine(v1,v2,Id=Id)
		elif typ=="Circle":
			vstart,idx=domain.getVertexById(vals[0])		
			vcenter,idx=domain.getVertexById(vals[1])
			vend,idx=domain.getVertexById(vals[2])
			domain.addArc(vstart,vcenter,vend,Id=Id)
		elif typ=="BSpline":
			vstart,idx=domain.getVertexById(vals[0])		
			vcenter,idx=domain.getVertexById(vals[1])
			vend,idx=domain.getVertexById(vals[2])
			domain.addBSpline(vals,Id=Id)
		elif typ=="Line Loop":
			domain.addLineLoop(Id=Id,edgeIDs=vals)
		elif typ=="Ruled Surface":
			domain.addRuledSurface(Id=Id,lineLoopID=vals[0])
		elif typ=="Surface Loop":
			domain.addSurfaceLoop(Id=Id,surfaceIDs=vals)
		elif typ=="Volume":
			domain.addVolume(Id=Id,surfaceLoopID=vals[0])	
		else:
			#This is a object like surface or volume which we don't care about
			pass
		
		return parmDic, typ, Id , vals,domain
	
	
	else:
		#Check if line is empty
		if len(line.strip())>0:
			
			#This is a parameter line, append new parameter to parmDic
			var,val=readParameter(line,parmDic)
			
			parmDic[var]=val
			return parmDic, "parameter", -2, [],domain
		
		else:
			return parmDic, "empty", -3, [], domain

def initBkgdField(line,domain):
	
	"""Initiates background field when reading a .geo file.
	
	Args:
		line (str): Line in geo file.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object.
	
	Returns:	
		pyfrp.modules.pyfrp_gmsh_geometry.domain: Updated domain.
	
	"""
	
	var,val = splitLine(line)
	field=domain.getFieldById(int(val))[0]
	field.setAsBkgdField()
	
	return domain
	
def initField(val,domain,Id):
	
	"""Adds the right type of field object to domain.
	
	Args:
		val (str): Value string.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object.
		Id (int): Id of new field.
		
	Returns:	
		pyfrp.modules.pyfrp_gmsh_geometry.domain: Updated domain.
		
	"""
	
	if 'Attractor' in val:
		domain.addAttractorField(Id=Id)
	elif 'Threshold' in val:
		domain.addThresholdField(Id=Id)
	elif 'Min' in val:
		domain.addMinField(Id=Id)
	elif 'Box' in val:	
		domain.addBoxField(Id=Id)
	elif 'Boundary'	in val:
		domain.addBoundaryLayerField(Id=Id)
	return domain
	
def readFieldLine(line,domain,parmDic):
	
	"""Reads line that belongs to field definition in .geo file.
	
	If line defines new field, will create new field using 
	:py:func:`initField`.
	
	Otherwise will try to find field in domain and set new 
	property value.
	
	Args:
		line (str): Line to read.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object.
		parmDic (dict): Parameter dictionary.
	
	Returns:	
		pyfrp.modules.pyfrp_gmsh_geometry.domain: Updated domain.
		
	"""
	
	var,val = splitLine(line)
	typ,Id = getId(var,delimOpen='[',delimClose=']')
	
	if "." in var:
		vals=getVals(val,parmDic)
		temp,prop=var.split(".")
		prop=prop.strip()

		domain.getFieldById(Id)[0].setFieldAttr(prop,vals)
		
	else:
		domain=initField(val,domain,Id)
	
	return domain
	
def readGeoFile(fn):
		
	"""Reads in .geo file and tries to extract geometry defined in .geo file
	into a :py:class:`pyfrp.modules.pyfrp_gmsh_geometry.domain`.
	
	Args:
		fn (str): Filename of .geo file.
	
	Returns:
		tuple: Tuple containing:
		
			* parmDic (dict): Updated parameter dictionary.
			* domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): Domain object.
			
	"""
	
	#new parameter dictionary
	parmDic={}
	
	#New domain
	domain=pyfrp_gmsh_geometry.domain()
	
	#Read file
	with open(fn,'r') as f:
		for line in f:
			parmDic,typ,Id,vals,domain=readGeoLine(line,parmDic,domain)
			
	return domain,parmDic

def readStlFile(fn,domain=None,volSizePx=20.):
	
	"""Reads stl file to domain.
	
	.. note:: Uses numpy-stl package. You may need to install via 
	   ``pip install numpy-stl``
	
	If no domain is given, will create new one
	
	Args:
		fn (str): Path to stl file.
		
	Keyword Args:
		volSizePx (float): Mesh density assigned at vertices.
		domain (pyfrp.modules.pyfrp_gmsh_geometry.domain): A domain object.
		
	Returns:
		pyfrp.modules.pyfrp_gmsh_geometry.domain: A domain object.
	
	"""
	
	
	
	#Load file
	mesh=meshstl.Mesh.from_file(fn,speedups=False)
	
	#New domain
	if domain==None:
		domain=pyfrp_gmsh_geometry.domain()
	
	#Loop through all surface triangles
	for triang in mesh.data:
		
		#Loop through all vertices of triangle
		vertices=[]
		edges=[]
		for i,x in enumerate(triang[1]):
			
			#Create vertex if it doesn't exist yet
			v,ind=domain.getVertexByX(x)
			if v==False:
				v=domain.addVertex(x,volSize=volSizePx)
			vertices.append(v)
			
			#Create edge if non-existent
			if i>0:	
				edge,ind=domain.getEdgeByVertices(vertices[i-1],vertices[i])
				if edge==False:
					edges.append(domain.addLine(vertices[i-1],vertices[i]))
				else:
					edges.append(edge)
					
			if i==len(triang[1])-1:
				edge,ind=domain.getEdgeByVertices(vertices[i],vertices[0])
				if edge==False:	
					edges.append(domain.addLine(vertices[i],vertices[0]))
				else:
					edges.append(edge)
			
		#Add line loop	
		loop=domain.addLineLoop(edgeIDs=pyfrp_misc_module.objAttrToList(edges,"Id"))
		loop.fix()
		
		#Add surface
		surface=domain.addRuledSurface(lineLoopID=loop.Id)
		surface.normal=triang[0]/np.linalg.norm(triang[0])
		
	return domain		
				
	

def updateParmGeoFile(fn,name,val):
	
	"""Updates parameter in .geo file.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		fn (str): Filename of .geo file.
		name (str): Name of parameter.
		val (float): Value of parameter.
			
	"""
		
	
	substr=name+"="+str(val)+";"+'\n'
	pyfrp_misc_module.txtLineReplace(fn,name,substr)
	
	return

def getAllIDsOfType(fn,elementType):
	
	"""Finds all IDs of a specific .geo element type in a .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementType (str): Type of parameter, for example ``"Point"``.
		
	Returns:
		list: List of IDs.
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	Ids=[]
	
	f = open (fn,'rb')
	for line in f:
		
		if line.strip().startswith(elementType):
			
			var,val = splitLine(line)
			if elementType=="Field":
				typ,Id = getId(var,delimOpen="[",delimClose="]")
			else:
				typ,Id = getId(var)
			
			Ids.append(Id)
	Ids=list(unique(Ids))			
	
	f.close()
	return Ids

def findComment(fn,comment):
	
	"""Finds a specific comment in .geo file and returns
	line in which it appears, otherwise -1.
	
	.. note:: Will only look for an exact match.
	
	Args:
		fn (str): Filename of .geo file.
		comment (str): Comment to look for.
		
	Returns:
		int: Line number of appearance.
	
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return -1
	
	f = open (fn,'rb')
	for i,line in enumerate(f):
		if line.strip().startswith("//"):
			if line.strip().replace("//","")==comment:
				return i
	
	return -1

def getFieldByComment(fn,comment,lineDiff=3):
	
	"""Returns field that is preceeded by comment ``comment``.
	
	.. note:: Will only look for an exact match.
	
	Args:
		fn (str): Filename of .geo file.
		comment (str): Comment to look for.
	
	Keyword Args:
		lineDiff (int): Maximum allowed difference between line of comment and field.
	
	Returns:
		int: Line number of appearance.

	"""
	
	#Get all field IDs
	fieldIDs=getAllIDsOfType(fn,"Field")
	
	#Find line of comment
	lineComment=findComment(fn,comment)

	for fieldID in fieldIDs:
		lines=getLinesByID(fn,fieldID,"Field")
		
		if lines[0]-lineComment<lineDiff:
			return fieldID, lines
		
	return -1,[]	
	
	
def getLargestIDOfType(fn,elementType):
	
	"""Finds largest ID of a specific .geo element type in a .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementType (str): Type of parameter, for example ``"Point"``.
		
	Returns:
		int: Largest ID.
	"""
	
	return max(getAllIDsOfType(fn,elementType))

def getBkgdFieldID(fn):
	
	"""Finds ID of background field in .geo file.
	
	.. note:: Will return ``None`` if .geo file has no background
	   field specified.
	
	Args:
		fn (str): Filename of .geo file.
			
	Returns:
		int: ID of background field.
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	f = open (fn,'rb')
	for line in f:
		if line.strip().startswith("Background Field"):
			var,val = splitLine(line)
			f.close()
			return int(val.strip())
	
	
	f.close()	
	return 

def getLastNonEmptyLine(fn):
	
	"""Finds index of last non-empty line in .geo file.
	
	Args:
		fn (str): Filename of .geo file.
			
	Returns:
		int: Index of last non-empty line.
	"""
	
	idx=0
	with open(fn,'rb') as f:
		for i,line in enumerate(f):
			if len(line.strip()):
				idx=i
	return idx
	
def removeTailingLines(filePath,idx):
	
	"""Removes all empty lines at the end of a .geo file.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	Args:
		filePath (str): Filename of .geo file.
		idx (int): Index of last non-empty line
			
	"""
	
	#Create temp file
	fh, absPath = mkstemp()
	newFile = open(absPath,'w')
	oldFile = open(filePath)
	
	#Loop through file and only write until line idx
	for i,line in enumerate(oldFile):
		
		if i<=idx:
			newFile.write(line)
			
	#close temp file
	newFile.close()
	os.close(fh)
	oldFile.close()
		
	#Remove original file
	os.remove(filePath)
	
	#Move new file
	shutil.move(absPath, filePath)
	return		
	
def copyIntoTempFile(fn,close=True):
	
	"""Copies file into tempfile.
	
	.. note:: Will create temporary file using ``tempfile.mkstemp()``. You should have 
	   read/write access to whereever ``mkstemp`` is putting files.
	
	.. note:: If ``close==True``, will return ``fh=None`` and ``tempFile=None``.
	
	Args:
		fn (str): Filename of file.
	
	Keyword Args:
		close (bool): Close files after copying.
	
	Returns:
		tuple: Tuple containing:
		
			* tempFile (file): File handle to temp file.
			* fh (file): File handle to original file.
			* tempPath (tempPath): Path to temp file.
	
	"""
	
	oldFile = open(fn)
	fh, tempPath = mkstemp()
	tempFile = open(tempPath,'w')
	
	for line in oldFile:
		tempFile.write(line)

	oldFile.close()
	
	if close:
		tempFile.close()
		os.close(fh)
		tempFile=None
		fh=None
		
	return tempFile, fh,tempPath

def getLinesByID(fn,elementId,elementType=""):
	
	"""Finds all lines in .geo file that contain geometric entitity with ID ``elementId``.
	
	.. note:: IDs in geometric files can be given per entitity type. That is, one can have
	   for example a point with ID=1 (``Point(1)``) aswell as a line with ID=1 (``Line(1)``).
	   Thus one may want to use ``elementType`` to restrict the search for a specific element type.
	
	Args:
		fn (str): Filename of .geo file.
		elementId (int): ID to look for.
	
	Keyword Args:
		elementType (str): Type of element to restrict search on.
	
	Returns:
		list: Line numbers at which element appears.
	
	"""
	
	if not os.path.isfile(fn):
		printWarning(fn + " does not exist.")
		return
	
	f = open (fn,'rb')
	
	lineNumbers=[]
	
	for i,line in enumerate(f):
		if line.strip().startswith(elementType):
			
			var,val = splitLine(line)
			if elementType=="Field":
				typ,Id = getId(var,delimOpen="[",delimClose="]")
			else:
				typ,Id = getId(var)
			
			if elementId==Id:
				lineNumbers.append(i)
	
	f.close()
	return lineNumbers
				
def removeCommentFromFile(fn,comment):
	
	"""Removes comment ``comment`` from .geo file. 
	
	.. note:: Will remove all appearances of ``comment``.
	
	.. note:: Will also remove comments that only start with ``comment``.
	
	Args:
		fn (str): Filename of .geo file.
		comment (str): Comment to remove
	
	"""
	
	pyfrp_misc_module.txtLineReplace(fn,"//"+comment,"")
	return

def removeElementFromFile(fn,elementType,elementId,delimOpen="(",delimClose=")"):
	
	"""Removes element with type ``elementType`` and ID ``elementID`` from .geo file.
	
	Args:
		fn (str): Filename of .geo file.
		elementId (int): ID of element to remove.
		elementType (str): Type of element to remove.
	
	Keyword Args:
		delimOpen (str): Openening delimiter of ID.
		delimClose (str): Closing delimiter of ID.
	
	"""
	
	pyfrp_misc_module.txtLineReplace(fn,elementType+delimOpen+str(elementId)+delimClose,"")
	return

	
def addBoxField(fn,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ,comment="",fnOut="",overwrite=True,sameComment=True):		
	
	"""Adds box field to .geo file by doing the following:
		
		* Copies file into temp file for backup using :py:func:`copyIntoTempFile`.
		* Finds all IDs of previous ``Field`` entities using :py:func:`getAllIDsOfType`.
		* Finds current background field using :py:func:`getBkgdFieldID` .
		* If previous fields exist, removes them from file using :py:func:`removeElementFromFile` .
		* Finds last non-empty line using  :py:func:`getLastNonEmptyLine` .
		* Removes empty lines at end of file using :py:func:`removeTailingLines` .
		* Writes comment using  :py:func:`writeComment` .
		* Writes box field using  :py:func:`writeBoxField` .
		* Writes background field using  :py:func:`writeBackgroundField` .
		
	.. note:: Comment is useful to describe in .geo file what the the box field actually does.
	
	.. note:: Generally, background field will use ``volSizeIn`` as background mesh volume size.
	
	.. note:: Unit for parameter is pixels.
	
	.. note:: If ``fnOut`` is not specified, will overwrite input file.
	
	.. note:: Will always remove previous background fields. If ``overwrite=True``, will remove all fields.
	   If additionally ``sameComment=True``, will look for the field that has the same comment as ``comment``
	   and only remove this particular one.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		fn (str): Filename of .geo file.
		volSizeIn (float): Mesh element volume inside box.
		volSizeOut (float): Mesh element volume outside box.
		rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
		rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
		rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
		
	Keyword Args:
		comment (str): Comment to be added before box field.
		fnOut (str): Filepath for output.
		overwrite (bool): Overwrite previously exisiting box fields.
		sameComment (bool): Only remove box field with particular comment.

	"""
	
	#Copy everything into tempfile
	tempFile,fh,tempPath = copyIntoTempFile(fn,close=True)
	
	#Find if there are already Fields defined
	fieldIDs=getAllIDsOfType(tempPath,"Field")
	
	#Find background fields
	bkgdID=getBkgdFieldID(tempPath)
	
	#Remove bkgdID from field IDs
	if bkgdID!=None:
		fieldIDs.remove(bkgdID)
	
	#If there is already a background field, remove all lines containing it
	if bkgdID!=None:
		removeElementFromFile(tempPath,"Field",bkgdID,delimOpen="[",delimClose="]")
		removeElementFromFile(tempPath,"Mesh.","",delimOpen="",delimClose="")
		removeElementFromFile(tempPath,"Background Field =",bkgdID,delimOpen="",delimClose="")
		
	#Remove other fields if selected:
	if overwrite:
		
		if sameComment:
			sameID,sameLines=getFieldByComment(tempPath,comment)
			if sameID>-1:
				removeElementFromFile(tempPath,"Field",sameID,delimOpen="[",delimClose="]")
				removeCommentFromFile(tempPath,comment)
				fieldIDs.remove(sameID)
		else:	
			for fieldID in fieldIDs:
				removeElementFromFile(tempPath,"Field",fieldID,delimOpen="[",delimClose="]")
		
	#Get index of last non-empty line
	idxNonEmpty=getLastNonEmptyLine(tempPath)
	
	#remove tailing lines
	removeTailingLines(tempPath,idxNonEmpty)
	
	#Open file again for appending new lines
	with open(tempPath,'a') as f:
		
		#Write Empty line
		f.write('\n')
		
		#Write Comment
		f=writeComment(f,comment)
		
		#Get ID of new field
		if len(fieldIDs)>0:
			newFieldID=max(fieldIDs)+1
		else:
			newFieldID=1
		
		#Write new box field
		f=writeBoxField(f,newFieldID,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ)
		
		#Append to field ids
		fieldIDs.append(newFieldID)
		
		#Write minimum field
		f=writeMinField(f,max(fieldIDs)+1,fieldIDs)
		
		#Write Background field
		f=writeBackgroundField(f,max(fieldIDs)+1)
			
	#Move new file either to fnOut or to fn itself
	if fnOut!="":
		shutil.move(tempPath, fnOut)
	else:
		shutil.move(tempPath, fn)
	
	return	
	
def writeComment(f,comment):
	
	"""Writes comment line into file.
	
	Args:
		f (file): Filehandle.
		comment (str): Comment to be written.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("//"+comment+"\n")
	return f
	
def writeBackgroundField(f,fieldID):
	
	"""Writes background field into into file.
	
	.. note:: Will take finest mesh for background field.
	   See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new background field.
	
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Background Field ="+str(fieldID) +";"+"\n")
	f.write('\n')
	return f

def writeMinField(f,fieldID,ids,charExtendFromBoundary=True):	
	
	"""Writes minimum field into into file.
	
	.. note:: Useful to determine background mesh. It's often reasonable to take the finest mesh for background field.
	   See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new background field.
		ids (list): List of field IDs used for background mesh computation.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Mesh.CharacteristicLengthExtendFromBoundary = "+str(int(charExtendFromBoundary))+";"+"\n")
	
	f.write("Field["+str(fieldID)+"] = Min"+";"+"\n")
	f.write("Field["+str(fieldID)+"].FieldsList = {")
	for i,d in enumerate(ids):
		f.write(str(d))
		if i<len(ids)-1:
			f.write(",")
	f.write("}"+";"+"\n")
	f.write('\n')
	
	return f
	
def writeBoxField(f,fieldID,volSizeIn,volSizeOut,rangeX,rangeY,rangeZ):	
	
	"""Writes box field into into file.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new box field.
		volSizeIn (float): Mesh element volume inside box.
		volSizeOut (float): Mesh element volume outside box.
		rangeX (list): Range of box field in x-direction given as ``[minVal,maxVal]``.
		rangeY (list): Range of box field in y-direction given as ``[minVal,maxVal]``.
		rangeZ (list): Range of box field in z-direction given as ``[minVal,maxVal]``.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Field["+str(fieldID)+"] = Box"+";"+"\n")
	f.write("Field["+str(fieldID)+"].VIn = "+str(volSizeIn)+";"+"\n")
	f.write("Field["+str(fieldID)+"].VOut = "+str(volSizeOut)+";"+"\n")
	f.write("Field["+str(fieldID)+"].XMin = "+str(rangeX[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].XMax = "+str(rangeX[1])+";"+"\n")
	f.write("Field["+str(fieldID)+"].YMin = "+str(rangeY[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].YMax = "+str(rangeY[1])+";"+"\n")
	f.write("Field["+str(fieldID)+"].ZMin = "+str(rangeZ[0])+";"+"\n")
	f.write("Field["+str(fieldID)+"].ZMax = "+str(rangeZ[1])+";"+"\n")
	f.write('\n')
	return f

def writeAttractorField(f,fieldID,NodesList):
	
	"""Writes attractor field into into file.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes and
	:py:class:`:pyfrp.modules.pyfrp_gmsh_geometry.attractorField`.
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new box field.
		NodesList (list): List of vertex IDs at which attractor is placed.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Field["+str(fieldID)+"] = Attractor"+";"+"\n")
	writeFieldProp(f,fieldID,"NodesList",str(NodesList).replace("[","{").replace("]","}"))
	f.write("\n")
	
	return f
		
def writeThresholdField(f,fieldID,IField,LcMin,LcMax,DistMin,DistMax):		
	
	"""Writes threshold field into into file.
	
	See also: http://gmsh.info/doc/texinfo/gmsh.html#Specifying-mesh-element-sizes . 
	and :py:class:`:pyfrp.modules.pyfrp_gmsh_geometry.thresholdField`.
	
	Args:
		f (file): Filehandle.
		fieldID (int): ID of new box field.
		IField (int): ID of vertex that is center to threshold field.
		LcMin (float): Minimum volSize of threshold field.
		LcMax (float): Maximum volSize of threshold field.
		DistMin (float): Minimun density of field.
		DistMax (float): Maximum density of field.
		
	Returns:
		file: Filehandle.
	
	"""
	
	f.write("Field["+str(fieldID)+"] = Threshold"+";"+"\n")
	
	#writeFieldPropByDict(f,fieldID,IField=IField,LcMin=LcMin,LcMax=LcMax,DistMin=DistMin,DistMax=DistMax)
	
	writeFieldProp(f,fieldID,"IField",IField)
	writeFieldProp(f,fieldID,"LcMin",LcMin)
	writeFieldProp(f,fieldID,"LcMax",LcMax)
	writeFieldProp(f,fieldID,"DistMin",DistMin)
	writeFieldProp(f,fieldID,"DistMax",DistMax)
	f.write("\n")
	
	return f

def writeBoundaryLayerField(f,fieldID,elements,fieldOpts):
	
	"""Writes boundary layer mesh.
	
	"""
	
	f.write("Field["+str(fieldID)+"] = BoundaryLayer"+";"+"\n")
	
	f=writeFieldPropByDict(f,fieldID,fieldOpts)
	f=writeFieldPropByDict(f,fieldID,elements)
	
	return f

def writeFieldPropByDict(f,fieldID,dic):
	
	"""Writes dictionary of field properties to file.
	
	Args:
		f (file): File to write to.
		fieldID (int): ID of field.
		
	Keyword Args:
		dic (dict): Keyword Arguments.
	
	Returns:
		file: Filehandle.
	
	"""
	
	for key, value in dic.items():
		f=writeFieldProp(f,fieldID,key,value)
		
	return f
	
def writeFieldProp(f,fieldID,prop,val):

	"""Writes field property to file.
	
	Args:
		f (file): File to write to.
		fieldID (int): ID of field.
		prop (str): Name of property to write.
		val (str): Value to write.
	
	Returns:
		file: Filehandle.
	
	"""

	val=str(val).replace("[","{").replace("]","}")
	f.write("Field["+str(fieldID)+"]."+prop+" = "+ str(val) + ";\n")
	
	return f
	
def repairDefaultGeoFiles(debug=False):
	
	"""Copies default geometry files from backup folder to meshfile folder.
	Useful if geometry files got somehow overwritten or corrupted.
	
	Keyword Args:
		debug (bool): Print debugging messages.
	
	Returns:
		bool: ``True`` if success, ``False`` otherwise.
	
	"""
	
	fnBackup=pyfrp_misc_module.getMeshfilesDir()+"backup"
	
	if debug:
		cmd = "cp -v " + fnBackup + "*.geo " + pyfrp_misc_module.getMeshfilesDir()
		print("Executing command:")
		print(cmd)
	else:
		cmd = "cp " + fnBackup + "*.geo " + pyfrp_misc_module.getMeshfilesDir()
	
	ret=os.system(cmd)
	
	if ret>0:
		return False
	else:
		return True

def getCorrespondingGeoFile(fn,meshFileExt=".msh"):
	
	"""Returns the corresponding geometry file 
	to a mesh file.
	
	Assumes that meshfile has same name and lives in the same folder.
	
	Args:
		fn (str): Path to mesh file.
		
	Keyword Args:
		meshFileExt (str): Extension of meshfile.
	
	Returns:
		tuple: Tuple containing:
		
			* exists (bool): Flag if corresponding file exits.
			* fnGeo (str): Path to geometry file.
		
	"""
	
	fnGeo=fn.replace(meshFileExt,".geo")
	
	if os.path.isfile(fnGeo):
		return True,fnGeo
	else:
		return False,""
	

def isMergeFile(fn):
	
	"""Checks if geo file is a merge file and returns 
	and returns the meshes used inside the merge file.
	
	Args:
		fn (str): Path to geo file.
			
	Returns:
		tuple: Tuple containing:
		
			* isMerge (bool): Flag if corresponding file is merge file.
			* mergedFiles (list): List of mesh files used in merge file.
		
	"""
	
	isMerge=False
	mergedFiles=[]
	
	try:
		with open(fn,'rb') as f:
			
			for line in f:
				if line.startswith("Merge"):
					isMerge=True
					splitted=line.split('"')
					mergedFiles.append(splitted[1])
	except TypeError:
		with open(fn,'r') as f:
			
			for line in f:
				if line.startswith("Merge"):
					isMerge=True
					splitted=line.split('"')
					mergedFiles.append(splitted[1])
		
	return isMerge,mergedFiles	

def genMergeGeoFile(meshFiles,fnGeo):
	
	"""Generates merged .geo file.
	
	Args:
		meshFiles (list): List of meshfiles that will be included.
		fnGeo (str): Output geometry file.
		
	Returns:
		bool: True if no error/warning occured.
	"""
	
	b=True
	with open(fnGeo,'wb') as f:
		
		for fn in meshFiles:
			if os.path.isfile(fn):
				writeMergeLine(f,fn)
			else:
				b=False
				printWarning("genMergeGeoFile: Mesh file "+fn+" cannot be detected. This might lead to problems.")	
		f.write("Coherence Mesh;")
		
	return b	
		
def writeMergeLine(f,fnMsh):
	
	"""Adds merge line to geo file.
	
	Args:
		f (file): File handle.
		fnMsh (str): Mesh file to be added.
	
	Returns:
		file: File handle.
	
	"""
	
	f.write('Merge "' +  fnMsh +'";\n' )
	
	return f
		
def mergeMeshes(meshFiles,fn,run=True,debug=False,redirect=False,fnStout=None,fnSterr=None,volSizeMax=None):
	
	"""Generates meshfile merging all meshes in meshFiles.
	
	If one of the files that is supposed to be merged is already a merged file, then this
	function will try to find the original mesh files to write a complete merged file.
	See also :py:func:`isMergeFile` and :py:func:`genMergeGeoFile`.
	
	If ``run==True`` is selected, then gmsh will be run via :py:func:`pyfrp.modules.pyfrp_gmsh_module.runGmsh` and 
	generate the corresponding .msh file of the merged .geo file.
	
	Args:
		meshFiles (list): List of path to mesh files.
		fn (str): Name of output .geo file.
	
	Keyword Args:	
		run (bool): Run gmsh on merged .geo file.
		debug (bool): Print debugging messages.
		redirect (bool): Redirect gmsh stout/sterr into seperate files.
		fnStout (str): File for gmsh stout.
		fnSterr (str): File for gmsh sterr.
		volSizeMax (float): Maximum allowed mesh element size.
		
	Returns:
		tuple: Tuple containing:
		
			* fn (str): Path to generated .geo file.
			* fnOut (str): Path to generated .msh file.

	"""
	
	meshFilesToMerge=[]
	
	for meshFile in meshFiles:
		
		meshFile=os.path.abspath(meshFile)
		
		b,geoFile=getCorrespondingGeoFile(meshFile)
		
		isMerge,mergedFiles=isMergeFile(geoFile)
		
		if isMerge:
			meshFilesToMerge=meshFilesToMerge+mergedFiles
		else:
			meshFilesToMerge.append(meshFile)
			
	genMergeGeoFile(meshFilesToMerge,fn)
	
	if run:
		fnOut=pyfrp_gmsh_module.runGmsh(fn,debug=debug,redirect=redirect,fnStout=fnStout,fnSterr=fnSterr,volSizeMax=volSizeMax)
	else:
		fnOut=""
		
	return fn,fnOut


	
	
	


	
