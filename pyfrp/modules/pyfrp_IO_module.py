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

"""Input/Output module for PyFRAP toolbox. 

Handles saving/loading PyFRAP projects into pickled files and the memory handling that comes with it.
"""

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

import pickle
import platform
import gc
import sys
import os
import csv
import shutil

from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_gmsh_IO_module
from pyfrp.modules.pyfrp_term_module import *

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

def saveToPickle(obj,fn=None):
	
	"""Saves obj into pickled format.
	
	.. note:: If ``fn==Non``, will try to save to ``obj.name``, otherwise unnamed.pk
	
	Keyword Args:
		fn (str): Output file name.	
	
	Returns: 
		str: Output filename.
	
	"""
	
	cleanUp()
	if fn==None:
		if hasattr(obj,"name"):
				fn=obj.name+".pk"
		else:
				fn="unnamed"+".pk"
		
	with open(fn, 'wb') as output: 
		pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
	
	return fn

def loadFromPickle(fn):
	
	"""Loads obj from pickled format.
	
	Args:
		fn (str): Filename.	
	
	Returns: 
		str: Output filename.
	
	"""
	
	cleanUp()
	
	#Need to do append subclasses folder here. Sometimes pickle has problem finding the classes
	
	sys.path.append(pyfrp_misc_module.getSubclassesDir()+'/')

	if platform.system() in ["Darwin","Linux"]:
			filehandler=open(fn, 'rb')
	elif platform.system() in ["Windows"]:
			filehandler=open(fn, 'rb')
			
	loadedFile=pickle.load(filehandler)
	
	return loadedFile

def loadMolecule(fn,update=True):
	
	"""Loads molecule object from pickle file
	and brings it up-to-date.
	
	Args:
		fn (str): Filename.	
	
	Keyword Args: 
		update (bool): Update to current version.
	
	Returns: 
		pyfrp.subclasses.pyfrp_molecule: Molecule file.
	
	"""
	
	mol=loadFromPickle(fn)
	if update:
		mol.update_version()
	return mol

def loadEmbryo(fn,update=True):
	
	"""Loads embryo object from pickle file
	and brings it up-to-date.
	
	Args:
		fn (str): Filename.	
	
	Keyword Args: 
		update (bool): Update to current version.
	
	Returns: 
		pyfrp.subclasses.pyfrp_embryo: Embryo file.
	
	"""
	
	emb=loadFromPickle(fn)
	if update:
		emb.update_version()
	return emb

def cleanUp():
	"""Calls garbage collector to clean up.
	"""
	
	gc.collect()
	return None

def copyMeshFiles(fn,fnGeo,fnMsh,debug=False):
	
	"""Copies meshfiles to new location. 
	
	If ``fn`` does not end on ``meshfiles``, will create a folder ``meshfiles`` 
	where to dump new files.
	
	.. note:: If ``fnGeo`` is a merged file, will try to copy all used .geo and .msh files
	   and also update the merged file such that it refers to the new mesh files.
	
	Args:
		fn (str): Filepath or parent directory where to put meshfiles.
		fnGeo (str): Filepath of geo file.
		fnMsh (str): Filepath of msh file.
	
	Keyword Args: 
		debug (bool): Print out debugging messages.
	
	Returns:
		tuple: Tuple containing:
		
			* fnGeoNew (str): New geo file location.
			* fnMshNew (str): New msh file location.
			
	"""
	
	# Check if fn is dir or not, if not, grab parenting directory
	if not os.path.isdir(fn):
		fn=os.path.realpath(fn)
		fn=os.path.dirname(fn)
	
	# Create meshfile folder if necessary
	if 'meshfiles'==os.path.split(fn)[-1]:
		if debug:
			print("Folder is already called meshfiles, will leave it as it is.")
	else:
		fn=pyfrp_misc_module.slashToFn(fn)
		try:
			os.mkdir(fn+"meshfiles")
		except OSError:
			if debug:
				printWarning("Cannot create folder " + fn + ". Already exists")
			
		fn=pyfrp_misc_module.slashToFn(fn+"meshfiles")
		if debug:
			print("Created new folder " + fn + " .")
	
	# Get list of files we want to copy
	files=[fnGeo,fnMsh]
	
	# If it is a merged file, grab all included .msh files
	isMerge,mergedFiles=pyfrp_gmsh_IO_module.isMergeFile(fnGeo)
	files=files+mergedFiles
	
	# Also try to grab all corresponding .geo files
	geoFiles=[]
	for m in mergedFiles:
		b,geoFile=pyfrp_gmsh_IO_module.getCorrespondingGeoFile(m)
		if b:
			geoFiles.append(geoFile)
	files=files+geoFiles
	
	# Finally copy
	for f in files:
		try:
			shutil.copy(f,fn)
		except:
			printWarning("copyMeshFiles: Was not able to copy file " + f )
		
	# If it is a merged file, we have to update the merged geo file, since the file
	# location of the merged mesh files changed
	if isMerge:
		mergedFilesNew=[]
		for m in mergedFiles:
			mergedFilesNew.append(fn+os.path.split(m)[-1])
		fnGeoNew=fn+os.path.split(fnGeo)[-1]
		
		fnGeoNew,fnMshNew=pyfrp_gmsh_IO_module.mergeMeshes(mergedFilesNew,fnGeoNew)
	else:
		fnGeoNew=fn+os.path.split(fnGeo)[-1]
		fnMshNew=fn+os.path.split(fnMsh)[-1]
		
	return fnGeoNew,fnMshNew
	
def copyAndRenameFile(fn,fnNew,debug=False):
	
	"""Copies file ``fn`` into same directory as ``fn`` and 
	renames it ``fnNew``.
	
	.. note:: If copying fails, then function will return old filename.
	
	Args:
		fn (str): Filepath of original file.
		fnNew (str): New filename.
	
	Keyword Args: 
		debug (bool): Print out debugging messages.
	
	Returns:
		str: Path to new file.
	
	"""
	
	head,tail=os.path.split(fn)
	ext=os.path.splitext(tail)[-1]
	fnNew=pyfrp_misc_module.slashToFn(head)+fnNew+ext
	
	if debug:
		ret=os.system("cp -v " + fn + " " + fnNew)
	else:
		ret=os.system("cp " + fn + " " + fnNew)
	
	if ret==0:
		return fnNew
	else:
		return fn
	
def writeTableToCSV(l,header,fn,col=False):
	
	"""Writes table to csv file.
	
	If ``col=True``, columns are given via ``l``, otherwise rows are given.
	
	Args:
		l (list): List of rows or columns.
		header (list): List of headers.
		col (bool): Flag on how rows/columns are given.
	
	Returns:
		tuple: Tuple containing:
			
			* header (list): Header of table.
			* table (list): Table as a list of rows.
	
	"""
	
	if col:
		
		lengths=[]
		for i in range(len(l)):
			lengths.append(len(l[i]))
		maxL=max(lengths)
		
		table=[]
		for i in range(maxL):
			row=[]
			for j in range(len(l)):
				try:
					row.append(l[j][i])
				except IndexError:
					row.append("")
			table.append(row)
	else:
		table=l
	
	with open(fn,'wb') as csvFile:
		fcsv = csv.writer(csvFile,delimiter=',')
		
		fcsv.writerow(header)
	
		for row in table:
			fcsv.writerow(row)
	
	return header,table
			
def readCSV(fn,hasHeader=True,delimiter=',',dtype=str):
	
	"""Reads csv sheet.
	
	Args:
		fn (str): Path to file.
		
	Keyword Args:	
		hasHeader (bool): File has header.
		delimiter (str): Delimiter of csv sheet.
		dtype (dtype): Datatype to convert array to.
	
	Returns:
		tuple: Tuple containing:
		
			* header (list): List of header entries.
			* rows (list): Read rows.
	
	"""
	
	rows=[]
	
	with open(fn,'rb') as f:
		
		fcsv=csv.reader(f,delimiter=delimiter)
		
		for i,row in enumerate(fcsv):
			
			
			
			if i==0:
				if hasHeader:
					header=row
				else:
					header=[]
					rows.append(pyfrp_misc_module.listToDtype(row,dtype))
			else:
				rows.append(pyfrp_misc_module.listToDtype(row,dtype))
			
	return header,rows	
	

