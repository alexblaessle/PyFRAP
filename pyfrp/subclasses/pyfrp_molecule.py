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

#Main molecule class saving all important data, parameters and results for particular molecule.
	
#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================

#Numpy
import numpy as np

#PyFRAP Modules
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules import pyfrp_stats_module

#PyFRAP Classes
from . import pyfrp_embryo

#Standard packages
import os


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Molecule object

class molecule:
	
	"""Molecule class, collecting information about a series of FRAP experiments.
	
	The main purpose of the molecule class is to gather and summarize multiple FRAP experiments.
	Embryo objects are stored in a ``embryos`` list. From those embryo objects, fit objects can be
	added to the ``selFits`` list to then be summarized to calculate measurement statistics.
	Fits can be forced to overlap in a set of parameters defined in ``crucialParameters``.
	
	"""
	
	
	#Creates new molecule object
	def __init__(self,name):
		
		#Data
		self.name=name
		self.embryos=[]
		
		#Results
		self.selFits=[]
		self.DOptMu=None
		self.DOptPx=None
		self.DOptPxStd=None
		self.DOptPxSterr=None
		self.prodOpt=None
		self.degrOpt=None
		self.MeanRsq=None
		self.Rsq=None
		self.DMuStd=None
		self.prodStd=None
		self.degrStd=None
		self.DMuStErr=None
		self.prodStErr=None
		self.degrStErr=None
		
		self.crucialParameters=["equOn","fitPinned","fitProd","fitDegr","LBD","LBProd","LBDegr","UBD","UBProd","UBDegr"]
					
	def addEmbryo(self,embryo):
		
		"""Appends embryo object to ``embryos`` list.
		
		Args:
			embryo (pyfrp.subclasses.pyfrp_embryo): Embryo to append.
		
		Returns:
			list: Updated pyfrp.subclasses.pyfrp_molecule.embryos list
		
		"""
		
		self.embryos.append(embryo)
		return self.embryos
	
	def replaceEmbryo(self,embryo,name=""):
		
		"""Replaces embryo with same name that ``embryo`` in molecule's ``embryos`` list.
		
		If no embryo with the same name exists, will simply add embryo. If ``name`` is given,
		will try to replace embryo with name ``name`` in ``embryos`` list.
		
		Args:
			embryo (pyfrp.subclasses.pyfrp_embryo): New embryo to be inserted.
			
		Keyword Args:
			name (str): Optional name of embryo embryo to be replaced.

		Returns:
			bool: True if anything was replaced, False if anything was appended.
			
		"""
		
		for i,emb in enumerate(self.embryos):
			if name!="":
				if emb.name==embryo.name:
					self.embryos[i]=embryo
					return True
			else:
				if emb.name==name:
					self.embryos[i]=embryo
					return True
			
		self.addEmbryo(embryo)
		
		return False
		
	def newEmbryo(self,name):
		
		"""Creates new embryo and appends it to ``embryos`` list.
		
		Args:
			name (str): Name of new embryo.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: New embryo object.
		
		"""
		
		emb=pyfrp_embryo.embryo(name)
		self.embryos.append(emb)
		return self.embryos[-1]
		
	def removeEmbryo(self,embryo):
			
		"""Removes embryo object from ``embryos`` list.
		
		Args:
			embryo (pyfrp.subclasses.pyfrp_embryo): Embryo object.
		
		Returns:
			list: Updated pyfrp.subclasses.pyfrp_molecule.embryos list
		
		"""
		
		self.embryos.remove(embryo)
		return self.embryos		
	
	def save(self,fn=None):
		
		"""Saves molecule to pickle file.
		
		.. note:: If ``fn`` is not specified, will assume ``fn=self.name``.
		
		Keyword Args:
			fn (str): Molecule file name.
		
		Returns:
			str: Filename of molecule file.
		
		"""
		
		if fn==None:
			fn=self.name+".mol"
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)	
		
		return fn
	
	def saveExtract(self,fn=None,copyMeshFiles=True,debug=False):
		
		"""Saves molecule to pickle file in compressed version by doing:
		
			* Extracts embryos in ``embryos`` list into seperate pickled files.
			* Clears all attributes of embryo objects
			* Saves molecule file
			
		This function is really useful if molecule file size gets out-of-hand.	
		
		.. note:: Embryo files will be saved in path/to/moculefile/moleculename/ .
		
		.. note:: If ``fn`` is not specified, will assume ``fn=self.name``.
		
		Keyword Args:
			fn (str): Molecule file name.
			copyMeshFiles (bool): Copy meshfiles to embryo file destination.
			debug (bool): Print out debugging messages.
			
		Returns:
			bool: True if success, False else.
		
		"""
		
		if not self.checkEmbryoNames():
			printError("Embryo names are not distinct, will not save.")
			return False
		
		if fn==None:
			fn=self.name+".mol"
		
		
		b=self.extractEmbryos2Files(fn=fn.replace('.mol','')+"_embryos",copyMeshFiles=copyMeshFiles,debug=debug)
		
		if not b:
			printError("Something went wrong extracting embryos, will not continue saving.")
			return False
		
		self.clearAllEmbryos()
		self.save(fn=fn)
		
		return True
		
		
	
	def extractEmbryos2Files(self,fn="",copyMeshFiles=True,debug=False):
		
		"""Extracts embryos in ``embryos`` list into seperate pickled files.
		
		.. note:: Will create folder ``fn`` if non-existent. If ``fn`` is not specified, will assume ``fn='embryoFiles/'`` .
		
		Keyword Args:
			fn (str): Path of folder where to save embryo files.
			copyMeshFiles (bool): Copy meshfiles to embryo file destination.
			debug (bool): Print out debugging messages.
			
		Returns:
			bool: True if success, False else.
		
		"""
		
		if fn=="":
			try:
				printWarning('fn not specified, will assume fn=embryoFiles/')
				os.system("mkdir -v embryoFiles/")
				fn="embryoFiles/"
			except OSError:
				pass
		
		if not os.path.isdir(fn):
			printWarning(fn+" does not  exist, will try to create it.")
			try:
				os.system("mkdir -v " +fn)
			except OSError:
				printError("Could not create " + fn+ ". Will not extract." )
				return False
		
		print("in mol.extract:", fn)
		
		for embryo in self.embryos:
			embryo.save(fn=fn+embryo.getName()+".emb",copyMeshFiles=copyMeshFiles,debug=debug)
		
		return True
		
	def clearAllEmbryos(self):
		
		"""Replaces all attribute values of each embryo in `embryos` list with ``None``, except ``name``.
		
		Useful if embryos are seperated and molecule file needs to be compressed.
		
		.. note:: Embryos should have all different names, so there will not be any missassignment when reimporting embryo files.
		
		Returns:
			bool: True if success, False else.
		
		"""
		
		b=True
		for embryo in self.embryos:
			bnew=embryo.clearAllAttributes()
			b=b+bnew
		
		return b
	
	def checkEmbryoNames(self):
		
		"""Check if all embryos in ``embryos`` list have different names.
		
		Returns:
			bool: True if all different, False else.	
		"""
		
		names=pyfrp_misc_module.objAttrToList(self.embryos,'name')
		b=True
		for name in names:
			if names.count(name)>1:
				printWarning('There is more than one embryo called '+ name + ' in molecule ' + self.name)
				b=False
		return b		
		
	def updateVersion(self):
		
		"""Updates molecule file to current version, making sure that it possesses
		all attributes.
		
		Creates a new molecule object and compares ``self`` with the new molecule file.
		If the new molecule object has a attribute that ``self`` does not have, will
		add attribute with default value from the new molecle file.
		
		.. note:: Will also update all subobject, making sure that embryo and fit objects are up-to-date.
		
		Returns:
			pyfrp.subclasses.pyfrp_molecule: ``self``
			
		"""
		
		#Create temporarly a blank molecule file
		moltemp=molecule("temp")
		
		#Update molecule file
		pyfrp_misc_module.updateObj(moltemp,self)
		
		#Update all embryos and fits
		for emb in self.embryos:
			
			emb.updateVersion()
			
			for fit in emb.fits:
				fit.updateVersion()
		
		return self
		
	def getEmbryoByName(self,s):
		
		"""Returns embryo with name ``s`` from ``embryos`` list, otherwise ``False``. 
		"""
		
		for emb in self.embryos:
			if str(s) == emb.name:
				return emb
		return False
	
	def setName(self,n):
		
		"""Sets molecule's name.
		
		Args:
			n (str): New name.
		
		Returns:
			str: Name.
			
		"""
		
		self.name=n
		return self.name
	
	def getName(self):
	
		"""Returns molecule's name.
		
		Returns:
			str: Name.
			
		"""
		
		return self.name
	
	
	
	def sumUpResults(self,sameSettings=False):
		
		"""Sums up results from all fits in ``selFits`` list.
		
		Keyword Args:
			sameSettings (bool): Fits must overlap in parameters defined in ``crucialParameters``.
		
		Returns:
			bool: True if success, False else.	
		"""
		
		if sameSettings:
			
			lastFit=self.selFits[0]
			
			#Check that all fits have roughly the same parameters
			for fit in self.selFits:
				
				#Check if fit is fitted
				if not fit.isFitted():
					printError("Cannot average fits, since fits " + fit.name + "has not been fitted yet.")
					return False
				
				same,different,notInBoth=pyfrp_misc_module.compareObjAttr(lastFit,fit)
				
				for item in self.crucialParameters:   
					if item in list(different.keys()):
						printError("Cannot average fits, since fits " + lastFit.name + " and " + fit.name + "do not have the same value for " + item +". However, this parameter is marked as crucial.")
						return False
					elif item in notInBoth:
						printError("Cannot average fits, since one of the fits " + lastFit.name + " and " + fit.name + " lacks the parameter " + item +".")
						return False
		
		#Compute Statistics
		self.DOptMu,self.DOptMuStd,self.DOptMuSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"DOptMu"))
		self.DOptPx,self.DOptPxStd,self.DOptPxSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"DOptPx"))
		self.prodOpt,self.prodOptStd,self.prodOptSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"prodOpt"))
		self.degrOpt,self.degrOptStd,self.degrOptSterr=pyfrp_stats_module.parameterStats(pyfrp_misc_module.objAttrToList(self.selFits,"degrOpt"))
		
		self.Rsq=np.mean(pyfrp_misc_module.objAttrToList(self.selFits,"Rsq"))
		self.MeanRsq=np.mean(pyfrp_misc_module.objAttrToList(self.selFits,"MeanRsq"))
			
		return	True
		
	def printResults(self):
		
		"""Prints results summarized in :py:func:`pyfrp.subclasses.pyfrp_molecule.sumUpResults`."""
		
		printObjAttr('name',self)
		printObjAttr('DOptMu',self)
		printObjAttr('DOptMuStd',self)
		printObjAttr('DOptMuSterr',self)
		printObjAttr('DOptPx',self)
		printObjAttr('DOptPxStd',self)
		printObjAttr('DOptPxSterr',self)
		printObjAttr('prodOpt',self)
		printObjAttr('prodOpt',self)
		printObjAttr('prodOptStd',self)
		printObjAttr('prodOptSterr',self)
		printObjAttr('degrOpt',self)
		printObjAttr('degrOptStd',self)
		printObjAttr('degrOptSterr',self)
		printObjAttr('Rsq',self)
		printObjAttr('MeanRsq',self)
		
	def getFitParm(self,parm):
		
		"""Collects all values of parameter from all selected fits of molecule.
		
		Args:
			parm (str): Name of parameter.
			
		Returns:
			list: List of values.
			
		"""	
		
		ps=[]
		for fit in self.selFits:
			
			try:
				ps.append(getattr(fit,parm))
			except AttributeError:
				printError("Fit " +fit.name + "does not have an attribute with name = ", parm )
			
		return ps	
		
	def getDOptMus(self):	
		
		"""Collects all values of DOptMu from all selected fits of molecule.
	
		Returns:
			list: List of diffusion rates.
			
		"""	
		
		return self.getFitParm("DOptMu")
	
	