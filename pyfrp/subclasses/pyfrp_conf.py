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

#! /usr/bin/python

#=====================================================================================================================================
#Description
#=====================================================================================================================================

#configuration is a simple class to save PyFRAP GUI configurations.

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#PyFRAP modules
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_misc_module

#OS
import os
import sys
import shutil

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Config object

class configuration:
	
	#Creates new molecule object
	def __init__(self):
		
		#Recently open files
		self.recentFiles=[]
		
		#Last view
		self.plotHidden=False
		self.termHidden=False
		self.propHidden=False
		self.backup2File=False
		self.backup2Memory=True
		self.pathFile=pyfrp_misc_module.getConfDir()+"paths"
		
		#Console History
		self.consoleHistory=[]
		
	def save(self,fn=None):
		
		if fn==None:
			fn=pyfrp_misc_module.getConfDir()+"lastConfiguration"+".conf"
		
		pyfrp_IO_module.saveToPickle(self,fn=fn)	
		
	def updateVersion(self):
		
		tempConf=configuration()
		pyfrp_misc_module.updateObj(tempConf,self)
		
		return self
	
	def setPlotHidden(self,h):
		self.plotHidden=h
		return self.plotHidden
	
	def setTermHidden(self,h):
		self.termHidden=h
		return self.termHidden
	
	def setPropHidden(self,h):
		self.propHidden=h
		return self.propHidden
	
	def getPlotHidden(self,h):
		return self.plotHidden
	
	def getTermHidden(self,h):
		return self.termHidden
	
	def getPropHidden(self,h):
		return self.propHidden
	
	def setBackup2File(self,h):
		self.backup2File=h
		return self.backup2File
	
	def getBackup2File(self,h):
		return self.backup2File
	
	def setBackup2Memory(self,h):
		self.backup2Memory=h
		return self.backup2Memory
	
	def getBackup2Memory(self,h):
		return self.backup2Memory
	
	def setRecentFiles(self,r):
		self.recentFiles=r
		return self.recentFiles
	
	def getPathFile(self):
		return self.pathFile
	
	def setPathFile(self,fn):
		self.pathFile=fn
		self.copyPathFileToDefaultLocation()
		return self.pathFile
	
	def copyPathFileToDefaultLocation(self):
		
		#Backup old paths file
		defaultFile=pyfrp_misc_module.getConfDir()+"paths"
		
		shutil.copy(defaultFile,defaultFile+".backup")
		try:
			shutil.copy(self.pathFile,defaultFile)
		except:
			printWarning('Was not able to copy '+self.pathFile+' to '+ defaultFile)
			pass
			
	def backupPathFile(self):	
		
		defaultFile=pyfrp_misc_module.getConfDir()+"paths"
		
		try:
			shutil.copy(defaultFile,self.pathFile)
		except:
			printWarning('Was not able to copy '+self.pathFile+' to '+ defaultFile)
			pass
	
	def getRecentFiles(self,r):
		return self.recentFiles
	
	def printConfiguration(self):
		pyfrp_term_module.printAllObjAttr(self)
	
	def addRecentFile(self,fn):
	
		if fn in self.recentFiles:
			ind=self.recentFiles.index(fn)
			self.recentFiles.pop(ind)
			
		self.recentFiles.insert(0,str(fn))
		
		return self.recentFiles