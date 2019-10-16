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

#PyQT Dialogs for fit class
#(1) 

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#QT
from PyQt4 import QtGui, QtCore

#PyFRAP GUI classes
from . import pyfrp_gui_basics

#PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_misc_module

#Numpy/Scipy
import numpy as np

#Misc 
import os


#===================================================================================================================================
#Dialog for editing pinning settings
#===================================================================================================================================

class defaultPinningDialog(pyfrp_gui_basics.basicSettingsDialog):
	
	def __init__(self,embryo,parent):
		
		super(defaultPinningDialog,self).__init__(parent)
		
		self.embryo = embryo
		
		self.useMin=True
		self.useMax=True
		self.debug=False
		
		#Labels
		self.lblUseMin = QtGui.QLabel("Use Min. Value for Background:", self)
		self.lblUseMax = QtGui.QLabel("Use Max. Value for Norming:", self)
		self.lblDebug = QtGui.QLabel("Show Debugging Output:", self)
		
		#Checkboxes
		self.cbUseMax = QtGui.QCheckBox('', self)
		self.cbUseMin = QtGui.QCheckBox('', self)
		self.cbDebug = QtGui.QCheckBox('', self)
		
		self.updateCBs()
		
		self.connect(self.cbUseMax, QtCore.SIGNAL('stateChanged(int)'), self.checkUseMax)
		self.connect(self.cbUseMin, QtCore.SIGNAL('stateChanged(int)'), self.checkUseMin)
		self.connect(self.cbDebug, QtCore.SIGNAL('stateChanged(int)'), self.checkDebug)
		
		#layout
		self.grid.addWidget(self.lblUseMin,1,1)
		self.grid.addWidget(self.lblUseMax,2,1)
		self.grid.addWidget(self.lblDebug,3,1)
		
		self.grid.addWidget(self.cbUseMin,1,2)
		self.grid.addWidget(self.cbUseMax,2,2)
		self.grid.addWidget(self.cbDebug,3,2)
		
		
		self.setWindowTitle("Basic Pinning Dialog")
		self.show()
	
	def updateCBs(self):	
		self.cbUseMax.setCheckState(2*int(self.useMax))
		self.cbUseMin.setCheckState(2*int(self.useMin))
		self.cbDebug.setCheckState(2*int(self.debug))
		
	def checkDebug(self,val):
		self.debug=bool(2*val)
	
	def checkUseMax(self,val):
		self.useMax=bool(2*val)
		
	def checkUseMin(self,val):
		self.useMin=bool(2*val)
	
	def donePressed(self):
		bkgdVal,normVal=self.embryo.computePinVals(useMin=self.useMin,useMax=self.useMax,bkgdVal=None,debug=self.debug)
		self.embryo.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdVal,normValSim=normVal,debug=self.debug)
		self.done(1)
		return
		
#===================================================================================================================================
#Dialog for editing ideal pinning settings
#===================================================================================================================================

class idealPinningDialog(defaultPinningDialog):
	
	def __init__(self,embryo,parent):
		
		super(idealPinningDialog,self).__init__(embryo,parent)
		
		self.switchThresh=0.95
		self.bkgdName='Bleached Square'
		self.normName='Slice'
		self.sepSim=True
		
		#Labels
		self.lblBkgd = QtGui.QLabel("ROI for Background Calculation:", self)
		self.lblNorm = QtGui.QLabel("ROI for Norm Calculation:", self)
		
		self.lblBkgdVal = QtGui.QLabel("", self)
		self.lblNormVal = QtGui.QLabel("", self)
		
		self.lblSepSim = QtGui.QLabel("Pin simulation seperately:", self)
		self.lblSwitchThresh = QtGui.QLabel("Switch Threshold:", self)
		
		self.updateBkgdROILbl()
		self.updateNormROILbl()

		#LineEdits
		self.doubleValid=QtGui.QDoubleValidator()
		self.qleSwitchThresh = QtGui.QLineEdit(str(self.switchThresh))
		self.qleSwitchThresh.setValidator(self.doubleValid)
		self.qleSwitchThresh.editingFinished.connect(self.setSwitchThresh)
		
		#Buttons
		self.btnSetBkgdROI=QtGui.QPushButton('Change')
		self.btnSetBkgdROI.connect(self.btnSetBkgdROI, QtCore.SIGNAL('clicked()'), self.setBkgdROI)
		
		self.btnSetNormROI=QtGui.QPushButton('Change')
		self.btnSetNormROI.connect(self.btnSetNormROI, QtCore.SIGNAL('clicked()'), self.setNormROI)
		
		
		#Checkboxes
		self.cbSepSim = QtGui.QCheckBox('', self)
		self.cbSepSim.setCheckState(2*int(self.sepSim))
	
		self.connect(self.cbSepSim, QtCore.SIGNAL('stateChanged(int)'), self.checkSepSim)
		
		#Layout
		nRows=self.grid.rowCount()
		self.grid.addWidget(self.lblBkgd,nRows+1,1)
		self.grid.addWidget(self.lblNorm,nRows+2,1)
		self.grid.addWidget(self.lblSwitchThresh,nRows+3,1)
		self.grid.addWidget(self.lblSepSim,nRows+4,1)
		
		self.grid.addWidget(self.lblBkgdVal,nRows+1,2)
		self.grid.addWidget(self.lblNormVal,nRows+2,2)
		self.grid.addWidget(self.qleSwitchThresh,nRows+3,2)
		self.grid.addWidget(self.cbSepSim,nRows+4,2)
		
		self.grid.addWidget(self.btnSetBkgdROI,nRows+1,3)
		self.grid.addWidget(self.btnSetNormROI,nRows+2,3)
		
		self.setWindowTitle("Ideal Pinning Dialog")
		self.show()
		
	def setSwitchThresh(self):
		self.switchThresh=float(str(self.qleSwitchThresh.text()))
	
	def checkSepSim(self,val):
		self.sepSim=bool(2*val)	
		
	def setBkgdROI(self):
		
		nameList=pyfrp_misc_module.objAttrToList(self.embryo.ROIs,'name')
		
		selectorDialog = pyfrp_gui_basics.basicSelectorDialog(nameList,self)
		if selectorDialog.exec_():
			selectedROIName = selectorDialog.getItem()
			if selectedROIName==None:
				return
		
		self.bkgdName=selectedROIName
		self.updateBkgdROILbl()	
			
	def updateBkgdROILbl(self):
		nameList=pyfrp_misc_module.objAttrToList(self.embryo.ROIs,'name')
		if self.bkgdName in nameList:
			pass
		else:
			self.bkgdName=""
			
		self.lblBkgdVal.setText(self.bkgdName)
	
	def setNormROI(self):
		
		nameList=pyfrp_misc_module.objAttrToList(self.embryo.ROIs,'name')
		
		selectorDialog = pyfrp_gui_basics.basicSelectorDialog(nameList,self)
		if selectorDialog.exec_():
			selectedROIName = selectorDialog.getItem()
			if selectedROIName==None:
				return
		
		self.normName=selectedROIName
		self.updateNormROILbl()	
			
	def updateNormROILbl(self):
		nameList=pyfrp_misc_module.objAttrToList(self.embryo.ROIs,'name')
		if self.normName in nameList:
			pass
		else:
			self.normName=""
			
		self.lblNormVal.setText(self.normName)
		
	def donePressed(self):

		bkgdVal, normVal, bkgdValSim, normValSim=self.embryo.computeIdealFRAPPinVals(bkgdName=self.bkgdName,normName=self.normName,debug=self.debug,useMin=self.useMin,useMax=self.useMax,sepSim=self.sepSim,switchThresh=self.switchThresh)
		self.embryo.pinAllROIs(bkgdVal=bkgdVal,normVal=normVal,bkgdValSim=bkgdValSim,normValSim=normValSim,debug=self.debug)
		self.done(1)
		return	
		
		