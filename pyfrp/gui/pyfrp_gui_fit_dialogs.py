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
from PyQt5 import QtGui, QtCore, QtWidgets

#PyFRAP GUI classes
from . import pyfrp_gui_basics

#PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_img_module
from pyfrp.modules import pyfrp_misc_module

#Numpy/Scipy
import numpy as np

#Misc 
import os


#===================================================================================================================================
#Dialog for editing simulation settings
#===================================================================================================================================

class fitSettingsDialog(pyfrp_gui_basics.basicSettingsDialog):
	
	def __init__(self,fit,parent):
		
		super(fitSettingsDialog,self).__init__(parent)
		
		self.fit = fit
		
		#Labels
		self.lblName = QtWidgets.QLabel("Name:", self)
		
		self.lblOptMeth = QtWidgets.QLabel("Optimization Method:", self)
		self.lblMaxFun = QtWidgets.QLabel("Max. Function Calls:", self)
		self.lblOptTol = QtWidgets.QLabel("Precision goal:", self)
		
		self.lblX0D = QtWidgets.QLabel("x0(D):", self)
		self.lblX0Prod = QtWidgets.QLabel("x0(production):", self)
		self.lblX0Degr = QtWidgets.QLabel("x0(degration):", self)
		
		self.lblBoundsD = QtWidgets.QLabel("<= D <=", self)
		self.lblBoundsProd = QtWidgets.QLabel("<= production <=", self)
		self.lblBoundsDegr = QtWidgets.QLabel("<= degration <=", self)
		
		self.lblKinetic = QtWidgets.QLabel("Kinetic Timescale:", self)
		self.lblSaveTrack = QtWidgets.QLabel("Save Track:", self)
		
		self.lblEquOn = QtWidgets.QLabel("Equalization:", self)
		self.lblFitPinned = QtWidgets.QLabel("Fit Pinned:", self)
		
		self.lblFitProd = QtWidgets.QLabel("Fit Production:", self)
		self.lblFitDegr = QtWidgets.QLabel("Fit Degradation:", self)
		
		self.lblFitCutOffT = QtWidgets.QLabel("Fit Cut-Off:", self)
		self.lblCutOffT = QtWidgets.QLabel("t(Cut-Off):", self)
		
		boldfont = QtGui.QFont()
		boldfont.setBold(True)
		
		self.lblHeadBounds = QtWidgets.QLabel("Bounds:", self)
		self.lblHeadBounds.setFont(boldfont)
		
		self.lblHeadX0 = QtWidgets.QLabel("Initial Guess:", self)
		self.lblHeadX0.setFont(boldfont)
		
		self.lblHeadGeneral = QtWidgets.QLabel("General:", self)
		self.lblHeadGeneral.setFont(boldfont)
		
		self.lblHeadROI = QtWidgets.QLabel("ROIs Fitted:", self)
		self.lblHeadROI.setFont(boldfont)
		
		self.lblHeadOptions = QtWidgets.QLabel("Fit Options:", self)
		self.lblHeadOptions.setFont(boldfont)
		
		#LineEdits
		self.qleName = QtWidgets.QLineEdit(self.fit.name)
		
		self.qleMaxFun = QtWidgets.QLineEdit(str(self.fit.maxfun))
		self.qleOptTol = QtWidgets.QLineEdit(str(self.fit.optTol))
		
		self.qleX0D = QtWidgets.QLineEdit(str(self.fit.x0[0]))
		self.qleX0Prod = QtWidgets.QLineEdit(str(self.fit.x0[1]))
		self.qleX0Degr = QtWidgets.QLineEdit(str(self.fit.x0[2]))
		
		self.qleLBD = QtWidgets.QLineEdit(str(self.fit.LBD))
		self.qleLBProd = QtWidgets.QLineEdit(str(self.fit.LBProd))
		self.qleLBDegr = QtWidgets.QLineEdit(str(self.fit.LBDegr))
		
		self.qleUBD = QtWidgets.QLineEdit(str(self.fit.UBD))
		self.qleUBProd = QtWidgets.QLineEdit(str(self.fit.UBProd))
		self.qleUBDegr = QtWidgets.QLineEdit(str(self.fit.UBDegr))
		
		self.qleKinetic = QtWidgets.QLineEdit(str(self.fit.kineticTimeScale))

		self.qleCutOffT = QtWidgets.QLineEdit(str(self.fit.cutOffT))
		
		self.doubleValid=QtGui.QDoubleValidator()
		self.intValid=QtGui.QIntValidator()
	
		self.qleX0D.setValidator(self.doubleValid)
		self.qleX0Prod.setValidator(self.doubleValid)
		self.qleX0Degr.setValidator(self.doubleValid)
		
		self.qleKinetic.setValidator(self.intValid)

		self.qleCutOffT.setValidator(self.doubleValid)
			
		self.qleName.editingFinished.connect(self.setName)
		self.qleMaxFun.editingFinished.connect(self.setMaxfun)
		
		self.qleOptTol.editingFinished.connect(self.setOptTol)
		self.qleKinetic.editingFinished.connect(self.setKineticTimeScale)
		
		self.qleLBProd.editingFinished.connect(self.setLBProd)
		self.qleLBDegr.editingFinished.connect(self.setLBDegr)
		self.qleLBD.editingFinished.connect(self.setLBD)
		
		self.qleUBProd.editingFinished.connect(self.setUBProd)
		self.qleUBDegr.editingFinished.connect(self.setUBDegr)
		self.qleUBD.editingFinished.connect(self.setUBD)
		
		self.qleX0Prod.editingFinished.connect(self.setX0Prod)
		self.qleX0Degr.editingFinished.connect(self.setX0Degr)
		self.qleX0D.editingFinished.connect(self.setX0D)
		
		self.qleCutOffT.editingFinished.connect(self.setCutOffT)
		
		#ComboBox
		self.comboMeth = QtWidgets.QComboBox(self)
		self.comboMeth.addItem("Constrained Nelder-Mead")
		self.comboMeth.addItem("TNC")
		self.comboMeth.addItem("Nelder-Mead")
		self.comboMeth.addItem("L-BFGS-B")
		self.comboMeth.addItem("SLSQP")
		self.comboMeth.addItem("brute")
		self.comboMeth.addItem("BFGS")
		self.comboMeth.addItem("CG")
		
		self.initComboMeth()
		
		self.comboMeth.activated[str].connect(self.setOptMeth)   
		
		#Checkboxes
		self.cbFitProd = QtWidgets.QCheckBox('', self)
		self.cbFitDegr = QtWidgets.QCheckBox('', self)
		self.cbFitPinned = QtWidgets.QCheckBox('', self)
		self.cbEquOn = QtWidgets.QCheckBox('', self)
		self.cbFitCutOffT = QtWidgets.QCheckBox('', self)
		self.cbSaveTrack = QtWidgets.QCheckBox('', self)
	
		self.updateCBs()
		
		self.cbFitProd.stateChanged.connect(self.checkFitProd)
		self.cbFitDegr.stateChanged.connect(self.checkFitDegr)
		self.cbFitPinned.stateChanged.connect(self.checkFitPinned)
		self.cbFitCutOffT.stateChanged.connect(self.checkFitCutOffT)
		self.cbEquOn.stateChanged.connect(self.checkEquOn)
		self.cbSaveTrack.stateChanged.connect(self.checkSaveTrack)
		
		#TreeWigget
		self.ROIList=QtWidgets.QTreeWidget()
		self.ROIList.setHeaderLabels(["Name"])
		self.ROIList.setColumnWidth(0,100)
		self.updateROIList()
		
		#Buttons
		self.btnAddROI=QtWidgets.QPushButton('Add')
		self.btnAddROI.clicked.connect(self.addROI)
		
		self.btnRemoveROI=QtWidgets.QPushButton('Remove')
		self.btnRemoveROI.clicked.connect(self.removeROI)
		
		#Layout	
		self.grid.addWidget(self.lblHeadGeneral,0,1,1,2,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.lblName,1,1)
		self.grid.addWidget(self.lblOptMeth,2,1)
		self.grid.addWidget(self.lblMaxFun,3,1)
		self.grid.addWidget(self.lblOptTol,4,1)
		self.grid.addWidget(self.lblKinetic,5,1)
		self.grid.addWidget(self.lblSaveTrack,6,1)
		
		self.grid.addWidget(self.qleName,1,2)
		self.grid.addWidget(self.comboMeth,2,2)
		self.grid.addWidget(self.qleMaxFun,3,2)
		self.grid.addWidget(self.qleOptTol,4,2)
		self.grid.addWidget(self.qleKinetic,5,2)
		self.grid.addWidget(self.cbSaveTrack,6,2)
		
		self.grid.addWidget(self.lblHeadX0,8,1,1,2,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.lblX0D,9,1)
		self.grid.addWidget(self.lblX0Prod,10,1)
		self.grid.addWidget(self.lblX0Degr,11,1)
		
		self.grid.addWidget(self.qleX0D,9,2)
		self.grid.addWidget(self.qleX0Prod,10,2)
		self.grid.addWidget(self.qleX0Degr,11,2)
		
		self.grid.addWidget(self.lblHeadBounds,0,3,1,3,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.qleLBD,1,3)
		self.grid.addWidget(self.qleLBProd,2,3)
		self.grid.addWidget(self.qleLBDegr,3,3)
		
		self.grid.addWidget(self.lblBoundsD,1,4,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.lblBoundsProd,2,4,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.lblBoundsDegr,3,4,QtCore.Qt.AlignHCenter)
		
		self.grid.addWidget(self.qleUBD,1,5)
		self.grid.addWidget(self.qleUBProd,2,5)
		self.grid.addWidget(self.qleUBDegr,3,5)
		
		self.grid.addWidget(self.lblHeadROI,5,3,1,3,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.ROIList,6,3,5,3)
		self.grid.addWidget(self.btnAddROI,11,3)
		self.grid.addWidget(self.btnRemoveROI,11,5)
		
		self.grid.addWidget(self.lblHeadOptions,0,6,1,2,QtCore.Qt.AlignHCenter)
		self.grid.addWidget(self.lblEquOn,1,6)
		self.grid.addWidget(self.lblFitPinned,2,6)
		self.grid.addWidget(self.lblFitProd,3,6)
		self.grid.addWidget(self.lblFitDegr,4,6)
		self.grid.addWidget(self.lblFitCutOffT,5,6)
		self.grid.addWidget(self.lblCutOffT,6,6)
		
		self.grid.addWidget(self.cbEquOn,1,7)
		self.grid.addWidget(self.cbFitPinned,2,7)
		self.grid.addWidget(self.cbFitProd,3,7)
		self.grid.addWidget(self.cbFitDegr,4,7)
		self.grid.addWidget(self.cbFitCutOffT,5,7)
		self.grid.addWidget(self.qleCutOffT,6,7)
		
		self.resize(900,700)
		
		self.setWindowTitle("Fit Settings")
		
		self.show()
	
	def validateBound(self,b):
		if b==None or b=="":
			return True
		try:
			c=float(b)
			return True
		except:
			printError("Invalid boundary entered.")
			return False
			
	def updateROIList(self):
		
		self.ROIList.clear()
		for r in self.fit.ROIsFitted:
			QtWidgets.QTreeWidgetItem(self.ROIList,[r.name])
		return self.ROIList
	
	def initComboMeth(self):
		idx=self.comboMeth.findText(self.fit.optMeth,QtCore.Qt.MatchExactly)
		self.comboMeth.setCurrentIndex(idx)
		
	def updateCBs(self):	
		self.cbEquOn.setCheckState(2*int(self.fit.equOn))
		self.cbFitPinned.setCheckState(2*int(self.fit.fitPinned))
		self.cbFitProd.setCheckState(2*int(self.fit.fitProd))
		self.cbFitDegr.setCheckState(2*int(self.fit.fitDegr))
		self.cbFitCutOffT.setCheckState(2*int(self.fit.fitCutOffT))
		self.cbSaveTrack.setCheckState(2*int(self.fit.saveTrack))
		
	def checkEquOn(self,val):
		self.fit.setEqu(bool(2*val))
	
	def checkFitPinned(self,val):
		self.fit.setFitPinned(bool(2*val))
	
	def checkFitProd(self,val):
		self.fit.setFitProd(bool(2*val))
	
	def checkFitDegr(self,val):
		self.fit.setFitDegr(bool(2*val))
	
	def checkFitCutOffT(self,val):
		self.fit.setFitCutOffT(bool(2*val))
	
	def checkSaveTrack(self,val):
		self.fit.setSaveTrack(bool(2*val))
	
	def setLBD(self):
		text=str(self.qleLBD.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setLBD(text)
		else:
			if self.validateBound(text):
				self.fit.setLBD(float(text))
	
	def setLBProd(self):
		text=str(self.qleLBProd.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setLBProd(text)
		else:
			if self.validateBound(text):
				self.fit.setLBProd(float(text))
	
	def setLBDegr(self):
		text=str(self.qleLBDegr.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setLBDegr(text)
		else:
			if self.validateBound(text):
				self.fit.setLBDegr(float(text))
		
	def setUBD(self):
		text=str(self.qleUBD.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setUBD(text)
		else:
			if self.validateBound(text):
				self.fit.setUBD(float(text))
	
	def setUBProd(self):
		text=str(self.qleUBProd.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setUBProd(text)
		else:
			if self.validateBound(text):
				self.fit.setUBProd(float(text))
	
	def setUBDegr(self):
		text=str(self.qleUBDegr.text())
		if text=="" or text=='None':
			text=None	
			self.fit.setUBDegr(text)
		else:
			if self.validateBound(text):
				self.fit.setUBDegr(float(text))
		
	def setX0D(self):
		text=str(self.qleX0D.text())
		if text=="":
			printWarning("You must give a value for x0.")
			return
		self.fit.setX0D(float(text))
		
	def setX0Prod(self):
		text=str(self.qleX0Prod.text())
		if text=="":
			printWarning("You must give a value for x0.")
			return
		self.fit.setX0Prod(float(text))
		
	def setX0Degr(self):
		text=str(self.qleX0Degr.text())
		if text=="":
			printWarning("You must give a value for x0.")
			return
		self.fit.setX0Degr(float(text))
		
	def setName(self):
		self.fit.setName(str(self.qleName.text()))
		
	def setOptMeth(self):
		self.fit.setOptMeth(str(self.comboMeth.currentText()))
	
	def setKineticTimeScale(self):
		self.fit.setKineticTimeScale(float(str(self.qleKinetic.text())))
	
	def setMaxfun(self):
		self.fit.setMaxfun(int(str(self.qleMaxFun.text())))
		
	def setOptTol(self):
		self.fit.setOptTol(float(str(self.qleOptTol.text())))
		
	def setCutOffT(self):
		self.fit.setCutOffT(float(str(self.qleCutOffT.text())))
	
	def addROI(self):
		
		nameList=pyfrp_misc_module.objAttrToList(self.fit.embryo.ROIs,'name')
		
		selectorDialog = pyfrp_gui_basics.basicSelectorDialog(nameList,self)
		if selectorDialog.exec_():
			selectedROIName = selectorDialog.getItem()
			if selectedROIName==None:
				return
		
		selectedROI=self.fit.embryo.ROIs[nameList.index(selectedROIName)]
		
		if selectedROI not in self.fit.ROIsFitted:
			self.fit.addROI(selectedROI)
		
		self.updateROIList()
		
	def removeROI(self):
		idx=self.ROIList.indexFromItem(self.ROIList.currentItem()).row()
		self.fit.removeROI(self.fit.ROIsFitted[idx])
		self.updateROIList()
		
#===================================================================================================================================
#Dialogs for fitting progress
#===================================================================================================================================

class fittingProgressDialog(pyfrp_gui_basics.waitDialog):
	
	def __init__(self,parent):
		super(fittingProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("Fitting in progress...")
		
		#Window title
		self.setWindowTitle('Fitting progress')
		    
		self.show()

class fittingThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, fit=None, parent=None):
		
		super(fittingThread,self).__init__(parent)
		self.obj=fit
		self.fit=fit
		
			
	def runTask(self,debug=False):
		self.fit.run(debug=debug)
			
			
