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

#PyQT Dialogs for analysis class
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
#Dialog for editing embryo datasets
#===================================================================================================================================

class analysisDialog(pyfrp_gui_basics.basicSettingsDialog):
	
	def __init__(self,analysis,parent):
		
		super(analysisDialog,self).__init__(parent)
		
		self.analysis = analysis
		self.parent=parent
		self.nCharDisplayed=50
		
		#Labels
		self.lblFnPreimage = QtWidgets.QLabel("Pre-Image:", self)
		self.lblFnFlatten = QtWidgets.QLabel("Flattening Folder:", self)
		self.lblFnBkgd = QtWidgets.QLabel("Background Folder:", self)
		
		self.lblMedianRadius = QtWidgets.QLabel("Median Radius:", self)
		self.lblGaussianSigma = QtWidgets.QLabel("Gaussian Sigma:", self)
		
		self.lblFnPreimageValue = QtWidgets.QLabel("", self)
		self.lblFnFlattenValue = QtWidgets.QLabel("", self)
		self.lblFnBkgdValue = QtWidgets.QLabel("", self)
		
		self.lblNPre = QtWidgets.QLabel("Number of Images used: ", self)
		self.lblNFlatten = QtWidgets.QLabel("Number of Images used: ", self)
		self.lblNBkgd = QtWidgets.QLabel("Number of Images used: ", self)
		
		self.updateFlattenLbl()
		self.updatePreImageLbl()
		self.updateBkgdLbl()
		
		#LineEdits
		self.qleMedianRadius = QtWidgets.QLineEdit(str(self.analysis.medianRadius))
		self.qleGaussianSigma = QtWidgets.QLineEdit(str(self.analysis.gaussianSigma))
		
		self.qleNPre = QtWidgets.QLineEdit(str(self.analysis.nPre))
		self.qleNFlatten = QtWidgets.QLineEdit(str(self.analysis.nFlatten))
		self.qleNBkgd = QtWidgets.QLineEdit(str(self.analysis.nBkgd))
		
		self.doubleValid=QtGui.QDoubleValidator()
		self.intValid=QtGui.QIntValidator()
		
		self.qleMedianRadius.setValidator(self.doubleValid)
		self.qleGaussianSigma.setValidator(self.doubleValid)
		
		self.qleNPre.setValidator(self.intValid)
		self.qleNFlatten.setValidator(self.intValid)
		self.qleNBkgd.setValidator(self.intValid)
		
		self.qleMedianRadius.editingFinished.connect(self.setMedianRadius)
		self.qleGaussianSigma.editingFinished.connect(self.setGaussianSigma)
		
		self.qleNPre.editingFinished.connect(self.setNPre)
		self.qleNBkgd.editingFinished.connect(self.setNBkgd)
		self.qleNFlatten.editingFinished.connect(self.setNFlatten)
	
		#Checkboxes
		self.cbNorm = QtWidgets.QCheckBox('Norm by pre-image?', self)
		self.cbMedian = QtWidgets.QCheckBox('Apply median filter?', self)
		self.cbGaussian = QtWidgets.QCheckBox('Apply gaussian filter?', self)
		self.cbFlatten = QtWidgets.QCheckBox('Apply flattening mask?', self)
		self.cbBkgd = QtWidgets.QCheckBox('Substract Background mask?', self)
		
		self.cbQuad = QtWidgets.QCheckBox('Flip to quadrant?', self)
		self.cbFlip = QtWidgets.QCheckBox('Flip before process?', self)
		
		self.updateCBs()
		
		self.cbMedian.stateChanged.connect(self.checkMedian)
		self.cbGaussian.stateChanged.connect(self.checkGaussian)
		self.cbFlatten.stateChanged.connect(self.checkFlatten)
		self.cbNorm.stateChanged.connect(self.checkNorm)
		self.cbQuad.stateChanged.connect(self.checkQuad)
		self.cbFlip.stateChanged.connect(self.checkFlip)
		self.cbBkgd.stateChanged.connect(self.checkBkgd)
		
		
		#Buttons
		self.btnFnPreImage=QtWidgets.QPushButton('Change')
		self.btnFnFlatten=QtWidgets.QPushButton('Change')
		self.btnFnBkgd=QtWidgets.QPushButton('Change')
		
		self.btnFnPreImage.clicked.connect(self.setFnPreImage)
		self.btnFnFlatten.clicked.connect(self.setFnFlatten)
		self.btnFnBkgd.clicked.connect(self.setFnBkgd)
		
		#Layout
		self.preImageGrid = QtWidgets.QGridLayout()
		self.preImageGrid.addWidget(self.lblFnPreimageValue,1,1)
		self.preImageGrid.addWidget(self.btnFnPreImage,1,2)
		self.preImageGrid.setColumnMinimumWidth(1,150)
		
		self.flattenGrid = QtWidgets.QGridLayout()
		self.flattenGrid.addWidget(self.lblFnFlattenValue,1,1)
		self.flattenGrid.addWidget(self.btnFnFlatten,1,2)
		self.flattenGrid.setColumnMinimumWidth(1,150)
		
		self.bkgdGrid = QtWidgets.QGridLayout()
		self.bkgdGrid.addWidget(self.lblFnBkgdValue,1,1)
		self.bkgdGrid.addWidget(self.btnFnBkgd,1,2)
		self.bkgdGrid.setColumnMinimumWidth(1,150)
		
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.cbNorm,nRows+2,1)
		self.grid.addWidget(self.cbFlatten,nRows+3,1)
		self.grid.addWidget(self.cbBkgd,nRows+4,1)
		self.grid.addWidget(self.cbMedian,nRows+5,1)
		self.grid.addWidget(self.cbGaussian,nRows+6,1)
		self.grid.addWidget(self.cbQuad,nRows+7,1)
		self.grid.addWidget(self.cbFlip,nRows+8,1)
		
		self.grid.addWidget(self.lblFnPreimage,nRows+2,3)
		self.grid.addWidget(self.lblFnFlatten,nRows+3,3)
		self.grid.addWidget(self.lblFnBkgd,nRows+4,3)
		self.grid.addWidget(self.lblMedianRadius,nRows+5,3)
		self.grid.addWidget(self.lblGaussianSigma,nRows+6,3)
		
		self.grid.addLayout(self.preImageGrid,nRows+2,4)
		self.grid.addLayout(self.flattenGrid,nRows+3,4)
		self.grid.addLayout(self.bkgdGrid,nRows+4,4)
		self.grid.addWidget(self.qleMedianRadius,nRows+5,4)
		self.grid.addWidget(self.qleGaussianSigma,nRows+6,4)
		
		self.grid.addWidget(self.lblNPre,nRows+2,5)
		self.grid.addWidget(self.lblNFlatten,nRows+3,5)
		self.grid.addWidget(self.lblNBkgd,nRows+4,5)
		
		self.grid.addWidget(self.qleNPre,nRows+2,6)
		self.grid.addWidget(self.qleNFlatten,nRows+3,6)
		self.grid.addWidget(self.qleNBkgd,nRows+4,6)
		
		self.grid.setColumnMinimumWidth(4,200)
		
		self.setWindowTitle('Analysis Settings')    
		self.show()
		
	def updateCBs(self):
		self.cbFlip.setCheckState(2*int(self.inProcess('flipBeforeProcess')))
		self.cbFlatten.setCheckState(2*int(self.inProcess('flatten')))
		self.cbQuad.setCheckState(2*int(self.inProcess('quad')))
		self.cbMedian.setCheckState(2*int(self.inProcess('median')))
		self.cbGaussian.setCheckState(2*int(self.inProcess('gaussian')))
		self.cbNorm.setCheckState(2*int(self.inProcess('norm')))
		
	def inProcess(self,key):
		return key in list(self.analysis.process.keys())
	
	def checkMedian(self,val):
		self.analysis.setMedian(bool(2*val))
		
	def checkGaussian(self,val):
		self.analysis.setGaussian(bool(2*val))
	
	def checkFlatten(self,val):
		self.analysis.setFlatten(bool(2*val))
	
	def checkNorm(self,val):
		self.analysis.setNorm(bool(2*val))
	
	def checkBkgd(self,val):
		self.analysis.setBkgd(bool(2*val))

	def checkQuad(self,val):
		self.analysis.setQuad(bool(2*val))
	
	def checkFlip(self,val):
		self.analysis.setFlipBeforeProcess(bool(2*val))
		
	def setFnPreImage(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Preimage Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		folder=pyfrp_misc_module.slashToFn(folder)
		
		self.analysis.setFnPre(folder)
		
		self.parent.lastopen=folder
		self.updatePreImageLbl()
				
	def setFnFlatten(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Flatten Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.analysis.setFnFlatten(folder)
		
		self.parent.lastopen=folder
		
		self.updateFlattenLbl()
		
	def setFnBkgd(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Background Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.analysis.setFnBkgd(folder)
		
		self.parent.lastopen=folder
		
		self.updateBkgdLbl()	
	
	def setNPre(self):
		self.analysis.setNPre(int(str(self.qleNPre.text())))
		
	def setNFlatten(self):
		self.analysis.setNFlatten(int(str(self.qleNFlatten.text())))
	
	def setNBkgd(self):
		self.analysis.setNBkgd(int(str(self.qleNBkgd.text())))
	
	def updatePreImageLbl(self):
		self.lblFnPreimageValue.setText("..."+self.analysis.fnPreimage[-self.nCharDisplayed:])
		
	def updateFlattenLbl(self):
		self.lblFnFlattenValue.setText("..."+self.analysis.fnFlatten[-self.nCharDisplayed:])
	
	def updateBkgdLbl(self):
		self.lblFnBkgdValue.setText("..."+self.analysis.fnBkgd[-self.nCharDisplayed:])
	
	
	def setMedianRadius(self):
		self.analysis.setMedianRadius(float(str(self.qleMedianRadius.text())))
		
	def setGaussianSigma(self):
		self.analysis.setGaussianSigma(float(str(self.qleGaussianSigma.text())))
		
	
#===================================================================================================================================
#Dialogs for anaylze progress
#===================================================================================================================================

class analysisProgressDialog(pyfrp_gui_basics.progressDialog):
	
	def __init__(self,parent):
		super(analysisProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("Data analysis in progress...")
		
		#Window title
		self.setWindowTitle('Analysis progress')
		    
		self.show()	

class analysisThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, embryo=None, parent=None):
		#QtCore.QThread.__init__(self)
		super(analysisThread,self).__init__(parent)
		self.obj=embryo
		self.embryo=embryo
		
		#self.prog_signal.connect(self.print_prog)
			
	def runTask(self,debug=False):
		self.embryo.analysis.run(signal=self.progressSignal,embCount=None,debug=debug)
		
			
