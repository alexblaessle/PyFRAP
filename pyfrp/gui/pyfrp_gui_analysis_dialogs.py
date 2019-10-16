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
from PyQt4 import QtGui, QtCore

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
		self.lblFnPreimage = QtGui.QLabel("Pre-Image:", self)
		self.lblFnFlatten = QtGui.QLabel("Flattening Folder:", self)
		self.lblFnBkgd = QtGui.QLabel("Background Folder:", self)
		
		self.lblMedianRadius = QtGui.QLabel("Median Radius:", self)
		self.lblGaussianSigma = QtGui.QLabel("Gaussian Sigma:", self)
		
		self.lblFnPreimageValue = QtGui.QLabel("", self)
		self.lblFnFlattenValue = QtGui.QLabel("", self)
		self.lblFnBkgdValue = QtGui.QLabel("", self)
		
		self.lblNPre = QtGui.QLabel("Number of Images used: ", self)
		self.lblNFlatten = QtGui.QLabel("Number of Images used: ", self)
		self.lblNBkgd = QtGui.QLabel("Number of Images used: ", self)
		
		self.updateFlattenLbl()
		self.updatePreImageLbl()
		self.updateBkgdLbl()
		
		#LineEdits
		self.qleMedianRadius = QtGui.QLineEdit(str(self.analysis.medianRadius))
		self.qleGaussianSigma = QtGui.QLineEdit(str(self.analysis.gaussianSigma))
		
		self.qleNPre = QtGui.QLineEdit(str(self.analysis.nPre))
		self.qleNFlatten = QtGui.QLineEdit(str(self.analysis.nFlatten))
		self.qleNBkgd = QtGui.QLineEdit(str(self.analysis.nBkgd))
		
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
		self.cbNorm = QtGui.QCheckBox('Norm by pre-image?', self)
		self.cbMedian = QtGui.QCheckBox('Apply median filter?', self)
		self.cbGaussian = QtGui.QCheckBox('Apply gaussian filter?', self)
		self.cbFlatten = QtGui.QCheckBox('Apply flattening mask?', self)
		self.cbBkgd = QtGui.QCheckBox('Substract Background mask?', self)
		
		self.cbQuad = QtGui.QCheckBox('Flip to quadrant?', self)
		self.cbFlip = QtGui.QCheckBox('Flip before process?', self)
		
		self.updateCBs()
		
		self.connect(self.cbMedian, QtCore.SIGNAL('stateChanged(int)'), self.checkMedian)
		self.connect(self.cbGaussian, QtCore.SIGNAL('stateChanged(int)'), self.checkGaussian)
		self.connect(self.cbFlatten, QtCore.SIGNAL('stateChanged(int)'), self.checkFlatten)
		self.connect(self.cbNorm, QtCore.SIGNAL('stateChanged(int)'), self.checkNorm)
		self.connect(self.cbQuad, QtCore.SIGNAL('stateChanged(int)'), self.checkQuad)
		self.connect(self.cbFlip, QtCore.SIGNAL('stateChanged(int)'), self.checkFlip)
		self.connect(self.cbBkgd, QtCore.SIGNAL('stateChanged(int)'), self.checkBkgd)
		
		
		#Buttons
		self.btnFnPreImage=QtGui.QPushButton('Change')
		self.btnFnFlatten=QtGui.QPushButton('Change')
		self.btnFnBkgd=QtGui.QPushButton('Change')
		
		self.btnFnPreImage.connect(self.btnFnPreImage, QtCore.SIGNAL('clicked()'), self.setFnPreImage)
		self.btnFnFlatten.connect(self.btnFnFlatten, QtCore.SIGNAL('clicked()'), self.setFnFlatten)
		self.btnFnBkgd.connect(self.btnFnBkgd, QtCore.SIGNAL('clicked()'), self.setFnBkgd)
		
		#Layout
		self.preImageGrid = QtGui.QGridLayout()
		self.preImageGrid.addWidget(self.lblFnPreimageValue,1,1)
		self.preImageGrid.addWidget(self.btnFnPreImage,1,2)
		self.preImageGrid.setColumnMinimumWidth(1,150)
		
		self.flattenGrid = QtGui.QGridLayout()
		self.flattenGrid.addWidget(self.lblFnFlattenValue,1,1)
		self.flattenGrid.addWidget(self.btnFnFlatten,1,2)
		self.flattenGrid.setColumnMinimumWidth(1,150)
		
		self.bkgdGrid = QtGui.QGridLayout()
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
		
		folder = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Preimage Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		folder=pyfrp_misc_module.slashToFn(folder)
		
		self.analysis.setFnPre(folder)
		
		self.parent.lastopen=folder
		self.updatePreImageLbl()
				
	def setFnFlatten(self):
		
		folder = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Flatten Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.analysis.setFnFlatten(folder)
		
		self.parent.lastopen=folder
		
		self.updateFlattenLbl()
		
	def setFnBkgd(self):
		
		folder = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Background Directory",  self.parent.lastopen,))
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
		
			