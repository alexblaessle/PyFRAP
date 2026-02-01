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

#PyQT Dialogs for editing embryo datasets
#(1) embryoDialog

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

class embryoDialog(pyfrp_gui_basics.basicCanvasDialog):
	
	def __init__(self,embryo,parent):
		
		super(embryoDialog,self).__init__(parent)
		
		#Passing embryo
		self.embryo=embryo
		self.parent=parent
		self.nCharDisplayed=50
		
		#Labels
		self.lblName = QtWidgets.QLabel("Name:", self)
		
		self.lblDataFT = QtWidgets.QLabel("Data Filetype:", self)
		self.lblDataEnc = QtWidgets.QLabel("Data Encoding:", self)
		
		self.lblFnDataFolder = QtWidgets.QLabel("Data Folder:", self)
	
		
		self.lblDataResPx = QtWidgets.QLabel("Resolution (px):", self)
		self.lblDataResMu = QtWidgets.QLabel("Resulution (um):", self)
		
		self.lblSliceDepth = QtWidgets.QLabel("Imaging Depth (um):", self)
		
		self.lblFrameInterval = QtWidgets.QLabel("Frame Interval (s):", self)
		self.lblnFrames = QtWidgets.QLabel("number of frames:", self)
		self.lbltStart = QtWidgets.QLabel("tStart (s):", self)
		self.lbltEnd = QtWidgets.QLabel("tEnd (s):", self)
		
		self.lblFnDataFolderValue = QtWidgets.QLabel("", self)
		
		self.updateDataFolderLbl()
		
		#LineEdits
		self.qleName = QtWidgets.QLineEdit(self.embryo.name)
		self.qleDataResPx = QtWidgets.QLineEdit(str(self.embryo.dataResPx))
		self.qleDataResMu = QtWidgets.QLineEdit(str(self.embryo.dataResMu))
		
		self.qleSliceDepth = QtWidgets.QLineEdit(str(self.embryo.sliceDepthMu))
		
		self.qleFrameInterval = QtWidgets.QLineEdit(str(self.embryo.frameInterval))
		self.qlenFrames = QtWidgets.QLineEdit(str(self.embryo.nFrames))
		self.qletStart = QtWidgets.QLineEdit(str(self.embryo.tStart))
		self.qletEnd = QtWidgets.QLineEdit(str(self.embryo.tEnd))
		
		self.doubleValid=QtGui.QDoubleValidator()
		
		self.qleDataResPx.setValidator(self.doubleValid)
		self.qleDataResMu.setValidator(self.doubleValid)
		self.qleSliceDepth.setValidator(self.doubleValid)
		self.qleFrameInterval.setValidator(self.doubleValid)
		self.qlenFrames.setValidator(self.doubleValid)
		self.qletStart.setValidator(self.doubleValid)
		self.qletEnd.setValidator(self.doubleValid)
		
		self.qleName.editingFinished.connect(self.setName)
		self.qleDataResMu.editingFinished.connect(self.setDataResMu)
		self.qleDataResPx.editingFinished.connect(self.setDataResPx)
		self.qleSliceDepth.editingFinished.connect(self.setSliceDepth)
		self.qleFrameInterval.editingFinished.connect(self.setFrameInterval)
		self.qletStart.editingFinished.connect(self.settStart)
		
		self.qletEnd.setReadOnly(True)
		self.qlenFrames.setReadOnly(True)
		
		#Combo
		self.comboDataFt = QtWidgets.QComboBox(self)
		self.comboDataFt.addItem("tif")
		self.comboDataFt.setCurrentIndex(self.comboDataFt.findText(self.embryo.getDataFT()))
		
		self.comboDataEnc = QtWidgets.QComboBox(self)
		self.comboDataEnc.addItem("uint8")
		self.comboDataEnc.addItem("uint16")
		self.comboDataEnc.setCurrentIndex(self.comboDataEnc.findText(self.embryo.getDataEnc())) 
		
		self.comboDataFt.activated[str].connect(self.setDataFt)   
		self.comboDataEnc.activated[str].connect(self.setDataEnc)   
	
		#Buttons
		self.btnFnDatafolder=QtWidgets.QPushButton('Change')
		
		self.btnFnDatafolder.clicked.connect(self.setFnDatafolder)
			
		#Layout
		self.dataFolderGrid = QtWidgets.QGridLayout()
		self.dataFolderGrid.addWidget(self.lblFnDataFolderValue,1,1)
		self.dataFolderGrid.addWidget(self.btnFnDatafolder,1,2)
		
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblName,nRows+1,1)
		self.grid.addWidget(self.lblDataFT,nRows+3,1)
		self.grid.addWidget(self.lblDataEnc,nRows+4,1)
		self.grid.addWidget(self.lblFnDataFolder,nRows+5,1)
		
		self.grid.addWidget(self.lblDataResPx,nRows+8,1)
		self.grid.addWidget(self.lblDataResMu,nRows+9,1)
		
		self.grid.addWidget(self.lblSliceDepth,nRows+10,1)
		
		self.grid.addWidget(self.lblFrameInterval,nRows+12,1)
		self.grid.addWidget(self.lblnFrames,nRows+13,1)
		self.grid.addWidget(self.lbltStart,nRows+14,1)
		self.grid.addWidget(self.lbltEnd,nRows+15,1)
		
		self.grid.addWidget(self.qleName,nRows+1,2)
		self.grid.addWidget(self.comboDataFt,nRows+3,2)
		self.grid.addWidget(self.comboDataEnc,nRows+4,2)
		self.grid.addLayout(self.dataFolderGrid,nRows+5,2)
		
		self.grid.addWidget(self.qleDataResPx,nRows+8,2)
		self.grid.addWidget(self.qleDataResMu,nRows+9,2)
		
		self.grid.addWidget(self.qleSliceDepth,nRows+10,2)
		
		self.grid.addWidget(self.qleFrameInterval,nRows+12,2)
		self.grid.addWidget(self.qlenFrames,nRows+13,2)
		self.grid.addWidget(self.qletStart,nRows+14,2)
		self.grid.addWidget(self.qletEnd,nRows+15,2)
		
		self.showFirstDataImg()
		
		self.show()
		
	def setFnDatafolder(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Data Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.embryo.fnDatafolder=folder
		
		self.parent.lastopen=folder
		
		self.embryo.fnDatafolder=pyfrp_misc_module.slashToFn(self.embryo.fnDatafolder)
	
		self.updateEmbryo()
		self.updateTimeQles()
		
		self.updateDataFolderLbl()
		
		self.showFirstDataImg()
		
	def setFnPreImage(self):
		
		fn = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', self.parent.lastopen,))
		if fn=='':
			return
		
		self.parent.lastopen,temp=os.path.split(fn)
		
		self.embryo.fnPreimage=fn
		self.updatePreImageLbl()
		
	def updateDataFolderLbl(self):
		self.lblFnDataFolderValue.setText("..."+self.embryo.fnDatafolder[-self.nCharDisplayed:])
	
	def updatePreImageLbl(self):
		self.lblFnPreimageValue.setText("..."+self.embryo.fnPreimage[-self.nCharDisplayed:])
	
	def updateTimeQles(self):
		self.qlenFrames.setText(str(self.embryo.getNFrames()))
		self.qletStart.setText(str(self.embryo.getTStart()))
		self.qletEnd.setText(str(self.embryo.getTEnd()))
			
	def setDataFt(self,text):
		
		self.embryo.setDataFT(str(text))
		
		self.updateEmbryo()
		self.updateTimeQles()
		
		self.showFirstDataImg()
		
	def setDataEnc(self,text):
	
		self.embryo.setDataEnc("uint16")
		self.showFirstDataImg()
		
	def setName(self):
		self.embryo.setName(str(self.qleName.text()))
		return self.embryo.getName()

	def setFrameInterval(self):
		self.embryo.setFrameInterval(float(str(self.qleFrameInterval.text())))
		self.updateEmbryo()
		self.updateTimeQles()
		return self.embryo.getTvecData()
		
	def settStart(self):
		self.embryo.setTStart(float(str(self.qletStart.text())))
		self.updateEmbryo()
		self.updateTimeQles()
		return self.embryo.getTStart()
		
	def setDataResPx(self):
		self.embryo.setDataResPx(float(str(self.qleDataResPx.text())))
		
		self.ax.set_xlim([0, self.embryo.dataResPx])
		self.ax.set_ylim([0, self.embryo.dataResPx])
		self.canvas.draw()
		return self.embryo.getDataResMu()
	
	def setSliceDepth(self):
		self.embryo.setSliceDepthMu(float(str(self.qleSliceDepth.text())))
	
	def setDataResMu(self):
		self.embryo.setDataResMu(float(str(self.qleDataResMu.text())))
		return self.embryo.getDataResPx()
		
	def updateEmbryo(self):
		self.embryo.updateFileList()
		self.embryo.updateNFrames()
		self.embryo.updateTimeDimensions()
		
	def showFirstDataImg(self):
		
		self.embryo.updateFileList()
		if len(self.embryo.fileList)>0:
		
			fnImg=self.embryo.fnDatafolder+self.embryo.fileList[0]
			img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
			
			self.showImg(img)
	

class lsmWizard(pyfrp_gui_basics.basicSettingsDialog):
	
	"""Dialog for lsm2embryo wizard."""
	
	def __init__(self,parent):
		
		super(lsmWizard,self).__init__(parent)
		
		#Passing parent GUI
		self.parent=parent
		
		#Default options
		self.nCharDisplayed=50
		self.fnDatafolder=""
		self.fnDestfolder=None
		self.name=""
		self.createEmbryo=True
		self.dataFT="lsm"
		self.embryo=None
		self.nChannel=1
		
		self.recoverIdent=['recover', 'post']
		self.bleachIdent=['bleach']
		self.preIdent=['pre']
		self.colorPrefix='_c00'
		
		self.cleanUp=True
		
		#Labels
		self.lblName = QtWidgets.QLabel("Name:", self)
		
		self.lblDataFT = QtWidgets.QLabel("Data Filetype:", self)
		
		self.lblFnDataFolder = QtWidgets.QLabel("Data Folder:", self)
	
		self.lblFnDestFolder = QtWidgets.QLabel("Destination Folder:", self)

		self.lblNChannel = QtWidgets.QLabel("nChannel:", self)
		self.lblCreate = QtWidgets.QLabel("Create Embryo:", self)
		
		self.lblFnDataFolderValue = QtWidgets.QLabel("", self)
		self.lblFnDestFolderValue = QtWidgets.QLabel("", self)
		
		self.lblRecoverIdent= QtWidgets.QLabel("Recover Identifiers", self)
		self.lblPreIdent= QtWidgets.QLabel("Pre Identifiers", self)
		self.lblBleachIdent= QtWidgets.QLabel("Bleach Identifiers", self)
		
		self.lblColorPrefix = QtWidgets.QLabel("Color Prefix", self)
		self.lblCleanUp = QtWidgets.QLabel("Clean Up Afterwards?", self)	
			
		self.updateDataFolderLbl()
		
		#LineEdits
		self.qleName = QtWidgets.QLineEdit("")
		self.qleRecoverIdent = QtWidgets.QLineEdit(str(self.recoverIdent))
		self.qlePreIdent = QtWidgets.QLineEdit(str(self.preIdent))
		self.qleBleachIdent = QtWidgets.QLineEdit(str(self.bleachIdent))
		
		self.qleColorPrefix = QtWidgets.QLineEdit(self.colorPrefix)
		
		self.qleName.editingFinished.connect(self.setName)
		self.qleRecoverIdent.editingFinished.connect(self.setRecoverIdent)
		self.qlePreIdent.editingFinished.connect(self.setPreIdent)
		self.qleBleachIdent.editingFinished.connect(self.setBleachIdent)
		
		self.qleColorPrefix.editingFinished.connect(self.setColorPrefix)
		
		
		#Combo
		self.comboDataFt = QtWidgets.QComboBox(self)
		self.comboDataFt.addItem("lsm")
		self.comboDataFt.addItem("czi")
		
		self.comboNChannel = QtWidgets.QComboBox(self)
		self.comboNChannel.addItem("1")
		self.comboNChannel.addItem("2")
		self.comboNChannel.addItem("3")
		self.comboNChannel.addItem("4")
		self.comboNChannel.addItem("5")
		self.comboNChannel.addItem("6")
		
		self.comboDataFt.activated[str].connect(self.setDataFt)   
		self.comboNChannel.activated[str].connect(self.setNChannel)   
		
		#Checkboxes
		self.cbCreate = QtWidgets.QCheckBox('', self)
		self.cbCreate.stateChanged.connect(self.checkCreate)
		self.cbCreate.setCheckState(2)
		
		self.cbCleanUp = QtWidgets.QCheckBox('', self)
		self.cbCleanUp.stateChanged.connect(self.checkCleanUp)
		self.cbCleanUp.setCheckState(2)
		
		#Buttons
		self.btnFnDatafolder=QtWidgets.QPushButton('Change')
		self.btnFnDatafolder.clicked.connect(self.setFnDatafolder)
		
		self.btnFnDestfolder=QtWidgets.QPushButton('Change')
		self.btnFnDestfolder.clicked.connect(self.setFnDestfolder)
		
		self.btnBuild=QtWidgets.QPushButton('Build Embryo')
		self.btnBuild.clicked.connect(self.buildEmbryo)
		
		#Layout
		self.dataFolderGrid = QtWidgets.QGridLayout()
		self.dataFolderGrid.addWidget(self.lblFnDataFolderValue,1,1)
		self.dataFolderGrid.addWidget(self.btnFnDatafolder,1,2)
		
		self.destFolderGrid = QtWidgets.QGridLayout()
		self.destFolderGrid.addWidget(self.lblFnDestFolderValue,1,1)
		self.destFolderGrid.addWidget(self.btnFnDestfolder,1,2)
		
		nRows=self.grid.rowCount()
		self.grid.addWidget(self.lblName,nRows+1,1)
		self.grid.addWidget(self.lblFnDataFolder,nRows+3,1)
		self.grid.addWidget(self.lblFnDestFolder,nRows+4,1)
		
		self.grid.addWidget(self.lblDataFT,nRows+6,1)
		self.grid.addWidget(self.lblNChannel,nRows+7,1)
		
		self.grid.addWidget(self.lblCreate,nRows+8,1)
		
		self.grid.addWidget(self.lblRecoverIdent,nRows+10,1)
		self.grid.addWidget(self.lblBleachIdent,nRows+11,1)
		self.grid.addWidget(self.lblPreIdent,nRows+12,1)
		
		self.grid.addWidget(self.lblColorPrefix,nRows+14,1)
		self.grid.addWidget(self.lblCleanUp,nRows+15,1)
		
		
		self.grid.addWidget(self.qleName,nRows+1,2)
		self.grid.addLayout(self.dataFolderGrid,nRows+3,2)
		self.grid.addLayout(self.destFolderGrid,nRows+4,2)
		
		self.grid.addWidget(self.comboDataFt,nRows+6,2)
		self.grid.addWidget(self.comboNChannel,nRows+7,2)
		
		self.grid.addWidget(self.cbCreate,nRows+8,2)
		
		self.grid.addWidget(self.qleRecoverIdent,nRows+10,2)
		self.grid.addWidget(self.qleBleachIdent,nRows+11,2)
		self.grid.addWidget(self.qlePreIdent,nRows+12,2)
	
		self.grid.addWidget(self.qleColorPrefix,nRows+14,2)
		
		self.grid.addWidget(self.cbCleanUp,nRows+15,2)
	
		
		self.grid.addWidget(self.btnBuild,nRows+16,2)
	
	
		self.setWindowTitle('From Microscope to Embryo Wizard')
		
		self.show()
	
	def fixIdent(self,ident):
		if ident[0]!="[":
			ident="["+ident
		if ident[-1]!="]":
			ident=ident+"]"
		
		return ident
	
	def setRecoverIdent(self):
		ident=str(self.qleRecoverIdent.text())
		ident=self.fixIdent(ident)
		self.recoverIdent=pyfrp_misc_module.str2list(ident,dtype="str")[0]
		
	def setBleachIdent(self):
		ident=str(self.qleBleachIdent.text())
		ident=self.fixIdent(ident)
		self.bleachIdent=pyfrp_misc_module.str2list(ident,dtype="str")[0]
		
	def setPreIdent(self):
		ident=str(self.qlePreIdent.text())
		ident=self.fixIdent(ident)
		self.preIdent=pyfrp_misc_module.str2list(ident,dtype="str")[0]
				
	def setName(self):
		self.name=str(self.qleName.text())
		
	def setColorPrefix(self):
		self.colorPrefix=str(self.qleColorPrefix.text())
		
	def buildEmbryo(self):	
		ret=pyfrp_misc_module.buildEmbryoWizard(self.fnDatafolder, self.dataFT, self.name, nChannel=self.nChannel-1, fnDest=self.fnDestfolder, createEmbryo=self.createEmbryo, recoverIdent=self.recoverIdent,
					  bleachIdent=self.bleachIdent,preIdent=self.preIdent,colorPrefix=self.colorPrefix,cleanUp=self.cleanUp)
		
		if ret==-1:
			printError("Embryo creation was not successful.")
		else:
			self.embryo=ret
			
				
	def setFnDatafolder(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Data Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.fnDatafolder=folder
	
		self.parent.lastopen=folder
		
		self.fnDatafolder=pyfrp_misc_module.slashToFn(self.fnDatafolder)
		
		if self.fnDestfolder==None:
			self.fnDestfolder=self.fnDatafolder
			self.updateDestFolderLbl()
			
		if self.name=="":
			files=pyfrp_misc_module.getSortedFileList(self.fnDatafolder,self.dataFT)
			self.name=files[0].replace("."+self.dataFT,"")
			self.qleName.setText(self.name)
			
		self.updateDataFolderLbl()
			
	def updateDataFolderLbl(self):
		self.lblFnDataFolderValue.setText("..."+self.fnDatafolder[-self.nCharDisplayed:])

	def setFnDestfolder(self):
		
		folder = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Dest Directory",  self.parent.lastopen,))
		if folder=='':
			return
		
		self.fnDestfolder=folder
	
		self.parent.lastopen=folder
		
		self.fnDestfolder=pyfrp_misc_module.slashToFn(self.fnDestfolder)
		
		self.updateDestFolderLbl()
			
	def updateDestFolderLbl(self):
		self.lblFnDestFolderValue.setText("..."+self.fnDestfolder[-self.nCharDisplayed:])

	
	def checkCreate(self,val):
		self.createEmbryo=(bool(2*val))
		
	def checkCleanUp(self,val):
		self.cleanUp=(bool(2*val))	
		
	def setDataFt(self,text):
		
		self.dataFT=(str(text))
		
	def setNChannel(self,text):
		self.nChannel=int(str(text))
	
	def getEmbryo(self):
		return self.embryo
	
class wizardSelector(QtWidgets.QDialog):
	
	"""Dialog to select if either use lsmWizard or not."""
	
	def __init__(self,parent):
		
		super(wizardSelector,self).__init__(parent)
		
		#Passing parent GUI
		self.parent=parent
		
		self.useLSM=None
		
		#Buttons
		self.btnUseLSM=QtWidgets.QPushButton('Create Embryo from Microscope Data')
		self.btnUseLSM.clicked.connect(self.setUseLSM)
		
		self.btnUseSorted=QtWidgets.QPushButton('Create Embryo from already prepared Data')
		self.btnUseSorted.clicked.connect(self.setUseSorted)
			
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.btnUseLSM)
		self.vbox.addWidget(self.btnUseSorted)
		
		self.setLayout(self.vbox)    
			
		self.setWindowTitle('Select Wizard Steps')    
		self.show()
		
	def setUseLSM(self):
		self.useLSM=True
		self.done(1)
		
	def setUseSorted(self):
		self.useLSM=False
		self.done(1)
		
	def getUseLSM(self):
		return self.useLSM
	
	
		
		
		
	
