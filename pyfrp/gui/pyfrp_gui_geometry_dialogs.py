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

#PyQT Dialogs for editing geometries
#(1) geometryDialog

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

#Matplotlib
from mpl_toolkits.mplot3d import Axes3D

#===================================================================================================================================
#Dialog for editing embryo datasets
#===================================================================================================================================

class geometryDialog(pyfrp_gui_basics.basicCanvasDialog):
	
	def __init__(self,geometry,parent):
		
		super(geometryDialog,self).__init__(parent)
		
		#Passing geometry
		self.geometry=geometry
		self.parent=parent
		self.nCharDisplayed=50
		self.ann=0
		
		#Labels
		self.lblTyp = QtGui.QLabel("Typ:", self)
		self.lblCenter = QtGui.QLabel("Center (px):", self)
		self.lblFnGeo = QtGui.QLabel("Geo File:", self)
		
		self.lblTypVal = QtGui.QLabel(self.geometry.typ, self)
		self.lblFnGeoVal = QtGui.QLabel(self.geometry.fnGeo, self)
		
		#LineEdits
		self.doubleValid=QtGui.QDoubleValidator()
		
		self.qleCenterX = QtGui.QLineEdit(str(self.geometry.center[0]))
		self.qleCenterY = QtGui.QLineEdit(str(self.geometry.center[1]))
		
		self.qleCenterX.setValidator(self.doubleValid)
		self.qleCenterY.setValidator(self.doubleValid)
		
		self.centerGrid = QtGui.QGridLayout()
		
		self.centerGrid.addWidget(self.qleCenterX)
		self.centerGrid.addWidget(self.qleCenterY)
		
		self.qleCenterX.editingFinished.connect(self.setCenter)
		self.qleCenterY.editingFinished.connect(self.setCenter)
		
		#Buttons
		self.btnFnGeo=QtGui.QPushButton('Change')
		self.btnFnGeo.connect(self.btnFnGeo, QtCore.SIGNAL('clicked()'), self.setFnGeo)
		
		self.btnGetCenterFromROI=QtGui.QPushButton('Grab from ROI')
		self.btnGetCenterFromROI.connect(self.btnGetCenterFromROI, QtCore.SIGNAL('clicked()'), self.getCenterFromROI)
		
		self.btnUpdateAll=QtGui.QPushButton('Update optimal All ROI')
		self.btnUpdateAll.connect(self.btnUpdateAll, QtCore.SIGNAL('clicked()'), self.updateAllROI)
		
		#Checkboxes
		self.cbAnn = QtGui.QCheckBox('Annotate?', self)
		self.cbAnn.setCheckState(2*int(self.ann))
		
		self.connect(self.cbAnn, QtCore.SIGNAL('stateChanged(int)'), self.setAnnotate)
		
		#Layout
		self.fnGeoGrid = QtGui.QGridLayout()
		self.fnGeoGrid.addWidget(self.lblFnGeoVal,1,1)
		self.fnGeoGrid.addWidget(self.btnFnGeo,1,2)
		
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblTyp,nRows+1,1)
		self.grid.addWidget(self.lblCenter,nRows+2,1)
		self.grid.addWidget(self.lblFnGeo,nRows+3,1)
		
		self.grid.addWidget(self.lblTypVal,nRows+1,2)
		self.grid.addLayout(self.centerGrid,nRows+2,2)
		
		self.grid.addLayout(self.fnGeoGrid,nRows+3,2)
		
		self.grid.addWidget(self.cbAnn,nRows+4,2)
		self.grid.addWidget(self.btnUpdateAll,nRows+5,2)
		
		self.grid.addWidget(self.btnGetCenterFromROI,nRows+2,3)
		
		self.setAxes3D()
		
		self.geometry.updateGeoFile()
		
		self.drawGeometry()
		
		self.setWindowTitle('Edit Geometry')    
		
		
		self.show()
	
	def setCenter(self):
		center=[float(self.qleCenterX.text()),float(self.qleCenterY.text())]
		self.geometry.setCenter(center)
		self.drawGeometry()
		return center
		
	def setFnGeo(self):
		mdir=pyfrp_misc_module.getMeshfilesDir()
		fn = str(QtGui.QFileDialog.getOpenFileName(self, 'Open file', mdir,"*.geo",))
		if fn=='':
			return
				
		self.geometry.setFnGeo(fn)
		self.updateFnGeoLbl()
		
		self.drawGeometry()
		
		return fn
	
	def getCenterFromROI(self):
		
		"""First finds all ROIs in ``geometry.embryo.ROIs`` that
		have an attribute called ``center``, then lets the user
		select one of them and then sets centerof geometry to the same 
		center as the ROI.
		
		"""
		
		center=self.getAttrFromROI('center')
		
		self.geometry.setCenter(center)
		self.updateCenterQles()
	
	def getAttrFromROI(self,attr):
		
		"""First finds all ROIs in ``geometry.embryo.ROIs`` that
		have an attribute called ``attr``, then lets the user
		select one of them and then returns the value of 
		``attr`` of the selected ROI.
		
		"""
		
		possROIs=list(self.geometry.embryo.ROIs)
		possROIsNew=[]
		for r in possROIs:
			if hasattr(r,attr):
				possROIsNew.append(r)
		possROIs=list(possROIsNew)		
			
		if len(possROIs)<1:
			printWarning("Cannot create selection of ROIs, there is none with attribute " + attr + " .")
			return None
		
		nameList=pyfrp_misc_module.objAttrToList(possROIs,'name')
		
		selectorDialog = pyfrp_gui_basics.basicSelectorDialog(nameList,self)
		if selectorDialog.exec_():
			selectedROIName = selectorDialog.getItem()
			if selectedROIName==None:
				return None
		
		r=self.geometry.embryo.getROIByName(selectedROIName)
	
		attrVal=getattr(r,attr)
		
		return attrVal
		
	def updateCenterQles(self):
		
		"""Updates the two center QLEs with current value in geometry.center."""
		
		self.qleCenterX.setText(str(self.geometry.getCenter()[0]))
		self.qleCenterY.setText(str(self.geometry.getCenter()[1]))
		self.setCenter()
	
	def updateFnGeoLbl(self):
		self.lblFnGeoVal.setText("..."+self.geometry.fnGeo[-self.nCharDisplayed:])
		
	def drawGeometry(self):
		self.ax.clear()
		self.geometry.plotGeometry(ax=self.ax,ann=self.ann)
	
	def setAnnotate(self,val):
		self.ann=bool(2*val)
		self.drawGeometry()
		
	def setAxes3D(self):
		self.fig.delaxes(self.ax)
		self.ax=self.fig.add_subplot(111,projection='3d')
		
		self.canvas.draw()
		self.plotFrame.adjustSize()
		return self.fig
	
	def updateAllROI(self):
		self.geometry.setAllROI()
	
	#def cropCenterFromImage(self):
		#pyfrp_plot_module.
	
class zebrafishDomeStageDialog(geometryDialog):
	
	def __init__(self,geometry,parent):	
		super(zebrafishDomeStageDialog,self).__init__(geometry,parent)	
		
		#Labels
		self.lblRadius = QtGui.QLabel("Imaging Radius (px):", self)
		self.lblScale = QtGui.QLabel("Radius Scale:", self)
		self.lblHeight = QtGui.QLabel("Imaging Height (px):", self)
		
		#LineEdits
		self.qleRadius = QtGui.QLineEdit(str(self.geometry.imagingRadius))
		self.qleScale = QtGui.QLineEdit(str(self.geometry.radiusScale))
		self.qleHeight = QtGui.QLineEdit(str(self.geometry.imagingHeight))
		
		self.qleRadius.setValidator(self.doubleValid)
		self.qleScale.setValidator(self.doubleValid)
		self.qleHeight.setValidator(self.doubleValid)
		
		self.qleRadius.editingFinished.connect(self.setRadius)
		self.qleScale.editingFinished.connect(self.setScale)
		self.qleHeight.editingFinished.connect(self.setHeight)
		
		#Buttons
		self.btnRestoreDefaults=QtGui.QPushButton('Restore Defaults')
		self.btnRestoreDefaults.connect(self.btnRestoreDefaults, QtCore.SIGNAL('clicked()'), self.restoreDefaults)
		
		self.btnGetRadiusFromROI=QtGui.QPushButton('Grab from ROI')
		self.btnGetRadiusFromROI.connect(self.btnGetRadiusFromROI, QtCore.SIGNAL('clicked()'), self.getRadiusFromROI)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblRadius,nRows+1,1)
		self.grid.addWidget(self.lblScale,nRows+2,1)
		self.grid.addWidget(self.lblHeight,nRows+3,1)
		
		self.grid.addWidget(self.qleRadius,nRows+1,2)
		self.grid.addWidget(self.qleScale,nRows+2,2)
		self.grid.addWidget(self.qleHeight,nRows+3,2)

		self.grid.addWidget(self.btnRestoreDefaults,nRows+4,2)
				
		self.grid.addWidget(self.btnGetRadiusFromROI,nRows+1,3)
		
		self.show()
	
	def setRadius(self):
		self.geometry.setImagingRadius(float(str(self.qleRadius.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getImagingRadius()
	
	def setScale(self):
		self.geometry.setRadiusScale(float(str(self.qleScale.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getRadiusScale()
		
	def setHeight(self):
		self.geometry.setImagingHeight(float(str(self.qleHeight.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getImagingHeight()
	
	def restoreDefaults(self):
		self.geometry.restoreDefault()
		self.geometry.updateGeoFile()
		self.drawGeometry()

	def getRadiusFromROI(self):
		
		"""First finds all ROIs in ``geometry.embryo.ROIs`` that
		have an attribute called ``radius``, then lets the user
		select one of them and then sets radius of geometry to the same 
		radius as the ROI.
		
		"""
		
		radius=self.getAttrFromROI('radius')
		
		self.geometry.setImagingRadius(radius)
		self.updateRadiusQle()
		
	def updateRadiusQle(self):
		
		"""Updates the imagingRadius QLE with current value in geometry.radius."""
		
		self.qleRadius.setText(str(self.geometry.getImagingRadius()))
		self.setRadius()	

class cylinderDialog(geometryDialog):
	
	def __init__(self,geometry,parent):	
		super(cylinderDialog,self).__init__(geometry,parent)	
		
		#Labels
		self.lblRadius = QtGui.QLabel("Radius (px):", self)
		self.lblHeight = QtGui.QLabel("Height (px):", self)
		
		#LineEdits
		self.qleRadius = QtGui.QLineEdit(str(self.geometry.radius))
		self.qleHeight = QtGui.QLineEdit(str(self.geometry.height))
		
		self.qleRadius.setValidator(self.doubleValid)
		self.qleHeight.setValidator(self.doubleValid)
		
		self.qleRadius.editingFinished.connect(self.setRadius)
		self.qleHeight.editingFinished.connect(self.setHeight)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblRadius,nRows+1,1)
		self.grid.addWidget(self.lblHeight,nRows+3,1)
		
		self.grid.addWidget(self.qleRadius,nRows+1,2)
		self.grid.addWidget(self.qleHeight,nRows+3,2)
		
		self.show()
	
	def setRadius(self):
		self.geometry.setRadius(float(str(self.qleRadius.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getRadius()
	
	def setHeight(self):
		self.geometry.setHeight(float(str(self.qleHeight.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getHeight()
	
	def getRadiusFromROI(self):
		
		"""First finds all ROIs in ``geometry.embryo.ROIs`` that
		have an attribute called ``radius``, then lets the user
		select one of them and then sets radius of geometry to the same 
		radius as the ROI.
		
		"""
		
		radius=self.getAttrFromROI('radius')
		
		self.geometry.setRadius(radius)
		self.updateRadiusQle()
		
	def updateRadiusQle(self):
		
		"""Updates the radius QLE with current value in geometry.radius."""
		
		self.qleRadius.setText(str(self.geometry.getRadius()))
		self.setRadius()	
		
	
class coneDialog(geometryDialog):
	
	def __init__(self,geometry,parent):	
		super(coneDialog,self).__init__(geometry,parent)	
		
		#Labels
		self.lblUpperRadius = QtGui.QLabel("Upper Radius (px):", self)
		self.lblLowerRadius = QtGui.QLabel("Lower Radius (px):", self)
		
		self.lblHeight = QtGui.QLabel("Height (px):", self)
		
		#LineEdits
		self.qleUpperRadius = QtGui.QLineEdit(str(self.geometry.upperRadius))
		self.qleLowerRadius = QtGui.QLineEdit(str(self.geometry.lowerRadius))
		self.qleHeight = QtGui.QLineEdit(str(self.geometry.height))
		
		self.qleUpperRadius.setValidator(self.doubleValid)
		self.qleLowerRadius.setValidator(self.doubleValid)
		self.qleHeight.setValidator(self.doubleValid)
		
		self.qleUpperRadius.editingFinished.connect(self.setUpperRadius)
		self.qleLowerRadius.editingFinished.connect(self.setLowerRadius)
		self.qleHeight.editingFinished.connect(self.setHeight)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblUpperRadius,nRows+1,1)
		self.grid.addWidget(self.lblLowerRadius,nRows+2,1)
		self.grid.addWidget(self.lblHeight,nRows+3,1)
		
		self.grid.addWidget(self.qleUpperRadius,nRows+1,2)
		self.grid.addWidget(self.qleLowerRadius,nRows+2,2)
		self.grid.addWidget(self.qleHeight,nRows+3,2)
		
		self.show()
	
	def setUpperRadius(self):
		self.geometry.setUpperRadius(float(str(self.qleUpperRadius.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getUpperRadius()
	
	def setLowerRadius(self):
		self.geometry.setLowerRadius(float(str(self.qleLowerRadius.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getLowerRadius()
	
	
	def setHeight(self):
		self.geometry.setHeight(float(str(self.qleHeight.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getHeight()
	
	
class xenopusBallDialog(geometryDialog):
	
	def __init__(self,geometry,parent):	
		super(xenopusBallDialog,self).__init__(geometry,parent)	
		
		#Labels
		self.lblRadius = QtGui.QLabel("Imaging Radius (px):", self)
		self.lblHeight = QtGui.QLabel("Imaging Height (px):", self)
			
		#LineEdits
		self.qleRadius = QtGui.QLineEdit(str(self.geometry.imagingRadius))
		self.qleHeight = QtGui.QLineEdit(str(self.geometry.imagingHeight))
		
		self.qleRadius.setValidator(self.doubleValid)
		self.qleHeight.setValidator(self.doubleValid)
		
		self.qleRadius.editingFinished.connect(self.setRadius)
		self.qleHeight.editingFinished.connect(self.setHeight)
		
		#Buttons
		self.btnRestoreDefaults=QtGui.QPushButton('Restore Defaults')
		self.btnRestoreDefaults.connect(self.btnRestoreDefaults, QtCore.SIGNAL('clicked()'), self.restoreDefaults)
		
		self.btnGetRadiusFromROI=QtGui.QPushButton('Grab from ROI')
		self.btnGetRadiusFromROI.connect(self.btnGetRadiusFromROI, QtCore.SIGNAL('clicked()'), self.getRadiusFromROI)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblRadius,nRows+1,1)
		self.grid.addWidget(self.lblHeight,nRows+2,1)
		
		self.grid.addWidget(self.qleRadius,nRows+1,2)
		self.grid.addWidget(self.qleHeight,nRows+2,2)
		
		self.grid.addWidget(self.btnRestoreDefaults,nRows+4,2)
		
		self.grid.addWidget(self.btnGetRadiusFromROI,nRows+1,3)
		
		self.show()
		
	def setRadius(self):
		self.geometry.setImagingRadius(float(str(self.qleRadius.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getImagingRadius()
		
	def setHeight(self):
		self.geometry.setImagingHeight(float(str(self.qleHeight.text())))
		self.geometry.updateGeoFile()
		self.drawGeometry()
		return self.geometry.getImagingHeight()
	
	def restoreDefaults(self):
		self.geometry.restoreDefaults()
		self.geometry.updateGeoFile()
		self.drawGeometry()
	
	
	def getRadiusFromROI(self):
		
		"""First finds all ROIs in ``geometry.embryo.ROIs`` that
		have an attribute called ``radius``, then lets the user
		select one of them and then sets radius of geometry to the same 
		radius as the ROI.
		
		"""
		
		radius=self.getAttrFromROI('radius')
		
		self.geometry.setImagingRadius(radius)
		self.updateRadiusQle()
		
	def updateRadiusQle(self):
		
		"""Updates the imagingRadius QLE with current value in geometry.radius."""
		
		self.qleRadius.setText(str(self.geometry.getImagingRadius()))
		self.setRadius()	
	
class geometrySelectDialog(QtGui.QDialog):
	
	def __init__(self,embryo,parent):
		
		super(geometrySelectDialog,self).__init__(parent)
		
		self.geometry=embryo.geometry
		self.embryo=embryo
		
		
		#Labels
		self.lblTyp = QtGui.QLabel("Geometry Type:", self)
		
		#Combo
		self.comboGeometry = QtGui.QComboBox(self)
		self.comboGeometry.addItem("zebraFishDomeStage")
		self.comboGeometry.addItem("cylinder")
		self.comboGeometry.addItem("xenopusBall")
		self.comboGeometry.addItem("cone")
		self.comboGeometry.addItem("custom")
		self.comboGeometry.addItem("zebraFishDomeStageQuad (experimental)")
		self.comboGeometry.addItem("cylinderQuad (experimental)")
		self.comboGeometry.addItem("xenopusBallQuad (experimental)")
		
		self.initComboGeometry()
	
		self.comboGeometry.activated[str].connect(self.setGeometry)   
		
		#Buttons
		self.btnDone=QtGui.QPushButton('Done')	
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Layout
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.lblTyp,1,1)
		
		self.grid.addWidget(self.comboGeometry,1,2)
		self.grid.addWidget(self.btnDone,2,2)
		
		self.setLayout(self.grid)    
			
		self.setWindowTitle('Select Geometry')    
		self.show()
	
	def initComboGeometry(self):
		
		
		
		if self.geometry!=None:
			
			idxs=[]
			matchLenghts=[]
			for i in range(self.comboGeometry.count()):
				text=str(self.comboGeometry.itemText(i))
				
				if self.geometry.typ in text:
					idxs.append(i)
					matchLenghts.append(len(text.replace(self.geometry.typ,"")))
			
			if len(idxs)>0:	
				idx=idxs[matchLenghts.index(max(matchLenghts))]
			self.comboGeometry.setCurrentIndex(idx)
		else:
			idx=0
			self.comboGeometry.setCurrentIndex(idx)
			self.setGeometry(str(self.comboGeometry.currentText()))
		
	
	def getOldAttr(self):
		if self.geometry!=None:
			if self.geometry.typ!='custom':
			
				if hasattr(self.geometry,'radius'):
					radius=self.geometry.getRadius()
				elif hasattr(self.geometry,'imagingRadius'):
					radius=self.geometry.getImagingRadius()
				elif hasattr(self.geometry,'upperRadius'):
					radius=self.geometry.computeRadiusFromSliceHeight(self.embryo.sliceHeightPx)
				else:
					radius=300.
					
				height=-30.	
				if hasattr(self.geometry,'imagingHeight'):
					height=self.geometry.getImagingHeight()
					
				return self.geometry.getCenter(), radius, height
			
			else:
				self.geometry.getCenter(),300.,-30.
			
		return [256,256],300.,-30.
			
				
	def setGeometry(self,text):
		text=str(text)
		text=text.replace("(experimental)","").strip()
		
		center,radius,height=self.getOldAttr()
		
		if text=='zebraFishDomeStage':
			
			self.embryo.setGeometry2ZebraFishDomeStage(center,radius,radiusScale=1.1)
			
		elif text=='cylinder':
			self.embryo.setGeometry2Cylinder(center,radius,height)
		elif text=='xenopusBall':
			self.embryo.setGeometry2Ball(center,radius)
		elif text=='custom':
			self.embryo.setGeometry2Custom(center)
		elif text=='zebraFishDomeStageQuad':
			self.embryo.setGeometry2ZebraFishDomeStageQuad(center,radius,radiusScale=1.1)
		elif text=='cylinderQuad':		
			self.embryo.setGeometry2CylinderQuad(center,radius,height)
		elif text=='xenopusBallQuad':
			self.embryo.setGeometry2BallQuad(center,radius)
		elif text=='cone':
			self.embryo.setGeometry2Cone(center,radius,radius,height)
		else:
			printWarning("Unknown geometry type " + text)
	
	def donePressed(self):
		self.done(1)
		return 
	
		


		
		
		
		
	

		
	

		
	
	
	
	
		
		
		
		
		
		
			
		
		
		
		
		
		
		
		
		
		
		
		
		