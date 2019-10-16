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

#PyQT Dialogs for managing and editing ROIs
#(1)  ROImanager
#(2)  ROIDialog
#(3)  sliceROIDialog
#(4)  radialROIDialog
#(5)  radialSliceROIDialog
#(6)  squareROIDialog
#(7)  squareSliceROIDialog
#(8)  rectangleROIDialog
#(9)  rectangleSliceROIDialog
#(10) polygonROIDialog
#(11) polygonSliceROIDialog
#(12) customROIDialog

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

#Matplotlib
import matplotlib

#Numpy/Scipy
import numpy as np

#Matplotlib
import matplotlib.patches as ptc

#===================================================================================================================================
#Main ROI Manager
#===================================================================================================================================

class ROImanager(QtGui.QDialog):
	
	def __init__(self,embryo,parent):
		super(ROImanager,self).__init__(parent)
		
		self.embryo=embryo
		self.parent=parent
		
		
		
		#-------------------------------------------------------------------------------------------------------------------
		#TreeWidget
		#-------------------------------------------------------------------------------------------------------------------
		
		self.ROIList=QtGui.QTreeWidget()
		self.ROIList.setHeaderLabels(["Name","Typ"])
		self.ROIList.setColumnWidth(0,200)
		self.ROIList.itemClicked.connect(self.roiSelected)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		self.btnDone=QtGui.QPushButton('Done')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		self.btnNewROI=QtGui.QPushButton('New')
		self.btnNewROI.connect(self.btnNewROI, QtCore.SIGNAL('clicked()'), self.newROI)
		
		self.btnRemoveROI=QtGui.QPushButton('Remove')
		self.btnRemoveROI.connect(self.btnRemoveROI, QtCore.SIGNAL('clicked()'), self.removeROI)
		
		self.btnEditROI=QtGui.QPushButton('Edit')
		self.btnEditROI.connect(self.btnEditROI, QtCore.SIGNAL('clicked()'), self.editROI)
		
		self.btnComputeIdxs=QtGui.QPushButton('Compute Idxs')
		self.btnComputeIdxs.connect(self.btnComputeIdxs, QtCore.SIGNAL('clicked()'), self.computeIdxs)
		
		self.btnComputeAllIdxs=QtGui.QPushButton('Compute Idxs')
		self.btnComputeAllIdxs.connect(self.btnComputeAllIdxs, QtCore.SIGNAL('clicked()'), self.computeAllIdxs)
		
		self.btnShowIdxs=QtGui.QPushButton('Show Idxs')
		#self.btnShowIdxs.connect(self.btnShowIdxs, QtCore.SIGNAL('clicked()'), self.showIdxs)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Combo Box
		#-------------------------------------------------------------------------------------------------------------------
		
		self.comboType = QtGui.QComboBox(self)
		self.comboType.addItem("radial")
		self.comboType.addItem("slice")
		self.comboType.addItem("radialSlice")
		self.comboType.addItem("square")
		self.comboType.addItem("squareSlice")
		self.comboType.addItem("rectangle")
		self.comboType.addItem("rectangleSlice")
		self.comboType.addItem("polygon")
		self.comboType.addItem("polygonSlice")
		self.comboType.addItem("custom")
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.grid = QtGui.QGridLayout()
		
		self.grid.addWidget(self.btnNewROI,1,1)
		self.grid.addWidget(self.btnEditROI,2,1)
		self.grid.addWidget(self.btnRemoveROI,3,1)
		self.grid.addWidget(self.btnComputeIdxs,4,1)
		self.grid.addWidget(self.btnComputeIdxs,5,1)
		self.grid.addWidget(self.btnShowIdxs,6,1)
		
		self.grid.addWidget(self.comboType,2,2)
		
		
		self.hbox = QtGui.QHBoxLayout()
		self.hbox.addLayout(self.grid)
		self.hbox.addWidget(self.ROIList)
		
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addLayout(self.hbox)
		self.vbox.addWidget(self.btnDone)
		
		self.setLayout(self.vbox)    
		
		self.setWindowTitle('ROI Manager')    
		
		self.updateROIList()
		
		self.setMinimumSize(650,500) 
		self.resize(650,500)
		
		self.show()
		
	def updateROIList(self):
		
		self.ROIList.clear()
		for r in self.embryo.ROIs:
			QtGui.QTreeWidgetItem(self.ROIList,[r.name,r.getType()])
		return self.ROIList
	
	def newROI(self):
		newID=self.embryo.getFreeROIId()
		typ=self.comboType.currentText()
		
		if self.embryo.geometry!=None:
			center=self.embryo.geometry.getCenter()
		else:
			center=[self.embryo.dataResPx]*2
		
		sidelength=100.
		offset=[center[0]-sidelength/2.,center[1]-sidelength/2.]
				
		if typ=='slice':
			self.embryo.newSliceROI('newSliceROI',newID,50.,20.,False)	
		elif typ=='radial':
			self.embryo.newRadialROI('newRadialROI',newID,center,200.)
		elif typ=='radialSlice':
			self.embryo.newRadialSliceROI('newRadialSliceROI',newID,center,200.,50.,20.,False)
		elif typ=='square':
			self.embryo.newSquareROI('newSquareROI',newID,offset,sidelength,200.)
		elif typ=='squareSlice':
			self.embryo.newSquareSliceROI('newSquareSliceROI',newID,offset,sidelength,50.,20.,False)
		elif typ=='rectangle':
			self.embryo.newRectangleROI('newRectangleROI',newID,offset,sidelength,sidelength)
		elif typ=='rectangleSlice':
			self.embryo.newRectangleSliceROI('newRectangleSliceROI',newID,offset,sidelength,sidelength,50.,20.,False)
		elif typ=='polygon':
			self.embryo.newPolyROI('newPolyROI',newID,[])
		elif typ=='polygonSlice':
			self.embryo.newPolySliceROI('newPolySliceROI',newID,[],50.,20.,False)	
		elif typ=='custom':
			self.embryo.newCustomROI('newCustomROI',newID)
		else:
			printWarning('Unknown ROI Type' + typ)
	
		self.openDialog(self.embryo.ROIs[-1],typ)
		
		self.updateROIList()
		
	def editROI(self):
		typ=self.currROI.getType()
	
		self.openDialog(self.currROI,typ)
		
		self.updateROIList()
	
	def openDialog(self,ROI,typ):
		if typ=='slice':
			ret=sliceROIDialog(ROI,self).exec_()
		elif typ=='radial':
			ret=radialROIDialog(ROI,self).exec_()
		elif typ=='radialSlice':
			ret=radialSliceROIDialog(ROI,self).exec_()
		elif typ=='square':
			ret=squareROIDialog(ROI,self).exec_()
		elif typ=='squareSlice':
			ret=squareSliceROIDialog(ROI,self).exec_()
		elif typ=='rectangle':
			ret=rectangleROIDialog(ROI,self).exec_()
		elif typ=='rectangleSlice':
			ret=rectangleSliceROIDialog(ROI,self).exec_()
		elif typ=='polygon':
			ret=polygonROIDialog(ROI,self).exec_()
		elif typ=='polygonSlice':
			ret=polygonSliceROIDialog(ROI,self).exec_()
		elif typ=='custom':
			ret=customROIDialog(ROI,self).exec_()
		else:
			printWarning('Unknown ROI Type' + typ)
	
	def removeROI(self):
		if len(self.currROI.findIncluded())>0:
			printWarning("ROI " + self.currROI.name + " is used in custom ROI. This might lead to problems." )
		self.embryo.ROIs.remove(self.currROI)
		self.updateROIList()
	
	def computeIdxs(self):
		self.currROI.computeIdxs()
	
	def computeAllIdxs(self):
		self.embryo.computeROIIdxs()
	
	#def showIdxs(self): ###NOTE: Still needs to be done
		#self.parent.
	
	def roiSelected(self):
		self.currNode=self.ROIList.currentItem()
		ind=self.ROIList.indexOfTopLevelItem(self.currNode)
		self.currROI=self.embryo.ROIs[ind]
		return self.currROI
	
	
		
	def donePressed(self):
		self.done(1)
		return self.embryo

#===================================================================================================================================
#Basic ROI Dialog
#===================================================================================================================================
	
class ROIDialog(pyfrp_gui_basics.basicCanvasDialog):
	
	def __init__(self,ROI,parent):
		super(ROIDialog,self).__init__(parent)
		
		self.ROI=ROI
		
		#Labels
		self.lblName = QtGui.QLabel("Name:", self)
		self.lblZmin = QtGui.QLabel("zmin:", self)
		self.lblZmax = QtGui.QLabel("zmax:", self)
		self.lblColor = QtGui.QLabel("Color:", self)
	
		#LineEdits
		self.qleName = QtGui.QLineEdit(self.ROI.name)
		self.qleZmin = QtGui.QLineEdit(str(self.ROI.zmin))
		self.qleZmax = QtGui.QLineEdit(str(self.ROI.zmax))
		
		self.doubleValid=QtGui.QDoubleValidator()
		
		self.qleZmin.setValidator(self.doubleValid)
		self.qleZmax.setValidator(self.doubleValid)
		
		self.qleName.editingFinished.connect(self.setName)
		self.qleZmin.editingFinished.connect(self.setZmin)
		self.qleZmax.editingFinished.connect(self.setZmax)
		
		#Checkboxes
		self.cbUseForRim = QtGui.QCheckBox('Use for rim?', self)
		self.cbUseForRim.setCheckState(2*int(self.ROI.useForRim))
		
		self.connect(self.cbUseForRim, QtCore.SIGNAL('stateChanged(int)'), self.checkUseForRim)
		
		#Buttons
		self.initColorButton()
		self.btnColor.connect(self.btnColor, QtCore.SIGNAL('clicked()'), self.setColor)
		
		#Layout
		self.grid.addWidget(self.lblName,1,1)
		self.grid.addWidget(self.lblZmin,2,1)
		self.grid.addWidget(self.lblZmax,3,1)
		self.grid.addWidget(self.lblColor,4,1)
		
		
		self.grid.addWidget(self.qleName,1,2)
		self.grid.addWidget(self.qleZmin,2,2)
		self.grid.addWidget(self.qleZmax,3,2)
		self.grid.addWidget(self.btnColor,4,2)
		
		
		self.grid.addWidget(self.cbUseForRim,5,2)
		
		self.showFirstDataImg()
		
		self.setWindowTitle('Edit ROI ' + self.ROI.name)   
		
		
		self.show()
	
	def setName(self):
		self.ROI.setName(self.qleName.text())
		return self.ROI.getName()
	
	def setZmin(self):
		self.ROI.setZExtend(float(self.qleZmin.text()),self.ROI.zmax)
		return self.ROI.getZExtend()
	
	def setZmax(self):
		self.ROI.setZExtend(self.ROI.zmin,float(self.qleZmax.text()))
		return self.ROI.getZExtend()
	
	def updateZExtendQles(self):
		self.qleZmax.setText(str(self.ROI.getZExtend()[1]))
		self.qleZmin.setText(str(self.ROI.getZExtend()[0]))
		
	def checkUseForRim(self,val):
		self.ROI.setUseForRim(bool(2*val))
		return self.ROI.getUseForRim()
	
	def initColorButton(self):
		self.btnColor=QtGui.QPushButton('         ')
		self.updateColorButton()
		return
	
	def updateColorButton(self):
		colorRGB=matplotlib.colors.colorConverter.to_rgb(self.ROI.color)
		colorHex=matplotlib.colors.rgb2hex(colorRGB) 
		self.btnColor.setStyleSheet("background-color: "+colorHex)
		
	def setColor(self):
		col = QtGui.QColorDialog.getColor(parent=self)	
		col=tuple(np.asarray([col.red(),col.green(),col.blue()])/256.)
		self.ROI.setColor(col)
		self.updateColorButton()
		
		self.setArtistColor(col)
		
		return col
	
	def drawPoint(self,x,y,color='r',idx=0):
		pt=ptc.Circle([x,y],radius=3,fill=True,color=color)
		self.replaceArtist(idx,self.ax.add_patch(pt))
		self.canvas.draw()
		return pt
	
	def drawCircle(self,center,radius,color='r',idx=1):
		c=ptc.Circle(center,radius=radius,fill=False,color=color)
		self.replaceArtist(idx,self.ax.add_patch(c))
		self.canvas.draw()
		return c
	
	def drawRectangle(self,offset,sidelengthX,sidelengthY,color='r',idx=1):
		r=ptc.Rectangle(offset,sidelengthX,sidelengthY,fill=False,color=color,linewidth=3)
		self.replaceArtist(idx,self.ax.add_patch(r))
		self.canvas.draw()
		return r
	
	def showFirstDataImg(self):
		
		self.ROI.embryo.updateFileList()
		
		fnImg=self.ROI.embryo.fnDatafolder
		if fnImg[-1]!='/':
			fnImg=fnImg+'/'
		fnImg=fnImg+self.ROI.embryo.fileList[0]
		
		img=pyfrp_img_module.loadImg(fnImg,self.ROI.embryo.dataEnc)
	
		self.showImg(img)


#===================================================================================================================================
#ROI Dialogs for specific ROI types
#===================================================================================================================================
	
	
class sliceROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(sliceROIDialog,self).__init__(ROI,parent)
		
		#Labels
		self.lblHeight = QtGui.QLabel("Height:", self)
		self.lblWidth = QtGui.QLabel("Width:", self)
		
		#LineEdits
		self.qleHeight = QtGui.QLineEdit(str(self.ROI.height))
		self.qleWidth = QtGui.QLineEdit(str(self.ROI.width))
		
		self.qleHeight.setValidator(self.doubleValid)
		self.qleWidth.setValidator(self.doubleValid)
		
		self.qleHeight.editingFinished.connect(self.setHeight)
		self.qleWidth.editingFinished.connect(self.setWidth)
		
		#Checkboxes
		self.cbSliceBottom = QtGui.QCheckBox('SliceBottom?', self)
		self.cbSliceBottom.setCheckState(2*int(self.ROI.sliceBottom))
		
		self.connect(self.cbSliceBottom, QtCore.SIGNAL('stateChanged(int)'), self.checkSliceBottom)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblHeight,nRows+1,1)
		self.grid.addWidget(self.lblWidth,nRows+2,1)
		
		self.grid.addWidget(self.qleHeight,nRows+1,2)
		self.grid.addWidget(self.qleWidth,nRows+2,2)
		
		self.grid.addWidget(self.cbSliceBottom,nRows+3,2)
		
		self.updateZExtendQles()
		
		self.show()
		
	def setHeight(self):
		self.ROI.setHeight(float(self.qleHeight.text()))
		self.updateZExtendQles()
		return self.ROI.getHeight()
	
	def setWidth(self):
		self.ROI.setWidth(float(self.qleWidth.text()))
		self.updateZExtendQles()
		return self.ROI.getWidth()
		
	def checkSliceBottom(self,val):
		self.ROI.setSliceBottom(bool(2*val))
		self.updateZExtendQles()
		return self.ROI.getSliceBottom()
	
class radialROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(radialROIDialog,self).__init__(ROI,parent)	
		
		#Labels
		self.lblRadius = QtGui.QLabel("Radius:", self)
		self.lblCenter = QtGui.QLabel("Center:", self)
		
		#LineEdits
		self.qleRadius = QtGui.QLineEdit(str(self.ROI.radius))
		self.qleCenterX = QtGui.QLineEdit(str(self.ROI.center[0]))
		self.qleCenterY = QtGui.QLineEdit(str(self.ROI.center[1]))
		
		self.qleRadius.setValidator(self.doubleValid)
		self.qleCenterX.setValidator(self.doubleValid)
		self.qleCenterY.setValidator(self.doubleValid)
		
		self.centerGrid = QtGui.QGridLayout()
		
		self.centerGrid.addWidget(self.qleCenterX)
		self.centerGrid.addWidget(self.qleCenterY)
		
		self.qleCenterX.editingFinished.connect(self.setCenter)
		self.qleCenterY.editingFinished.connect(self.setCenter)
		self.qleRadius.editingFinished.connect(self.setRadius)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblRadius,nRows+1,1)
		self.grid.addWidget(self.lblCenter,nRows+2,1)
		
		self.grid.addWidget(self.qleRadius,nRows+1,2)
		self.grid.addLayout(self.centerGrid,nRows+2,2)
		
		#Connect Canvas
		self.connectCanvas()
		
		self.show()
	
	def setCenter(self):
		center=[float(self.qleCenterX.text()),float(self.qleCenterY.text())]
		self.ROI.setCenter(center)
		self.drawCenter()
		self.drawRadius()
		return center
	
	def setRadius(self):
		self.ROI.setRadius(float(self.qleRadius.text()))
		self.drawRadius()
		return self.ROI.getRadius()
	
	def computeRadiusFromCoordinate(self,x,y):
		return np.sqrt((x-self.ROI.getCenter()[0])**2+(y-self.ROI.getCenter()[1])**2)
	
	def drawCenter(self):
		pt=self.drawPoint(self.ROI.getCenter()[0],self.ROI.getCenter()[1],color=self.ROI.color)
		return pt
	
	def drawRadius(self):
		c=self.drawCircle(self.ROI.getCenter(),self.ROI.getRadius(),color=self.ROI.color)
		return c
			
	def getMouseCanvas(self,event):
		
		"""Overwrite function for mouse input. 
		
		Directs left clicks to creation of new artist, right click to deletion of artist.
		
		"""
		
		#Left click to define center and then radius
		if event.button==1:
			
			#Check if clicked withing axes
			if event.xdata==None:
				return
			
			#Check if this is center 
			if len(self.artists)==0:
				
				#Update qles
				self.updateCenterQles([event.xdata,event.ydata])
				
				#Set center (and draw)
				self.setCenter()
					
			#Check if this is radius
			elif len(self.artists)>1:
				
				#Update radius
				radius=self.computeRadiusFromCoordinate(event.xdata,event.ydata)
				
				#Update qle
				self.updateRadiusQle(radius)
				
				#Set radius (and draw)
				self.setRadius()
		
		#Remove last artist
		elif event.button==3:
			self.removeArtist()
	
	def updateCenterQles(self,center):
		
		"""Updates the two center Qles with new center."""
		
		self.qleCenterX.setText(str(center[0]))
		self.qleCenterY.setText(str(center[1]))
	
	def updateRadiusQle(self,radius):
		
		"""Updates radius Qle with given radius."""
		
		self.qleRadius.setText(str(radius))
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseRadius()
		elif event.key=='ctrl+down':
			self.decreaseRadius()
			
	def moveLeft(self):
		
		"""Moves center 1 px to the left."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.ROI.getCenter()[0]-1,self.ROI.getCenter()[1]])
			self.setCenter()
		
		
	def moveRight(self):
		
		"""Moves center 1 px to the right."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.ROI.getCenter()[0]+1,self.ROI.getCenter()[1]])
			self.setCenter()
	
	def moveUp(self):
		
		"""Moves center 1 px up."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.ROI.getCenter()[0],self.ROI.getCenter()[1]+1])
			self.setCenter()
	
	def moveDown(self):
		
		"""Moves center 1 px down."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.ROI.getCenter()[0],self.ROI.getCenter()[1]-1])
			self.setCenter()
			
	def increaseRadius(self):
		
		"""Increases radius by 1 px."""
		
		if len(self.artists)>1:
			self.updateRadiusQle(self.ROI.getRadius()+1)
			self.setRadius()
			
	def decreaseRadius(self):
		
		"""Decreases radius by 1 px."""
		
		if len(self.artists)>1:
			self.updateRadiusQle(self.ROI.getRadius()-1)
			self.setRadius()
			
class radialSliceROIDialog(radialROIDialog,sliceROIDialog):
	
	def __init__(self,ROI,parent):
			
		super(radialSliceROIDialog,self).__init__(ROI,parent)
		#sliceROIDialog.__init__(self,ROI,parent)
		#radialROIDialog.__init__(self,ROI,parent)
	
class squareROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(squareROIDialog,self).__init__(ROI,parent)	
		
		#Labels
		self.lblSidelength = QtGui.QLabel("Sidelength:", self)
		self.lblOffset = QtGui.QLabel("Offset:", self)
		
		#LineEdits
		self.qleSidelength = QtGui.QLineEdit(str(self.ROI.sidelength))
		self.qleOffsetX = QtGui.QLineEdit(str(self.ROI.offset[0]))
		self.qleOffsetY = QtGui.QLineEdit(str(self.ROI.offset[1]))
		
		self.qleSidelength.setValidator(self.doubleValid)
		self.qleOffsetX.setValidator(self.doubleValid)
		self.qleOffsetY.setValidator(self.doubleValid)
		
		self.offsetGrid = QtGui.QGridLayout()
		
		self.offsetGrid.addWidget(self.qleOffsetX)
		self.offsetGrid.addWidget(self.qleOffsetY)
		
		self.qleOffsetX.editingFinished.connect(self.setOffset)
		self.qleOffsetY.editingFinished.connect(self.setOffset)
		self.qleSidelength.editingFinished.connect(self.setSideLength)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblSidelength,nRows+1,1)
		self.grid.addWidget(self.lblOffset,nRows+2,1)
		
		self.grid.addWidget(self.qleSidelength,nRows+1,2)
		self.grid.addLayout(self.offsetGrid,nRows+2,2)
		
		#Connect Canvas
		self.connectCanvas()
	
	def setOffset(self):
		offset=[float(self.qleOffsetX.text()),float(self.qleOffsetY.text())]
		self.ROI.setOffset(offset)
		self.drawOffset()
		self.drawSidelength()
		return offset
	
	def setSideLength(self):
		self.ROI.setSideLength(float(self.qleSidelength.text()))
		self.drawSidelength()
		return self.ROI.getSideLength()
	
	def computeSidelengthFromCoordinate(self,x,y):
		return max([x,y])
	
	def drawOffset(self):
		pt=self.drawPoint(self.ROI.getOffset()[0],self.ROI.getOffset()[1],color=self.ROI.color)
		return pt
	
	def drawSidelength(self):
		r=self.drawRectangle(self.ROI.getOffset(),self.ROI.getSideLength(),self.ROI.getSideLength(),color=self.ROI.color)
		return r
	
	def getMouseCanvas(self,event):
		
		if event.button==1:
	
			if event.xdata==None:
				return
			
			if len(self.artists)==0:
			
				self.updateOffsetQles([event.xdata,event.ydata])
				self.setOffset()
					
			elif len(self.artists)>1:
				
		
				sidelength=self.computeSidelengthFromCoordinate(event.xdata,event.ydata)
				self.updateSidelengthQle(sidelength)
				self.setSideLength()
			
		#Remove last artist
		elif event.button==3:
			self.removeArtist()
	
	def updateOffsetQles(self,offset):
		
		"""Updates the two offset Qles with new offset."""
		
		self.qleOffsetX.setText(str(offset[0]))
		self.qleOffsetY.setText(str(offset[1]))
	
	def updateSidelengthQle(self,sidelength):
		
		"""Updates sidelength Qle with given sidelength."""
		
		self.qleSidelength.setText(str(sidelength))
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseSidelength()
		elif event.key=='ctrl+down':
			self.decreaseSidelength()
			
	def moveLeft(self):
		
		"""Moves offset 1 px to the left."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0]-1,self.ROI.getOffset()[1]])
			self.setOffset()
			
	def moveRight(self):
		
		"""Moves offset 1 px to the right."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0]+1,self.ROI.getOffset()[1]])
			self.setOffset()
	
	def moveUp(self):
		
		"""Moves offset 1 px up."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0],self.ROI.getOffset()[1]+1])
			self.setOffset()
	
	def moveDown(self):
		
		"""Moves offset 1 px down."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0],self.ROI.getOffset()[1]-1])
			self.setOffset()
			
	def increaseSidelength(self):
		
		"""Increases sidelength by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLength()+1)
			self.setSideLength()
			
	def decreaseSidelength(self):
		
		"""Decreases sidelength by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLength()-1)
			self.setSideLength()
	
	
class squareSliceROIDialog(squareROIDialog,sliceROIDialog):
	
	def __init__(self,ROI,parent):
			
		super(squareSliceROIDialog,self).__init__(ROI,parent)	
		#sliceROIDialog.__init__(self,ROI,parent)
		#squareROIDialog.__init__(self,ROI,parent)	

class rectangleROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(rectangleROIDialog,self).__init__(ROI,parent)	
		
		#Labels
		self.lblSidelengthX = QtGui.QLabel("SidelengthX:", self)
		self.lblSidelengthY = QtGui.QLabel("SidelengthY:", self)
		self.lblOffset = QtGui.QLabel("Offset:", self)
		
		#LineEdits
		self.qleSidelengthX = QtGui.QLineEdit(str(self.ROI.sidelengthX))
		self.qleSidelengthY = QtGui.QLineEdit(str(self.ROI.sidelengthY))
		self.qleOffsetX = QtGui.QLineEdit(str(self.ROI.offset[0]))
		self.qleOffsetY = QtGui.QLineEdit(str(self.ROI.offset[1]))
		
		self.qleSidelengthX.setValidator(self.doubleValid)
		self.qleSidelengthY.setValidator(self.doubleValid)
		self.qleOffsetX.setValidator(self.doubleValid)
		self.qleOffsetY.setValidator(self.doubleValid)
		
		self.offsetGrid = QtGui.QGridLayout()
		
		self.offsetGrid.addWidget(self.qleOffsetX)
		self.offsetGrid.addWidget(self.qleOffsetY)
		
		self.qleOffsetX.editingFinished.connect(self.setOffset)
		self.qleOffsetY.editingFinished.connect(self.setOffset)
		self.qleSidelengthX.editingFinished.connect(self.setSideLengthX)
		self.qleSidelengthY.editingFinished.connect(self.setSideLengthY)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblSidelengthX,nRows+1,1)
		self.grid.addWidget(self.lblSidelengthY,nRows+2,1)
		self.grid.addWidget(self.lblOffset,nRows+3,1)
		
		self.grid.addWidget(self.qleSidelengthX,nRows+1,2)
		self.grid.addWidget(self.qleSidelengthY,nRows+2,2)
	
		self.grid.addLayout(self.offsetGrid,nRows+3,2)
		
		#Connect Canvas
		self.connectCanvas()
		
		self.show()
	
	def setOffset(self):
		offset=[float(self.qleOffsetX.text()),float(self.qleOffsetY.text())]
		self.ROI.setOffset(offset)
		self.drawOffset()
		self.drawSidelength()
		return offset
	
	def setSideLengthX(self):
		self.ROI.setSideLengthX(float(self.qleSidelengthX.text()))
		self.drawSidelength()
		return self.ROI.getSideLengthX()
	
	def setSideLengthY(self):
		self.ROI.setSideLengthY(float(self.qleSidelengthY.text()))
		self.drawSidelength()
		return self.ROI.getSideLengthY()
	
	def computeSidelengthsFromCoordinate(self,x,y):
		return x-self.ROI.offset[0],y-self.ROI.offset[1]
	
	def drawOffset(self):
		pt=self.drawPoint(self.ROI.getOffset()[0],self.ROI.getOffset()[1],color=self.ROI.color)
		return pt
	
	def drawSidelength(self):
		r=self.drawRectangle(self.ROI.getOffset(),self.ROI.getSideLengthX(),self.ROI.getSideLengthY(),color=self.ROI.color)
		return r
	
	def getMouseCanvas(self,event):
		
		if event.button==1:
	
			if event.xdata==None:
				return
			
			if len(self.artists)==0:
			
				self.updateOffsetQles([event.xdata,event.ydata])
				self.setOffset()
					
			elif len(self.artists)>1:
				
				sidelengthX,sidelengthY=self.computeSidelengthsFromCoordinate(event.xdata,event.ydata)
		
				self.updateSidelengthQle(sidelengthX,sidelengthY)
				
				self.setSideLengthX()
				self.setSideLengthY()
				
		#Remove last artist
		elif event.button==3:
			self.removeArtist()
	
	def updateOffsetQles(self,offset):
		
		"""Updates the two offset Qles with new offset."""
		
		self.qleOffsetX.setText(str(offset[0]))
		self.qleOffsetY.setText(str(offset[1]))
	
	def updateSidelengthQle(self,sidelengthX,sidelengthY):
		
		"""Updates sidelength Qle with given sidelength."""
		
		self.qleSidelengthX.setText(str(sidelengthX))
		self.qleSidelengthY.setText(str(sidelengthY))
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseSidelengthY()
		elif event.key=='ctrl+down':
			self.decreaseSidelengthY()
		elif event.key=='ctrl+right':
			self.increaseSidelengthX()
		elif event.key=='ctrl+left':
			self.decreaseSidelengthX()
			
	def moveLeft(self):
		
		"""Moves offset 1 px to the left."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0]-1,self.ROI.getOffset()[1]])
			self.setOffset()
			
	def moveRight(self):
		
		"""Moves offset 1 px to the right."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0]+1,self.ROI.getOffset()[1]])
			self.setOffset()
	
	def moveUp(self):
		
		"""Moves offset 1 px up."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0],self.ROI.getOffset()[1]+1])
			self.setOffset()
	
	def moveDown(self):
		
		"""Moves offset 1 px down."""
		
		if len(self.artists)>0:
			self.updateOffsetQles([self.ROI.getOffset()[0],self.ROI.getOffset()[1]-1])
			self.setOffset()
			
	def increaseSidelengthY(self):
		
		"""Increases sidelengthY by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLengthX(),self.ROI.getSideLengthY()+1)
			self.setSideLengthY()
			
	def decreaseSidelengthY(self):
		
		"""Decreases sidelengthY by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLengthX(),self.ROI.getSideLengthY()-1)
			self.setSideLengthY()
			
	def increaseSidelengthX(self):
		
		"""Increases sidelengthX by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLengthX()+1,self.ROI.getSideLengthY())
			self.setSideLengthX()
			
	def decreaseSidelengthX(self):
		
		"""Decreases sidelengthX by 1 px."""
		
		if len(self.artists)>1:
			self.updateSidelengthQle(self.ROI.getSideLengthX()-1,self.ROI.getSideLengthY())
			self.setSideLengthX()
	
class rectangleSliceROIDialog(rectangleROIDialog,sliceROIDialog):
	
	def __init__(self,ROI,parent):
		
		super(rectangleSliceROIDialog,self).__init__(ROI,parent)
		
		#sliceROIDialog.__init__(self,ROI,parent)
		#rectangleROIDialog.__init__(self,ROI,parent)			
			
	
class polygonROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(polygonROIDialog,self).__init__(ROI,parent)	
		
		#bookkeeping
		self.currCorner=None
		self.edges=None
		self.highlighted=None
		self.currInd=None
		
		#Labels
		self.lblCorners = QtGui.QLabel("Corners:", self)
		
		#TreeWidget
		self.cornerList=QtGui.QTreeWidget()
		self.cornerList.setHeaderLabels(["x","y"])
		self.cornerList.itemClicked.connect(self.cornerSelected)
		
		#Buttons
		self.btnRemoveCorner=QtGui.QPushButton('Remove')
		self.btnRemoveCorner.connect(self.btnRemoveCorner, QtCore.SIGNAL('clicked()'), self.removeCornerPressed)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblCorners,nRows+1,2)
		self.grid.addWidget(self.cornerList,nRows+2,2)
		
		self.grid.addWidget(self.btnRemoveCorner,nRows+2,1)
		
		#Connect Canvas
		self.connectCanvas()
		
		
		self.show()	
	
	def updateCornerList(self):
		self.cornerList.clear()
		for corner in self.ROI.corners:
			QtGui.QTreeWidgetItem(self.cornerList,[str(corner[0]),str(corner[1])])
		return self.cornerList
	
	def cornerSelected(self):
		self.currNode=self.cornerList.currentItem()
		ind=self.cornerList.indexOfTopLevelItem(self.currNode)
		self.currCorner=self.ROI.corners[ind]
		self.currInd=ind
		
		self.removeHighlighted()
		
		self.highlightCorner(self.currCorner)
		
		return self.currCorner
	
	def newCorner(self,x,y):
		self.ROI.appendCorner([x,y])
		self.drawCorner([x,y])
		self.updateCornerList()
		
		self.drawPolygon()
		self.removeHighlighted()
		self.highlightCorner([x,y])
		self.currCorner=self.ROI.corners[-1]
		self.currInd=-1
		return [x,y]
	
	def removeCornerPressed(self):
		self.currNode=self.cornerList.currentItem()
		ind=self.cornerList.indexOfTopLevelItem(self.currNode)
		self.removeCorner(ind)
		self.removeHighlighted()
	
	def removeCorner(self,idx):
		if len(self.ROI.corners)<idx or len(self.ROI.corners)==0:
			return
		self.ROI.removeCorner(idx)
		self.updateCornerList()
		
		self.drawPolygon()
		
	def drawCorner(self,corner):
		pt=ptc.Circle([corner[0],corner[1]],radius=3,fill=True,color=self.ROI.color)
		self.artists.append(self.ax.add_patch(pt))
		return pt
	
	def highlightCorner(self,corner):
		c=ptc.Circle([corner[0],corner[1]],radius=6,fill=False,color=self.ROI.color)
		self.highlighted=self.ax.add_patch(c)
		self.canvas.draw()
		return c
	
	def removeHighlighted(self):
		if self.highlighted!=None:
			try:
				self.highlighted.remove()
			except ValueError:
				self.canvas.draw()
				#pass
				
			return self.highlighted
	
	def drawAllEdges(self,closed=True):
		if self.edges!=None:
			self.edges.remove()
			
		if len(self.ROI.corners)<2:
			self.edges=None
			return self.edges
		
		temp=list(self.ROI.corners)
		if closed:
			temp.append(self.ROI.corners[0])
		
		temp=np.asarray(temp)
		self.edges,=self.ax.plot(temp[:,0],temp[:,1],color=self.ROI.color)
		
		return self.edges
	
	def drawPolygon(self):
	
		self.removeAllCorners()
		
		for c in self.ROI.corners:
			self.drawCorner(c)
				
		self.drawAllEdges()
		
		self.canvas.draw()
				
	def removeAllCorners(self):
		for i in range(len(self.artists)):
			self.removeArtist()
		return self.artists	
					
		
	def getMouseCanvas(self,event):
		
		if event.button==1:
	
			if event.xdata==None:
				return
			
			self.newCorner(event.xdata,event.ydata)
						
		elif event.button==3:
			self.removeCorner(-1)
			
	def setColor(self):
		col = QtGui.QColorDialog.getColor(parent=self)	
		col=tuple(np.asarray([col.red(),col.green(),col.blue()])/256.)
		self.ROI.setColor(col)
		self.updateColorButton()
		
		self.drawPolygon()
		
		return col	

	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		
	def moveLeft(self):
		
		"""Moves corner 1 px to the left."""
		
		if len(self.artists)>0 and self.highlighted!=None:
			
			self.ROI.moveCorner(self.currInd,self.currCorner[0]-1,self.currCorner[1])
			print(self.ROI.getCorners())
			
			self.updateCornerList()
			self.drawPolygon()
			
	def moveRight(self):
		
		"""Moves corner 1 px to the right."""
		
		if len(self.artists)>0 and self.highlighted!=None:
			self.ROI.moveCorner(self.currInd,self.currCorner[0]+1,self.currCorner[1])
			self.updateCornerList()
			self.drawPolygon()
			
	
	def moveUp(self):
		
		"""Moves corner 1 px up."""
		
		if len(self.artists)>0 and self.highlighted!=None:
			self.ROI.moveCorner(self.currInd,self.currCorner[0],self.currCorner[1]+1)
			self.updateCornerList()
			self.drawPolygon()
	
	def moveDown(self):
		
		"""Moves corner 1 px down."""
		
		if len(self.artists)>0 and self.highlighted!=None:
			self.ROI.moveCorner(self.currInd,self.currCorner[0],self.currCorner[1]-1)
			self.updateCornerList()
			self.drawPolygon()
			
class polygonSliceROIDialog(polygonROIDialog,sliceROIDialog):
	
	def __init__(self,ROI,parent):
		
		super(polygonSliceROIDialog,self).__init__(ROI,parent)
		
		#sliceROIDialog.__init__(self,ROI,parent)
		#polygonROIDialog.__init__(self,ROI,parent)	

class customROIDialog(ROIDialog):
	
	def __init__(self,ROI,parent):	
		super(customROIDialog,self).__init__(ROI,parent)	
		
		#Labels
		self.lblAvailable = QtGui.QLabel("Available ROIs:", self)
		self.lblIncluded = QtGui.QLabel("Included ROIs:", self)
		
		#TreeWidget
		self.includedList=QtGui.QTreeWidget()
		self.includedList.setHeaderLabels(["name","operation"])
		self.includedList.itemClicked.connect(self.includedSelected)
		
		self.availableList=QtGui.QTreeWidget()
		self.availableList.setHeaderLabels(["name"])
		self.availableList.itemClicked.connect(self.availableSelected)
			
		#Buttons
		self.btnRemoveROI=QtGui.QPushButton('Remove')
		self.btnRemoveROI.connect(self.btnRemoveROI, QtCore.SIGNAL('clicked()'), self.removeROI)
		
		self.btnAddROI=QtGui.QPushButton('Add')
		self.btnAddROI.connect(self.btnAddROI, QtCore.SIGNAL('clicked()'), self.addROI)
		
		self.btnSubsctractROI=QtGui.QPushButton('Substract')
		self.btnSubsctractROI.connect(self.btnSubsctractROI, QtCore.SIGNAL('clicked()'), self.substractROI)
		
		self.btnMoveUp=QtGui.QPushButton('Move Up')
		self.btnMoveUp.connect(self.btnMoveUp, QtCore.SIGNAL('clicked()'), self.moveROIUp)
		
		self.btnMoveDown=QtGui.QPushButton('Move Down')
		self.btnMoveDown.connect(self.btnMoveDown, QtCore.SIGNAL('clicked()'), self.moveROIDown)
		
		#Layout
		nRows=self.grid.rowCount()
		
		self.grid.addWidget(self.lblAvailable,nRows+1,2)
		self.grid.addWidget(self.availableList,nRows+2,2,4,1)
		
		self.grid.addWidget(self.btnAddROI,nRows+6,1)
		self.grid.addWidget(self.btnSubsctractROI,nRows+6,2)
		
		self.grid.addWidget(self.lblIncluded,nRows+7,2)
		self.grid.addWidget(self.includedList,nRows+8,2,4,1)
		
		self.grid.addWidget(self.btnRemoveROI,nRows+8,1)
		self.grid.addWidget(self.btnMoveUp,nRows+9,1)
		self.grid.addWidget(self.btnMoveDown,nRows+10,1)
		
		self.updateAvailableList()
		self.updateIncludedList()
		
		self.show()	

	def updateAvailableList(self):
		self.availableList.clear()
		for roi in self.ROI.embryo.ROIs:
			if roi not in self.ROI.ROIsIncluded and roi!=self:
				item=QtGui.QTreeWidgetItem(self.availableList,[roi.name])
				item.setBackgroundColor(0,self.getQColor(roi.color))
				a=self.getQColor(roi.color)
				
		return self.availableList
	
	def updateIncludedList(self):
		self.includedList.clear()
		for i,roi in enumerate(self.ROI.ROIsIncluded):
			
			if self.ROI.procedures[i]==1:
				op='+'
			elif self.ROI.procedures[i]==-1:
				op='-'
			else:
				op=' '
				printWarning('Unknown procedure ' + str(roi.procedures[i]))
				
			item=QtGui.QTreeWidgetItem(self.includedList,[roi.name,op])
			item.setBackgroundColor(0,self.getQColor(roi.color))
			item.setBackgroundColor(1,self.getQColor(roi.color))
			
		self.updateAxes()		
				
		return self.includedList
	
	def availableSelected(self):
		ind,node=self.getCurrentAvailableNode()
		self.currAvailable=self.ROI.embryo.ROIs[ind]
		return self.currAvailable
	
	def includedSelected(self):
		ind,node=self.getCurrentIncludedNode()
		self.currIncluded=self.ROI.embryo.ROIs[ind]
		return self.currIncluded
	
	def getQColor(self,color):
		colorRGB=np.asarray(matplotlib.colors.colorConverter.to_rgb(color))
		
		colorRGB=colorRGB*255
		
		return QtGui.QColor(int(colorRGB[0]),int(colorRGB[1]),int(colorRGB[2]))
	
	def addROI(self):
		self.ROI.addROI(self.currAvailable,1)
		self.updateAvailableList()
		self.updateIncludedList()
		
	def substractROI(self):
		self.ROI.addROI(self.currAvailable,-1)
		self.updateAvailableList()
		self.updateIncludedList()
		
	def updateAxes(self):
		#Will need something that removes all children of a certain type
		self.ROI.showBoundary(ax=self.ax,color='each')
		self.canvas.draw()
	
	def getCurrentIncludedNode(self):
		self.currNode=self.includedList.currentItem()
		if self.currNode==None:
			return None,None
		
		ind=self.includedList.indexOfTopLevelItem(self.currNode)
		return ind,self.currNode
	
	def getCurrentAvailableNode(self):
		self.currNode=self.availableList.currentItem()
		if self.currNode==None:
			return None,None
		
		ind=self.availableList.indexOfTopLevelItem(self.currNode)
		return ind,self.currNode
		
	def moveROIUp(self):
		ind,node=self.getCurrentIncludedNode()
		if ind==0 or ind==None:
			return 
		
		self.ROI.ROIsIncluded[ind-1] , self.ROI.ROIsIncluded[ind]  = self.ROI.ROIsIncluded[ind], self.ROI.ROIsIncluded[ind-1]
		self.ROI.procedures[ind-1] , self.ROI.procedures[ind]  = self.ROI.procedures[ind], self.ROI.procedures[ind-1]
	
		self.updateIncludedList()
		
	def moveROIDown(self):
		
		ind,node=self.getCurrentIncludedNode()
		if ind==len(self.ROI.ROIsIncluded)-1 or ind==None:
			return 
		
		self.ROI.ROIsIncluded[ind+1] , self.ROI.ROIsIncluded[ind]  = self.ROI.ROIsIncluded[ind], self.ROI.ROIsIncluded[ind+1]
		self.ROI.procedures[ind+1] , self.ROI.procedures[ind]  = self.ROI.procedures[ind], self.ROI.procedures[ind+1]
		self.updateIncludedList()
		
	def removeROI(self):
		ind,node=self.getCurrentIncludedNode()
		if ind==None:
			return 
		
		self.ROI.ROIsIncluded.pop(ind)
		self.ROI.procedures.pop(ind)
		
		self.updateIncludedList()
		self.updateAvailableList()
		
#===================================================================================================================================
#Default ROI Dialog
#===================================================================================================================================		
		
class defaultROIsDialog(pyfrp_gui_basics.basicCanvasDialog):
	
	def __init__(self,embryo,parent):
		
		super(defaultROIsDialog,self).__init__(parent)
		
		self.embryo=embryo
		self.radius=300
		
		#Labels
		self.lblRadius = QtGui.QLabel("Radius of ROIs:", self)
		self.lblCenter = QtGui.QLabel("Center of ROIs:", self)
		
		#LineEdits
		self.qleRadius = QtGui.QLineEdit(str(self.radius))
		
		self.initCenterQle()
		
		self.doubleValid=QtGui.QDoubleValidator()
		
		self.qleRadius.setValidator(self.doubleValid)
		self.qleCenterX.setValidator(self.doubleValid)
		self.qleCenterY.setValidator(self.doubleValid)
		
		self.qleCenterX.editingFinished.connect(self.setCenter)
		self.qleCenterY.editingFinished.connect(self.setCenter)
		self.qleRadius.editingFinished.connect(self.setRadius)
		
		#Layout
		nRows=self.grid.rowCount()
		self.grid.addWidget(self.lblRadius,nRows+1,1)
		self.grid.addWidget(self.lblCenter,nRows+2,1)
		
		self.centerGrid = QtGui.QGridLayout()
		
		self.centerGrid.addWidget(self.qleCenterX)
		self.centerGrid.addWidget(self.qleCenterY)
		
		self.grid.addWidget(self.qleRadius,nRows+1,2)
		self.grid.addLayout(self.centerGrid,nRows+2,2)
		
		
		self.showFirstDataImg()
		
		#Connect Canvas
		self.connectCanvas()
		
		self.setWindowTitle('Create default ROIs')
		
		self.show()
			
	def initCenterQle(self):
		
		if self.embryo.geometry==None:
			self.center=[256,256]	
		else:
			self.center=list(self.embryo.geometry.center)
		self.qleCenterX = QtGui.QLineEdit(str(self.center[0]))
		self.qleCenterY = QtGui.QLineEdit(str(self.center[1]))
	
	def setCenter(self):
		self.center=[float(self.qleCenterX.text()),float(self.qleCenterY.text())]
		self.drawCenter()
		self.drawRadius()
		return self.center
	
	def setRadius(self):
		self.radius=float(self.qleRadius.text())
		self.drawRadius()
		return self.radius
	
	def computeRadiusFromCoordinate(self,x,y):
		return np.sqrt((x-self.center[0])**2+(y-self.center[1])**2)
	
	def drawCenter(self):
		pt=self.drawPoint(self.center[0],self.center[1],color='r')
		return pt
	
	def drawRadius(self):
		c=self.drawCircle(self.center,self.radius,color='r')
		return c
	
	def drawPoint(self,x,y,color='r',idx=0):
		pt=ptc.Circle([x,y],radius=3,fill=True,color=color)
		self.replaceArtist(idx,self.ax.add_patch(pt))
		self.canvas.draw()
		return pt
	
	def drawCircle(self,center,radius,color='r',idx=1):
		c=ptc.Circle(center,radius=radius,fill=False,color=color)
		self.replaceArtist(idx,self.ax.add_patch(c))
		self.canvas.draw()
		return c
	
	def showFirstDataImg(self):
		
		self.embryo.updateFileList()
		
		fnImg=self.embryo.fnDatafolder
		if fnImg[-1]!='/':
			fnImg=fnImg+'/'
		fnImg=fnImg+self.embryo.fileList[0]
		
		img=pyfrp_img_module.loadImg(fnImg,self.embryo.dataEnc)
	
		self.showImg(img)
	
	def getMouseCanvas(self,event):
		
		#Left click to define center and then radius
		if event.button==1:
			
			#Check if clicked withing axes
			if event.xdata==None:
				return
			
			#Check if this is center 
			if len(self.artists)==0:
				
				#Update qles
				self.qleCenterX.setText(str(event.xdata))
				self.qleCenterY.setText(str(event.ydata))
				
				#Set center (and draw)
				self.setCenter()
					
			#Check if this is radius
			elif len(self.artists)>1:
				
				#Update radius
				radius=self.computeRadiusFromCoordinate(event.xdata,event.ydata)
				
				#Update qle
				self.qleRadius.setText(str(radius))
				
				#Set radius (and draw)
				self.setRadius()
			
		#Remove last artist
		elif event.button==3:
			self.removeArtist()
	
	def updateCenterQles(self,center):
		
		"""Updates the two center Qles with new center."""
		
		self.qleCenterX.setText(str(center[0]))
		self.qleCenterY.setText(str(center[1]))
	
	def updateRadiusQle(self,radius):
		
		"""Updates radius Qle with given radius."""
		
		self.qleRadius.setText(str(radius))
	
	def keyPressed(self,event):
		
		"""Directs all key press events to the respective functions."""
		
		if event.key=='left':
			self.moveLeft()
		elif event.key=='right':
			self.moveRight()
		elif event.key=='up':
			self.moveUp()
		elif event.key=='down':
			self.moveDown()
		elif event.key=='ctrl+up':
			self.increaseRadius()
		elif event.key=='ctrl+down':
			self.decreaseRadius()
			
	def moveLeft(self):
		
		"""Moves center 1 px to the left."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.center[0]-1,self.center[1]])
			self.setCenter()
		
		
	def moveRight(self):
		
		"""Moves center 1 px to the right."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.center[0]+1,self.center[1]])
			self.setCenter()
	
	def moveUp(self):
		
		"""Moves center 1 px up."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.center[0],self.center[1]+1])
			self.setCenter()
	
	def moveDown(self):
		
		"""Moves center 1 px down."""
		
		if len(self.artists)>0:
			self.updateCenterQles([self.center[0],self.center[1]-1])
			self.setCenter()
			
	def increaseRadius(self):
		
		"""Increases radius by 1 px."""
		
		if len(self.artists)>1:
			self.updateRadiusQle(self.radius+1)
			self.setRadius()
			
	def decreaseRadius(self):
		
		"""Decreases radius by 1 px."""
		
		if len(self.artists)>1:
			self.updateRadiusQle(self.radius-1)
			self.setRadius()
	
	def donePressed(self):
		self.embryo.genDefaultROIs(self.center,self.radius)
		self.done(1)
		return
	
#===================================================================================================================================
#Default ROI Dialog
#===================================================================================================================================		
		
class ROISelector(QtGui.QDialog):
	
	"""Simple dialog allowing to select between different ROI types.
	
	""" 
	
	def __init__(self,embryo,initType,parent,**kwargs):
		super(ROISelector,self).__init__(parent)
		
		#Parsing in args
		self.parseInInput(kwargs,embryo)
		
		self.ROI=None	
		
		#Labels
		self.lblMsg=QtGui.QLabel(self.msg, self)
		
		#Buttons	
		self.btnDefault=QtGui.QPushButton('Use Default')
		self.btnDefault.connect(self.btnDefault, QtCore.SIGNAL('clicked()'), self.defaultPressed)
		
		self.btnDone=QtGui.QPushButton('Create ROI')
		self.btnDone.connect(self.btnDone, QtCore.SIGNAL('clicked()'), self.donePressed)
		
		#Combo box
		self.comboType = QtGui.QComboBox(self)
		self.comboType.addItem("radial")
		self.comboType.addItem("slice")
		self.comboType.addItem("radialSlice")
		self.comboType.addItem("square")
		self.comboType.addItem("squareSlice")
		self.comboType.addItem("rectangle")
		self.comboType.addItem("rectangleSlice")
		self.comboType.addItem("polygon")
		self.comboType.addItem("polygonSlice")
		self.comboType.addItem("custom")
		
		self.initComboType(initType)
		
		#Layout
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.lblMsg)
		self.vbox.addWidget(self.comboType)
		self.vbox.addWidget(self.btnDefault)
		self.vbox.addWidget(self.btnDone)
		
		self.setLayout(self.vbox)   
		
		self.setWindowTitle(self.title)    
		
		self.show()
		
	def initComboType(self,initType):
		
		"""Sets combo box for type selection to entry with label ``initType``.
		"""
	
		idx=self.comboType.findText(initType,QtCore.Qt.MatchExactly)
		self.comboType.setCurrentIndex(idx)
	
	
	def parseInInput(self,kwargs,embryo):
		
		"""Parses in input from dialog and saves it into dialog attributes.
		
		If attributes get not parsed in, will use default values instead.

		"""
		
		#Parse embryo
		self.embryo=embryo
		
		#Grab ROI name
		if "name" in list(kwargs.keys()):
			self.name=str(kwargs["name"])
		else:
			self.name="newROI"
		
		#Grab center
		if "center" in list(kwargs.keys()):
			self.center=kwargs["center"]
		else:
			if self.embryo.geometry!=None:
				self.center=self.embryo.geometry.getCenter()
			else:
				self.center=[self.embryo.dataResPx/2.]*2
		
		#Grab sidelength
		if "sidelength" in list(kwargs.keys()):
			self.sidelength=kwargs["sidelength"]
		else:
			self.sidelength=100
		
		#Grab offset
		if "offset" in list(kwargs.keys()):
			self.offset=kwargs["offset"]
		else:
			self.offset=np.asarray(self.center)-self.sidelength/2.
			
		#Grab radius
		if "radius" in list(kwargs.keys()):
			self.radius=kwargs["radius"]
		else:
			self.radius=200.
		
		#Grab color
		if "color" in list(kwargs.keys()):
			self.color=kwargs["color"]
		else:
			self.color='r'
	
		#Grab slice height
		if "sliceHeight" in list(kwargs.keys()):
			self.sliceHeight=kwargs["sliceHeight"]
		else:
			self.sliceHeight=self.embryo.sliceHeightPx
		
		#Grab slice width
		if "sliceWidth" in list(kwargs.keys()):
			self.sliceWidth=kwargs["sliceWidth"]
		else:
			self.sliceWidth=self.embryo.sliceWidthPx
		
		#Grab title
		if "title" in list(kwargs.keys()):
			self.title=kwargs["title"]
		else:
			self.title="ROI Selector"
		
		#Grab Message
		if "msg" in list(kwargs.keys()):
			self.msg=kwargs["msg"]
		else:
			self.msg="Select ROI type"
			
		#Grab asMaster
		if "asMaster" in list(kwargs.keys()):
			self.asMaster=kwargs["asMaster"]
		else:
			self.asMaster=False
			
				
	def newROI(self):
		newID=self.embryo.getFreeROIId()
		typ=self.comboType.currentText()
		
		if typ=='slice':
			self.embryo.newSliceROI(self.name,newID,self.sliceHeight,self.sliceWidth,False,color=self.color,asMaster=self.asMaster)	
		elif typ=='radial':
			self.embryo.newRadialROI(self.name,newID,self.center,self.radius,color=self.color,asMaster=self.asMaster)
		elif typ=='radialSlice':
			self.embryo.newRadialSliceROI(self.name,newID,self.center,self.radius,self.sliceHeight,self.sliceWidth,False,color=self.color,asMaster=self.asMaster)
		elif typ=='square':
			self.embryo.newSquareROI(self.name,newID,self.offset,self.sidelength,self.radius,color=self.color,asMaster=self.asMaster)
		elif typ=='squareSlice':
			self.embryo.newSquareSliceROI(self.name,newID,self.offset,self.sidelength,self.sliceHeight,self.sliceWidth,False,color=self.color,asMaster=self.asMaster)
		elif typ=='rectangle':
			self.embryo.newRectangleROI(self.name,newID,self.offset,self.sidelength,self.sidelength,color=self.color)
		elif typ=='rectangleSlice':
			self.embryo.newRectangleSliceROI(self.name,newID,self.offset,self.sidelength,self.sidelength,self.sliceHeight,self.sliceWidth,False,color=self.color,asMaster=self.asMaster)
		elif typ=='polygon':
			self.embryo.newPolyROI(self.name,newID,[],color=self.color,asMaster=self.asMaster)
		elif typ=='polygonSlice':
			self.embryo.newPolySliceROI(self.name,newID,[],self.sliceHeight,self.sliceWidth,False,color=self.color,asMaster=self.asMaster)	
		elif typ=='custom':
			self.embryo.newCustomROI(self.name,newID,color=self.color,asMaster=self.asMaster)
		else:
			printWarning('Unknown ROI Type' + typ)
	
		self.openDialog(self.embryo.ROIs[-1],typ)

		return self.embryo.ROIs[-1]
		
	def openDialog(self,ROI,typ):
		if typ=='slice':
			ret=sliceROIDialog(ROI,self).exec_()
		elif typ=='radial':
			ret=radialROIDialog(ROI,self).exec_()
		elif typ=='radialSlice':
			ret=radialSliceROIDialog(ROI,self).exec_()
		elif typ=='square':
			ret=squareROIDialog(ROI,self).exec_()
		elif typ=='squareSlice':
			ret=squareSliceROIDialog(ROI,self).exec_()
		elif typ=='rectangle':
			ret=rectangleROIDialog(ROI,self).exec_()
		elif typ=='rectangleSlice':
			ret=rectangleSliceROIDialog(ROI,self).exec_()
		elif typ=='polygon':
			ret=polygonROIDialog(ROI,self).exec_()
		elif typ=='polygonSlice':
			ret=polygonSliceROIDialog(ROI,self).exec_()
		elif typ=='custom':
			ret=customROIDialog(ROI,self).exec_()
		else:
			printWarning('Unknown ROI Type' + typ)	
		
	def getROI(self):
		return self.ROI
	
	def defaultPressed(self):
		self.ROI=None
		self.done(1)
		return
	
	def donePressed(self):
		self.ROI=self.newROI()
		self.done(1)
		return
	
#===================================================================================================================================
#Dialogs for indexing progress
#===================================================================================================================================

class indexProgressDialog(pyfrp_gui_basics.progressDialog):
	
	def __init__(self,parent):
		super(indexProgressDialog,self).__init__(parent)
		
		#Labels
		self.lblName.setText("ROI index computation in progress...")
		
		#Window title
		self.setWindowTitle('Indexing progress')
		    
		self.show()	

class indexThread(pyfrp_gui_basics.pyfrpThread):
	
	def __init__(self, embryo=None, parent=None):
		
		super(indexThread,self).__init__(parent)
		self.obj=embryo
		self.embryo=embryo
		
		#self.prog_signal.connect(self.print_prog)
			
	def runTask(self,debug=False):
		self.embryo.computeROIIdxs(signal=self.progressSignal,debug=True)

#===================================================================================================================================
#Dialogs for wizard selection
#===================================================================================================================================
	
class wizardSelector(QtGui.QDialog):
	
	"""Dialog to select if either use default ROIs, Wizard, or ROI manager."""
	
	def __init__(self,parent):
		
		super(wizardSelector,self).__init__(parent)
		
		#Passing parent GUI
		self.parent=parent
		
		self.mode=None
		
		#Buttons
		self.btnUseDefault=QtGui.QPushButton('Create default ROIs')
		self.btnUseDefault.connect(self.btnUseDefault, QtCore.SIGNAL('clicked()'), self.setUseDefault)
		
		self.btnUseWizard=QtGui.QPushButton('Use ROI Wizard')
		self.btnUseWizard.connect(self.btnUseWizard, QtCore.SIGNAL('clicked()'), self.setUseWizard)
		
		self.btnUseManager=QtGui.QPushButton('Use ROI Manager')
		self.btnUseManager.connect(self.btnUseManager, QtCore.SIGNAL('clicked()'), self.setUseManager)
			
		self.vbox = QtGui.QVBoxLayout()
		self.vbox.addWidget(self.btnUseDefault)
		self.vbox.addWidget(self.btnUseWizard)
		self.vbox.addWidget(self.btnUseManager)
		
	
		self.setLayout(self.vbox)    
			
		self.setWindowTitle('Select Wizard Steps')    
		self.show()
		
	def setUseDefault(self):
		self.mode=0
		self.done(1)
		
	def setUseWizard(self):
		self.mode=1
		self.done(1)
	
	def setUseManager(self):
		self.mode=2
		self.done(1)
	
	def getMode(self):
		return self.mode