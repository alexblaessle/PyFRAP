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

#Module containing basic custom PyQT classes:
#1) basicCanvasDialog



#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Misc
import sys
import time
import os, os.path

#PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *

#PyQT
from PyQt5 import QtGui, QtCore, QtWidgets

#matplotlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

"""
Apparently the NavigationToolbar naming has changed in newer matplotlib versions, thus
we need to test out some cases.
"""

try:
	from matplotlib.backends.backend_qt5agg import NavigationToolbar2QTAgg as  NavigationToolbar
except ImportError:
	try:
		from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as  NavigationToolbar
	except ImportError:
		printWarning("Cannot import NavigationToolbar.")
		
from matplotlib.figure import Figure
	
#===================================================================================================================================
#Basic Dialog with space for QLEs/Btns on left hand side and canvas on right hand side
#===================================================================================================================================

class basicCanvasDialog(QtWidgets.QDialog):
	
	def __init__(self,parent,xlim=[0,512],ylim=[0,512]):
		
		super(basicCanvasDialog,self).__init__(parent)
				
		self.dpi = 100
		self.setMinimumSize(1000,500) 
		self.resize(1300,500)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Bookkeeping variables
		#-------------------------------------------------------------------------------------------------------------------
		
		self.artists=[]
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		#Done button
		self.btnDone=QtWidgets.QPushButton('Done')
		self.btnDone.clicked.connect(self.donePressed)
	
		#-------------------------------------------------------------------------------------------------------------------
		#Plot frame
		#-------------------------------------------------------------------------------------------------------------------
		
		self.plotFrame = QtWidgets.QWidget()
		self.plotFrame.setMaximumWidth(1)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.grid = QtWidgets.QGridLayout()		
		self.grid.setColumnMinimumWidth(2,200) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Create Canvas
		#-------------------------------------------------------------------------------------------------------------------
		
		self.createCanvas(xlim=xlim,ylim=ylim)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.canvas)
		self.vbox.addWidget(self.btnDone)
		
		#Add everything to Horizontal Box
		self.hbox = QtWidgets.QHBoxLayout()
		self.hbox.addLayout(self.grid)
		self.hbox.addLayout(self.vbox)
		
		self.setLayout(self.hbox)    
			
		self.setWindowTitle('basicCanvasDialog')    
		self.show()
	
	def createCanvas(self,xlim=None,ylim=None):
			
		h=500/self.dpi
		v=500/self.dpi
		
		self.fig = Figure( dpi=self.dpi)
		self.fig.set_size_inches(h,v,forward=True)
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.plotFrame)
		
		self.ax = self.fig.add_subplot(111)
		
		if xlim!=None:
			self.ax.set_xlim(xlim)
		if ylim!=None:
			self.ax.set_ylim(ylim)
		
		self.canvas.draw()
		#self.plotFrame.adjustSize()
		
		return 
		
	def showImg(self,img):
		
		self.ax.imshow(img)
		self.ax.set_xlim([1,img.shape[0]])
		self.ax.set_ylim([1,img.shape[1]])
		
		self.canvas.draw()
		
		return self.canvas
	
	def connectCanvas(self):
		self.canvas.mpl_connect('button_press_event', self.getMouseCanvas)
		self.canvas.mpl_connect('key_press_event', self.keyPressed)
		
		self.canvas.setFocusPolicy( QtCore.Qt.ClickFocus )
		self.canvas.setFocus()
		
		
	def keyPressed(self):
		printWarning("No Key-Press Action defined.")
	
	def removeArtist(self,idx=-1):
		if len(self.artists)>0:
			self.artists[idx].remove()
			self.artists.pop(idx)
			self.canvas.draw()
		return self.artists
	
	def replaceArtist(self,idx,newArtist):
		if len(self.artists)>idx:
			self.removeArtist(idx=idx)
			self.artists.insert(idx,newArtist)
			self.canvas.draw()
		else:
			self.artists.append(newArtist)
		
		return self.artists
	
	def setArtistColor(self,color):
		for artist in self.artists:
			
			artist.set_color(color)
		self.canvas.draw()	
		return self.artists
	
	def donePressed(self):
		self.done(1)
		return 
	
#===================================================================================================================================
#Basic Dialog for settings of any kind.
#===================================================================================================================================

class basicSettingsDialog(QtWidgets.QDialog):
	
	def __init__(self,parent):
		
		super(basicSettingsDialog,self).__init__(parent)
				
		self.setMinimumSize(500,500) 
		self.resize(700,500)
		
		#-------------------------------------------------------------------------------------------------------------------
		#Buttons
		#-------------------------------------------------------------------------------------------------------------------
		
		#Done button
		self.btnDone=QtWidgets.QPushButton('Done')
		self.btnDone.clicked.connect(self.donePressed)
	
		#-------------------------------------------------------------------------------------------------------------------
		#Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.grid = QtWidgets.QGridLayout()		
		self.grid.setColumnMinimumWidth(2,20) 
		
		#-------------------------------------------------------------------------------------------------------------------
		#Validators
		#-------------------------------------------------------------------------------------------------------------------
		
		self.doubleValid=QtGui.QDoubleValidator()
		self.intValid=QtGui.QIntValidator()
		
		#-------------------------------------------------------------------------------------------------------------------
		#Final Layout
		#-------------------------------------------------------------------------------------------------------------------
		
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addLayout(self.grid)
		
		self.hbox = QtWidgets.QHBoxLayout()
		self.hbox.addWidget(self.btnDone,stretch=0,alignment=QtCore.Qt.AlignRight)	
		self.vbox.addLayout(self.hbox)
			
		self.setLayout(self.vbox)    
			
		self.setWindowTitle('basicSettingsDialog')    

	def donePressed(self):
		self.done(1)
		return

#===================================================================================================================================
#Basic Selector for a single item out of a list
#===================================================================================================================================

class basicSelectorDialog(QtWidgets.QDialog):
	
	def __init__(self,List,parent):
		
		super(basicSelectorDialog,self).__init__(parent)
		
		self.item=None
		self.List=List
		
		#Done button
		self.btnDone=QtWidgets.QPushButton('Done')
		self.btnDone.clicked.connect(self.donePressed)
			
		#QTreeWidget
		self.ListWidget=QtWidgets.QTreeWidget()
		self.ListWidget.setColumnWidth(0,100)
		self.ListWidget.itemClicked.connect(self.itemClicked)
		self.updateList()
		
		#Layout
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.ListWidget)
		
		self.hbox = QtWidgets.QHBoxLayout()
		self.hbox.addWidget(self.btnDone,stretch=0,alignment=QtCore.Qt.AlignRight)	
		self.vbox.addLayout(self.hbox)
			
		self.setLayout(self.vbox)    
		
		self.setWindowTitle('basicSelectorDialog')    
		self.show()
		
	def updateList(self):
		
		self.ListWidget.clear()
		for r in self.List:
			QtWidgets.QTreeWidgetItem(self.ListWidget,[r])
		return 
	
	def itemClicked(self):
		idx=self.ListWidget.indexFromItem(self.ListWidget.currentItem()).row()
		self.item=self.List[idx]
		
	def getItem(self):
		return self.item
			
	def donePressed(self):
		self.done(1)
		return		

#===================================================================================================================================
#Basic Selector for a list of items out of a list
#===================================================================================================================================

class listSelectorDialog(QtWidgets.QDialog):
	
	def __init__(self,parent,List,leftTitle="",rightTitle="",itemsRight=[]):
		super(listSelectorDialog,self).__init__(parent)
		#print type(self), type(parent)
		
		#QtWidgets.QDialog.__init__()
		
		
		self.itemsRight=itemsRight
		self.itemsLeft=list(List)
		
		self.List=List
		
		#Buttons
		self.btnAdd=QtWidgets.QToolButton()
		self.btnAdd.clicked.connect(self.addItem)
		self.btnAdd.setArrowType(QtCore.Qt.RightArrow)
		
		self.btnRemove=QtWidgets.QToolButton()
		self.btnRemove.clicked.connect(self.removeItem)
		self.btnRemove.setArrowType(QtCore.Qt.LeftArrow)
		
		self.btnDone=QtWidgets.QPushButton('Done')
		self.btnDone.clicked.connect(self.donePressed)
		
		#Left QtreeWidgetItem
		self.leftList=QtWidgets.QTreeWidget()
		self.leftList.setHeaderLabels([leftTitle])
		self.leftList.setColumnWidth(0,200)
		self.leftList.setColumnWidth(1,75)
		self.leftList.itemDoubleClicked.connect(self.addItem)
		
		#right QtreeWidgetItem
		self.rightList=QtWidgets.QTreeWidget()
		self.rightList.setHeaderLabels([rightTitle])
		self.rightList.setColumnWidth(0,200)
		self.rightList.setColumnWidth(1,75)
		self.rightList.itemDoubleClicked.connect(self.removeItem)
		
		#Layout
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.btnAdd)
		self.vbox.addWidget(self.btnRemove)
		
		self.hbox = QtWidgets.QHBoxLayout()
		self.hbox.addWidget(self.leftList)
		self.hbox.addLayout(self.vbox)
		self.hbox.addWidget(self.rightList)
		
		self.vbox2 = QtWidgets.QVBoxLayout()
		self.vbox2.addLayout(self.hbox)
		self.vbox2.addWidget(self.btnDone)
		
		#Init lists
		self.initLeftList()
		self.initRightList()
		
		self.resize(400,500)
		self.setLayout(self.vbox2)
		self.setWindowTitle("list Selector Dialog")
		
		self.show()
	
	def getListDifference(self):
		
		for item in self.itemsLeft:
			if item in self.itemsRight:
				self.itemsLeft.remove(item)
		
		return self.itemsLeft
				
	
	def initLeftList(self):
		
		self.getListDifference()
		
		for item in self.itemsLeft:
			QtWidgets.QTreeWidgetItem(self.leftList,[item])
			
	def initRightList(self):
		
		for item in self.itemsRight:
			QtWidgets.QTreeWidgetItem(self.rightList,[item])
				
	def addItem(self):

		#Determine selected item
		self.currentItem=str(self.leftList.currentItem().data(0,0).toString())
		
		#Insert new node in right list
		newNode=QtWidgets.QTreeWidgetItem(self.rightList,[self.currentItem])
		
		#Remove node in left list
		self.currLeftInd=self.leftList.indexFromItem(self.leftList.currentItem()).row()
		self.leftList.takeTopLevelItem(self.currLeftInd)
		
		self.itemsRight.append(self.currentItem)
		self.itemsLeft.remove(self.currentItem)
		
		
	def removeItem(self):
	
		#Determine selected item
		self.currentItem=str(self.rightList.currentItem().data(0,0).toString())
		
		#Insert new node in left list
		newNode=QtWidgets.QTreeWidgetItem(self.leftList,[self.currentItem])
		
		#Remove node in right list
		self.currRightInd=self.rightList.indexFromItem(self.rightList.currentItem()).row()
		self.rightList.takeTopLevelItem(self.currRightInd)
		
		self.itemsRight.remove(self.currentItem)
		self.itemsLeft.append(self.currentItem)

	def getSelection(self):
		return self.itemsRight
	
	def donePressed(self):
		self.done(1)
		return


#===================================================================================================================================
#Basic Selector for a list of items out of a list
#===================================================================================================================================

class advancedListSelectorDialog(listSelectorDialog):
	
	def __init__(self,parent,List,leftTitle="",rightTitle="",itemsRight=[]):
		
		super(advancedListSelectorDialog,self).__init__(parent,List,leftTitle=leftTitle,rightTitle=rightTitle)
		#print type(self),type(parent)
		#raw_input()
		
		#listSelectorDialog.__init__(parent,List,leftTitle=leftTitle,rightTitle=rightTitle,itemsRight=itemsRight)
		
		self.btnUp=QtWidgets.QToolButton()
		self.btnUp.clicked.connect(self.upItem)
		self.btnUp.setArrowType(QtCore.Qt.UpArrow)
		
		self.btnDown=QtWidgets.QToolButton()
		self.btnDown.clicked.connect(self.downItem)
		self.btnDown.setArrowType(QtCore.Qt.DownArrow)
		
		self.vbox.addWidget(self.btnUp)
		self.vbox.addWidget(self.btnDown)
		
		self.setWindowTitle("advanced List Selector Dialog")
		
	def upItem(self):
		
		#Determine selected item
		self.currentItem=str(self.rightList.currentItem().data(0,0).toString())
		
		#Determine index in sel prop list
		ind=self.itemsRight.index(self.currentItem)
		
		#Swap in list
		if ind>0:
			self.itemsRight[ind-1], self.itemsRight[ind] = self.itemsRight[ind], self.itemsRight[ind-1]
		
		#Clear list and recreate it
		self.rightList.clear()
		self.initRightList()
		self.currentRightItem=self.rightList.topLevelItem(ind-1)
		self.rightList.setCurrentItem(self.currentRightItem)
		
	def downItem(self):
		
		#Determine selected item
		self.currentItem=str(self.rightList.currentItem().data(0,0).toString())
		
		#Determine index in sel prop list
		ind=self.itemsRight.index(self.currentItem)
		
		#Swap in list
		if ind<len(self.itemsRight)-1:
			self.itemsRight[ind+1], self.itemsRight[ind] = self.itemsRight[ind], self.itemsRight[ind+1]
		
		#Clear list and recreate it
		self.rightList.clear()
		self.initRightList()
		self.currentRightItem=self.rightList.topLevelItem(ind+1)
		self.rightList.setCurrentItem(self.currentRightItem)	
		
#===================================================================================================================================
#Basic Progress Dialog
#===================================================================================================================================
		
class progressDialog(QtWidgets.QDialog):
	
	def __init__(self,parent):
		super(progressDialog,self).__init__(parent)
		
		#Labels
		self.lblName = QtWidgets.QLabel("Something in progress...", self)
		
		#Buttons
		self.btnCancel=QtWidgets.QPushButton('Cancel')
		self.btnCancel.clicked.connect(self.cancel)	
		
		#ProgressBar
		self.progressbar = QtWidgets.QProgressBar()
		self.progressbar.setMinimum(1)
		self.progressbar.setMaximum(100)
		
		#Layout
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.lblName)
		self.vbox.addWidget(self.progressbar)
		self.vbox.addWidget(self.btnCancel)
		
		self.setLayout(self.vbox)
		self.setWindowTitle('Progress Dialog')
		
		self.show()	
	
	def cancel(self):
	
		self.accepted.emit()

#===================================================================================================================================
#Basic PyFRAP Thread
#===================================================================================================================================

#Simple worker class        
class pyfrpWorker(QtCore.QObject):
	
	taskFinished = QtCore.pyqtSignal()
	start = QtCore.pyqtSignal()
	
	def __init__(self, function, *args, **kwargs):
		super(pyfrpWorker, self).__init__()
		
		self.function = function
		self.args = args
		self.kwargs = kwargs
		self.start.connect(self.run)
		
	#@QtCore.pyqtSlot() (Interesting: the pyqtslot decorator seems to block the start signal...)
	def run(self):
		self.function(*self.args, **self.kwargs)
		self.taskFinished.emit()
		
		
class pyfrpThread(QtCore.QThread):
	
	progressSignal = QtCore.pyqtSignal(int)
	
	def __init__(self, parent=None):
		QtCore.QThread.__init__(self)
		
	def __del__(self):
		self.wait()
    
		
#===================================================================================================================================
#Basic Wait Dialog
#===================================================================================================================================
		
class waitDialog(QtWidgets.QDialog):
	
	def __init__(self,parent):
		super(waitDialog,self).__init__(parent)
		
		#Labels
		self.lblName = QtWidgets.QLabel("Something in progress...", self)
		
		#Buttons
		self.btnCancel=QtWidgets.QPushButton('Cancel')
		self.btnCancel.clicked.connect(self.cancel)	
		
		#Layout
		self.vbox = QtWidgets.QVBoxLayout()
		self.vbox.addWidget(self.lblName)
		self.vbox.addWidget(self.btnCancel)
		
		self.setLayout(self.vbox)
		self.setWindowTitle('Progress Dialog')
		
		self.show()	
	
	def cancel(self):
		self.accepted.emit()
	
	
