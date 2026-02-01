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

#PyQT Dialogs for editing PyFRAP settings
#(1) moleculeDialog

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#QT
from PyQt5 import QtGui, QtCore, QtWidgets

#===================================================================================================================================
#Dialog for select/edit molecule
#===================================================================================================================================
	
class pathDialog(QtWidgets.QDialog):
	
	"""Dialog to modify path settings.
	
	Gives user the ability to easily set a path for the paths file.
	
	"""
	
	def __init__(self,identifier,path, parent):
	
		QtWidgets.QDialog.__init__(self, parent)

		self.parent=parent
		self.identifier=identifier
		self.path=path
		
		#Labels
		self.lblIdentifier=QtWidgets.QLabel("Identifier:", self)
		self.lblPath=QtWidgets.QLabel("Path:", self)
		self.lblPathVal=QtWidgets.QLabel(self.path, self)
		
		#QLEs
		self.qleIdentifier = QtWidgets.QLineEdit(str(self.identifier))
		self.qleIdentifier.editingFinished.connect(self.setIdentifier)
		
		#Button
		self.btnPath=QtWidgets.QPushButton('Change')
		self.btnPath.clicked.connect(self.setPath)
		
		self.btnDone=QtWidgets.QPushButton('Done')
		self.btnDone.clicked.connect(self.donePressed)
		
		#Layout
		self.grid = QtWidgets.QGridLayout()		
		self.grid.setColumnMinimumWidth(2,200) 
		
		self.grid.addWidget(self.lblIdentifier,1,1)
		self.grid.addWidget(self.lblPath,2,1)
		self.grid.addWidget(self.lblPathVal,2,2)
		self.grid.addWidget(self.qleIdentifier,1,2)
		self.grid.addWidget(self.btnPath,2,3)
		self.grid.addWidget(self.btnDone,3,3)
		
		self.setLayout(self.grid)    
			
		self.setWindowTitle('Path Dialog')   
		
		self.show()
			
	def setIdentifier(self):
		self.identifier=str(self.qleIdentifier.text())
	
	def updateLblPathVal(self,n=50):	
		self.lblPathVal.setText("..."+self.path[-n:])		
		
	def getPath(self):
		return self.identifier,self.path
	
	def setPath(self):
		
		fn = str(QtWidgets.QFileDialog.getOpenFileName(self, 'Open file',self.parent.lastopen,))
		if fn=='':
			return
		
		self.path=fn
		self.updateLblPathVal()
		
	def donePressed(self):
		self.done(1)
		return			
		
