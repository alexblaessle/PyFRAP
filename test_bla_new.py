self.btnDone=QtWidgets.QPushButton('Done')
self.btnDone.clicked.connect(self.donePressed)

#Labels
self.lblName = QtWidgets.QLabel("Name:", self)

#LineEdit
self.qleName = QtWidgets.QLineEdit(self.molecule.name)
self.qleName.editingFinished.connect(self.setName)
	
