# =====================================================================================================================================
# Copyright
# =====================================================================================================================================

# Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
# This software is distributed under the terms of the GNU General Public License.

# This file is part of PyFRAP.

# PyFRAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# ===========================================================================================================================================================================
# Module Description
# ===========================================================================================================================================================================

# PyQT Dialogs for mesh class
# (1)

# ===========================================================================================================================================================================
# Importing necessary modules
# ===========================================================================================================================================================================

# QT
from PyQt5 import QtGui, QtCore, QtWidgets

# PyFRAP GUI classes
from . import pyfrp_gui_basics

# PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_misc_module

# Numpy/Scipy
import numpy as np

# Misc
import os

# ===================================================================================================================================
# Dialog for selecting fits
# ===================================================================================================================================


class fitSelector(pyfrp_gui_basics.listSelectorDialog):
    def __init__(self, molecule, singleFit, parent):

        self.molecule = molecule
        self.embryosInRightList = []
        self.fitsInRightList = []

        self.singleFit = singleFit

        pyfrp_gui_basics.listSelectorDialog.__init__(
            self, parent, [], leftTitle="Available Fits", rightTitle="Selected Fits"
        )

        # self.initLeftList()
        self.resize(400, 500)

        self.setWindowTitle("Fit Selector")
        self.show()

    def initLeftList(self):

        """Sets up left list of selector."""

        for emb in self.molecule.embryos:

            self.currEmbryoNode = QtWidgets.QTreeWidgetItem(self.leftList, [emb.name])

            for fit in emb.fits:
                if fit.isFitted():
                    if fit not in self.molecule.selFits:
                        QtWidgets.QTreeWidgetItem(self.currEmbryoNode, [fit.name])

            self.leftList.expandItem(self.currEmbryoNode)

    def initRightList(self):

        """Sets up right list of selector."""

        for fit in self.molecule.selFits:

            if fit.embryo not in self.embryosInRightList:
                self.currEmbryoNode = QtWidgets.QTreeWidgetItem(
                    self.rightList, [fit.embryo.name]
                )
                self.embryosInRightList.append(fit.embryo)

                QtWidgets.QTreeWidgetItem(self.currEmbryoNode, [fit.name])
                self.fitsInRightList.append(fit)

    def addItem(self):

        self.getLeftSelections()

        if (
            self.leftList.currentItem() == None
            or self.leftList.currentItem().parent() == None
        ):
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return

        else:

            # Check if current embryo is already in right list
            if self.currEmbryo in self.embryosInRightList:

                # Check if fit is already in right list
                if self.currFit in self.fitsInRightList:
                    pass
                else:

                    self.getRightInd(embryoName=self.currEmbryo.name)

                    if self.singleFit:
                        # Check if target embr already has a fit, if so, remove it
                        if self.currTargetEmbryoNode.childCount() > 0:
                            self.getRightInd(
                                fitName=self.currTargetEmbryoNode.child(0)
                                .data(0, 0)
                                .toString()
                            )
                            self.currTargetEmbryoNode.takeChild(0)
                            self.fitsInRightList.remove(self.currTargetFit)

                    # Add to list of fits on right side
                    self.fitsInRightList.append(self.currFit)
                    newNode = QtWidgets.QTreeWidgetItem(
                        self.currTargetEmbryoNode, [self.currFit.name]
                    )
                    self.rightList.expandItem(self.currTargetEmbryoNode)

                    # Remove in left list
                    # self.leftList.removeItemWidget(self.currFitNode,0)
                    ind = self.currEmbryoNode.indexOfChild(self.leftList.currentItem())
                    self.currEmbryoNode.takeChild(ind)
                    # self.currEmbryoNode.takeChild()

            else:

                # Add currEmbryo to list
                self.embryosInRightList.append(self.currEmbryo)
                newNode = QtWidgets.QTreeWidgetItem(
                    self.rightList, [self.currEmbryo.name]
                )

                # Add to list of fits on right side
                self.fitsInRightList.append(self.currFit)
                QtWidgets.QTreeWidgetItem(newNode, [self.currFit.name])
                self.rightList.expandItem(newNode)

    def removeItem(self):

        """Removes item from right list."""

        # Get selection
        self.getRightSelections()

        # Check selection
        if (
            self.rightList.currentItem() == None
            or self.rightList.currentItem().parent() == None
        ):
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return
        else:

            # Remove
            ind = self.currEmbryoNode.indexOfChild(self.rightList.currentItem())
            self.currEmbryoNode.takeChild(ind)
            self.fitsInRightList.remove(self.currFit)

            # Add to left list
            self.getLeftInd(embryoName=self.currFit.embryo.name)
            QtWidgets.QTreeWidgetItem(self.currTargetEmbryoNode, [self.currFit.name])

    def getLeftSelections(self):

        # Get current embryo selection
        for emb in self.molecule.embryos:
            if self.leftList.currentItem().parent() == None:

                # This is a embryo
                if self.leftList.currentItem().data(0, 0).toString() == emb.name:
                    self.currEmbryo = emb
                    self.currEmbryoNode = self.leftList.currentItem()
                    self.currFit = None
                    self.currFitNode = None
                    break
            else:
                # This is a fit
                if (
                    self.leftList.currentItem().parent().data(0, 0).toString()
                    == emb.name
                ):
                    self.currEmbryo = emb
                    self.currEmbryoNode = self.leftList.currentItem().parent()
                    for fit in self.currEmbryo.fits:
                        if (
                            fit.name
                            == self.leftList.currentItem().data(0, 0).toString()
                        ):
                            self.currFit = fit
                            self.currFitNode = self.leftList.currentItem()
                            break

    def getRightSelections(self):

        # Get current embryo selection
        for emb in self.molecule.embryos:
            if self.rightList.currentItem().parent() == None:

                # This is a embryo
                if self.rightList.currentItem().data(0, 0).toString() == emb.name:
                    self.currEmbryo = emb
                    self.currEmbryoNode = self.rightList.currentItem()
                    self.currFit = None
                    self.currFitNode = None
                    break
            else:
                # This is a fit
                if (
                    self.rightList.currentItem().parent().data(0, 0).toString()
                    == emb.name
                ):
                    self.currEmbryo = emb
                    self.currEmbryoNode = self.rightList.currentItem().parent()
                    for fit in self.currEmbryo.fits:
                        if (
                            fit.name
                            == self.rightList.currentItem().data(0, 0).toString()
                        ):
                            self.currFit = fit
                            self.currFitNode = self.rightList.currentItem()
                            break

    def getLeftInd(self, embryoName=None, fitName=None):
        if embryoName != None:
            self.currTargetEmbryoNode = self.leftList.findItems(
                embryoName, QtCore.Qt.MatchExactly, 0
            )
            self.currTargetEmbryoNode = self.currTargetEmbryoNode[0]
        if fitName != None:
            self.currTargetFit_node = self.leftList.findItems(
                fitName, QtCore.Qt.MatchExactly, 0
            )
            for fit in self.currEmbryo.fits:
                if fit.name == fitName:
                    self.currTargetFit = fit

    def getRightInd(self, embryoName=None, fitName=None):
        if embryoName != None:
            self.currTargetEmbryoNode = self.rightList.findItems(
                embryoName, QtCore.Qt.MatchExactly, 0
            )
            self.currTargetEmbryoNode = self.currTargetEmbryoNode[0]
        if fitName != None:
            self.currTargetFit_node = self.rightList.findItems(
                fitName, QtCore.Qt.MatchExactly, 0
            )
            for fit in self.currEmbryo.fits:
                if fit.name == fitName:
                    self.currTargetFit = fit

    def getSelectedFits(self):
        return self.fitsInRightList

    def donePressed(self):
        self.molecule.selFits = self.fitsInRightList
        self.done(1)


# ===================================================================================================================================
# Dialog for selecting crucial Parameters
# ===================================================================================================================================


class crucialParameterSelector(pyfrp_gui_basics.listSelectorDialog):
    def __init__(self, molecule, parent):

        self.molecule = molecule

        fittingParms = self.getFittingParms()

        pyfrp_gui_basics.listSelectorDialog.__init__(
            self,
            parent,
            fittingParms,
            leftTitle="Available Parameters",
            rightTitle="Selected Parameters",
            itemsRight=self.molecule.crucialParameters,
        )

        self.setWindowTitle("Select crucial parameters")
        self.show()

    def getFittingParms(self):

        fit = None

        for embryo in self.molecule.embryos:
            if embryo.isFitted():
                try:
                    fit = embryo.fits[0]
                except IndexError:
                    pass

        if fit == None:
            printError(
                "No suitable fit found to select fitting parameters from, will not be able to iniate widget."
            )
            self.done(1)

        return list(vars(fit).keys())

    def donePressed(self):
        self.molecule.crucialParameters = self.itemsRight
        self.done(1)
