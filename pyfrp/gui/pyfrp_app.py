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

#! /usr/bin/python

# =====================================================================================================================================
# Description
# =====================================================================================================================================

# PyFRAP is a free Python based software to analyze FRAP measurements. This file contains the main GUI helping to access the different
# modules and classes:

# =====================================================================================================================================
# Importing Packages and Modules
# =====================================================================================================================================

# Standard packages
import sys
import os, os.path
import time
import copy as cpy
import functools
import code

# Numpy/Scipy
import numpy as np

# PyFRAP Modules
from .pyfrp_term import *
from pyfrp.modules.pyfrp_term_module import *
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_IO_module
from pyfrp.modules import pyfrp_plot_module
from pyfrp.modules import pyfrp_stats_module


# PyFRAP GUIs
from pyfrp.gui import pyfrp_gui_molecule_dialogs
from pyfrp.gui import pyfrp_gui_embryo_dialogs
from pyfrp.gui import pyfrp_gui_ROI_manager
from pyfrp.gui import pyfrp_gui_geometry_dialogs
from pyfrp.gui import pyfrp_gui_gmsh_editor
from pyfrp.gui import pyfrp_gui_analysis_dialogs
from pyfrp.gui import pyfrp_gui_simulation_dialogs
from pyfrp.gui import pyfrp_gui_mesh_dialogs
from pyfrp.gui import pyfrp_gui_fit_dialogs
from pyfrp.gui import pyfrp_gui_pinning_dialogs
from pyfrp.gui import pyfrp_gui_statistics_dialogs
from pyfrp.gui import pyfrp_gui_basics
from pyfrp.gui import pyfrp_gui_settings_dialogs

# PyFRAP Classes
from pyfrp.subclasses import pyfrp_conf
from pyfrp.subclasses import pyfrp_molecule

# QT
from PyQt5 import QtGui, QtCore, QtWidgets

# Matplotlib
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas

"""
Apparently the NavigationToolbar naming has changed in newer matplotlib versions, thus
we need to test out some cases.
"""

try:
    from matplotlib.backends.backend_qt4agg import (
        NavigationToolbar2QTAgg as NavigationToolbar,
    )
except ImportError:
    try:
        from matplotlib.backends.backend_qt4agg import (
            NavigationToolbar2QT as NavigationToolbar,
        )
    except ImportError:
        printWarning("Cannot import NavigationToolbar.")

from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

# VTK
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

# bioformats
import javabridge
import bioformats

# =====================================================================================================================================
# Main Simulation window
# =====================================================================================================================================


class pyfrp(QtWidgets.QMainWindow):

    # Initializes main GUI window
    def __init__(self, parent=None, redirect=False, overwrite_conf=False):
        QtWidgets.QWidget.__init__(self, parent)

        # -------------------------------------------
        # Version
        # -------------------------------------------

        from pyfrp import __version__ as version

        # -------------------------------------------
        # Initializing window
        # -------------------------------------------

        self.setWindowTitle("PyFRAP")
        self.setMinimumSize(400, 300)
        self.resize(1500, 1000)
        self.dpi = 100
        self.version = version
        self.website = "http://www.fml.tuebingen.mpg.de/de/mueller-group/software.html"
        self.pyfrpDir = pyfrp_misc_module.getConfDir().replace("configurations/", "")

        # -------------------------------------------
        # Statusbar
        # -------------------------------------------

        self.statusBar().showMessage("Idle")

        # -------------------------------------------
        # Some variables
        # -------------------------------------------

        # Keep track of plotting tabs
        self.tabAxes = []
        self.tabFigs = []

        # Keep track of all molecules currently open
        self.molecules = []

        # Keep track of which objects are currently used
        self.currMolecule = None
        self.currNode = None
        self.currMoleculeNode = None
        self.currObj = None

        # Keep track of which folder was recently used
        self.lastopen = os.getcwd()

        # -------------------------------------------
        # Menubar
        # -------------------------------------------

        self.menubar = self.menuBar()

        # Init Menubars
        self.initMainMenubar()
        self.initMenubars()

        # Setting
        # -------------------------------------------
        # Embryo list
        # -------------------------------------------

        self.objectBar = QtWidgets.QTreeWidget()
        self.objectBar.setHeaderLabels(["Object", "Analyzed", "Simulated", "Fitted"])
        self.objectBar.setColumnWidth(0, 200)
        self.objectBar.setColumnWidth(1, 75)
        self.objectBar.setColumnWidth(2, 75)
        self.objectBar.itemClicked.connect(self.updatePropBar)

        # -------------------------------------------
        # Property bar
        # -------------------------------------------

        self.propBar = QtWidgets.QListWidget()
        # self.propBar.itemDoubleClicked.connect(self.edit_prop)

        # -------------------------------------------
        # Console
        # -------------------------------------------

        self.console = PyInterp(self, redirect=redirect)
        self.console.initInterpreter(locals())

        # -------------------------------------------
        # Splitter
        # -------------------------------------------

        # Creating splitters
        self.horizontalSplitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.verticalSplitter = QtWidgets.QSplitter(QtCore.Qt.Vertical)

        # -------------------------------------------
        # Frames
        # -------------------------------------------

        # Creating frames for widgets
        self.objectBarFrame = QtWidgets.QFrame()
        self.objectBarFrame = self.createEmptyFrame(self.objectBarFrame)

        self.propBarFrame = QtWidgets.QFrame()
        self.propBarFrame = self.createEmptyFrame(self.propBarFrame)

        self.terminalFrame = QtWidgets.QFrame()
        self.terminalFrame = self.createEmptyFrame(self.terminalFrame)

        # -------------------------------------------
        # Plotting tabs
        # -------------------------------------------

        self.plotTabs = QtWidgets.QTabWidget()
        self.plotTabs.setTabsClosable(False)
        self.plotTabs.currentChanged.connect(self.currentTabChanged)
        self.plotTabs.tabCloseRequested.connect(self.currentTabClosed)

        # -------------------------------------------
        # Final Layout
        # -------------------------------------------

        ##Adding widgets to splitters
        self.horizontalSplitter.addWidget(self.objectBar)
        self.horizontalSplitter.addWidget(self.plotTabs)
        self.horizontalSplitter.addWidget(self.propBar)
        self.verticalSplitter.addWidget(self.horizontalSplitter)
        self.verticalSplitter.addWidget(self.console)

        # Setting default sizes of splitters
        self.verticalSplitter.setSizes([750, 250])
        self.horizontalSplitter.setSizes([350, 900, 250])

        # Connecting splitter movement to figure size adjustment
        self.horizontalSplitter.splitterMoved.connect(self.adjustCanvas)
        self.verticalSplitter.splitterMoved.connect(self.adjustCanvas)

        # Create first plot tab
        self.createDummpyTab()

        # Load config file
        self.initConfiguration(overwrite=overwrite_conf)

        self.setCentralWidget(self.verticalSplitter)
        QtWidgets.QApplication.setStyle(QtWidgets.QStyleFactory.create("Cleanlooks"))

        self.show()

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # List of Methods
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecInit: Initialization Methods
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def initMainMenubar(self):

        """Initiates main menubar."""

        self.mbFile = self.menubar.addMenu("&File")
        self.mbEdit = self.menubar.addMenu("&Edit")
        self.mbEmbryo = self.menubar.addMenu("&Embryo")
        self.mbAnalysis = self.menubar.addMenu("&Data Analysis")
        self.mbSimulation = self.menubar.addMenu("&Simulation")
        self.mbPinning = self.menubar.addMenu("&Pinning")
        self.mbFitting = self.menubar.addMenu("&Fitting")
        self.mbStatistics = self.menubar.addMenu("&Statistics")
        self.mbSettings = self.menubar.addMenu("&Settings")

        return

    def initMenubars(self):

        """Initiates all menubars."""

        self.initFileMenubar()
        self.initEditMenubar()
        self.initEmbryoMenubar()
        self.initAnalysisMenubar()
        self.initPinningMenubar()
        self.initSimulationMenubar()
        self.initFittingMenubar()
        self.initStatsMenubar()
        # self.initPlottingMenubar()
        self.initSettingsMenubar()

    def initFileMenubar(self):

        """Creates entries of file menubar and connects actions with gui methods.
		"""

        newMoleculeButton = QtWidgets.QAction("New Molecule", self)
        newMoleculeButton.triggered.connect(self.newMolecule)

        loadMoleculeButton = QtWidgets.QAction("Open Molecule", self)
        loadMoleculeButton.triggered.connect(self.loadMolecule)

        saveMoleculeButton = QtWidgets.QAction("Save Molecule", self)
        saveMoleculeButton.triggered.connect(self.saveMolecule)

        exitButton = QtWidgets.QAction("Exit", self)
        exitButton.setShortcut("Ctrl+Q")
        exitButton.triggered.connect(self.close)

        self.mbFile.addAction(newMoleculeButton)
        self.mbFile.addAction(loadMoleculeButton)
        self.recentMB = self.mbFile.addMenu("&Open recent")
        self.mbFile.addAction(saveMoleculeButton)

        self.mbFile.addAction(exitButton)

        return

    def initEditMenubar(self):

        """Creates entries of edit menubar and connects actions with gui methods.
		"""

        editMoleculeButton = QtWidgets.QAction("Edit Molecule", self)
        editMoleculeButton.triggered.connect(self.editMolecule)

        removeMoleculeButton = QtWidgets.QAction("Remove Molecule", self)
        removeMoleculeButton.triggered.connect(self.removeMolecule)

        openWizardButton = QtWidgets.QAction("PyFRAP Wizard", self)
        openWizardButton.triggered.connect(self.openWizard)

        self.mbEdit.addAction(editMoleculeButton)
        self.mbEdit.addAction(removeMoleculeButton)
        self.mbEdit.addAction(openWizardButton)

    def initEmbryoMenubar(self):

        """Creates entries of embryo menubar and connects actions with gui methods.
		"""

        newEmbryoButton = QtWidgets.QAction("New Embryo", self)
        newEmbryoButton.triggered.connect(self.newEmbryo)

        removeEmbryoButton = QtWidgets.QAction("Remove Embryo", self)
        removeEmbryoButton.triggered.connect(self.removeEmbryo)

        editEmbryoButton = QtWidgets.QAction("Edit Embryo", self)
        editEmbryoButton.triggered.connect(self.editEmbryo)

        loadEmbryoButton = QtWidgets.QAction("Load Embryo", self)
        loadEmbryoButton.triggered.connect(self.loadEmbryo)

        selectGeometryButton = QtWidgets.QAction("Select Geometry", self)
        selectGeometryButton.triggered.connect(self.selectGeometry)

        editGeometryButton = QtWidgets.QAction("Edit Geometry", self)
        editGeometryButton.triggered.connect(self.editGeometry)

        editGeoFileButton = QtWidgets.QAction("Edit GeoFile", self)
        editGeoFileButton.triggered.connect(self.editGeoFile)

        drawGeometryButton = QtWidgets.QAction("Plot Geometry", self)
        drawGeometryButton.triggered.connect(self.drawGeometry)

        roiManagerButton = QtWidgets.QAction("ROI Manager", self)
        roiManagerButton.triggered.connect(self.openROIManager)

        defaultROIButton = QtWidgets.QAction("Create Default ROIs", self)
        defaultROIButton.triggered.connect(self.createDefaultROIs)

        defaultROIWizardButton = QtWidgets.QAction("Default ROIs Wizard", self)
        defaultROIWizardButton.triggered.connect(self.defaultROIsWizard)

        indexROIButton = QtWidgets.QAction("Update ROIs indices", self)
        indexROIButton.triggered.connect(self.updateROIIdxs)

        self.mbEmbryo.addAction(newEmbryoButton)
        self.mbEmbryo.addAction(editEmbryoButton)
        self.mbEmbryo.addAction(loadEmbryoButton)
        self.mbEmbryo.addAction(removeEmbryoButton)

        self.geometryMB = self.mbEmbryo.addMenu("&Geometry")
        self.geometryMB.addAction(selectGeometryButton)
        self.geometryMB.addAction(editGeometryButton)
        self.geometryMB.addAction(editGeoFileButton)
        self.geometryMB.addAction(drawGeometryButton)

        self.roiMB = self.mbEmbryo.addMenu("&ROIs")
        self.roiMB.addAction(roiManagerButton)
        self.roiMB.addAction(defaultROIButton)
        self.roiMB.addAction(defaultROIWizardButton)

        self.roiMB.addAction(indexROIButton)

        return

    def initAnalysisMenubar(self):

        """Creates entries of analysis menubar and connects actions with gui methods.
		"""

        editAnalysisButton = QtWidgets.QAction("Analysis Settings", self)
        editAnalysisButton.triggered.connect(self.editAnalysis)

        analyzeEmbryoButton = QtWidgets.QAction("Analyze Embryo", self)
        analyzeEmbryoButton.triggered.connect(self.analyzeEmbryo)

        plotEmbryoAnalysisButton = QtWidgets.QAction(
            "Plot Analysis Result of Embryo", self
        )
        plotEmbryoAnalysisButton.triggered.connect(self.plotAllDataTSOfEmbryo)

        self.mbAnalysis.addAction(editAnalysisButton)
        self.runAnalysisMB = self.mbAnalysis.addMenu("&Run Analysis")
        self.runAnalysisMB.addAction(analyzeEmbryoButton)

        self.plotAnalysisMB = self.mbAnalysis.addMenu("&Plotting")
        self.plotAnalysisMB.addAction(plotEmbryoAnalysisButton)

        return

    def initSimulationMenubar(self):

        """Creates entries of simulation menubar and connects actions with gui methods.
		"""

        editSimulationButton = QtWidgets.QAction("Simulation Settings", self)
        editSimulationButton.triggered.connect(self.editSimulation)

        simulateEmbryoButton = QtWidgets.QAction("Simulate Embryo", self)
        simulateEmbryoButton.triggered.connect(self.simulateEmbryo)

        editMeshButton = QtWidgets.QAction("Mesh Settings", self)
        editMeshButton.triggered.connect(self.editMesh)

        genMeshButton = QtWidgets.QAction("Generate Mesh", self)
        genMeshButton.triggered.connect(self.generateMesh)

        refineMeshButton = QtWidgets.QAction("Refine Mesh", self)
        refineMeshButton.triggered.connect(self.refineMesh)

        addBoundaryLayerMeshAroundROIButton = QtWidgets.QAction(
            "Add boundary layer around ROI", self
        )
        addBoundaryLayerMeshAroundROIButton.triggered.connect(
            self.addBoundaryLayerMeshAroundROI
        )

        forceMeshButton = QtWidgets.QAction("Force Mesh Density", self)
        forceMeshButton.triggered.connect(self.forceMeshDensity)

        forceROIMeshButton = QtWidgets.QAction("Refine Mesh in ROI", self)
        forceROIMeshButton.triggered.connect(self.forceROIMeshDensity)

        printMeshButton = QtWidgets.QAction("Mesh Stats", self)
        printMeshButton.triggered.connect(self.printMeshStats)

        plotMeshButton = QtWidgets.QAction("Plot Mesh", self)
        plotMeshButton.triggered.connect(self.plotMesh)

        plotMeshDensityButton = QtWidgets.QAction("Plot Mesh Density", self)
        plotMeshDensityButton.triggered.connect(self.plotMeshDensity)

        plotSimButton = QtWidgets.QAction("Plot Simulation", self)
        plotSimButton.triggered.connect(self.plotAllSimTSOfEmbryo)

        plotSimAndDataButton = QtWidgets.QAction("Plot Simulaten & Data", self)
        plotSimAndDataButton.triggered.connect(self.plotAllSimAndDataTSOfEmbryo)

        self.mbSimulation.addAction(editSimulationButton)
        self.runSimulationMB = self.mbSimulation.addMenu("&Run Simulation")
        self.runSimulationMB.addAction(simulateEmbryoButton)

        self.meshSimulationMB = self.mbSimulation.addMenu("&Mesh")
        self.meshSimulationMB.addAction(editMeshButton)
        self.meshSimulationMB.addAction(genMeshButton)
        self.meshSimulationMB.addAction(refineMeshButton)
        self.meshSimulationMB.addAction(forceMeshButton)
        self.meshSimulationMB.addAction(forceROIMeshButton)
        self.meshSimulationMB.addAction(addBoundaryLayerMeshAroundROIButton)
        self.meshSimulationMB.addAction(printMeshButton)
        self.meshSimulationMB.addAction(plotMeshButton)
        self.meshSimulationMB.addAction(plotMeshDensityButton)

        self.plotSimulationMB = self.mbSimulation.addMenu("&Plotting")
        self.plotSimulationMB.addAction(plotSimButton)
        self.plotSimulationMB.addAction(plotSimAndDataButton)

    def initPinningMenubar(self):

        """Creates entries of pinning menubar and connects actions with gui methods.
		"""

        defaultPinButton = QtWidgets.QAction("Default Pinning", self)
        defaultPinButton.triggered.connect(self.defaultPinEmbryo)

        idealPinButton = QtWidgets.QAction("Ideal Pinning", self)
        idealPinButton.triggered.connect(self.idealPinEmbryo)

        self.mbPinning.addAction(defaultPinButton)
        self.mbPinning.addAction(idealPinButton)

    def initFittingMenubar(self):

        """Creates entries of fitting menubar and connects actions with gui methods.
		"""

        newFitButton = QtWidgets.QAction("New fit", self)
        newFitButton.triggered.connect(self.newFit)

        removeFitButton = QtWidgets.QAction("Remove fit", self)
        removeFitButton.triggered.connect(self.removeFit)

        editFitButton = QtWidgets.QAction("Edit fit", self)
        editFitButton.triggered.connect(self.editFit)

        # editmultfit = QtWidgets.QAction('Edit multiple fits', self)
        # editmultfit.triggered.connect(self.edit_mult_fit)

        # copyfit = QtWidgets.QAction('Copy fit', self)
        # copyfit.triggered.connect(self.copy_fit)

        # copyfitforall = QtWidgets.QAction('Copy fit into all embryos', self)
        # copyfitforall.triggered.connect(self.copy_fit_to_all)

        performFitButton = QtWidgets.QAction("Perform fit", self)
        performFitButton.triggered.connect(self.performFit)

        # fitall = QtWidgets.QAction('Perform all fits in molecule', self)
        # fitall.triggered.connect(self.perform_fits_molecule)

        plotFitButton = QtWidgets.QAction("Plot fit", self)
        plotFitButton.triggered.connect(self.plotFit)

        # plottrackfit = QtWidgets.QAction('Plot fitting progress', self)
        # plottrackfit.triggered.connect(self.plot_track_fit)

        printFitButton = QtWidgets.QAction("Print fit results", self)
        printFitButton.triggered.connect(self.printFitResults)

        self.mbFitting.addAction(newFitButton)
        self.mbFitting.addAction(editFitButton)
        self.mbFitting.addAction(removeFitButton)
        self.mbFitting.addAction(performFitButton)
        self.mbFitting.addAction(printFitButton)

        self.plotFittingMB = self.mbFitting.addMenu("&Plotting")
        self.plotFittingMB.addAction(plotFitButton)

        return

    def initStatsMenubar(self):

        """Creates entries of statistics menubar and connects actions with gui methods.
		"""

        selectFitsButton = QtWidgets.QAction("Select Fits", self)
        selectFitsButton.triggered.connect(self.selectFits)

        selectCrucialParametersButton = QtWidgets.QAction(
            "Select CrucialParameters", self
        )
        selectCrucialParametersButton.triggered.connect(self.selectCrucialParameters)

        summarizeMoleculeButton = QtWidgets.QAction("Summarize Molecule", self)
        summarizeMoleculeButton.triggered.connect(self.summarizeMolecule)

        tTest = QtWidgets.QAction("Perform standard t-test", self)
        tTest.triggered.connect(self.performtTest)

        tTestWelch = QtWidgets.QAction("Perform Welchs t-test", self)
        tTestWelch.triggered.connect(self.performtTestWelch)

        wilcoxonTest = QtWidgets.QAction("Perform Wilcoxon test", self)
        wilcoxonTest.triggered.connect(self.performWilcoxon)

        mannWhitneyUTest = QtWidgets.QAction("Perform Mann-Whitney-U test", self)
        mannWhitneyUTest.triggered.connect(self.performMannWhitneyUTest)

        shapiroTest = QtWidgets.QAction("Perform Shaprio test", self)
        shapiroTest.triggered.connect(self.performShapiroTest)

        AIC = QtWidgets.QAction("Perform Akaike Information Criterion", self)
        AIC.triggered.connect(self.performAIC)

        self.mbStatistics.addAction(selectFitsButton)
        self.mbStatistics.addAction(selectCrucialParametersButton)
        self.mbStatistics.addAction(summarizeMoleculeButton)

        self.statisticsTestMB = self.mbStatistics.addMenu("&Tests")
        self.statisticsTestMB.addAction(tTest)
        self.statisticsTestMB.addAction(tTestWelch)
        self.statisticsTestMB.addAction(wilcoxonTest)
        self.statisticsTestMB.addAction(mannWhitneyUTest)
        self.statisticsTestMB.addAction(shapiroTest)

        self.statisticsComparisonMB = self.mbStatistics.addMenu("&Model comparison")
        self.statisticsComparisonMB.addAction(AIC)

    # def initPlottingMenubar(self):

    # """Creates entries of statistics menubar and connects actions with gui methods.
    # """

    # newMPLTabButton = QtWidgets.QAction('New Matplotlib tab', self)
    # newMPLTabButton.triggered.connect(self.createPlotTab)

    # newVTKTabButton = QtWidgets.QAction('New Matplotlib tab', self)
    # newVTKTabButton.triggered.connect(self.createVtkTab)

    # self.newTabMB=self.mbPlotting.addMenu('&New Tab')
    # self.newTabMB.addAction(newMPLTabButton)
    # self.newTabMB.addAction(newVTKTabButton)

    # showImgIndButton = QtWidgets.QAction('Show image indices', self)
    # showImgIndButton.triggered.connect(self.showROIImgIdx)

    # showMeshIndButton = QtWidgets.QAction('Show mesh indices', self)
    # showImgIndButton.triggered.connect(self.showROIMeshIdx)

    # showExtIndButton = QtWidgets.QAction('Show extended indices', self)
    # showExtIndButton.triggered.connect(self.showROIExtIdx)

    # showROIIndButton = QtWidgets.QAction('Show all indices', self)
    # showROIIndButton.triggered.connect(self.showROIIdx)

    # showROIBoundButton = QtWidgets.QAction('Show ROI boundaries', self)
    # showROIBoundButton.triggered.connect(self.showROIBoundaries)

    # showROIBoundButton = QtWidgets.QAction('Show all ROI boundaries', self)
    # showROIBoundButton.triggered.connect(self.showAllROIBoundaries)

    # self.embryoPlotMB=self.mbPlotting.addMenu('&Embryo')
    # self.ROIPlotMB=self.embryoPlotMB.addMenu('&ROI')

    # self.geomPlotMB=self.embryoPlotMB.addMenu('&Geometry')
    # self.analysisPlotMB=self.embryoPlotMB.addMenu('&Analysis')
    # self.simulationPlotMB=self.embryoPlotMB.addMenu('&Simulation')
    # self.fittingPlotMB=self.embryoPlotMB.addMenu('&Fitting')
    # self.statisticsPlotMB=self.embryoPlotMB.addMenu('&Statistics')

    def initSettingsMenubar(self):

        """Creates entries of settings menubar and connects actions with gui methods.
		"""

        setPathFileButton = QtWidgets.QAction("Set Path File", self)
        setPathFileButton.triggered.connect(self.setPathFile)

        setGmshButton = QtWidgets.QAction("Set Gmsh Path", self)
        setGmshButton.triggered.connect(self.setGmshPath)

        setFijiButton = QtWidgets.QAction("Set Fiji Path", self)
        setFijiButton.triggered.connect(self.setFijiPath)

        setOpenscadButton = QtWidgets.QAction("Set Openscad Path", self)
        setOpenscadButton.triggered.connect(self.setOpenscadPath)

        printPathsButton = QtWidgets.QAction("Print Paths", self)
        printPathsButton.triggered.connect(self.printPaths)

        printPathFileButton = QtWidgets.QAction("Print Path File Location", self)
        printPathFileButton.triggered.connect(self.printPathFile)

        checkPathsButton = QtWidgets.QAction("Check paths in path file", self)
        checkPathsButton.triggered.connect(self.checkPaths)

        self.mbSettings.addAction(setPathFileButton)
        self.mbSettings.addAction(setGmshButton)
        self.mbSettings.addAction(setFijiButton)
        self.mbSettings.addAction(setOpenscadButton)
        self.mbSettings.addAction(printPathsButton)
        self.mbSettings.addAction(printPathFileButton)
        self.mbSettings.addAction(checkPathsButton)

    def closeEvent(self, event):

        """Closes GUI and saves configuration.
		
		Args:
			event (QtCore.closeEvent): Close event triggered by close signal.
		
		"""

        reply = QtWidgets.QMessageBox.question(
            self,
            "Message",
            "Are you sure you want to quit?",
            QtWidgets.QMessageBox.Yes,
            QtWidgets.QMessageBox.No,
        )

        if reply == QtWidgets.QMessageBox.Yes:
            fn = pyfrp_misc_module.getConfDir() + "lastConfiguration.conf"
            self.config.consoleHistory = self.console.history
            self.config.backupPathFile()
            self.config.save(fn=fn)

            javabridge.kill_vm()

            event.accept()
        else:
            event.ignore()

        return

    # def show_about(self):

    # ret=pyfrp_subwin.about_dialog(self).exec_()

    def createEmptyFrame(self, frame):

        """Creates frame around widget.
		
		Args:
			frame (QtWidgets.QWidget): Widget to be framed
		
		Returns: 
			QtWidgets.QWidget: Framed Widget
		
		"""

        frame.setFrameStyle(QtWidgets.QFrame.StyledPanel)
        frame.setBackgroundRole(QtGui.QPalette.Light)
        frame.setAutoFillBackground(True)
        frame.setLineWidth(1)
        frame.setFrameShadow(frame.Sunken)

        return frame

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecConfig: Configuration handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def initConfiguration(self, overwrite=False):

        """Initialize Configuration from file.
		
		Does this by:
		
			* Loading config file.
			* Updating history
			* Calling :py:func:`updateConfig`.
		
		"""

        fn = pyfrp_misc_module.getConfDir() + "lastConfiguration.conf"

        if os.path.isfile(fn) and not overwrite:
            self.config = pyfrp_IO_module.loadFromPickle(fn)
        else:
            self.config = pyfrp_conf.configuration()

        self.console.history = list(self.config.consoleHistory)

        self.updateConfig()

        return self.config

    def updateConfig(self):

        """Update configuration of GUI to match with config.
		
		That is:
			
			* hide/show parts of GUI.
			* Update recently opened files menubar.
			* Update Path file.
			
		"""

        self.console.setHidden(self.config.termHidden)
        self.propBar.setHidden(self.config.propHidden)
        self.plotTabs.setHidden(self.config.plotHidden)

        self.verticalSplitter.refresh()
        self.horizontalSplitter.refresh()
        # self.adjustCanvas()

        self.updateRecentMBs()

        self.config.copyPathFileToDefaultLocation()

        return self.config

    def appendRecent(self, fn):

        """Makes filename the most recently opened file, adds it to config.recentFiles and updates menubar.  
		
		Args:
			fn (str): Some filename.
		"""

        self.config.addRecentFile(fn)
        self.recentMB.clear()
        self.updateRecentMBs()

        return self.recentMB

    def updateRecentMBs(self):

        """Updates recently opened menubar.
		
		Loops through current loaded configuration ``config.recentFiles`` and appends the 5 newest entries to 
		``self.recentActions``.
		
		"""

        self.recentActions = []

        for i in range(len(self.config.recentFiles)):
            if i > 5:
                self.config.recentFiles.pop(i)
            else:
                self.recentActions.append(
                    QtWidgets.QAction(self.config.recentFiles[i], self)
                )
                item = self.recentActions[i]
                item.triggered.connect(
                    functools.partial(
                        self.openMolecule,
                        pyfrp_misc_module.fixPath(self.config.recentFiles[i]),
                    )
                )

                self.recentMB.addAction(item)

        return self.recentMB

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # secObjectBar: Object Bar handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def highlightCurrentObjectNode(self):

        """Highlists the currently used Node of the objectBar.
		"""

        self.objectBar.setCurrentItem(self.currNode)
        return self.objectBar

    def checkSelectedNode(self, typ="any"):

        """Checks if the current selected node is of certain type and returns respective error popups
		
		Keyword Args:
			typ (str): Desired Node type. Use ``"any"`` if type doesn't matter.
		
		Returns:
			bool: True if proper node is selected
		
		"""

        if self.objectBar.currentItem() == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "Nothing selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return False
        if typ == "any":
            return True
        if typ == self.currNodeType:
            return True

        return False

    def embryo2ObjectBar(self, embryo, parentNode):

        """Adds Embryo to ObjectBar. 
		
	
		
		"""

        # Embryo Status
        analyzed = str(int(embryo.isAnalyzed()))
        simulated = str(int(embryo.isSimulated()))
        fitted = str(int(embryo.isFitted()))

        # Create new node
        newEmbryoNode = QtWidgets.QTreeWidgetItem(
            parentNode, [embryo.name, analyzed, simulated, fitted]
        )

        # Geometry node
        newGeometryNode = QtWidgets.QTreeWidgetItem(
            newEmbryoNode, ["Geometry", "", "", ""]
        )

        # ROIs node
        newROIsNode = QtWidgets.QTreeWidgetItem(newEmbryoNode, ["ROIs", "", "", ""])

        # Add ROIs
        for roi in embryo.ROIs:
            newROINode = QtWidgets.QTreeWidgetItem(
                newROIsNode,
                [
                    roi.name,
                    str(int(roi.isAnalyzed())),
                    str(int(roi.isSimulated())),
                    str(int(roi.isFitted())),
                ],
            )

        # Analysis node
        newAnalysisNode = QtWidgets.QTreeWidgetItem(
            newEmbryoNode, ["Analysis", analyzed, "", ""]
        )

        # Simulation node
        newSimulationNode = QtWidgets.QTreeWidgetItem(
            newEmbryoNode, ["Simulation", "", simulated, ""]
        )

        # if embryo.simulation!=None:
        newMeshNode = QtWidgets.QTreeWidgetItem(newSimulationNode, ["Mesh", "", "", ""])

        # Fits node
        newFitsNode = QtWidgets.QTreeWidgetItem(newEmbryoNode, ["Fits", "", "", ""])

        # Add Fits
        for fit in embryo.fits:
            newFitNode = QtWidgets.QTreeWidgetItem(
                newFitsNode, [fit.name, "", "", str(int(fit.isFitted()))]
            )

        # Expand Fits node
        self.objectBar.expandItem(newFitsNode)

        return newEmbryoNode

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Update ObjectBar Properties of embryo Node

    def updateEmbryoNodeProps(self, embryoNode):
        currEmbryoNode = self.getCurrentEmbryoNode()
        currEmbryo = self.getCurrentEmbryo()

        currEmbryoNode.setText(1, str(int(currEmbryo.isAnalyzed())))
        currEmbryoNode.setText(2, str(int(currEmbryo.isSimulated())))
        currEmbryoNode.setText(3, str(int(currEmbryo.isFitted())))

        return

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Update ObjectBar Properties of embryo Node

    def updateEmbryoNodeChildren(self):
        self.updateFitsNodeChildren()
        self.updateROIsNodeChildren()
        return

    def getChildByName(self, node, name):
        for i in range(node.childCount()):
            try:
                curr_name = str(node.child(i).data(0, 0).toString())
            except AttributeError:
                curr_name = str(node.child(i).data(0, 0))

            if name == str(curr_name):
                return node.child(i)
        return None

    def removeAllNodeChildren(self, node):
        nChild = node.childCount()
        for i in range(nChild):
            node.takeChild(0)

        return node

    def updateFitsNodeChildren(self):

        currEmbryo = self.getCurrentEmbryo()
        currEmbryoNode = self.getCurrentEmbryoNode()
        if currEmbryoNode != None:
            fitsNode = self.getChildByName(currEmbryoNode, "Fits")
            self.removeAllNodeChildren(fitsNode)
            for fit in currEmbryo.fits:
                newFitNode = QtWidgets.QTreeWidgetItem(
                    fitsNode, [fit.name, "", "", str(int(fit.isFitted()))]
                )
        self.getCurrentObjNode()
        return newFitNode

    def updateROIsNodeChildren(self):

        currEmbryo = self.getCurrentEmbryo()
        currEmbryoNode = self.getCurrentEmbryoNode()
        if currEmbryoNode != None:
            roisNode = self.getChildByName(currEmbryoNode, "ROIs")
            self.removeAllNodeChildren(roisNode)
            for roi in currEmbryo.ROIs:
                newROINode = QtWidgets.QTreeWidgetItem(
                    roisNode,
                    [
                        roi.name,
                        str(int(roi.isAnalyzed())),
                        str(int(roi.isSimulated())),
                        str(int(roi.isFitted())),
                    ],
                )
        self.getCurrentObjNode()

    def getCurrentObjNode(self):

        """Returns Current Object Node from objectBar
		"""

        self.currNode = self.objectBar.currentItem()

        if self.currNode.parent() == None:
            self.currNodeType = "molecule"
        else:
            if self.currNode.parent().parent() == None:
                self.currNodeType = "embryo"
            elif self.currNode.data(0, 0).toString() == "Fits":
                self.currNodeType = "fits"
            elif self.currNode.parent().data(0, 0).toString() == "Fits":
                self.currNodeType = "fit"
            elif self.currNode.data(0, 0).toString() == "Simulation":
                self.currNodeType = "simulation"
            elif self.currNode.data(0, 0).toString() == "Analysis":
                self.currNodeType = "analysis"
            elif self.currNode.data(0, 0).toString() == "Mesh":
                self.currNodeType = "mesh"
            elif self.currNode.data(0, 0).toString() == "Geometry":
                self.currNodeType = "geometry"
            elif self.currNode.data(0, 0).toString() == "ROIs":
                self.currNodeType = "rois"
            elif self.currNode.parent().data(0, 0).toString() == "ROIs":
                self.currNodeType = "roi"

        return self.currNode, self.currNodeType

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Returns Current Object selected in objectBar

    def getCurrentObj(self):

        # Find out current node and which type it is
        self.getCurrentObjNode()

        # Find current molecule node
        self.getCurrentMoleculeNode()

        # Find corresponding molecule object
        self.getCurrentMolecule()

        if self.currNodeType == "molecule":
            self.currObj = self.currMolecule
        else:
            currEmbryo = self.getCurrentEmbryo()

            if self.currNodeType == "embryo":
                self.currObj = currEmbryo
            elif self.currNodeType == "analysis":
                self.currObj = currEmbryo.analysis
            elif self.currNodeType == "geometry":
                self.currObj = currEmbryo.geometry
            elif self.currNodeType == "simulation":
                self.currObj = currEmbryo.simulation
            elif self.currNodeType == "mesh":
                self.currObj = currEmbryo.simulation.mesh
            elif self.currNodeType == "fit":
                self.currObj = currEmbryo.getFitByName(
                    self.currNode.data(0, 0).toString()
                )
            elif self.currNodeType == "roi":
                self.currObj = currEmbryo.getROIByName(
                    self.currNode.data(0, 0).toString()
                )
            else:
                self.currObj = None

        return self.currObj

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Returns Current Molecule Node

    def getCurrentMoleculeNode(self):
        if self.currNodeType == "molecule":
            self.currMoleculeNode = self.currNode
        elif self.currNodeType == "embryo":
            self.currMoleculeNode = self.currNode.parent()
        elif self.currNodeType in [
            "analysis",
            "simulation",
            "fits",
            "geometry",
            "rois",
        ]:
            self.currMoleculeNode = self.currNode.parent().parent()
        elif self.currNodeType in ["mesh", "roi", "fit"]:
            self.currMoleculeNode = self.currNode.parent().parent().parent()

        ###NOTE: Try later if this is the smarter way
        # ind=self.objectBar.indexOfTopLevelItem(self.currNode)
        # self.currMoleculeNode=self.objectBar.topLevelItem(ind)

        return self.currMoleculeNode

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Returns Current Embryo Node

    def getCurrentEmbryoNode(self):
        if self.currNodeType == "molecule":
            currEmbryoNode = None
        elif self.currNodeType == "embryo":
            currEmbryoNode = self.currNode
        elif self.currNodeType in [
            "analysis",
            "simulation",
            "fits",
            "geometry",
            "rois",
        ]:
            currEmbryoNode = self.currNode.parent()
        elif self.currNodeType in ["mesh", "roi", "fit"]:
            currEmbryoNode = self.currNode.parent().parent()

        return currEmbryoNode

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Returns Current Molecule

    def getCurrentMolecule(self):
        molNames = pyfrp_misc_module.objAttrToList(self.molecules, "name")
        try:
            self.currMolecule = self.molecules[
                molNames.index(self.currMoleculeNode.data(0, 0).toString())
            ]
        except AttributeError:
            self.currMolecule = self.molecules[
                molNames.index(self.currMoleculeNode.data(0, 0))
            ]
        return self.currMolecule

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Returns Current Embryo

    def getCurrentEmbryo(self):
        currEmbryoNode = self.getCurrentEmbryoNode()
        if currEmbryoNode == None:
            return None

        embryoNames = pyfrp_misc_module.objAttrToList(self.currMolecule.embryos, "name")
        try:
            currEmbryo = self.currMolecule.embryos[
                embryoNames.index(currEmbryoNode.data(0, 0).toString())
            ]
        except AttributeError:
            currEmbryo = self.currMolecule.embryos[
                embryoNames.index(currEmbryoNode.data(0, 0))
            ]

        return currEmbryo

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecPropertyBar: Porperty Bar handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def showObjPropsInBar(self, obj, maxArraySize=3):

        """Shows property bar.
		
		Args:
			obj (object): Object selected in objectBar.
			
		Keyword Args:
			maxArraySize (int): Maximal length of array to still be displayed.
				
		"""

        if obj == None:
            return self.propBar

        for item in vars(obj):

            ###NOTE: Creates deprecation warning on Windows

            # First check if numpy array or list
            if isinstance(vars(obj)[str(item)], (np.ndarray, list)):

                # Check if array is not too long
                if (
                    len(np.shape(vars(obj)[str(item)])) > 0
                    and np.shape(vars(obj)[str(item)])[0] < maxArraySize
                ):
                    pass
                else:
                    continue

            # All other variables, check if type is OK
            if (
                isinstance(vars(obj)[str(item)], (int, float, str))
                or vars(obj)[str(item)] == None
            ):
                pass
            else:
                continue

            self.propBar.addItem(item + "=" + str(vars(obj)[str(item)]))

        return self.propBar

    def updatePropBar(self):

        """Updates Property bar."""

        # Clearing property list
        self.propBar.clear()

        # Get Current Object
        self.getCurrentObj()

        # Show properties
        self.showObjPropsInBar(self.currObj)

        # Sort prop_list
        self.propBar.sortItems()

        return self.propBar

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecMolecule: Molecule handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Creates new molecule object and asks for name, then automatically adds it to objectBar

    def newMolecule(self):

        # Generate name
        moleculeNames = pyfrp_misc_module.objAttrToList(self.molecules, "name")
        newName = pyfrp_misc_module.enumeratedName(
            "newMolecule", moleculeNames, sep="_"
        )

        # Create new molecule object
        mol = pyfrp_molecule.molecule(newName)
        pyfrp_gui_molecule_dialogs.moleculeDialog(mol, self).exec_()

        # Make current molecule
        self.currMolecule = mol

        # Add new molecule to list of embryos
        self.molecules.append(mol)

        # Add to objectBar
        self.currNode = QtWidgets.QTreeWidgetItem(
            self.objectBar, [self.currMolecule.name]
        )

        # Memorize molecule node
        self.currMoleculeNode = self.currNode

        # Highligth new node
        self.highlightCurrentObjectNode()

        # Show molecule properties
        self.updatePropBar()

        # Check if names are alright
        self.checkObjNames()

        return self.currMolecule

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Edit Molecule

    def editMolecule(self):

        self.getCurrentMolecule()
        ret = pyfrp_gui_molecule_dialogs.moleculeDialog(self.currMolecule, self).exec_()

        self.currMoleculeNode.setText(0, self.currMolecule.name)

        self.updatePropBar()

        return self.currMolecule

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Removes Molecule from PyFRAP

    def removeMolecule(self):

        # Check if molecule or subnode selected
        if not self.checkSelectedNode():
            return

        # Remove from list of molecules used
        self.molecules.remove(self.currMolecule)

        # Remove from objectBar
        ind = self.objectBar.indexOfTopLevelItem(self.currMoleculeNode)
        self.objectBar.takeTopLevelItem(ind)

        # Update current Node
        self.currNode = self.objectBar.currentItem()

        # Update PropBar
        if self.currNode != None:
            self.updatePropBar()
        else:
            self.propBar.clear()

        return self.molecules

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Saves molecule object into pickle file

    def saveMolecule(self):

        # Check if molecule or subnode selected
        if not self.checkSelectedNode():
            return

        # Filename dialog
        fn = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save file",
            self.lastopen + "/" + self.currMolecule.name + ".mol",
            "*.mol",
        )
        fn = str(fn)

        # Check if folder exists
        if not os.path.exists(os.path.dirname(fn)):
            return

        # Remember folder
        self.lastopen = os.path.dirname(str(fn))

        # Ask if extract save
        # reply = QtWidgets.QMessageBox.question(self, 'Message',"Do you want to save  (Will only keep essential results and delete everything else)?", QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        # No slim save
        # if reply == QtWidgets.QMessageBox.No:
        self.currMolecule.save(fn=fn)
        # else:

        ##Temporarily save molecule
        # if self.config.backup2File:
        # fn_backup=self.lastopen+"/"+self.currMolecule.name+"_backup.mol"
        # self.currMolecule.save(fn=fn_backup)
        # tempMolecule=self.currMolecule

        ##Temporarily make a copy of molecule into memory
        # if self.config.backup2Memory:
        # tempMolecule=cpy.deepcopy(self.currMolecule)

        ####NOTE: Need to specify what needs to be deleted for slim save
        ##Make a function in embryo object for this
        ##Delete big stuff
        ##for temp_emb in temp_mol.embryos:

        ###Deleting all image data to reduce file size
        ##temp_emb.masks_embryo=[]

        # if self.config.backup2File:
        # self.currMolecule=pyfrp_IO_module.loadFromPickle(fn_backup)
        # os.remove(fn_backup)
        # if self.config.backup2Memory:
        # self.currMolecule=cpy.deepcopy(tempMolecule)
        # tempMolecule=None

        self.appendRecent(fn)

        return fn

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Load molecule

    def loadMolecule(self):

        fnLoad = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open file", self.lastopen, "*.mol"
        )
        if fnLoad == "":
            return

        self.currMolecule = self.openMolecule(fnLoad)

        return self.currMolecule

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Open molecule

    def openMolecule(self, fnLoad):

        # Memorize last path and append to recently opened
        self.lastopen = os.path.dirname(str(fnLoad))
        self.appendRecent(fnLoad)

        # Load molecule object
        self.currMolecule = pyfrp_IO_module.loadFromPickle(fnLoad)

        # Update Version
        self.currMolecule.updateVersion()

        # Add molecule to list of molecules
        self.molecules.append(self.currMolecule)

        # Adding molecule to sidebar
        self.currNode = QtWidgets.QTreeWidgetItem(
            self.objectBar, [self.currMolecule.name, "", "", ""]
        )
        self.currMoleculeNode = self.currNode

        for embryo in self.currMolecule.embryos:
            self.embryo2ObjectBar(embryo, self.currMoleculeNode)

        self.objectBar.expandItem(self.currNode)

        self.checkObjNames()

        return self.currMolecule

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Embryo handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def newEmbryo(self):

        """Creates new embryo object and adds it to molecule.
		
		Will first run :py:func:`selectEmbryoGenMode` to choose way of embryo generation.
		Will create embryo instance via :py:func:`LSM2EmbryoWizard` if selected, ohterwise
		will create default embryo.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: New embryo object.
		
		"""

        # Check if molecule or subnode selected
        if not self.checkSelectedNode():
            return

        ret = self.selectEmbryoGenMode()
        if ret == -1:
            return 0
        elif ret == 0:

            # Generate name
            embryoNames = pyfrp_misc_module.objAttrToList(
                self.currMolecule.embryos, "name"
            )
            newName = pyfrp_misc_module.enumeratedName(
                "newEmbryo", embryoNames, sep="_"
            )

            # Create new embryo
            newEmbryo = self.currMolecule.newEmbryo(newName)

            # Append to object bar
            newNode = self.embryo2ObjectBar(newEmbryo, self.currMoleculeNode)
            self.objectBar.setCurrentItem(newNode)
            self.getCurrentObj()

        # Call embryo editor
        ret = self.editEmbryo()
        if ret == 0:
            return 0

        # Geometry Editor
        ret = self.selectGeometry()
        if ret == 0:
            return 0

        ret = self.editGeometry()
        if ret == 0:
            return 0

        return 1

    def loadEmbryo(self):

        """Loads saved embryo object and adds it to currently selected molecule file.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: Loaded embryo object.
		
		"""

        # Check if molecule or subnode selected
        if not self.checkSelectedNode():
            return

        # Get filename
        fnLoad = QtWidgets.QFileDialog.getOpenFileName(
            self, "Open file", self.lastopen, "*.emb"
        )

        # Load and add to molecule
        newEmbryo = pyfrp_misc_module.loadFromPickle(fnLoad)
        self.currMolecule.addEmbryo(newEmbryo)

        # Append to object bar
        newNode = self.embryo2ObjectBar(newEmbryo, self.currMoleculeNode)
        self.objectBar.setCurrentItem(newNode)
        self.getCurrentObj()

        return newEmbryo

    def editEmbryo(self):

        """Opens main embryo editing dialog.
		
		Returns:
			pyfrp.subclasses.pyfrp_embryo: Edited embryo object.
		
		"""

        currEmbryoNode = self.getCurrentEmbryoNode()
        currEmbryo = self.getCurrentEmbryo()

        if (
            self.checkSelectedNode(typ="molecule")
            or self.objectBar.currentItem() == None
        ):
            return 0

        ret = pyfrp_gui_embryo_dialogs.embryoDialog(currEmbryo, self).exec_()

        if ret == 0:
            return 0

        currEmbryoNode.setText(0, currEmbryo.name)

        self.updatePropBar()

        return 1

    def removeEmbryo(self):

        """Removes selected embryo object both from molecule object as well as 
		main GUI.
		
		Returns:
			list: Updated pyfrp.subclasses.pyfrp_molecule.embryos list
		
		"""

        currEmbryoNode = self.getCurrentEmbryoNode()
        currEmbryo = self.getCurrentEmbryo()

        if (
            self.checkSelectedNode(typ="molecule")
            or self.objectBar.currentItem() == None
        ):
            return

        # Remove from list of embryo objects
        self.currMolecule.embryos.remove(currEmbryo)

        # Remove from sidebar
        ind = self.currMoleculeNode.indexOfChild(currEmbryoNode)
        self.currMoleculeNode.takeChild(ind)

        self.objectBar.setCurrentItem(self.currMoleculeNode)

        self.updatePropBar()

        return self.currMolecule.embryos

    def checkObjNames(self):

        """Checks if molecules or embryos have the same name.
		
		Returns:
			bool: True if there is a molecule/embryo with the same name.
		"""

        b = False

        b = self.checkMolNames()

        for mol in self.molecules:
            b = self.checkEmbryoNames(mol=mol)

        return b

    def checkMolNames(self):

        """Checks if molecules have the same name.
		
		Returns:
			bool: True if there is a molecule the same name.
		"""

        b = False

        molNames = pyfrp_misc_module.objAttrToList(self.molecules, "name")
        for name in molNames:
            if molNames.count(name) > 1:
                b = True
                printWarning(
                    "Molecule with name "
                    + name
                    + " exists twice. This can lead to problems"
                )

        return b

    def checkEmbryoNames(self, mol=None):

        """#Checks if embryos have the same name.
		
		Returns:
			bool: True if there is a embryo the same name.
		"""

        mol = pyfrp_misc_module.assignIfVal(mol, self.currMolecule, None)

        b = False

        embryoNames = pyfrp_misc_module.objAttrToList(self.currMolecule.embryos, "name")
        for name in embryoNames:
            if embryoNames.count(name) > 1:
                b = True
                printWarning(
                    "Embryo with name "
                    + name
                    + " exists twice in molecule "
                    + mol.name
                    + ". This can lead to problems"
                )

        return b

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecWizards: Wizards
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def LSM2EmbryoWizard(self):

        """Opens :py:class:`pyfrp.gui.pyfrp_gui_embryo_dialogs.lsmWizard` which lets
		one to extract and sort microscope files and then creates new 
		:py:class:`pyfrp.subclasses.pyfrp_embryo.embryo` including new data.
		
		.. note:: Uses Fiji. Fiji path needs to be set properly for success.
		
		"""

        lsmWizard = pyfrp_gui_embryo_dialogs.lsmWizard(self)
        if lsmWizard.exec_():
            embryo = lsmWizard.getEmbryo()
            if embryo == None:
                return

            self.currMolecule.addEmbryo(embryo)

            newNode = self.embryo2ObjectBar(embryo, self.currMoleculeNode)
            self.objectBar.setCurrentItem(newNode)
            self.getCurrentObj()

    def openWizard(self):

        """Wizard for PyFRAP analysis.
		
		Leads user through following steps:
		
			* Embryo creation :py:func:`newEmbryo`.
			* ROI creation :py:func:`createDefaultROIs`.
			* ROI handling :py:func:`openROIManager`.
			* Mesh generation :py:func:`generateMesh`.
			* Data analysis :py:func:`analyzeEmbryo`.
			* FRAP simulation :py:func:`simulateEmbryo`.
			* Fit creation :py:func:`newFit`.
			* Fit plotting :py:func:`plotFit`.
			
		"""

        ret = self.newEmbryo()
        if ret == 0:
            return 0

        ret = self.ROIWizardSelector()
        if ret == 0:
            return 0

        ret = self.editGeometry()
        if ret == 0:
            return 0

        self.getCurrentEmbryo().newSimulation()
        ret = self.generateMesh()
        if ret == 0:
            return 0

        self.updateROIIdxs()
        ret = self.analyzeEmbryo()
        if ret == 0:
            return 0

        ret = self.simulateEmbryo()
        if ret == 0:
            return 0

        ret = self.newFit()
        if ret == 0:
            return 0

        self.plotFit()

        return 1

    def selectEmbryoGenMode(self):

        """Opens :py:class:`pyfrp.gui.pyfrp_gui_embryo_dialogs.wizardSelector`
		and lets user select embryo generation mode. 
		
		Returns:
			int: 1 if user chose from microscope, 0 if not, -1 else.
		
		"""

        selectWizard = pyfrp_gui_embryo_dialogs.wizardSelector(self)
        if selectWizard.exec_():
            useLSM = selectWizard.getUseLSM()
            if useLSM == None:
                return -1

            if useLSM:
                ret = self.LSM2EmbryoWizard()
                if ret == 0:
                    return -1

                return 1
            return 0

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecEmbryo: Embryo Menubar Methods
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def openROIManager(self):

        """Open ROI Manager for current embryo.
		
		Returns:
			int: Dialog Code (0 when cancelled, 1 when accepted.)
		"""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            ret = pyfrp_gui_ROI_manager.ROImanager(currEmbryo, self).exec_()

            if ret == 0:
                return ret
            elif ret == 1:
                self.updateROIsNodeChildren()
                return ret
        return 0

    def createDefaultROIs(self):

        """Create default ROIs for current embryo.
		
		Returns:
			int: Dialog Code (0 when cancelled, 1 when accepted.)
		"""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            ret = pyfrp_gui_ROI_manager.defaultROIsDialog(currEmbryo, self).exec_()
            if ret == 0:
                return ret
            elif ret == 1:
                self.updateROIsNodeChildren()
                return ret
        return 0

    def ROIWizardSelector(self):

        """Lets the user select between the createDefaultROIs, defaultROIsWizard or ROImanager.
		
		Opens :py:class:`pyfrp.gui.pyfrp_gui_ROI_manager.wizardSelector`
		and lets user select ROI generation mode. 
		
		Launches the right wizard afterwards.
		
		Returns:
			int: 0 when default, 1 for wizard and 2 for manager
		
		"""

        selectWizard = pyfrp_gui_ROI_manager.wizardSelector(self)
        if selectWizard.exec_():
            mode = selectWizard.getMode()

            print("mode is", mode)

            if mode == None:
                return 0

            if mode == 0:
                ret = self.createDefaultROIs()
            elif mode == 1:
                ret = self.defaultROIsWizard()
            elif mode == 2:
                ret = self.openROIManager()

            return ret

        return 0

    def defaultROIsWizard(self):

        """Small Wizard that helps creating a set of default ROIs for current embryo.
		
		Returns:
			int: Dialog Code (0 when cancelled, 1 when accepted.)
		"""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Select Slice ROI
        sliceSelector = pyfrp_gui_ROI_manager.ROISelector(
            currEmbryo,
            "radialSlice",
            self,
            msg="Select Slice ROI",
            sliceHeight=currEmbryo.sliceHeightPx,
            sliceWidth=currEmbryo.sliceWidthPx,
            color="g",
            name="Slice",
        )
        if sliceSelector.exec_():
            sliceROI = sliceSelector.getROI()

        # Select Bleached ROI
        squSelector = pyfrp_gui_ROI_manager.ROISelector(
            currEmbryo,
            "squareSlice",
            self,
            msg="Select Bleached ROI",
            sliceHeight=currEmbryo.sliceHeightPx,
            sliceWidth=currEmbryo.sliceWidthPx,
            color="b",
            name="Bleached Square",
        )
        if squSelector.exec_():
            bleachedROI = squSelector.getROI()

        # Select Rim ROI
        if hasattr(sliceROI, "radius"):
            radius = sliceROI.radius * 0.66
            center = sliceROI.center
            rimSelector = pyfrp_gui_ROI_manager.ROISelector(
                currEmbryo,
                "radialSlice",
                self,
                msg="Select Rim ROI",
                radius=radius,
                center=center,
                sliceHeight=currEmbryo.sliceHeightPx,
                sliceWidth=currEmbryo.sliceWidthPx,
                color="y",
                name="Slice rim",
            )
        else:
            radius = 133.0
            rimSelector = pyfrp_gui_ROI_manager.ROISelector(
                currEmbryo,
                "radialSlice",
                self,
                msg="Select Rim ROI",
                radius=radius,
                sliceHeight=currEmbryo.sliceHeightPx,
                sliceWidth=currEmbryo.sliceWidthPx,
                color="y",
                name="Slice rim",
            )

        if rimSelector.exec_():
            rimROI = rimSelector.getROI()

        # Build default ROIs
        radius = 300.0
        center = [256.0, 256.0]
        if hasattr(sliceROI, "radius"):
            radius = sliceROI.radius
            center = sliceROI.center
        else:
            if currEmbryo.geometry != None:
                if hasattr(currEmbryo.geometry, "imagingRadius"):
                    radius = currEmbryo.geometry.imagingRadius
                else:
                    if hasattr(currEmbryo.geometry, "radius"):
                        radius = currEmbryo.geometry.radius

        currEmbryo.genDefaultROIs(
            center,
            radius,
            rimFactor=0.66,
            masterROI=sliceROI,
            bleachedROI=bleachedROI,
            rimROI=rimROI,
            clean=True,
        )

        self.updateROIsNodeChildren()

        return 1

    def updateROIIdxs(self):

        """Updates all ROI Idxs."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo != None:

            # Genereate wait popup
            self.progressDialog = pyfrp_gui_ROI_manager.indexProgressDialog(None)

            # Make backup copy of embryo
            self.originalObj = currEmbryo
            self.backupObj = cpy.deepcopy(currEmbryo)

            self.statusBar().showMessage("Indexing ROIs of embryo  " + currEmbryo.name)

            # Generate worker and generate Qthread
            self.task = pyfrp_gui_basics.pyfrpThread()
            self.worker = pyfrp_gui_basics.pyfrpWorker(
                currEmbryo.computeROIIdxs, signal=self.task.progressSignal
            )

            # Start
            self.initTask()
        else:
            return 0

        return 1

    def selectGeometry(self):

        """Open Geometry Editor for current embryo."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo != None:
            ret = pyfrp_gui_geometry_dialogs.geometrySelectDialog(
                currEmbryo, self
            ).exec_()
        else:
            return 0

        return ret

    def editGeometry(self):

        """Open Geometry Editor for current embryo."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo != None:
            if "zebrafishDomeStage" in currEmbryo.geometry.typ:
                ret = pyfrp_gui_geometry_dialogs.zebrafishDomeStageDialog(
                    currEmbryo.geometry, self
                ).exec_()
            elif "cylinder" in currEmbryo.geometry.typ:
                ret = pyfrp_gui_geometry_dialogs.cylinderDialog(
                    currEmbryo.geometry, self
                ).exec_()
            elif "xenopusBall" in currEmbryo.geometry.typ:
                ret = pyfrp_gui_geometry_dialogs.xenopusBallDialog(
                    currEmbryo.geometry, self
                ).exec_()
            elif "cone" in currEmbryo.geometry.typ:
                ret = pyfrp_gui_geometry_dialogs.coneDialog(
                    currEmbryo.geometry, self
                ).exec_()
            else:
                ret = pyfrp_gui_geometry_dialogs.geometryDialog(
                    currEmbryo.geometry, self
                ).exec_()

            return ret

        return 0

    def editGeoFile(self):

        """Open Geometry Editor for current embryo."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            ret = pyfrp_gui_gmsh_editor.gmshFileEditor(
                currEmbryo.geometry, self
            ).exec_()

        return ret

    def drawGeometry(self, ann=False):

        """Draw Geometry in plot tab."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            self.createPlotTab("xyz", plotName="Geometry", proj=["3d"])
            currEmbryo.geometry.plotGeometry(ax=self.ax, ann=ann)
            self.adjustCanvas()
        return

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecAnalysis: Analysis handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def editAnalysis(self):

        """Edit Analysis Settings."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            if currEmbryo.analysis == None:
                currEmbryo.newAnalysis()

            ret = pyfrp_gui_analysis_dialogs.analysisDialog(
                currEmbryo.analysis, self
            ).exec_()
        else:
            return 0

        return ret

    def analyzeEmbryo(self):

        """Analyze embryo task/progressbar."""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Launching Dialog
        ret = self.editAnalysis()
        if ret == 0:
            return 0

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_analysis_dialogs.analysisProgressDialog(None)

        # Make backup copy of embryo
        self.originalObj = currEmbryo
        self.backupObj = cpy.deepcopy(currEmbryo)

        self.statusBar().showMessage("Analyzing Dataset " + currEmbryo.name)

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()
        self.worker = pyfrp_gui_basics.pyfrpWorker(
            currEmbryo.analysis.run, signal=self.task.progressSignal
        )

        # Init and start
        self.initTask()

        return 1

    def plotAllDataTSOfEmbryo(self):

        """Plots data analysis results of all ROIs of embryo."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            return

        self.createPlotTab(
            "intensityTS", plotName=currEmbryo.name + " data", size=[1, 1]
        )
        currEmbryo.plotAllData(ax=self.ax)

        self.adjustCanvas()

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecSimulation: Simulation handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def editSimulation(self):

        """Edit Simulation Settings."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            if currEmbryo.simulation == None:
                currEmbryo.newSimulation()

            ret = pyfrp_gui_simulation_dialogs.simulationSettingsDialog(
                currEmbryo.simulation, self
            ).exec_()

        else:
            return 0

        return ret

    def simulateEmbryo(self):

        """Simulate embryo task/progressbar."""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Launching Dialog
        ret = self.editSimulation()
        if ret == 0:
            return 0

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_simulation_dialogs.simulationProgressDialog(
            None
        )

        # Make backup copy of embryo
        self.originalObj = currEmbryo
        self.backupObj = cpy.deepcopy(currEmbryo)

        self.statusBar().showMessage("Simulating Dataset " + currEmbryo.name)

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()
        self.worker = pyfrp_gui_basics.pyfrpWorker(
            currEmbryo.simulation.run, signal=self.task.progressSignal
        )

        # Init and start
        self.initTask()

        return 1

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecMesh: Mesh handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def editMesh(self):

        """Edit Mesh Settings.
		
		Launches :py:class:`pyfrp.gui.pyfrp_gui_mesh_dialogs.meshSettingsDialog` and lets user define standard mesh properties.
		
		"""
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo != None:
            if currEmbryo.simulation == None:
                currEmbryo.newSimulation()
            ret = pyfrp_gui_mesh_dialogs.meshSettingsDialog(
                currEmbryo.simulation.mesh, self
            ).exec_()
            return ret

        else:
            return 0

    def generateMesh(self):

        """Generate Mesh.
		
		Does this via submitting mesh generation to a seperate `QThread`. Meanwhile displays
		:py:class:`pyfrp.gui.pyfrp_gui_mesh_dialogs.genMeshProgressDialog`.
		
		Prints out mesh statistics at the end.
		
		"""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Launching Dialog
        ret = self.editMesh()
        if ret == 0:
            return 0

        # Genereate wait popup
        # self.progressDialog=pyfrp_gui_mesh_dialogs.genMeshProgressDialog(None)

        ##Make backup copy of embryo
        # self.originalObj=currEmbryo
        # self.backupObj=cpy.deepcopy(currEmbryo)

        # self.statusBar().showMessage("Generating Mesh " + currEmbryo.name)

        ##Generate Qthread and pass analysis there
        # self.task=pyfrp_gui_basics.pyfrpThread()
        # self.worker=pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.mesh.genMesh)
        printNote("Generating Mesh... This can take a few seconds...")
        currEmbryo.simulation.mesh.genMesh()

        # Print out mesh stats
        self.printMeshStats()

        # Init and start
        # self.initTask()

        return 1

    def refineMesh(self):

        """Refine Mesh by tetrahedron splitting.
		
		Submits py:func:`pyfrp.subclasses.pyfrp_mesh.mesh.refine` to a seperate QThread. 
		Meanwhile displays :py:class:`pyfrp.gui.pyfrp_gui_mesh_dialogs.refineMeshProgressDialog`.
		
		"""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_mesh_dialogs.refineMeshProgressDialog(None)

        # Make backup copy of embryo
        self.originalObj = currEmbryo
        self.backupObj = cpy.deepcopy(currEmbryo)

        self.statusBar().showMessage("Refining Mesh " + currEmbryo.name)

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()
        self.worker = pyfrp_gui_basics.pyfrpWorker(currEmbryo.simulation.mesh.refine)

        # Init and start
        self.initTask()

        return 1

    def forceMeshDensity(self):

        """Force Mesh to have given density globally.
		
		Submits py:func:`pyfrp.subclasses.pyfrp_mesh.mesh.refine` to a seperate QThread. 
		Meanwhile displays :py:class:`pyfrp.gui.pyfrp_gui_mesh_dialogs.refineMeshProgressDialog`.
		
		"""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Get options
        forceMeshDialog = pyfrp_gui_mesh_dialogs.forceMeshSettingsDialog(
            currEmbryo.simulation.mesh, self
        )
        if forceMeshDialog.exec_():
            roiUsed, density, stepPercentage, debug, findIdxs, method, maxCells = (
                forceMeshDialog.getVals()
            )

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_mesh_dialogs.forceMeshProgressDialog(None)

        # Make backup copy of embryo
        self.originalObj = currEmbryo
        self.backupObj = cpy.deepcopy(currEmbryo)

        self.statusBar().showMessage(
            "Computing new mesh with required density " + currEmbryo.name
        )

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()
        self.worker = pyfrp_gui_basics.pyfrpWorker(
            currEmbryo.simulation.mesh.forceMinMeshDensityInROI,
            roiUsed,
            density,
            stepPercentage,
            debug=debug,
            findIdxs=findIdxs,
            method=method,
            maxCells=maxCells,
        )

        # Init and start
        self.initTask()

        return 1

    def forceROIMeshDensity(self):

        """Force Mesh to have given density in ROI.
		
		See also :py:class:`pyfrp.gui.pyfrp_gui_mesh_dialogs.refineROIMeshSettingsDialog`.
		
		"""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        # Get options
        ret = pyfrp_gui_mesh_dialogs.refineROIMeshSettingsDialog(
            currEmbryo.simulation.mesh, self
        ).exec_()

        return ret

    def printMeshStats(self):

        """Prints out mesh statistics."""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return 0

        currEmbryo.simulation.mesh.printStats()

    def plotMesh(self):

        """Plots mesh of selected embryo.
		"""

        # Grab embryo
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No embryo selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return

        # try:
        self.renderer = currEmbryo.simulation.mesh.plotMesh(renderWindow=1)
        # except AttributeError:
        # QtWidgets.QMessageBox.critical(None, "Error","Embryo has no simulation or mesh yet.",QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default)

        self.createVtkTab()

        self.vtkCanvas.GetRenderWindow().AddRenderer(self.renderer)
        self.vtkCanvas.GetRenderWindow().GetInteractor().Initialize()

    def plotMeshDensity(self):

        """Plots mesh density in x/y/z direction."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return

        self.createPlotTab(
            "meshDensity", plotName=currEmbryo.name + " mesh density", size=[2, 2]
        )

        currEmbryo.simulation.mesh.plotDensity(axes=self.axes)

    def addBoundaryLayerMeshAroundROI(self):

        """Adds boundary layer around ROI.
		
		"""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return

        dialog = pyfrp_gui_mesh_dialogs.boundaryLayerAroundROISettingsDialog(
            currEmbryo, self
        )
        if dialog.exec_():
            fnOut, roiUsed, segments, simplify, iterations, triangIterations, fixSurfaces, debug, volSizePx, thickness, volSizeLayer, cleanUp, approxBySpline, angleThresh, faces, onlyAbs = (
                dialog.getVals()
            )

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_mesh_dialogs.boundaryLayerProgressDialog(None)

        # Make backup copy of embryo
        self.originalObj = currEmbryo
        self.backupObj = cpy.deepcopy(currEmbryo)

        self.statusBar().showMessage(
            "Computing new mesh with required density " + currEmbryo.name
        )

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()

        # For now don't submit to QThread and see how it works out.
        self.worker = pyfrp_gui_basics.pyfrpWorker(
            currEmbryo.simulation.mesh.addBoundaryLayerAroundROI,
            roiUsed,
            fnOut=fnOut,
            segments=segments,
            simplify=simplify,
            iterations=iterations,
            triangIterations=triangIterations,
            fixSurfaces=fixSurfaces,
            debug=debug,
            volSizePx=volSizePx,
            volSizeLayer=volSizeLayer,
            thickness=thickness,
            cleanUp=cleanUp,
            approxBySpline=approxBySpline,
            angleThresh=angleThresh,
            faces=faces,
            onlyAbs=onlyAbs,
        )

        # Init and start
        self.initTask()

    def plotAllSimTSOfEmbryo(self):

        """Plots simulation results of all ROIs of embryo."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            return

        self.createPlotTab(
            "intensityTS", plotName=currEmbryo.name + " data", size=[1, 1]
        )
        currEmbryo.plotAllSim(ax=self.ax)

        self.adjustCanvas()

    def plotAllSimAndDataTSOfEmbryo(self):

        """Plots simulation results of all ROIs of embryo."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            return

        self.createPlotTab(
            "intensityTS", plotName=currEmbryo.name + " data", size=[1, 1]
        )
        currEmbryo.plotAllData(ax=self.ax)
        currEmbryo.plotAllSim(ax=self.ax)

        self.adjustCanvas()

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecPinning: Pinning handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def defaultPinEmbryo(self):

        """Default Pinning.
		
		Lets the user selected his/her pinning options themself by launching 
		:py:class:`pyfrp.gui.pyfrp_gui_pinning_dialogs.defaultPinningDialog`.
		
		
		"""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No embryo selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Launch dialog
        ret = pyfrp_gui_pinning_dialogs.defaultPinningDialog(currEmbryo, self).exec_()

        return ret

    def idealPinEmbryo(self):

        """Ideal Pinning.
		
		Lets the user select some options for ideal pinning via 
		:py:class:`pyfrp.gui.pyfrp_gui_pinning_dialogs.idealPinningDialog`.
		
		
		"""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No embryo selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Launch dialog
        ret = pyfrp_gui_pinning_dialogs.idealPinningDialog(currEmbryo, self).exec_()

        return ret

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecFit: Fit handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def newFit(self):

        """Creates new fit."""

        currEmbryo = self.getCurrentEmbryo()

        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No embryo selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Generate name
        fitNames = pyfrp_misc_module.objAttrToList(currEmbryo.fits, "name")
        newName = pyfrp_misc_module.enumeratedName("newFit", fitNames, sep="_")

        # Generate new fit
        newFit = currEmbryo.newFit(newName)

        # Launch editor
        ret = pyfrp_gui_fit_dialogs.fitSettingsDialog(newFit, self).exec_()

        if ret == 0:
            return 0

        # Update Object Bar
        newNode = self.updateFitsNodeChildren()
        self.objectBar.setCurrentItem(newNode)

        # Perform the new fit
        self.performFit()

        return 1

    def editFit(self):

        """Opens edit fit dialog."""

        self.getCurrentObj()

        if self.currNodeType != "fit":
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Launch editor
        ret = pyfrp_gui_fit_dialogs.fitSettingsDialog(self.currObj, self).exec_()

        if ret == 0:
            return 0

        # Update Object Bar
        self.updateFitsNodeChildren()

        return 1

    def removeFit(self):

        """Removes fit."""

        self.getCurrentObj()

        if self.currNodeType != "fit":
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        currEmbryo = self.getCurrentEmbryo()

        currEmbryo.deleteFit(currEmbryo.fits.index(self.currObj))

        # Update Object Bar
        self.updateFitsNodeChildren()

        return 1

    def performFit(self):

        """Runs fitting progess for selected fit."""

        self.getCurrentObj()

        if self.currNodeType != "fit":
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Check if pinned
        if self.currObj.fitPinned:
            if not self.currObj.checkPinned():
                self.idealPinEmbryo()

        # Genereate wait popup
        self.progressDialog = pyfrp_gui_fit_dialogs.fittingProgressDialog(None)

        # Make backup copy of embryo
        self.originalObj = self.currObj
        self.backupObj = cpy.deepcopy(self.currObj)

        self.statusBar().showMessage("Performing fit " + self.currObj.name)

        # Generate Qthread and pass analysis there
        self.task = pyfrp_gui_basics.pyfrpThread()
        self.worker = pyfrp_gui_basics.pyfrpWorker(self.currObj.run)

        # Init and start
        self.initTask()

        return 1

    def printFitResults(self):

        """Print fit results."""

        self.getCurrentObj()

        if self.currNodeType != "fit":
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        self.currObj.printResults()

        return 1

    def plotFit(self):

        """Plots currently selected fit."""

        self.getCurrentObj()

        if self.currNodeType != "fit":
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No fit selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        if self.currObj.isFitted():

            self.createPlotTab(
                "intensityTS", plotName=self.currObj.name + " ", size=[1, 1]
            )
            self.currObj.plotFit(ax=self.ax)

        self.adjustCanvas()

        return 1

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecStatistics: Statistics handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def selectFits(self):

        """Opens fit selector for selected molecule."""

        self.getCurrentObj()

        if self.currMoleculeNode == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No molecule selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        if len(self.currMolecule.embryos) == 0:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "Molecule does not have any embryos to average.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Open Fit selector dialog
        ret = pyfrp_gui_statistics_dialogs.fitSelector(
            self.currMolecule, True, self
        ).exec_()

        return ret

    def selectCrucialParameters(self):

        """Opens crucial parameter selector."""

        self.getCurrentObj()

        if self.currMoleculeNode == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No molecule selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        if len(self.currMolecule.embryos) == 0:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "Molecule does not have any embryos to average.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # Open Fit selector dialog
        ret = pyfrp_gui_statistics_dialogs.crucialParameterSelector(
            self.currMolecule, self
        ).exec_()

        return ret

    def summarizeMolecule(self):

        """Summarizes currently selected molecule."""

        self.getCurrentObj()

        if self.currMoleculeNode == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No molecule selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return 0

        # if len(self.currMolecule.selFits)==0:
        ret = self.selectFits()
        if ret == 0:
            return 0

        ret = self.selectCrucialParameters()
        if ret == 0:
            return 0

        self.currMolecule.sumUpResults()

        return 1

    def selMolecules(self, n=2, nmin=2):

        """Opens listSelectorDialog and lets user select n molecules.
		
		Keyword Args:
			n (int): Number of molecules that should be selected
			
		
		Returns:
			list: List of molecules.
		
		"""

        names = []
        for mol in self.molecules:
            names.append(mol.name)

        selDialog = pyfrp_gui_basics.listSelectorDialog(
            self, names, leftTitle="Available", rightTitle="Selected", itemsRight=[]
        )

        selectedMols = []

        if selDialog.exec_():
            selected = selDialog.getSelection()

            for mol in self.molecules:
                if mol.name in selected:
                    selectedMols.append(mol)

        if len(selected) > n:
            printWarning(
                "More than "
                + str(n)
                + " molecules selected. Will only use the first and second molecule."
            )
            selectedMols = selectedMols[0:n]

        if len(selected) < nmin:
            printError("Less than " + str(nmin) + " molecules selected.")
            return []

        return selectedMols

    def performtTest(self):

        """Lets user select two molecules and then performs standard t-test on both of them."""

        selected = self.selMolecules()

        opt1 = selected[0].getDOptMus()
        opt2 = selected[1].getDOptMus()

        stat, pval = pyfrp_stats_module.tTestStandard(opt1, opt2)

    def performtTestWelch(self):

        """Lets user select two molecules and then performs Welch's t-test on both of them."""

        selected = self.selMolecules()

        opt1 = selected[0].getDOptMus()
        opt2 = selected[1].getDOptMus()

        stat, pval = pyfrp_stats_module.tTestWelch(opt1, opt2)

    def performWilcoxon(self):

        """Lets user select two molecules and then performs Wilcoxon on both of them."""

        selected = self.selMolecules()

        opt1 = selected[0].getDOptMus()
        opt2 = selected[1].getDOptMus()

        stat, pval = pyfrp_stats_module.wilcoxonTest(opt1, opt2)

    def performMannWhitneyUTest(self):

        """Lets user select two molecules and then performs Mann-Whitney-U test on both of them."""

        selected = self.selMolecules()

        opt1 = selected[0].getDOptMus()
        opt2 = selected[1].getDOptMus()

        stat, pval = pyfrp_stats_module.mannWhitneyUTest(opt1, opt2)

    def performShapiroTest(self):

        """Lets user select multiple molecules and then performs shapiro test on both of them."""

        selected = self.selMolecules(n=1, nmin=1)

        opt1 = selected[0].getDOptMus()

        stat, pval = pyfrp_stats_module.shapiroTest(opt1)

    def performAIC(self):

        """Peforms Akaike-Information-Criterion for model comparison."""

        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "Nothing selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return

        AICs, deltaAICs, weights, acc, ks, ns = currEmbryo.compareFitsByCorrAIC(
            thresh=0.2
        )

        # fitNames = pyfrp_misc_module.objAttrToList(currEmbryo.fits,'name')

        # header, table = printTable([fitNames,AICs,deltaAICs,weights,ks,ns],["fitNames","AICs","deltaAICs","weights","k","n"],col=True)

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecTask: Task handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def initTask(self):

        """Initializes task by connecting task and progressDialog to signals."""

        # Set Main GUI set disabled
        self.setDisabled(True)

        # Connect signals
        self.worker.taskFinished.connect(self.taskFinished)
        self.task.progressSignal.connect(self.updateProgressDialog)
        self.progressDialog.accepted.connect(self.taskCanceled)

        # Move worker and start
        self.worker.moveToThread(self.task)
        self.worker.start.emit()

    def updateProgressDialog(self, n):

        """Updates ProgressBar."""

        self.progressDialog.progressbar.setValue(n)

    def taskFinished(self):

        """Launched when task is finished. Sets Main GUI available again, closes all dialogs, updates ObjectBar."""

        # Quit thread
        self.task.quit()

        # Close dialog and make GUI available again
        self.progressDialog.close()
        self.statusBar().showMessage("Idle")
        self.setEnabled(True)

        # Updating analyzed column in ObjectBar
        currNode = self.getCurrentEmbryoNode()
        self.updateEmbryoNodeProps(currNode)

    def taskCanceled(self):

        """Launched when task is canceled. Sets Main GUI available again, closes all dialogs, updates ObjectBar."""

        self.statusBar().showMessage("Idle")
        self.setEnabled(True)

        self.task.terminate()

        # Map back backup
        self.originalObj = cpy.deepcopy(self.backupObj)
        self.backupObj = None

        self.progressDialog.close()

        # Updating analyzed column in ObjectBar
        currNode = self.getCurrentEmbryoNode()
        self.updateEmbryoNodeProps(currNode)

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecTab: PlotTab handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def createCanvas(
        self, parent, size=[1, 1], titles=None, sup=None, proj=None, tight=False
    ):

        """Create plot canvas."""

        h = 10000 / self.dpi
        v = 10000 / self.dpi
        self.fig = Figure(dpi=self.dpi)
        self.fig.set_size_inches(h, v, forward=True)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.currTab)

        self.fig, self.axes = pyfrp_plot_module.makeSubplot(
            size, titles=titles, tight=tight, sup=sup, proj=proj, fig=self.fig
        )

        self.ax = self.axes[0]

        return self.fig, self.canvas, self.ax

    def adjustAxes(self, plotType):

        """Adjusts axes spacing and labels."""

        if plotType == "bar":
            self.fig.subplots_adjust(bottom=0.3)

        elif plotType == "xyz":
            for ax in self.axes:
                ax.set_xlabel("x (um)")
                ax.set_ylabel("y (um)")
                ax.set_zlabel("z (um)")

        elif plotType == "intensityTS":
            self.fig.subplots_adjust(right=0.75)
            for ax in self.axes:
                ax.set_xlabel("Time (s)")
                ax.set_ylabel("Intensity (AU)")

        elif plotType == "meshDensity":
            c = ["x", "y", "z", ""]
            for i, ax in enumerate(self.axes):
                ax.set_xlabel(c[i])
                ax.set_ylabel("cell Volume")

        return

    def createPlotTab(
        self,
        plotType,
        plotName="",
        size=[1, 1],
        titles=None,
        sup=None,
        proj=None,
        tight=False,
    ):

        """Creates new plotting tab for matplotlib plots.
		"""

        # Grab embryo (might need to change this if whole molecule plots)
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            return

        # Create tab name
        tabname = currEmbryo.name + "/" + plotName + "#1"

        # Increment tab counter in name
        for i in range(self.plotTabs.count()):
            if self.plotTabs.tabText(i) == tabname:
                nme, nmbr = tabname.split("#")
                nmbr = str(int(nmbr) + 1)
                tabname = nme + "#" + nmbr

        # Add new Tab
        self.currTab = QtWidgets.QWidget()
        self.plotTabs.addTab(self.currTab, tabname)

        self.createCanvas(
            self.currTab, size=size, titles=titles, sup=sup, proj=proj, tight=tight
        )

        self.currTab.typ = plotType

        self.adjustAxes(plotType)

        # Append for bookkeeping
        self.tabAxes.append(self.ax)
        self.tabFigs.append(self.fig)

        # Update Canvas
        self.adjustCanvas()

        # Check for dummy plot tab
        if self.plotTabs.tabText(0) == "PlotTab":

            self.plotTabs.removeTab(self.plotTabs.currentIndex())
            self.plotTabs.setTabsClosable(True)

        self.plotTabs.setCurrentWidget(self.currTab)

    def currentTabChanged(self, value):

        """Callback if current tab is changed.
		
		Makes sure that widgets still have proper size by calling :py:func:`adjustCanvas`.
		
		Args:
			value (int): Index of new tab.
		
		"""

        self.currTab = self.plotTabs.widget(value)
        if len(self.tabAxes) > 0:

            self.ax = self.tabAxes[value]
            self.fig = self.tabFigs[value]

        if self.currTab != None:
            self.verticalSplitter.refresh()
            self.horizontalSplitter.refresh()
            self.currTab.setHidden(True)
            self.currTab.setVisible(True)

            self.adjustCanvas()

    def currentTabClosed(self, value):

        """Callback if the current tab is closed.
		
		Removes widgets from bookkeeping lists and creates new dummy tab if there is no tab left.
		"""

        self.currTab = self.plotTabs.widget(value)
        if len(self.tabAxes) > 0:
            self.tabAxes.pop(value)
            self.tabFigs.pop(value)

        self.plotTabs.removeTab(value)

        if self.plotTabs.count() == 0:
            self.createDummpyTab()

    def createDummpyTab(self):

        """Creates dummy tab and adds it to plot tabs."""

        self.currTab = QtWidgets.QWidget()
        self.firstTab = self.plotTabs.addTab(self.currTab, "PlotTab")
        self.plotTabs.setTabsClosable(False)

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecVTK: vtk handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def createVtkTab(self, plotName="", addRenderer=False):

        """Creates VTK tab for plotting.
		
		Creates new tab and adds it to plotTabs.
		
		Keyword Args:
			plotName (str): Name of plot displayed on tab.
			addRenderer (bool): Adds renderer.
			
		Returns:
			vtk.QVTKRenderWindowInteractor
		
		"""

        # Grab embryo (might need to change this if whole molecule plots)
        currEmbryo = self.getCurrentEmbryo()
        if currEmbryo == None:
            QtWidgets.QMessageBox.critical(
                None,
                "Error",
                "No embryo selected.",
                QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Default,
            )
            return

        # Create tab name
        tabname = currEmbryo.name + "/" + plotName + "#1"

        # Increment tab counter in name
        for i in range(self.plotTabs.count()):
            if self.plotTabs.tabText(i) == tabname:
                nme, nmbr = tabname.split("#")
                nmbr = str(int(nmbr) + 1)
                tabname = nme + "#" + nmbr

        # New empty widget
        self.currTab = QtWidgets.QWidget()

        # Create canvas
        self.createVtkCanvas(addRenderer=addRenderer)

        # Add to bookkeeping lists
        self.tabAxes.append(self.vtkCanvas)
        self.tabFigs.append(None)

        # Add tab
        self.plotTabs.addTab(self.currTab, tabname)

        # Check for dummy plot tab
        if self.plotTabs.tabText(0) == "PlotTab":
            self.plotTabs.removeTab(self.plotTabs.currentIndex())
            self.plotTabs.setTabsClosable(True)

        # Set as current tab
        self.plotTabs.setCurrentWidget(self.currTab)

        # Adjust size
        self.adjustCanvas()

        return self.vtkCanvas

    def createVtkCanvas(self, addRenderer=False):

        """Create vtk canvas.
		
		Keyword Args:
			addRenderer (bool): Adds renderer.
		
		Returns:
			vtk.QVTKRenderWindowInteractor
		
		"""

        self.vtkCanvas = QVTKRenderWindowInteractor(self.currTab)
        self.vtkCanvas.UpdateSize(int(self.currTab.width()), int(self.currTab.height()))

        self.vboxVTKCanvas = QtWidgets.QVBoxLayout()
        self.vboxVTKCanvas.addWidget(self.vtkCanvas, stretch=1)

        self.currTab.setLayout(self.vboxVTKCanvas)

        if addRenderer:
            self.renderer = vtk.vtkRenderer()
            self.vtkCanvas.GetRenderWindow().AddRenderer(self.renderer)

        return self.vtkCanvas

    def adjustCanvas(self):

        """Adjust canvas if slider/splitter changes to currently maximum size. 
		"""

        if hasattr(self, "currTab"):
            h = int(self.horizontalSplitter.sizes()[1])
            v = int(self.verticalSplitter.sizes()[0])
            self.currTab.resize(h, v)
        else:
            return

        if not hasattr(self, "fig"):
            return

        if self.fig != None:

            h = float(self.horizontalSplitter.sizes()[1]) / float(self.dpi)
            v = float(self.verticalSplitter.sizes()[0]) / float(self.dpi)

            if hasattr(self.currTab, "currSlider"):

                hSlider = float(self.currTab.currSlider.size().width()) / float(
                    self.dpi
                )
                vSlider = float(self.currTab.currSlider.size().height()) / float(
                    self.dpi
                )

                v = v - 5 * vSlider

            self.fig.set_size_inches(h, v, forward=False)

            self.canvas.draw()
            return

        elif self.fig == None:

            h = int(self.horizontalSplitter.sizes()[1])
            v = int(self.verticalSplitter.sizes()[0])

            self.vtkCanvas.UpdateSize(h, v)
            return

    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # SecSettings: Settings handling
    # ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    def setPath(self, identifier="", path=""):

        """Opens path setting dialog and afterwards sets path in path file.
		
		Keyword Args:
			identifier (str): Identifier displayed at dialog startup.
			path (str): Path displayed at dialog startup.
				
		"""

        pathDialog = pyfrp_gui_settings_dialogs.pathDialog(identifier, path, self)
        if pathDialog.exec_():
            identifier, path = pathDialog.getPath()

        pyfrp_misc_module.setPath(identifier, path)

        return 1

    def setGmshPath(self):

        """Opens path setting dialog and lets user set path to gmsh binary.
		"""

        path = pyfrp_misc_module.getPath("gmshBin")
        self.setPath(identifier="gmshBin", path=path)

        return 1

    def setFijiPath(self):

        """Opens path setting dialog and lets user set path to gmsh binary.
		"""

        path = pyfrp_misc_module.getPath("fijiBin")
        self.setPath(identifier="fijiBin", path=path)

        return 1

    def setOpenscadPath(self):

        """Opens path setting dialog and lets user set path to gmsh binary.
		"""

        path = pyfrp_misc_module.getPath("openscadBin")
        self.setPath(identifier="openscadBin", path=path)

        return 1

    def printPathFile(self):

        """Prints out path file."""

        print(pyfrp_misc_module.getPathFile())

        return 1

    def printPaths(self):

        """Prints out all content of path file."""

        pyfrp_misc_module.printPaths()

        return 1

    def checkPaths(self):

        """Checks if all paths in path file exist."""

        pyfrp_misc_module.checkPaths()

        return 1

    def setPathFile(self):

        """Allows setting the path file."""

        fn = str(
            QtWidgets.QFileDialog.getOpenFileName(
                self, "Choose path file", pyfrp_misc_module.getConfDir()
            )
        )
        if fn == "":
            return 0

        self.config.setPathFile(fn)

        return 1


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Main function calling PyFRAP main GUI
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


def main():

    # Creating application
    # font=QtGui.QFont()
    app = QtWidgets.QApplication(sys.argv)
    font = app.font()
    font.setPointSize(12)
    app.setFont(font)

    # Check if stout/sterr should be redirected
    try:
        print("redirect = %s" % sys.argv[1])
        redirect = bool(int(sys.argv[1]))
        # print(redirect)
    except IndexError:
        redirect = True

    try:
        print("overwrite_conf = %s" % sys.argv[2])
        overwrite_conf = bool(int(sys.argv[2]))
        # print(redirect)
    except IndexError:
        overwrite_conf = False

    # Start javabridge
    javabridge.start_vm(class_path=bioformats.JARS)

    mainWin = pyfrp(redirect=redirect, overwrite_conf=overwrite_conf)
    mainWin.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
