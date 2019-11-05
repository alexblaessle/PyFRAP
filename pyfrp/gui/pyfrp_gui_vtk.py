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

# Module containing vtk visualizing classes:


# ===========================================================================================================================================================================
# Importing necessary modules
# ===========================================================================================================================================================================

# Misc
import sys
import time
import os, os.path

# PyFRAP modules
from pyfrp.modules.pyfrp_term_module import *

# PyQT
from PyQt5 import QtGui, QtCore, QtWidgets

# VTK
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor


class vtkSimVisualizer(QtWidgets.QMainWindow):

    """Basic simulation visualizing class.
	
	.. note:: Only works if the simulation results 
	   are saved in embryo via ``embryo.simulation.saveSim=True``.
	
	.. image:: ../imgs/pyfrp_gui_vtk/vtkSimVisualizer.png
		
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): An embryo
	
	"""

    def __init__(self, embryo, parent=None):

        super(vtkSimVisualizer, self).__init__(parent)

        # Load embryo
        self.embryo = embryo

        # Check if simulation is saved
        if len(self.embryo.simulation.vals) == 0:
            printWarning(
                "Embryo does not have saved simulation. Rerun simulation with saveSim=True."
            )
            return

        # To numpy array
        self.embryo.simulation.vals = np.array(self.embryo.simulation.vals)

        # Frame
        self.frame = QtWidgets.QFrame()

        # Vtk widgets
        self.vtkWidget, self.ren, self.iren = self.initVTKWidget(self.frame)

        # Add slider for 3D
        self.slider, self.label = self.initSlider(
            0, len(self.embryo.simulation.tvecSim) - 1, self.sliderCallback
        )

        # Show simulation visualization
        self.initSimVis()

        # Add Layout
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.vtkWidget, 1, 1)
        self.layout.addWidget(self.slider, 2, 1)
        self.layout.addWidget(self.label, 3, 1)

        # Set central widget
        self.frame.setLayout(self.layout)
        self.setCentralWidget(self.frame)

        # Window title
        self.setWindowTitle("vtkSimVisualizer of embryo " + self.embryo.name)

        # Show everything
        self.show()
        self.iren.Initialize()

    def initSlider(self, sMin, sMax, callback):

        """Sets up slider and corresponding label with given limits and callback.
		
		Args:
			sMin (int): Minimum index of slider.
			sMax (int): Maximum index of slider.
			callback (function): Callback function of slider.
			
		Returns:
			tuple: Tuple containing:
			
				* slider (QtWidgets.QSlider)
				* lbl (QtWidgets.QLabel)
		"""

        slider = QtWidgets.QSlider()
        slider.setRange(sMin, sMax)
        slider.valueChanged.connect(callback)
        slider.setOrientation(QtCore.Qt.Horizontal)
        lbl = QtWidgets.QLabel("")
        return slider, lbl

    def initVTKWidget(self, frame):

        """Sets up vtkWidget inside frame.
		
		Also sets up corresponding renderer and interactor.
		
		Args:
			frame (QtWidgets.QFrame): Parenting frame.
			
		Returns:
			tuple: Tuple containing:
			
				* vtkWidget (vtk.qt.QVTKRenderWindowInteractor.QVTKRenderWindowInteractor): Qt vtk window interactor.
				* ren (vtk.vtkRenderer)
				* iren (vtk.vtkRenderWindowInteractor)
				
		"""

        vtkWidget = QVTKRenderWindowInteractor(frame)
        ren = vtk.vtkRenderer()
        vtkWidget.GetRenderWindow().AddRenderer(ren)
        iren = vtkWidget.GetRenderWindow().GetInteractor()
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        ren.ResetCamera()

        return vtkWidget, ren, iren

    def showVals(self, idx):

        """Updates simulation values in renderer.
		
		Args:
			idx (int): Index of timepoint.
		
		"""

        # Convert to right cmap
        vals = self.embryo.simulation.vals[idx] / self.embryo.simulation.vals.max()
        vals = 255 * self.cm(vals)
        vals = vals[:, :3]

        self.initScalarArray()
        for i, val in enumerate(vals):
            self.scalar.InsertNextTupleValue(val)
        self.grid.GetCellData().SetScalars(self.scalar)
        self.grid.Modified()
        self.iren.Render()

        try:
            self.CutMesh()
        except AttributeError:
            pass

    def initSimVis(self):

        """Inits the simulation renderer.
		
		This contains the following steps:
		
			* Initializing points.
			* Initializing grid ontop of points.
			* Initialzing tetrahedras.
			* Initializing a scalar array.
			* Initializing a mapper and actor.
			* Initializing a color table.
			* Shows initial time point in renderer.
			
		"""

        self.initPoints()
        self.initGrid()
        if self.embryo.geometry.dim == 2:
            self.initTriangles()
        else:
            self.initTetras()
        self.initScalarArray()
        self.initMapper()
        self.initActor()
        self.initColorTable()
        self.showVals(0)

    def initActor(self):

        """Sets up simulation actor."""

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.ren.AddActor(self.actor)

    def initMapper(self):

        """Sets up simulation mapper."""

        self.mapper = vtk.vtkDataSetMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            self.mapper.SetInputConnection(self.grid.GetProducerPort())
        else:
            self.mapper.SetInputData(self.grid)

    def initGrid(self):

        """Sets up simulation grid."""

        self.grid = vtk.vtkUnstructuredGrid()
        self.grid.SetPoints(self.points)

        return self.grid

    def initTetras(self):

        """Sets up tetrahedras describing mesh."""

        verts = self.embryo.simulation.mesh.mesh._getOrderedCellVertexIDs().T

        cellArray = vtk.vtkCellArray()
        for j, vert in enumerate(verts):

            tetra = vtk.vtkTetra()

            for i, v in enumerate(vert):
                tetra.GetPointIds().SetId(i, v)
            cellArray.InsertNextCell(tetra)

        self.grid.SetCells(vtk.VTK_TETRA, cellArray)

    def initTriangles(self):

        """Sets up triangles describing mesh."""

        verts = self.embryo.simulation.mesh.mesh._getOrderedCellVertexIDs().T

        cellArray = vtk.vtkCellArray()
        for j, vert in enumerate(verts):

            tetra = vtk.vtkTriangle()

            for i, v in enumerate(vert):
                tetra.GetPointIds().SetId(i, v)
            cellArray.InsertNextCell(tetra)

        self.grid.SetCells(vtk.VTK_TRIANGLE, cellArray)

    def initScalarArray(self):

        """Sets up scalar array describing tetrahedra values."""

        self.scalar = vtk.vtkUnsignedCharArray()
        self.scalar.SetNumberOfComponents(3)
        self.scalar.SetName("Colors")

        return self.scalar

    def initPoints(self):

        """Sets up point array."""

        if self.embryo.geometry.dim == 2:
            x, y = self.embryo.simulation.mesh.mesh.vertexCoords
        else:
            x, y, z = self.embryo.simulation.mesh.mesh.vertexCoords

        self.points = vtk.vtkPoints()
        for i in range(len(x)):
            if self.embryo.geometry.dim == 2:
                self.points.InsertNextPoint(x[i], y[i], 0)
            else:
                self.points.InsertNextPoint(x[i], y[i], z[i])

        return self.points

    def getRenderer(self):

        """Returns renderer."""

        return self.ren

    def sliderCallback(self):

        """Call back function for slider movement."""

        index = self.sender().value()
        self.label.setText("t = " + str(self.embryo.simulation.tvecSim[index]))
        self.showVals(index)
        return

    def initColorTable(self):

        """Sets up color table."""

        from matplotlib import cm

        self.cm = cm.jet


class vtkSimVisualizerCutter(QtWidgets.QMainWindow):

    """Simulation visualizing class that allows to draw the simulation in 3D and slice a plane 
	out of the 3D simulation.
	
	.. note:: Only works for 3D embryos that already have been simulated and the simulation results 
	   are saved in embryo via ``embryo.simulation.saveSim=True``.
	
	.. image:: ../imgs/pyfrp_gui_vtk/vtkSimVisualizerCutter.png
		
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): An embryo
	
	"""

    def __init__(self, embryo, parent=None):

        super(vtkSimVisualizerCutter, self).__init__(parent)

        # Load embryo
        self.embryo = embryo

        # Check if simulation is saved
        if len(self.embryo.simulation.vals) == 0:
            printWarning(
                "Embryo does not have saved simulation. Rerun simulation with saveSim=True."
            )

        # Check if 3D embryo
        if self.embryo.geometry.dim == 2:
            printError(
                "vtkSimVisualizerCutter does only work and make sense for 3D geometries. Use vtkSimVisualizer instead."
            )

        # To numpy array
        self.embryo.simulation.vals = np.array(self.embryo.simulation.vals)

        # Frame
        self.frame = QtWidgets.QFrame()

        # Vtk widgets
        self.vtkWidget, self.ren, self.iren = self.initVTKWidget(self.frame)
        self.vtkWidgetCut, self.renCut, self.irenCut = self.initVTKWidget(self.frame)

        # Add slider for 3D
        self.slider, self.label = self.initSlider(
            0, len(self.embryo.simulation.tvecSim) - 1, self.sliderCallback
        )

        # Setup 3D simulation plot
        self.init3D()

        # Plane used for cutting
        self.plane = vtk.vtkPlane()

        # Cutter
        self.cutter = vtk.vtkCutter()

        # Plane widget
        self.initPlaneWidget()

        # Add Layout
        self.layout = QtWidgets.QGridLayout()
        self.layout.addWidget(self.vtkWidget, 1, 1)
        self.layout.addWidget(self.vtkWidgetCut, 1, 2)
        self.layout.addWidget(self.slider, 2, 1)
        self.layout.addWidget(self.label, 3, 1)

        # Set central widget
        self.frame.setLayout(self.layout)
        self.setCentralWidget(self.frame)

        self.setWindowTitle("vtkSimVisualizerCutter of embryo " + self.embryo.name)

        # Run cut mesh once
        self.CutMesh()

        # Show everything
        self.show()
        self.iren.Initialize()
        self.irenCut.Initialize()

    def closeIt(self):
        self.close()

    def initSlider(self, sMin, sMax, callback):

        """Sets up slider and corresponding label with given limits and callback.
		
		Args:
			sMin (int): Minimum index of slider.
			sMax (int): Maximum index of slider.
			callback (function): Callback function of slider.
			
		Returns:
			tuple: Tuple containing:
			
				* slider (QtWidgets.QSlider)
				* lbl (QtWidgets.QLabel)
		"""

        slider = QtWidgets.QSlider()
        slider.setRange(sMin, sMax)
        slider.valueChanged.connect(callback)
        slider.setOrientation(QtCore.Qt.Horizontal)
        lbl = QtWidgets.QLabel("")
        return slider, lbl

    def initVTKWidget(self, frame):

        """Sets up vtkWidget inside frame.
		
		Also sets up corresponding renderer and interactor.
		
		Args:
			frame (QtWidgets.QFrame): Parenting frame.
			
		Returns:
			tuple: Tuple containing:
			
				* vtkWidget (vtk.qt.QVTKRenderWindowInteractor.QVTKRenderWindowInteractor): Qt vtk window interactor.
				* ren (vtk.vtkRenderer)
				* iren (vtk.vtkRenderWindowInteractor)
				
		"""

        vtkWidget = QVTKRenderWindowInteractor(frame)
        ren = vtk.vtkRenderer()
        vtkWidget.GetRenderWindow().AddRenderer(ren)
        iren = vtkWidget.GetRenderWindow().GetInteractor()
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        ren.ResetCamera()

        return vtkWidget, ren, iren

    def initPlaneWidget(self):

        """Sets up vtkImplicitPlaneWidget for plane selection.
		
		Sets also up necessary actor and mapper.
		
		Returns:
			vtk.vtkImplicitPlaneWidget: Widget.
		
		"""

        # Initialize a plane widget
        self.planeWidget = vtk.vtkImplicitPlaneWidget()
        self.planeWidget.SetInteractor(self.iren)
        self.planeWidget.SetPlaceFactor(1.25)
        self.planeWidget.SetInput(self.grid)
        self.planeWidget.PlaceWidget()

        # Callback connection of interaction event
        self.planeWidget.AddObserver("InteractionEvent", self.StartPlaneCallback)
        self.planeWidget.AddObserver("EndInteractionEvent", self.EndPlaneCallback)

        # Mapper
        self.selectMapper = vtk.vtkDataSetMapper()
        self.selectMapper.SetInput(self.grid)

        # Actor
        self.selectActor = vtk.vtkLODActor()
        self.selectActor.SetMapper(self.selectMapper)
        self.selectActor.GetProperty().SetColor(0, 1, 0)
        self.selectActor.VisibilityOff()
        self.selectActor.SetScale(1.01, 1.01, 1.01)

        self.ren.AddActor(self.selectActor)

        # Place widget at right position
        self.planeWidget.SetNormal(0, 0, 1)
        self.planeWidget.SetOrigin(self.grid.GetCenter())

        return self.planeWidget

    def StartPlaneCallback(self, obj, event):

        """Callback for interaction event.
		"""

        obj.GetPlane(self.plane)
        self.selectActor.VisibilityOn()

    def EndPlaneCallback(self, widget, event_string):

        """Callback for end of interaction event.
		"""

        self.CutMesh()

    def CutMesh(self):

        """Cuts plane out of 3D simulation object.
		
		Grabs plane and creates and actor and mapper to
		display cut out plane in second renderer.
		
		"""

        self.cutter.SetInput(self.grid)
        self.planeWidget.GetPlane(self.plane)
        self.cutter.SetCutFunction(self.plane)
        self.cutter.Update()

        self.Actor = vtk.vtkActor()

        mapper = vtk.vtkPolyDataMapper()

        mapper.SetInputConnection(self.cutter.GetOutputPort())
        self.Actor.SetMapper(mapper)
        self.renCut.AddActor(self.Actor)

        self.cutter.GetOutput().Modified()
        self.irenCut.Render()

    def showVals(self, idx):

        """Updates simulation values in renderer.
		
		Args:
			idx (int): Index of timepoint.
		
		"""

        # Convert to right cmap
        vals = self.embryo.simulation.vals[idx] / self.embryo.simulation.vals.max()
        vals = 255 * self.cm(vals)
        vals = vals[:, :3]

        self.initScalarArray()
        for i, val in enumerate(vals):
            self.scalar.InsertNextTupleValue(val)
        self.grid.GetCellData().SetScalars(self.scalar)
        self.grid.Modified()
        self.iren.Render()

        try:
            self.CutMesh()
        except AttributeError:
            pass

    def init3D(self):

        """Inits the simulation renderer.
		
		This contains the following steps:
		
			* Initializing points.
			* Initializing grid ontop of points.
			* Initialzing tetrahedras.
			* Initializing a scalar array.
			* Initializing a mapper and actor.
			* Initializing a color table.
			* Shows initial time point in renderer.
			
		"""

        self.initPoints()
        self.initGrid()
        self.initTetras()
        self.initScalarArray()
        self.initMapper()
        self.initActor()
        self.initColorTable()
        self.showVals(0)

    def initActor(self):

        """Sets up simulation actor."""

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        self.ren.AddActor(self.actor)

    def initMapper(self):

        """Sets up simulation mapper."""

        self.mapper = vtk.vtkDataSetMapper()
        if vtk.VTK_MAJOR_VERSION <= 5:
            self.mapper.SetInputConnection(self.grid.GetProducerPort())
        else:
            self.mapper.SetInputData(self.grid)

    def initGrid(self):

        """Sets up simulation grid."""

        self.grid = vtk.vtkUnstructuredGrid()
        self.grid.SetPoints(self.points)

        return self.grid

    def initTetras(self):

        """Sets up tetrahedras describing mesh."""

        verts = self.embryo.simulation.mesh.mesh._getOrderedCellVertexIDs().T

        cellArray = vtk.vtkCellArray()
        for j, vert in enumerate(verts):

            tetra = vtk.vtkTetra()

            for i, v in enumerate(vert):
                tetra.GetPointIds().SetId(i, v)
            cellArray.InsertNextCell(tetra)

        self.grid.SetCells(vtk.VTK_TETRA, cellArray)

    def initScalarArray(self):

        """Sets up scalar array describing tetrahedra values."""

        self.scalar = vtk.vtkUnsignedCharArray()
        self.scalar.SetNumberOfComponents(3)
        self.scalar.SetName("Colors")

        return self.scalar

    def initPoints(self):

        """Sets up point array."""

        x, y, z = self.embryo.simulation.mesh.mesh.vertexCoords

        self.points = vtk.vtkPoints()
        for i in range(len(x)):
            self.points.InsertNextPoint(x[i], y[i], z[i])

        return self.points

    def getRenderer(self):

        """Returns renderer."""

        return self.ren

    def sliderCallback(self):

        """Call back function for slider movement."""

        index = self.sender().value()
        self.label.setText("t = " + str(self.embryo.simulation.tvecSim[index]))
        self.showVals(index)
        return

    def initColorTable(self):

        """Sets up color table."""

        from matplotlib import cm

        self.cm = cm.jet


class vtkSimVisualizerCutterBySub(vtkSimVisualizer):

    """Simulation visualizing class that allows to draw the simulation in 3D and slice a plane 
	out of the 3D simulation.
	
	This one subclasses from vtkSimVisualizer. 
	
	.. warning:: Currently not working. Crashes at ``self.irenCut.Initialize()`` with ``segFault``. Currently no way to fix this, 
	   hence vtkSimVisualizerCutter is own instance instead of subclassing.
	
	Args:
		embryo (pyfrp.subclasses.pyfrp_embryo.embryo): An embryo
	
	"""

    def __init__(self, embryo, parent=None):

        # super(vtkSimVisualizerCutter, self).__init__(embryo,parent)
        vtkSimVisualizer.__init__(self, embryo, parent=parent)

        # Check if 3D embryo
        if self.embryo.geometry.dim == 2:
            printWarning(
                "vtkSimVisualizerCutter does only work and make sense for 3D geometries. Will not initialize cutter and fall back tostandad vtkSimVisualizer."
            )
            return

        self.vtkWidgetCut, self.renCut, self.irenCut = self.initVTKWidget(self.frame)

        # Plane used for cutting
        self.plane = vtk.vtkPlane()

        # Cutter
        self.cutter = vtk.vtkCutter()

        # Plane widget
        self.initPlaneWidget()

        # Add Layout
        self.layout.addWidget(self.vtkWidgetCut, 1, 2)

        # Run cut mesh once
        self.CutMesh()

        self.irenCut.Initialize()

    def initPlaneWidget(self):

        """Sets up vtkImplicitPlaneWidget for plane selection.
		
		Sets also up necessary actor and mapper.
		
		Returns:
			vtk.vtkImplicitPlaneWidget: Widget.
		
		"""

        # Initialize a plane widget
        self.planeWidget = vtk.vtkImplicitPlaneWidget()
        self.planeWidget.SetInteractor(self.iren)
        self.planeWidget.SetPlaceFactor(1.25)
        self.planeWidget.SetInput(self.grid)
        self.planeWidget.PlaceWidget()

        # Callback connection of interaction event
        self.planeWidget.AddObserver("InteractionEvent", self.StartPlaneCallback)
        self.planeWidget.AddObserver("EndInteractionEvent", self.EndPlaneCallback)

        # Mapper
        self.selectMapper = vtk.vtkDataSetMapper()
        self.selectMapper.SetInput(self.grid)

        # Actor
        self.selectActor = vtk.vtkLODActor()
        self.selectActor.SetMapper(self.selectMapper)
        self.selectActor.GetProperty().SetColor(0, 1, 0)
        self.selectActor.VisibilityOff()
        self.selectActor.SetScale(1.01, 1.01, 1.01)

        self.ren.AddActor(self.selectActor)

        # Place widget at right position
        self.planeWidget.SetNormal(0, 0, 1)
        self.planeWidget.SetOrigin(self.grid.GetCenter())

        return self.planeWidget

    def StartPlaneCallback(self, obj, event):

        """Callback for interaction event.
		"""

        obj.GetPlane(self.plane)
        self.selectActor.VisibilityOn()

    def EndPlaneCallback(self, widget, event_string):

        """Callback for end of interaction event.
		"""

        self.CutMesh()

    def CutMesh(self):

        """Cuts plane out of 3D simulation object.
		
		Grabs plane and creates and actor and mapper to
		display cut out plane in second renderer.
		
		"""

        self.cutter.SetInput(self.grid)
        self.planeWidget.GetPlane(self.plane)
        self.cutter.SetCutFunction(self.plane)
        self.cutter.Update()

        self.Actor = vtk.vtkActor()

        mapper = vtk.vtkPolyDataMapper()

        mapper.SetInputConnection(self.cutter.GetOutputPort())
        self.Actor.SetMapper(mapper)
        self.renCut.AddActor(self.Actor)

        self.cutter.GetOutput().Modified()
        self.irenCut.Render()
