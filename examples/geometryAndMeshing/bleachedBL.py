"""Example for PyFRAP on how generate BL meshes.

Draws mesh and saves it if filename is given.

USAGE: python bleachedBL.py outputFilePath

"""

#Import necessary modules
from pyfrp.subclasses import pyfrp_embryo
from pyfrp.modules import pyfrp_misc_module
from pyfrp.modules import pyfrp_gmsh_IO_module

from pyfrp.modules.pyfrp_term_module import *

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import cm

import csv
import sys
import os

#--------------------------------------------------------
# Define parameters
#--------------------------------------------------------

# Imaging depth
sliceDepth=36.
sliceWidth=5. 

# Resolution
dataResMu=566.79
dataResPx=512

# Define bleached square (assuming its centered)
sidelength=2*73.71

#Try to make rim 100px wide
#rimFactor=0.8
rimFactor=1.-101./316.74

# Define cone geometry
cylinderHeight=90.3332804037
coneUpperRadius=635.3/2.
coneLowerRadius=448.5/2.

# Slice
center=[249.00918273645539, 263.82920110192839] 
imagingRadius=216.742098118
 
#--------------------------------------------------------
# Embryo setup
#--------------------------------------------------------

# Create embryo
emb=pyfrp_embryo.embryo("Test")

# Set experimental details
emb.setDataResMu(dataResMu)
emb.setSliceDepthMu(sliceDepth)
emb.offsetBleachedPx=[center[0]-sidelength/2.,center[1]-sidelength/2.]
	
# Geometry
emb.setGeometry2Cone(center,coneUpperRadius,coneLowerRadius,cylinderHeight)
		
# Update geometry properties in geo-file
emb.geometry.updateGeoFile()

# Create default ROIs
emb.genDefaultROIs(emb.geometry.getCenter(),imagingRadius,rimFactor=rimFactor,sliceHeightPx=-sliceDepth)
emb.newAnalysis()
	
# Add simulation and mesh
sim=emb.newSimulation()

# Generate original mesh
sim.mesh.genMesh()

# Genereate BL mesh around All Square
r=emb.getROIByName("All Square")	
fnOut=os.path.dirname(emb.geometry.fnGeo)+"/BL.geo"

sim.mesh.addBoundaryLayerAroundROI(r,fnOut=fnOut,simplify=True,iterations=3,triangIterations=0,
			fixSurfaces=False,debug=False,volSizePx=35.,volSizeLayer=10.,thickness=20.,cleanUp=True,faces=['x','y'])

# Show mesh
sim.mesh.plotMesh()

# Save mesh
if len(sys.argv)>1:
	sim.mesh.saveMeshToImg(sys.argv[1])

raw_input("Press to quit")
