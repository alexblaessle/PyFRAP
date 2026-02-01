"""Example for PyFRAP on how generate merged meshes.

Draws mesh and saves it if filename is given.

USAGE: python mergedROIMesh.py outputFilePath

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

# Grab slice and refine
fnOut=emb.geometry.fnGeo
fnMeshSlice=emb.getROIByName("Slice").genMeshFile(volSizePx=10)

print fnOut,fnMeshSlice

# Merge
fnGeoNew,fnMshNew=pyfrp_gmsh_IO_module.mergeMeshes([sim.mesh.fnMesh,fnMeshSlice],fnOut.replace(".geo","_merged.geo"))
sim.mesh.setFnMesh(fnMshNew)

# Show mesh
sim.mesh.plotMesh()

# Save mesh
if len(sys.argv)>0:
	sim.mesh.saveMeshToImg(sys.argv[1])

raw_input("Press to quit")


