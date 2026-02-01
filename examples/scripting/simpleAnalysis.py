"""Example script that outlines how to analyze a FRAP dataset with PyFRAP
using scripts.

(1) Creates embryo object.
(2) Analyzes image data.
(3) Runs simulation.
(4) Fits data.
(5) Shows result output.
(6) If an output filename is given, will save embryo file to it.

Run as follows:

python simpleAnalysis.py path/to/recovery/data/ outputFilePath

"""

# Import modules
import os
import sys
import time
from pyfrp.subclasses import pyfrp_embryo
from pyfrp.modules import pyfrp_gmsh_geometry

# Data parameters
frameinterval=1. 
dataResMu=512.
sliceRadius=300.
sidelength=140.
sliceDepth=15.
sliceWidth=5 
center=[256,256]
rimFactor=0.8

# Simulation parameters
D=50
steps=1001
tEnd=1500

# Start timer
starTime=time.time()

# Create embryo
emb=pyfrp_embryo.embryo("Demo")

# Define where data is
emb.setDataFolder(sys.argv[1])




# Set data parameters
emb.setDataResMu(dataResMu)
emb.setFrameInterval(frameinterval)

# Set experimental details
emb.setDataResMu(dataResMu)
emb.setSideLengthBleachedMu(sidelength)
emb.setSliceDepthMu(sliceDepth)

# Create a .geo file of a 2 dimensional domain
d=pyfrp_gmsh_geometry.domain()
d.addCircleByParameters(center,sliceRadius,-sliceDepth,10,genLoop=True,genSurface=True)
d.writeToFile("demo.geo")

# Define geometry
emb.setGeometry2Custom(center,fnGeo='demo.geo',dim=2)

#Create default ROIs
emb.genDefaultROIs(emb.geometry.getCenter(),sliceRadius,rimFactor=rimFactor)

# Create analysis
emb.newAnalysis()

# Create Simulation
sim=emb.newSimulation()

# Generate mesh
sim.mesh.genMesh(debug=True)

# Simulation settings
emb.simulation.setD(D)
emb.simulation.setTimesteps(steps)
emb.simulation.setTEnd(tEnd)
emb.simulation.toLogTimeScale()

# Create a fit
fit=emb.newFit('Testfit')

# Fit settings
fit.addROIByName('Bleached Square')
fit.addROIByName('Slice')

fit.setEqu(True)
fit.setFitPinned(True)
fit.setOptMeth('Constrained Nelder-Mead')

# Run everything
emb.quickAnalysis()

# Show results
fit.plotFit()
fit.printResults()

stopTime=time.time()

# Print time it took
print("Analysis done in ", stopTime-starTime , " seconds.")

try:
	emb.save(sys.argv[2])
except IndexError:
	print("Did not save embryo, no output given.")

input()
