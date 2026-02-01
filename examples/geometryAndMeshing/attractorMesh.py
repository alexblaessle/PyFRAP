"""Example for PyFRAP on how generate an attractor mesh in 2D.

USAGE: python attractorMesh.py outputFilePath

"""

# Import some modules
from pyfrp.modules import pyfrp_gmsh_geometry
from pyfrp.subclasses import pyfrp_embryo
import sys

# Create domain
d = pyfrp_gmsh_geometry.domain()
volSize=30
	
# Add circular domain
d.addCircleByParameters([256,256],300,0,volSize,genLoop=True,genSurface=True)

# Add attractor mesh around center
v=d.vertices[0]
v.addToAttractor(attrField=None,LcMin=5.,LcMax=20.,DistMin=30.,DistMax=60.)

# Write to file
d.writeToFile(sys.argv[1])

# Create embryo
emb=pyfrp_embryo.embryo("test")
		
# Set newly generated geofile as geometry
emb.setGeometry2Custom([256,256],fnGeo=sys.argv[1],dim=2)

# Create mesh
sim=emb.newSimulation()
sim.mesh.genMesh()

# Plot mesh
sim.mesh.plotMesh()

raw_input()