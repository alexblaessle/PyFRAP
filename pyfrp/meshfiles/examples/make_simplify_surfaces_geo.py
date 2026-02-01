"""Script to build simple geometry to test simplify surfaces.

This file can be used to test the simplifySurfaces method 
of the pyfrp.modules.pyfrp_gmsh_geometry.domain module. 

See also pyfrp/scripts/geomeotry/simplifySurfaces.py

"""

# Import modules
from pyfrp.modules import pyfrp_gmsh_geometry

# Create geometry
d=pyfrp_gmsh_geometry.domain()

# Add cuboid
d.addCuboidByParameters([0,0,0],100,150,50,30.,plane="z",genLoops=False,genSurfaces=False,genVol=False)

# Add 4 lines cutting sides into triangles
d.addLine(d.getVertexById(5)[0],d.getVertexById(2)[0])
d.addLine(d.getVertexById(6)[0],d.getVertexById(3)[0])
d.addLine(d.getVertexById(7)[0],d.getVertexById(4)[0])
d.addLine(d.getVertexById(8)[0],d.getVertexById(1)[0])

# Add loops
d.addLineLoop(edgeIDs=[5,10,13])
d.addLineLoop(edgeIDs=[1,9,13])

d.addLineLoop(edgeIDs=[6,11,14])
d.addLineLoop(edgeIDs=[2,10,14])

d.addLineLoop(edgeIDs=[7,12,15])
d.addLineLoop(edgeIDs=[3,11,15])

d.addLineLoop(edgeIDs=[8,9,16])
d.addLineLoop(edgeIDs=[4,12,16])

d.addLineLoop(edgeIDs=[8,7,6,5])
d.addLineLoop(edgeIDs=[4,3,2,1])

# Make sure loops are proper
for loop in d.lineLoops:
	loop.fix()

# Add Surfaces
for loop in d.lineLoops:
	d.addRuledSurface(lineLoopID=loop.Id)

# Write to file
d.writeToFile("simplifySurfaces.geo")

