"""Example for PyFRAP on how generate circular geometry.

Draws geometry and saves it if filename is given.

USAGE: python createCircularGeometry.py outputFilePath

"""


from pyfrp.modules import pyfrp_gmsh_geometry
import sys

# Create domain
d = pyfrp_gmsh_geometry.domain()
volSize=30.
	
# Add some vertices
vcenter=d.addVertex([0,0,0],volSize=volSize)
v1=d.addVertex([1,0,0],volSize=volSize)
v2=d.addVertex([0,1,0],volSize=volSize)
v3=d.addVertex([-1,0,0],volSize=volSize)
v4=d.addVertex([0,-1,0],volSize=volSize)

# Add Arcs
a1=d.addArc(v1,vcenter,v2)
a2=d.addArc(v2,vcenter,v3)
a3=d.addArc(v3,vcenter,v4)
a4=d.addArc(v4,vcenter,v1)

# Draw
ax=d.draw()
if len(sys.argv)>0:
	ax.get_figure().savefig(sys.argv[1])

raw_input("Press to quit")
