"""Example for PyFRAP on how generate cylindrical geometry.

Draws geometry and saves it if filename is given.

USAGE: python createCylinderGeometry.py outputFilePath

"""


from pyfrp.modules import pyfrp_gmsh_geometry
import sys

# Create domain
d = pyfrp_gmsh_geometry.domain()

	
# Some parameters
center=[0,0]
radius=1.
height=50.
zOffset=100.
volSize=30.

# Add cylinder
d.addCylinderByParameters(center,radius,zOffset,height,volSize,plane="z",genLoops=True,genSurfaces=True,genVol=True)

# Draw
d.draw()

# Draw
ax=d.draw()
if len(sys.argv)>0:
	ax.get_figure().savefig(sys.argv[1])

raw_input("Press to quit")
