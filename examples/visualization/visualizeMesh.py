"""This script is an example how to load an embryo file and use PyFRAP's functions to visualize the mesh of the embryo.

Run as follows:

python visualizeMesh.py path/To/Embryo/File fromFile

fromFile can be 0 and 1.

fromFile controls if mesh plotting does extra step through writing into vtk file or not.

"""

from pyfrp.modules import pyfrp_IO_module
import sys

# Load embryo
emb=pyfrp_IO_module.loadFromPickle(sys.argv[1])

# Plot mesh
mesh=emb.simulation.mesh.plotMesh(fromFile=bool(int(sys.argv[2])))


