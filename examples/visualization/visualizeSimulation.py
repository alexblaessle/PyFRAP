"""This script is an example how to load an embryo file and use PyFRAP's functions to visualize the simulation.

Run as follows:

python visualizeSimulation.py path/To/Embryo/File cut

cut can be 0 and 1.

"""

from pyfrp.modules import pyfrp_IO_module
import sys

emb=pyfrp_IO_module.loadFromPickle(sys.argv[1])
emb.simulation.visualize(cut=bool(int(sys.argv[2])))
