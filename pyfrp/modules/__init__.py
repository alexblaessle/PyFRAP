#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFRAP: A Python based FRAP analysis tool box. Basic modules.
"""

import os
import sys
import platform

# Here we check if OS is OSX. If som we change matplotlib backend
# (We have to do this to make sure that pyplot import works in non-framework
# installations)
# This might cause warnings if pyplot is already imported.
	
if platform.system() in ["Darwin"]:
	import matplotlib 
	matplotlib.use('qt5agg')

#Basic PyFRAP modules
from . import pyfrp_term_module
from . import pyfrp_IO_module 
from . import pyfrp_misc_module
from . import pyfrp_plot_module
from . import pyfrp_vtk_module
from . import pyfrp_img_module
from . import pyfrp_optimization_module
from . import pyfrp_stats_module
from . import pyfrp_fit_module
from . import pyfrp_gmsh_module
from . import pyfrp_gmsh_IO_module
from . import pyfrp_sim_module
from . import pyfrp_integration_module
from . import pyfrp_idx_module
from . import pyfrp_geometry_module
from . import pyfrp_gmsh_geometry
from . import pyfrp_openscad_module

#Obsolete/Not-integrated modules  
#from . import pyfrp_zstack_module

