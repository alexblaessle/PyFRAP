#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFRAP: A Python based FRAP analysis tool box. Subclass module.
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

#PyFRAP classes
from . import pyfrp_analysis
from . import pyfrp_conf
from . import pyfrp_fit
from . import pyfrp_mesh
from . import pyfrp_simulation
from . import pyfrp_geometry
from . import pyfrp_ROI
from . import pyfrp_embryo
from . import pyfrp_molecule
