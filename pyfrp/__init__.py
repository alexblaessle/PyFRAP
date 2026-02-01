#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFRAP: A Python based FRAP analysis tool box
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
from . import modules

#PyFRAP classes
from . import subclasses

#PyFRAP GUI classes

#Only import if not RTD or not selected. RTD is currently having problems 
#with importing GUI classes. Will need to have fix for this 
#at some point

if os.environ.get('READTHEDOCS', None) == 'True':
	pass
else:
	#from . import gui
	from .gui.pyfrp_app import main
	
	
		
	
__version__ = '1.1'
__author__ = u"Alexander Blaessle"
__license__ = "GNU GPL v3"