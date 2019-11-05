#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PyFRAP: A Python based FRAP analysis tool box. Gui package.
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

    matplotlib.use("qt5agg")

# PyFRAP main GUI
from . import pyfrp_gui_builder
from . import pyfrp_gui_vtk
from . import pyfrp_app

# from . import pyfrp_term
