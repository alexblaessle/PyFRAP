.. PyFRAP documentation master file, created by
   sphinx-quickstart on Fri Mar 25 21:59:01 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyFRAP's API!
==================================

PyFRAP is an extensive Python based FRAP (Fluorescence Recovery After Photobleaching) analysis software, featuring various tools that help
analyzing FRAP datasets, such as

* Import FRAP datasets from timelapse experiments and analyze image data with various options such as
	* various filters
	* background substraction
	* illumination correction
* Simulate the FRAP experiment with exact interpolated initial conditions
* Fit simulated experiment to analyzed data and extract diffusion coefficient
* Statistical analysis of fitting results
* Hierarchical data structure making data exchange/sharing easy
* Comprehensive GUI, making almost all PyFRAP tools available

I have tried to keep the API short but clear. If it is unclear, don't hesitate to mail.

More information about PyFRAP can be found on  `Github <https://github.com/alexblaessle/PyFRAP>`_ .


Installation (setup.py)
-----------------------

.. toctree::
   :maxdepth: 2
   
   setup.rst

The modules package
-------------------

.. toctree::
   :maxdepth: 2
   
   pyfrp.modules.rst
   

The subclasses package
----------------------
   
.. toctree::
   :maxdepth: 2
   
   pyfrp.subclasses.rst
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

