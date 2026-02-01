# PyFRAP: A Python based FRAP analysis tool box. 

Fluorescence Recovery After Photobleaching (FRAP) is a common way to assess molecular diffusion. PyFRAP is a novel simulation-based analysis software
that makes use of PDE simulations to analyze FRAP experiments in 3D geometries. It uses the first post-bleach image as initial condition, making assumptions 
about the underlying initial conditions obsolete. PyFRAP can fit different reaction-diffusion models to FRAP data, providing quantitative information about
effective diffusion.

## Features

PyFRAP comes with a full image analysis and simulation toolbox. In particular, PyFRAP can

- Import FRAP datasets from timelapse experiments and analyze image data with various options such as
	+ image filters
	+ background substraction
	+ illumination correction
- Simulate the FRAP experiment with exact interpolated initial conditions
- Fit simulated experiment to analyzed data and extract diffusion coefficient
- Statistical analysis of fitting results
- Hierarchical data structure making data exchange/sharing easy
- Comprehensive GUI, making almost all PyFRAP tools available

## Installation

PyFRAP can be installed in different ways. We provide installation scripts that allow an easy installation of PyFRAP in combination with Anaconda, including all necessary
Python packages and external softwares. Installation
instructions can be found [here](https://github.com/alexblaessle/PyFRAP/wiki/Installation#short).

If you are familiar with Python and git, you can install PyFRAP via 

	git clone https://github.com/alexblaessle/PyFRAP
	
and:

	python setup.py install --user
	
We highly recommend installing with the *--user* option, since PyFRAP needs to read/write data files in the installation folder. In some cases, this
might lead to file permission issues.
For a full installation documentation, have a look at the [wiki](https://github.com/alexblaessle/PyFRAP/wiki/Installation).

### Requirements

PyFRAP depends on 

- numpy>=1.8.2
- matplotlib>=1.4.3
- scipy>=0.13.3
- scikit-image>=0.11.3
- fipy>=3.1
- PyQT4>4.10.4
- vtk>=5.8.0
- colorama>=0.2.5
- wget>=3.2
- gmsh (compiled with TetGen Algorithm) MUST BE Version 2.14.0!

Note that the installation described [here](https://github.com/alexblaessle/PyFRAP/wiki/Installation#short) installs all necessary packages.

## Getting Started

### Running the PyFRAP GUI

PyFRAP comes with a comprehensive GUI. It can easily be started by clicking on *runPyFRAP.bat* (Windows), *runPyFRAP.command* (OSX) or *runPyFRAP.sh* (Linux). 

If you are using a terminal, *cd* to your PyFRAP directory and type

	python pyfrp/PyFRAP.py

If you are already in a python session, you can simply run 

	import pyfrp
	pyfrp.main()

Note that in the latter method PyFRAP's stdout might get redirected to the python shell you are executing it from.

### Running PyFRAP from the command line

PyFRAP is a complete python package and can be imported via

```python
	import pyfrp
```

Note that PyFRAP has three main submodules: *pyfrp.modules*, *pyfrp.subclasses* and *pyfrp.gui*. If you do not want to use any GUI elements, we recommend only importing the modules you need.

### Using the PyFRAP GUI to analyze a FRAP experiment

Check out the PyFRAP wiki's [First Steps Section](https://github.com/alexblaessle/PyFRAP/wiki/FirstSteps).

## API

PyFRAP is fully documented, allowing easy creation of scripts and extensions of the PyFRAP toolbox. The API of PyFRAP can be found [here](http://pyfrap.readthedocs.org/en/latest/) .

## Documentation

To learn more about PyFRAP, check out the PyFRAP [wiki](https://github.com/alexblaessle/PyFRAP/wiki)





