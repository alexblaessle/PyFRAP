REM This script installs all dependencies of PyFRAP for OSX, assuming Anaconda is installed.

REM Downgrade to PyQt4
conda remove -y qt pyqt
conda install -y qt=4 pyqt=4.10.4

REM Reinstall MPL
conda remove -y pywavelets
conda install -y matplotlib=1.4.3

REM Reinstall skimage
conda install -y scikit-image=0.11.3

REM Install extra packages
conda install -y vtk
conda install -y colorama

REM Install extra packages using 
pip install numpy-stl
pip install solidpython
pip install wget
pip install pysparse
pip install fipy
pip install python-bioformats

REM Install PyFRAP
python setup.py install 
