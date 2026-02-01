#! /bin/bash

# Install everything that is in apt
sudo apt-get install -y git
sudo apt-get install default-jdk
sudo apt-get install default-jre
sudo apt-get install -y python3-numpy
sudo apt-get install -y python3-scipy
sudo apt-get install -y python3-matplotlib
sudo apt-get install -y python3-pyqt5
sudo apt-get install -y python3-skimage
sudo apt-get install -y python3-vtk

# Install everything that is in pip
sudo pip3 install fipy
sudo pip3 install solidpython
sudo pip3 install numpy-stl
sudo pip3 install wget
sudo pip3 install python-bioformats

# Install openscad
sudo apt-get install openscad

# Clone PyFRAP and install
# cd PyFRAP-master
# sudo python setup.py install

 


