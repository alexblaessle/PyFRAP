#=====================================================================================================================================
#Copyright
#=====================================================================================================================================

#Copyright (C) 2014 Alexander Blaessle, Patrick Mueller and the Friedrich Miescher Laboratory of the Max Planck Society
#This software is distributed under the terms of the GNU General Public License.

#This file is part of PyFRAP.

#PyFRAP is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

"""Z-stack module for PyFRAP toolbox. Used for reading z-stack datasets and converting them into meshes
that then can be used for a FRAP simulation.

.. warning:: This module is still experimental and unfinished and should not be used. 
   Will be added as stable in further versions.

"""

#===========================================================================================================================================================================
#Importing necessary modules
#===========================================================================================================================================================================

#Numpy/Scipy
import numpy as np

#Image processing
import skimage.measure as skimsr
import cv2

#matplotlib
import matplotlib.pyplot as plt

#Misc
import sys

#PyFRAP
from . import pyfrp_integration_module as pyfrp_integr

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Deletes small patches via opening, then closes holes via closing and then returns contours and areas of contours

def getContours(img,kernel=(10,10)):

	#Define kernel
	kernel = np.ones(kernel, np.uint8)
	
	#Open to erode small patches
	thresh = cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel)
	
	#Close little holes
	thresh = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE,kernel, iterations=4)
	
	#Find contours
	#contours=skimsr.find_contours(thresh,0)

	thresh=thresh.astype('uint8')
	contours, hierarchy = cv2.findContours(thresh, cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
	
	areas=[]
	for c in contours:
		areas.append(cv2.contourArea(c))
	
	return contours,thresh,areas

def fillEndStacks(n,debug=False):

	inpoints_all=[]

	if debug==1:
		fig=plt.figure()
		fig.show()
		j=1

	for i in range(len(files)):
		
		if i==0 or i==len(files)-1:
			inpoints=[]
			if fill_mode=="random":
			
				npts=fill_samples
				
				while len(inpoints)<npts:
				
					rx=random.random() 	
					ry=random.random() 	
					
					x=mpolys[i][:,1].min()+rx*(mpolys[i][:,1].max()-mpolys[i][:,1].min())
					y=mpolys[i][:,0].min()+ry*(mpolys[i][:,0].max()-mpolys[i][:,0].min())
					
					poly=list(zip(list(mpolys[i][:,1]),list(mpolys[i][:,0])))
					
					if point_inside_polygon(x,y,poly):
						
						inpoints.append([x,y])
			
			elif fill_mode=="regular":
				
				dx=abs(mpolys[i][:,1].min()-mpolys[i][:,1].max())
				dy=abs(mpolys[i][:,0].min()-mpolys[i][:,0].max())
				
				d=sqrt((dx*dy)/fill_samples)
					
				xvec=arange(mpolys[i][:,1].min(),mpolys[i][:,1].max(),d)
				yvec=arange(mpolys[i][:,0].min(),mpolys[i][:,0].max(),d)
				
				poly=list(zip(list(mpolys[i][:,1]),list(mpolys[i][:,0])))
					
				for x in xvec:
					for y in yvec:
						if point_inside_polygon(x,y,poly):
							
							inpoints.append([x,y])
				
			inpoints_all.append(inpoints)
			
			if debug==1:
			
				ax=fig.add_subplot(2,1,j)
				ax.plot(mpolys[i][:, 1], mpolys[i][:, 0],'r-', linewidth=2)
				
				for pt in inpoints:
					ax.plot(pt[0],pt[1],'g*')
				j=j+1
		
	if debug==1:
		plt.draw()
		input("Done filling first and last zstack, press ENTER to continue")


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Tries to interpolate gaps between contours (experimental/still not working)

def interpolateContourGaps(contours):

		
		#Take first contour to start with
		contour=contours[0]
		
		#Grab first and last point of contour
		pStart=contour[0,:]
		pEnd=contour[-1,:]
		
		#Build List of contours that are left
		contoursLeft=list(contours)
		contoursLeft.pop(0)
		
		#Remember contour
		contourLast=contour
		
		#Remember idx last
		idxLast=1
		
		#Create result array
		ptsFinal=contours[0].copy()
		
		#Create figure and plot contour and image
		fig=plt.figure()
		fig.show()
		ax=fig.add_subplot(111)
		ax.imshow(img,cmap='Greys')
		ax.plot(contour[:,0],contour[:,1],'r')
		plt.draw()
			
		while len(contoursLeft)>0:
			
			#Bookkeeping list
			distPEnd=[]
			
			#Now loop through contoursLeft and find the one that is closest
			for contour in contoursLeft:
				
				#Get Start/Endpoint
				pStartNew=contour[0,:]
				pEndNew=contour[-1,:]
				
				#Compute distance to both end and start point
				distPEnd.append(np.linalg.norm(pEnd-pEndNew))
				distPEnd.append(np.linalg.norm(pEnd-pStartNew))
				
				#Get the closest of the two
				idx=distPEnd.index(min(distPEnd))
				
				#Choose points for interpolation
				if idx==0:
					ptsNew=contour[-3:,:]
				elif idx==1:
					ptsNew=contour[:3,:]
					
				if idxLast==0:
					ptsOld=contourLast[:3,:]
				elif idxLast==1:
					ptsOld=contourLast[-3:,:]
				
				#Concatenate interpolation points
				ptsInterp=np.concatenate((ptsOld,ptsNew),axis=0)
				
				
				ax.plot(contour[:,0],contour[:,1],'r')
				ax.plot(ptsInterp[:,0],ptsInterp[:,1],'b*')
				plt.draw()
				
				
				
				#Generate Interpolation function
				tck,u=interp.splprep(ptsInterp.T,s=0.0)
				
				#Get interpolation range
				minx=min(ptsInterp[:,0])
				miny=min(ptsInterp[:,1])
				maxx=max(ptsInterp[:,0])
				maxy=max(ptsInterp[:,1])
				
				
				
				#Evaluate over ppoints
				x_i,y_i= interp.splev(np.linspace(minx,maxx,100),tck)
				
			
				
				ptsResult=np.zeros((len(x_i),2))
				ptsResult[:,0]=x_i
				ptsResult[:,1]=y_i
				
				
				
				#Append result
				ptsFinal=np.concatenate((ptsFinal,ptsResult),axis=0)
				ptsFinal=np.concatenate((ptsFinal,contour),axis=0)
				
				#Remeber idxLast
				idxLast=idx
			
				contourLast=contour
				
				
				#ax.plot(ptsResult[:,0],ptsResult[:,1],'g')
			
				plt.draw()
				
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Checks if a point is within a polygon via crossing number algorithm

def point_inside_polygon(x,y,poly):
	#Taken from http://www.ariel.com.au/a/python-point-int-poly.html
	n = len(poly)
	inside =False

	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y

	return inside