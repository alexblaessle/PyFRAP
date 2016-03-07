#===========================================================================================================================================================================
#Module Description
#===========================================================================================================================================================================

#Integration module for PyFRAP toolbox, including following functions:

#(1)  

#===========================================================================================================================================================================
#Improting necessary modules
#===========================================================================================================================================================================


from numpy import linalg as LA

#===========================================================================================================================================================================
#Module Functions
#===========================================================================================================================================================================

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Returns 

def getAvgConc(val,cvs,ind):
	if len(ind)>0:
		return sum(val.value[ind]*cvs[ind])/sum(cvs[ind])
	else:
		return 0.



#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Calculates sidelengths of tetrahedron given by 4 points

def calcTetSidelengths(point0,point1,point2,point3):

	#Taking point0 as base point, calculating vectors
	vec1=point1-point0
	vec2=point2-point0
	vec3=point3-point0
	
	norms=[LA.norm(vec1),LA.norm(vec2), LA.norm(vec3)]
	
	return norms
