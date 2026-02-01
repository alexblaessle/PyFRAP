import sys
import argparse


# Create parser
parser = argparse.ArgumentParser(description='PyFRAP command line.')

parser.add_argument('-f','--fIn', dest='fIn',action='append',help='Input file names. Can be either embryo files or paths to image data. Note that all inputs should be of the same type.')

# Options for project creation

# Options for image analysis
parser.add_argument('-G','--gaussian', dest='gaussian', action='const_store',help='Perform gaussian blur. Specified as tuple flag,sigma',default="0,5")
parser.add_argument('-M','--median', dest='median', action='const_store',help='Perform median filter. Specified as tuple flag,radius',default="0,2")
parser.add_argument('-F','--flatten', dest='flatten', action='const_store',help='Perform flattening. Specified as tuple flag,pathToFlattenData',default="0,")
parser.add_argument('-B','--bkgd', dest='bkgd', action='const_store',help='Perform background subtraction. Specified as tuple flag,pathToBkgdData',default="0,")
parser.add_argument('-N','--norm', dest='norm', action='const_store',help='Perform normalization. Specified as tuple flag,pathToPreData',default="0,")

# Options for meshing



# Options for simulations
parser.add_argument('--IC', dest='IC', action='const_store',help='Choose initial conditions',default=4)
parser.add_argument('-t','--tsteps', dest='tsteps', action='const_store',help='Number of simulations timesteps',default=3000)
parser.add_argument('-s','--stepping', dest='stepping', action='const_store',help='Time stepping (0=linear,1=log)',default=1)
parser.add_argument('--IC', dest='IC', action='const_store',help='Choose initial conditions',default=4)

# Options for fitting








args = parser.parse_args()


print(args)

# Choose what kind of data is given and what to do with it.
#inputTypes=[]
#for inp in args.fIn:
	#if os.path.isdir(inp):
		#fns=os.listdir(inp)
		
		
		
		#for f in fns:
			#if os.path.isdir(os.path.join(inp,f)):
				#pass
			#else:
				#path,ext=os.path.splitext(inp)
				
				#if ext in ['.czi','.lsm','.tif','']
				
				
	#else:
		#path,ext=os.path.splitext(inp)
		
		#if ext=='.emb':
			#inputTypes.append('embryo')
		#elif ext=='.mol':
			#inputTypes.append('molecule')
		#else:
			#printWarning("Unknown input type " + inp)
			
			
		
		


