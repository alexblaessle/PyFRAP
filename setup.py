#Import setuptools
from setuptools import setup

#We need those two for doing some file permission magic later
import os
import platform
import shutil

#Overwrite setuptools install to be able to set file permissions
from setuptools.command.install import install
from distutils import log 
from setuptools.command.install_scripts import install_scripts

#Import option parser so we can parse in options
import sys

# Define version
version='1.1.13'

def getOptions():
	
	"""Checks options given to script.
	
	If --fiji is in sys.argv, will set dFiji=1. \n
	If --gmsh is in sys.argv, will set dGmsh=1.
	If --silent is in sys.argv, will set silent=1.
	
	
	
	Note: Makes dGmsh and dFiji global: Not nice but seems easiest 
	way to get options into OverrideInstall.
	"""
	
	global dGmsh
	global dFiji
	global silent
	
	dFiji=getOpt("--fiji")
	dGmsh=getOpt("--noGmsh",default=0)
	silent=getOpt("--silent")
	
def getOpt(optStr,default=0):	
	
	"""Checks if optStr is in sys.argv. If this is the case,
	returns 1 and removes it form sys.argv so setup.py will not crash,
	otherwise returns default.
	"""
	
	if optStr in sys.argv:
		opt=1
		sys.argv.remove(optStr)
	else:
		opt=default
	return opt	

#Get Options
if __name__ == '__main__':
	getOptions()

class OverrideInstall(install):
	
	"""Override class subclassing install class from setuptools.
	
	The Main purpose of this class is to give more possibilities when installing PyFRAP, such as:
	
		* Download Gmsh and enter it automatically into path spec file
		* Download Fiji and enter it automatically into path spec file
		* Set ownership of data files so that even PyFRAP gets installed as superuser,
		users will be able to use its full capacities.
	
	Idea taken from http://stackoverflow.com/questions/5932804/set-file-permission-in-setup-py-file (thanks a bunch!)
	
	"""
	
	def initOptions(self):
		"""Parses options into override class.
		"""
		self.dFiji=bool(dFiji)
		self.dGmsh=not bool(dGmsh)
		self.silent=bool(silent)
		
		log.info("Install options are:")
		log.info("Install Fiji: " + str(self.dFiji))
		log.info("Install gmsh: " + str(self.dGmsh))
		log.info("Silent Mode: " + str(self.silent))
		
		#Define pathFile
		self.pathFile='paths'
			
	def run(self):
		
		"""Runs install. 
		
		"""
		
		self.initOptions()
		
		#Try to download gmsh
		if self.dGmsh:
			self.downloadGmsh()
		else:
			self.gmshDownloaded=False
		
		#Try to download fiji
		if self.dFiji:
			self.downloadFiji()
		else:
			self.fijiDownloaded=False
		
		#Run setuptools install
		install.run(self) 
		
		#Print log info
		if not self.silent:
			log.info("Overriding setuptools mode of scripts ...")
		
		#Add Data and edit file permissions
		self.addData()
			
	def addData(self):
		
		"""Adds Datafiles to PyFRAP installation. 
		
		Makes sure that $USER has proper read/write/execute rights. Note that for Windows it will change rights,
		since it is not necesary. \n 
		Also makes sure that gmsh/Fiji bin ins properly linked.
		"""
		
		if not self.silent:
			log.info("in add data")
		
		uid,gid,mode=self.getPermDetails()
		
		#Overwrite file permissions
		for filepath in self.get_outputs():
			
			#log.info("Copying files.")
			
			if platform.system() not in ["Windows"]:
			
				if "meshfiles" in filepath or "configurations" in filepath or "executables" in filepath:
					
					#Change permissions for folder containing files
					folderpath=os.path.dirname(os.path.realpath(filepath))
					self.changePermissions(folderpath,uid,gid,mode)
					
					#Change permissions of file
					self.changePermissions(filepath,uid,gid,mode)
					
					#Make some more necessary data folders
					if folderpath.endswith("meshfiles"):
						self.makeAdditionalDataFolders(folderpath,"field",uid,gid,mode)
						self.makeAdditionalDataFolders(folderpath,"field/custom",uid,gid,mode)
						
					if folderpath.endswith("configurations"):
						self.makeAdditionalDataFolders(folderpath,"macros",uid,gid,mode)
			
			#log.info("Adding executables to path file")
			
			#Add gmsh into paths.default if download was successful
			if self.pathFile == os.path.basename(filepath):
				
				if self.gmshDownloaded:
					self.setGmshPath(filepath)
				if self.fijiDownloaded:
					self.setFijiPath(filepath)
				if platform.system() not in ["Windows"]:
					folderpath=os.path.dirname(os.path.realpath(filepath))
					self.changePermissions(folderpath,uid,gid,mode)
					self.changePermissions(filepath,uid,gid,mode)

	
	def getPermDetails(self):
		
		"""Returns the permission details used to change permissions.
		"""
		
		if platform.system() not in ["Windows"]:
			import pwd
			
			#Grab user ID and group ID of actual use
			try:
				uid=pwd.getpwnam(os.getlogin())[2]
				gid=pwd.getpwnam(os.getlogin())[3]
			
			except: 
				if not self.silent:
					log.info("Was not able to retrieve UID via os.getlogin, using os.getuid instead.")
				uid=os.getuid()
				gid=os.getgid()
				
			#Mode for files (everyone can read/write/execute. This is somewhat an overkill, but 0666 seems somehow not to work.)
			#mode="0777"
			mode="0o777"
			
			
			return uid,gid,mode
			
		return	0,0,0
	
	def cleanUpExe(self,fnDL,folderFn,filesBefore,exePath):		
		
		"""Moves it to executables directory and cleans up afterwards. 
		
		"""
		
		#Copy file to pyfrp/executables/
		try:
			shutil.rmtree(exePath)
		except:
			pass
		
		
		
		if os.path.isdir(folderFn):
			shutil.copytree(folderFn+"/",exePath)
		else:
			try:
				os.mkdir(os.path.dirname(exePath))
			except:
				log.info("Was not able to create folder " + os.path.dirname(exePath))
			shutil.copy(folderFn,exePath)
				
		#Remove downloaded files
		os.remove(fnDL)
		
		#Get fileList before
		filesAfter=os.listdir('.')
		
		self.cleanDiff(filesBefore,filesAfter)
		
	def cleanDiff(self,filesBefore,filesAfter):
		
		#Difference between files
		filesDiff=list(set(filesAfter)-set(filesBefore))
		
		for fn in filesDiff:
				
			try:
				if os.path.isdir(fn):
					shutil.rmtree(fn)
				else:
					os.remove(fn)
			except:
				log.info("cleanDiff report: Was not able to delete file:" + fn)
				
	def downloadGmsh(self):
		
		"""Downloads Gmsh, moves it to executables directory and cleans up afterwards. 
		
		Note that this will only work if *wget* is installed. 
		"""
		
		#Define gmshVersion (might need to update this line once in a while)
		gmshVersion='2.14.0'
		
		#Flag to see if gmsh DL went through
		self.gmshDownloaded=False
		
		self.makeExeFolder()
		
		#Get fileList before
		filesBefore=os.listdir('.')
		
		#Try to import wget
		try: 
			import wget
			
			#Get Architecture
			arch=platform.architecture()[0].replace('bit','')
			
			if platform.system() in ["Windows"]:
				fnDL,folderFn=self.downloadGmshWin(arch,gmshVersion)
				
			elif platform.system() in ["Linux"]:
				fnDL,folderFn=self.downloadGmshLinux(arch,gmshVersion)
				
			elif platform.system() in ["Darwin"]:
				fnDL,folderFn=self.downloadGmshOSX(arch,gmshVersion)
			
			#Remove files
			self.cleanUpExe(fnDL,folderFn,filesBefore,'pyfrp/executables/gmsh/')
			
			uid,gid,mode=self.getPermDetails()
			if platform.system() not in ["Windows"]:
				self.changePermissions(self.gmshPath,uid,gid,mode)
			
			#self.addPathToWinPATHs(self.gmshPath)
			
			log.info("Installed gmsh to "+ self.gmshPath)
			
			#Set Flag=True
			self.gmshDownloaded=True
			
		except ImportError:
			log.info("Cannot find wget, will not be downloading gmsh. You will need to install it later manually")	
	
	def downloadFileIfNotExist(self,url):
		
		"""Downloads URL if file does not already exist.
		
		Args:
			url (str): URL to download.
			
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		import wget
		
		cwd=os.getcwd()+"/"
		
		if not os.path.exists(cwd+os.path.basename(url)):
			log.info(cwd+os.path.basename(url) +" does not exist, will download it.")
			folderFn=wget.download(url)
		else:
			log.info(cwd+os.path.basename(url) +" alreay exists, will not download.")
			folderFn=os.path.basename(url)
		fnDL=str(folderFn)
		print
		
		return folderFn, fnDL
			
	def downloadGmshWin(self,arch,gmshVersion):
		
		"""Downloads Gmsh from Gmsh website for Windows
		
		Args:
			arch (str): System architecture, e.g. 64/32.
			gmshVersion (str): gmshVersion String, e.g. 2.12.0 .
			
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		#Download Gmsh
		url='http://gmsh.info/bin/Windows/gmsh-'+gmshVersion+'-Windows'+arch+'.zip'
		folderFn, fnDL=self.downloadFileIfNotExist(url)
		
		#Decompress
		import zipfile 
		with zipfile.ZipFile(folderFn) as zf:
			zf.extractall()
			
		folderFn='gmsh-'+gmshVersion+'-Windows'	
		
		self.gmshPath='executables/gmsh/gmsh.exe'
		
		return fnDL,folderFn
		
	def downloadGmshOSX(self,arch,gmshVersion):
		
		"""Downloads Gmsh from Gmsh website for OSX.
		
		Args:
			arch (str): System architecture, e.g. 64/32.
			gmshVersion (str): gmshVersion String, e.g. 2.12.0 .
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		
		
		#Download Gmsh (if file isn't there yet)
		url='http://gmsh.info/bin/MacOSX/gmsh-'+gmshVersion+'-MacOSX'+'.dmg'
		folderFn, fnDL=self.downloadFileIfNotExist(url)
		
		#Mount dmg file (Here the user need to read through LICENSE, don't know how to fix this)
		log.info("executing: "+ 'hdiutil attach '+folderFn)
		os.system('hdiutil attach '+folderFn)
		folderFn=folderFn.replace('.dmg','')
		
		#try:
			#os.mkdir(folderFn)
		#except OSError:
			#pass
		
		cwd=os.getcwd()
		
		#Copy gmsh executable to cwd 
		#Note: It seems to vary where gmsh executable is in mounted dmg file, hence we 
		#just have to try out, take the one that actually worked and remember it
		rets=[]
		possFiles=["bin/","share/","gmsh"]
		
		rets.append(os.system('cp -rv /Volumes/'+folderFn+'/Gmsh.app/Contents/MacOS/bin/ '+ cwd))
		rets.append(os.system('cp -rv /Volumes/'+folderFn+'/Gmsh.app/Contents/MacOS/share/ '+ cwd))
		rets.append(os.system('cp -rv /Volumes/'+folderFn+'/Gmsh.app/Contents/MacOS/gmsh '+ cwd))
		
		fnWorked=possFiles[rets.index(0)]
		
		#Unmount gmsh
		os.system('hdiutil detach /Volumes/'+folderFn+'/')
		
		#Build filename of acutally copied file
		folderFn=cwd+"/"+fnWorked
		
		self.gmshPath='executables/gmsh/./gmsh'
		
		return fnDL,folderFn
		
	def downloadGmshLinux(self,arch,gmshVersion):
		
		"""Downloads Gmsh from Gmsh website for Linux.
		
		Args:
			arch (str): System architecture, e.g. 64/32.
			gmshVersion (str): gmshVersion String, e.g. 2.12.0 .
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
	
		#Download Gmsh
		url='http://gmsh.info/bin/Linux/gmsh-'+gmshVersion+'-Linux'+arch+'.tgz'
		folderFn, fnDL=self.downloadFileIfNotExist(url)
		
		#Decompress
		import tarfile
		with tarfile.open(folderFn,mode='r:gz') as zf:
			zf.extractall()
		
		folderFn='gmsh-'+gmshVersion+'-Linux'
		
		self.gmshPath='executables/gmsh/bin/./gmsh'
		
		return fnDL,folderFn	
	
	def makeExeFolder(self):
		
		#Make executables folder if it doesn't exist yet
		try:
			os.mkdir('pyfrp/executables')	
		except OSError:
			log.info('Was not able to create directory pyfrp/executables')
	
	def addPathToWinPATHs(self,path):
	
		"""Adds a path to Windows' PATH list. 
		
		.. note:: Only adds path if file exits.
		
		.. note:: You will need to restart the terminal to 
		be sure that the change has any effect.
		
		Args:
			path (str): Path to be added.
			
		Returns:
			bool: True if successful.
		"""
		
		if platform.system() not in ["Windows"]:
			log.info("OS is not Windows, won't set path")
			return False
		
		if path in os.environ['PATH']:
			log.info("Path is already in PATH, won't set path")
			return False
		
		if os.path.exists(path):
			os.system("set PATH=%PATH%;"+path)
			return True
		else:
			log.info(path + " does not exist, won't set path")
			return False
	
	def downloadFiji(self):
		
		"""Downloads Gmsh, moves it to executables directory and cleans up afterwards. 
		
		Note that this will only work if *wget* is installed. 
		"""
		
		#Flag to see if gmsh DL went through
		self.fijiDownloaded=False
		
		self.makeExeFolder()
		
		#Get fileList before
		filesBefore=os.listdir('.')
		
		#Try to import wget
		try: 
			import wget
			
			#Get Architecture
			arch=platform.architecture()[0].replace('bit','')
			
			if platform.system() in ["Windows"]:
				fnDL,folderFn=self.downloadFijiWin(arch)
				
			elif platform.system() in ["Linux"]:
				fnDL,folderFn=self.downloadFijiLinux(arch)
				
			elif platform.system() in ["Darwin"]:
				fnDL,folderFn=self.downloadFijiOSX(arch)
			
			#Remove files
			self.cleanUpExe(fnDL,folderFn,filesBefore,'pyfrp/executables/Fiji.app/')
			
			uid,gid,mode=self.getPermDetails()
			if platform.system() not in ["Windows"]:
				self.changePermissions(self.fijiPath,uid,gid,mode)
			
			#self.addPathToWinPATHs(self.fijiPath)
			
			log.info("Installed Fiji to "+ self.fijiPath)
			
			#Set Flag=True
			self.fijiDownloaded=True
			
		except ImportError:
			log.info("Cannot find wget, will not be downloading fiji. You will need to install it later manually")	
	
	def downloadFijiLinux(self,arch):
		
		"""Downloads Fiji from Fiji website for Linux.
		
		Args:
			arch (str): System architecture, e.g. 64/32.
		
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		import wget
		
		#Download Fiji
		url='http://downloads.imagej.net/fiji/latest/fiji-linux'+arch+'.zip'
		
		folderFn=wget.download(url)
		fnDL=str(folderFn)
		print
		
		#Decompress
		import zipfile 
		with zipfile.ZipFile(folderFn) as zf:
			zf.extractall()
		
		
		folderFn='Fiji.app'
		
		self.fijiPath='executables/Fiji.app/./ImageJ-linux64'
		
		return fnDL,folderFn	

	def downloadFijiWin(self,arch):
		
		"""Downloads Fiji from Fiji website for Windows.
		
		Args:
			arch (str): System architecture, e.g. 64/32.
		
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		import wget
		
		#Download fiji
		url='http://downloads.imagej.net/fiji/latest/fiji-win'+arch+'.zip'
		
		folderFn=wget.download(url)
		fnDL=str(folderFn)
		print
		
		#Decompress
		import zipfile 
		with zipfile.ZipFile(folderFn) as zf:
			zf.extractall()
		
		folderFn='Fiji.app'
		
		self.fijiPath='executables/Fiji.app/ImageJ-linux64.exe'
		
		return fnDL,folderFn	
	
	def downloadFijiOSX(self,arch):
		
		"""Downloads Fiji from Fiji website for OSX.
		
		Returns:
			tuple: Tuple containing:
			
				* fnDL (str): Donwload filename
				* folderFn (str): Filename of extracted download files
			
		"""
		
		import wget
		
		#Download fiji
		url='http://downloads.imagej.net/fiji/latest/fiji-macosx.dmg'
		
		folderFn=wget.download(url)
		fnDL=str(folderFn)
		print
	
		
		#Mount dmg file 
		os.system('hdiutil attach '+folderFn)
		
		cwd=os.getcwd()
		
		#Copy fiji executable to cwd
		os.system('cp -rv /Volumes/Fiji/Fiji.app '+ cwd)
		
		#Unmount gmsh
		os.system('hdiutil detach /Volumes/Fiji')
	
		folderFn='Fiji.app'
		
		self.fijiPath='executables/Fiji.app/Contents/MacOS/./ImageJ-macosx'
		
		return fnDL,folderFn	
	
	def winToLinPath(self,path):
		
		return path.replace('\\','/')
		
	def setExePath(self,fn,identifier,exePath):
		
		"""Enters executable path into path spec file.
		
		Args:
			fn (str): Path to gmsh executable.
			identifier (str): Identifier in spec file.
			exePath (str): Path to exe file
			
		"""
		
		#Make backup of default path file
		shutil.copy(fn,fn+'_backup')
		
		#Get filepath to PyFRAP
		fnPyfrp=fn.split('configurations')[0]
		
		# Remove possible backslashes from path
		fn=self.winToLinPath(fn)
		
		#Open file and enter new gmsh bin
		with open(fn,'rb') as fPath:
			with open(fn+"_new",'wb') as fPathNew:
				for line in fPath:
					if line.strip().startswith(identifier):
						ident,path=line.split('=')
						path=path.strip()
						lineNew=ident+"="+fnPyfrp+exePath
						fPathNew.write(lineNew+'\n')
					else:
						fPathNew.write(line)
			
		#Rename file
		shutil.move(fn+'_new',fn)
		
	def setGmshPath(self,fn):
		
		"""Enters gmsh executable path into path spec file.
		
		Args:
			fn (str): Path to gmsh executable.
			
		"""
		
		self.setExePath(fn,'gmshBin',self.gmshPath)
	
	def setFijiPath(self,fn):
		
		"""Enters fiji executable path into path spec file.
		
		Args:
			fn (str): Path to fiji executable.
			
		"""
		
		self.setExePath(fn,'fijiBin',self.fijiPath)
	
	def changePermissions(self,filepath,uid,gid,mode):
		
		"""Sets File Permissions.
		
		Args:
			filepath (str): Path to file.
			uid (int): user ID.
			gid (int): group ID.
			mode (int): Permission mode.
		
		Returns:
			bool: True if success
		
		"""
		
		ret=True
		try:
			os.chown(filepath, uid, gid)
			if not self.silent:
				log.info("Changing ownership of %s to uid:%s gid %s" %(filepath, uid, gid))
		except:
			if not self.silent:
				log.info("Was not able to change ownership of file %s" %(filepath))
			ret=False

		try:
			if not self.silent:
				log.info("Changing permissions of %s to %s" %(filepath, mode))
			os.chmod(filepath, mode)
		except:	
			if not self.silent:
				log.info("Was not able to change file permissions of file %s" %(filepath))
			ret=False
		return ret
		
	def makeAdditionalDataFolders(self,folder,fn,uid,gid,mode):
		
		"""Tries to generate additional data folders.
		
		Args:
			folder (str): Path to containing folder.
			fn (str): New folder name
			uid (int): user ID.
			gid (int): group ID.
			mode (int): Permission mode.
		
		Returns:
			bool: True if success
		
		"""
		
		if not folder.endswith("/"):
			folder=folder+"/"
			
		if os.path.isdir(folder+fn):
			return False
		else:
			try:
				os.mkdir(folder+fn)
				self.changePermissions(folder+fn,uid,gid,mode)
				return True
			except:
				log.info("Unable to create folder %s" %(folder+fn))
				return False
			
#Define setup

#Check if setup.py is used to build RTD, then don't overwrite install command
if os.environ.get('READTHEDOCS', None) == 'True':
	
	print("Installing on RTD, will not overwrite install command.")
	
	setup(name='pyfrp',
		version=version,
		description='PyFRAP: A Python based FRAP analysis tool box',
		url='https://github.com/alexblaessle/PyFRAP',
		author='Alexander Blaessle',
		author_email='alexander.blaessle@tuebingen.mpg.de',
		license='GNU GPL Version 3',
		packages=['pyfrp','pyfrp.modules','pyfrp.subclasses','pyfrp.gui'],
		package_dir={'pyfrp': 'pyfrp',
				'pyfrp.modules': 'pyfrp/modules',
				'pyfrp.gui' : 'pyfrp/gui'
				},
		#package_data = {'pyfrp':['meshfiles','configurations']},
		include_package_data=True,
		classifiers= [
			'Development Status :: 4 - Beta',
			'Operating System :: OS Independent',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',\
			'Programming Language :: Python :: 2.7',
			],
		install_requires=['pyopenssl','ndg-httpsclient','pyasn1','ez_setup','numpy','scipy','matplotlib','scikit-image','FiPy','colorama','numpy-stl','solidpython','wget','python-bioformats'],	
		platforms=['ALL'],
		keywords=["FRAP", "fluorescence",'recovery','after','photobleaching','reaction-diffusion','fitting'
			],
		zip_safe=False
		)


else:

	setup(name='pyfrp',
		version=version,
		description='PyFRAP: A Python based FRAP analysis tool box',
		url='https://github.com/alexblaessle/PyFRAP',
		author='Alexander Blaessle',
		author_email='alexander.blaessle@tuebingen.mpg.de',
		license='GNU GPL Version 3',
		packages=['pyfrp','pyfrp.modules','pyfrp.subclasses','pyfrp.gui'],
		package_dir={'pyfrp': 'pyfrp',
				'pyfrp.modules': 'pyfrp/modules',
				'pyfrp.gui' : 'pyfrp/gui'
				},
		#package_data = {'pyfrp':['meshfiles','configurations']},
		include_package_data=True,
		classifiers= [
			'Development Status :: 4 - Beta',
			'Operating System :: OS Independent',
			'Programming Language :: Python :: 2.7',
			'Topic :: Scientific/Engineering :: Bio-Informatics',
			'Intended Audience :: Science/Research',
			'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',\
			'Programming Language :: Python :: 2.7',
			],
		install_requires=['pyopenssl','ndg-httpsclient','pyasn1','ez_setup','numpy','scipy','matplotlib','scikit-image','FiPy','colorama','numpy-stl','solidpython','wget','python-bioformats'],	
		platforms=['ALL'],
		keywords=["FRAP", "fluorescence",'recovery','after','photobleaching','reaction-diffusion','fitting'
			],
		zip_safe=False,
		cmdclass={'install': OverrideInstall} #Need this here to overwrite our install
		)
