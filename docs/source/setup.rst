Installing PyFRAP via setup.py
==============================

PyFRAP can conveniently be installed via::

	python setup.py install
	
Since PyFRAP requires on I/O of data files that come with it we recommend using::
	
	python setup.py install --user
	
*setuptools* comes with multiple installation options, to check them out, type::

	python setup.py install --help
	
Note that PyFRAP will automatically install *Gmsh*. 	
	
PyFRAP has some additional installation options:

+------------------------+------------------------------------------------------+
| Option                 | Effect                                               |
+========================+======================================================+
| ``--fiji``             | Download and install Fiji during PyFRAP installation.| 
|                        | Link Fiji with PyFRAP                                | 
+------------------------+------------------------------------------------------+
| ``--noGmsh``           | Deactivate the download and installation of Gmsh     | 
|                        | during PyFRAP installation.                          | 
+------------------------+------------------------------------------------------+
| ``--silent``           | Print out less log messages.                         | 
|                        |                                                      | 
+------------------------+------------------------------------------------------+


For ``--fiji`` and the Gmsh installation to work, you need to install `wget <https://pypi.python.org/pypi/wget>`_. To install *wget*, type::

	pip install wget
	
or if you use the *Ananconda* distribution::

	conda install pywget

PyFRAP setup.py API
===================
	
.. py:function:: getOptions()	
	
	Checks options given to script:
	
		* If ``--fiji`` is in sys.argv, will set dFiji=1. \n
		* If ``--noGmsh`` is in sys.argv, will set dGmsh=0.
		* If ``--silent`` is in sys.argv, will set silent=1.
		
	
	.. note:: Makes dGmsh and dFiji global: Not nice but seems easiest way to get options into OverrideInstall.
	
.. py:function:: getOpt(optStr)	
	
	Checks if optStr is in ``sys.argv``. If this is the case,
	returns *1* and removes it form ``sys.argv`` so ``setup.py`` will not crash,
	otherwise returns *0*.
	
.. py:class:: OverrideInstall(install)
	
	Override class subclassing ``install`` class from ``setuptools``.
	
	The Main purpose of this class is to give more possibilities when installing PyFRAP, such as:
	
		* Download Gmsh and enter it automatically into path spec file
		* Download Fiji and enter it automatically into path spec file
		* Set ownership of data files so that even PyFRAP gets installed as superuser,
		  users will be able to use its full capacities.
	
	Idea taken from `here <http://stackoverflow.com/questions/5932804/set-file-permission-in-setup-py-file>`_ (thanks a bunch!)
	
	.. py:method:: initOptions()

		Parses options into override class.
		
	.. py:method::  run()
		
		Runs install. 
		
			
	.. py:method::  addData()
		
		Adds Datafiles to PyFRAP installation. 
		
		Makes sure that $USER has proper read/write/execute rights. Note that for Windows it will change rights,
		since it is not necesary. \n 
		Also makes sure that gmsh/Fiji bin ins properly linked.
		
		
	.. py:method:: cleanUpExe(fnDL,folderFn,filesBefore,exePath):	
		
		Moves it to executables directory and cleans up afterwards. 
		
		:param str fnDL: Filename of downloaded file.
		:param str folderFn: Filename of folder containing extracted files.
		:param list filesBefore: Snapshot of ``cwd``.
		:param str exePath: Path where executables go.
		
	.. py:method:: downloadGmsh()
		
		Downloads Gmsh, moves it to executables directory and cleans up afterwards. 
		
		.. note::  Only works if ``wget`` is installed. 
		
	.. py:method:: downloadGmshWin(arch,gmshVersion)
		
		Downloads Gmsh from Gmsh website for Windows
		
		:param str arch: System architecture, e.g. 64/32.
		:param str gmshVersion: gmshVersion String, e.g. 2.12.0 .
	
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)
					
	.. py:method:: downloadGmshOSX(arch,gmshVersion)
		
		Downloads Gmsh from Gmsh website for OSX.
		
		:param str arch: System architecture, e.g. 64/32.
		:param str gmshVersion: gmshVersion String, e.g. 2.12.0 .
	
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)

	.. py:method:: downloadGmshLinux(arch,gmshVersion)
		
		Downloads Gmsh from Gmsh website for Linux.
			
		:param str arch: System architecture, e.g. 64/32.
		:param str gmshVersion: gmshVersion String, e.g. 2.12.0 .
	
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)
		
		
	.. py:method:: makeExeFolder()
		
		Make executables folder if it doesn't exist yet
		
	
	.. py:method:: downloadFiji()
		
		Downloads Gmsh, moves it to executables directory and cleans up afterwards. 
		
		.. note::  Only works if ``wget`` is installed. 
		
	.. py:method:: downloadFijiLinux(arch)
		
		Downloads Fiji from Fiji website for Linux.
		
		:param str arch: System architecture, e.g. 64/32.
		
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)
			
	.. py:method:: downloadFijiWin(arch)
		
		Downloads Fiji from Fiji website for Windows.
		
		:param str arch: System architecture, e.g. 64/32.
		
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)
			

	.. py:method:: downloadFijiOSX()
		
		Downloads Fiji from Fiji website for OSX.
		
		:return: (Donwload filename, Filename of extracted download files)
		:rtype: (str ,str)
			
	
	.. py:method:: setExePath(fn,identifier,exePath)
		
		Enters executable path into path spec file.
		
		:param str fn: Path to gmsh executable.
		:param str identifier: Identifier in spec file.
		:param str exePath: Path to exe file
			
		
	.. py:method:: setGmshPath(fn)
		
		Enters gmsh executable path into path spec file.
		
		:param str fn: Path to gmsh executable.
			
		
	.. py:method:: setFijiPath(fn)
		
		Enters fiji executable path into path spec file.
		
		:param str fn: Path to fiji executable.
			
	.. py:method:: changePermissions(filepath,uid,gid,mode)
		
		Sets File Permissions.
		
		:param str filepath: Path to file.
		:param int uid: user ID.
		:param int gid: group ID.
		:param int mode: Permission mode.
		
		:return: True if success
		:rtype: bool
		
		
		
	.. py:method:: makeAdditionalDataFolders(folder,fn,uid,gid,mode)
		
		Tries to generate additional data folders.
		
		:param str folder: Path to containing folder.
		:param str fn: New folder name
		:param int uid: user ID.
		:param int gid: group ID.
		:param int mode: Permission mode.
		
		:return: True if success
		:rtype: bool
		
		
		
		
	
	