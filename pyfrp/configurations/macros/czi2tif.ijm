// This macro opens a list of lsm files and saves them into tif

//Read parsed arguments
arg = getArgument()

//Split and assign arguments
args=split(arg);
dir=args[0];

print("Going to convert images in folder" + dir);

//Retrive list of files
list2 = getFileList(dir);

// print file list
for (j=0;j<list2.length;j++) {
	print(list2[j]);
	
	path = dir+list2[j];
		
	if (File.isDirectory(path)==0){
		
		//Check to only open lsm files
		if (endsWith(path,".czi")){
			
			
			//Print Current file
			print("Converting file " +path);
				
			//Open
			run("Bio-Formats", "open="+path+" autoscale color_mode=Default view=Hyperstack stack_order=XYCZT");
			
			//Create path/filename for output
			fileName=replace(list2[j],".czi","");
			
			outpath = dir+fileName+"_c001_t001";
			
			//Save as image sequence
			run("Image Sequence... ", "format=TIFF name="+fileName+" start=0 digits=3 save="+outpath+".tif");
			
			close();
			}
		else {
			print("Skipping file ", path);
			}	
			
		}
	}

print ("Done.")
