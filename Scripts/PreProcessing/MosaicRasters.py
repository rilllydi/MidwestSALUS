# Mosiac multiple rasters together 
import arcpy
from arcpy import env
from arcpy.sa import *
import glob

# Set environment settings
env.workspace = r"C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles"
env.overwriteOutput = True


#states = ['MI','OH','IN','IL','IA','SD','MN','MO']
states = ['WI']
for state in states:

	# get all the shapefiles for the state
	pth = "C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles/soils_SSURGSDM_" + state + "/"
	files = glob.glob(pth + "*NAD.tif")
	num = len(files)
	print num
	chunk = round(int(num)/4)
	print chunk

	for i in range(0,len(files),chunk):
			j = i + chunk
			if j > len(files):
				j = len(files)
				
			print i, j
			target = "C:/Users/Lydia Rill/Documents/SSURGO_Soils/Raster_pieces/" + state + "_" + str(i)
			print files[i]
			print "FIRST"
			print files[i:j]
			
			# copy an existing raster (f0) to the output directory
			#arcpy.CopyRaster_management(files[i],target)

			##Background value: 0
			##Nodata value: 0
			##Merge to the exisiting raster (f0) which has been renamed
			#arcpy.Mosaic_management(files[i:j],target,"LAST","FIRST","0", "0", "", "", "")




