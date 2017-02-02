#PyQGIS script to merge all the individual county soil shapefiles into one large state shapefile
#This script is NOT standalone. It needs to be run through the python application in QGIS.
#Note the state parameter may need to be changed by the user
import glob
import processing
#from processing.core.Processing import Processing
#Processing.initialize()
#from processing.tools import *
#general.runalg('saga:mergelayers',files,True,True,out)
#Processing.updateAlgsList()
#states = ["WI"]
states = ['MI','OH','IN','IL','IA','SD','MN','MO']
for state in states:

	# get all the shapefiles for the state
	pth = "C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles/soils_SSURGSDM_" + state + "/"
	files = glob.glob(pth + "*NAD.shp")

	# save the shapefiles as NAD 83 Conus Albers
	crs = 'EPSG:5070'  
	#files = ["C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles/soils_SSURGSDM_WI\soilmu_a_wi139_NAD.shp"]
	for file in files:
        
		output = file.replace(".shp","_NAD.shp")
		#output = file
		print output
		processing.runalg("qgis:reprojectlayer", file, crs, output)