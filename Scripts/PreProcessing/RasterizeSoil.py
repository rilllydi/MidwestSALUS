import arcpy
from arcpy import env
import glob

# Set environment settings
env.workspace = "C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles"
env.overwriteOutput = True

# Set local variables
valField = "MUKEY"
outRaster = "c:/output/ca_counties"
assignmentType = "MAXIMUM_AREA"
#priorityField = ""
cellSize = 30 #meters

states = ['MI','OH','IN','IL','IA','SD','MN','MO']
#states = ['WI']
for state in states:

	# get all the shapefiles for the state
	pth = "C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles/soils_SSURGSDM_" + state + "/"
	files = glob.glob(pth + "*NAD.shp")
	
	#files = ["C:/Users/Lydia Rill/Documents/SSURGO_Soils/Original_County_Shapefiles/soils_SSURGSDM_WI\soilmu_a_wi139_NAD.shp"]
	for file in files:
			print file
			output = file.replace(".shp",".tif")

			# Execute PolygonToRaster
			arcpy.PolygonToRaster_conversion(file, valField, output, assignmentType, "NONE", cellSize)