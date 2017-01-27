# Extract by attribute
# find a crop type from the regions CDL layers with all crops for each year we want
# Description: Extracts the cells of a raster based on a logical query. 

# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *

# Set environment settings
env.workspace = r"Z:/Users/rilllydi/MidwestSALUS/Scratch"
env.scratchWorkspace = r"Z:/Users/rilllydi/MidwestSALUS/scratch.gdb"
env.overwriteOutput = True

# Check out the ArcGIS Spatial Analyst extension license
#arcpy.CheckOutExtension("Spatial")

# Define crop type for naming convention
Cr = "Corn"

##loop through each year (2010 to 2015)
for Year in range(2010, 2016):
    if Year == 2014 or Year == 2015:
        inRaster = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/AllCrops/cdl_" + str(Year) + ".tif"
    else:
        inRaster = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/AllCrops/cdl_" + str(Year)
    if Year == 2010 or Year == 2011:
        inSQLClause = "CLASS_NAMES = 'Corn' OR CLASS_NAMES = 'Dbl Crop Barley/Corn' OR CLASS_NAMES = 'Dbl Crop Corn/Soybeans' OR CLASS_NAMES = 'Dbl Crop Oats/Corn' OR CLASS_NAMES = 'Dbl Crop WinWht/Corn'"
    else:
        inSQLClause = "CLASS_NAME = 'Corn' OR CLASS_NAME = 'Dbl Crop Barley/Corn' OR CLASS_NAME = 'Dbl Crop Corn/Soybeans' OR CLASS_NAME = 'Dbl Crop Oats/Corn' OR CLASS_NAME = 'Dbl Crop WinWht/Corn'" 

    print "extracting"
    attExtract = ExtractByAttributes(inRaster, inSQLClause) 
    print "saving the output"
    attExtract.save("Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/Corn/cdl_" + str(Year) + "_" + Cr)