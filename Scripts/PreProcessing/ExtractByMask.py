# Name: ExtractByMask_Ex_02.py
# Description: Extracts the cells of a raster that correspond with the areas
#    defined by a mask.
# Requirements: Spatial Analyst Extension

# in this case the raster is the soils and the mask is the crop (corn) raster

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

Cr = "CornOnly"
sr = arcpy.SpatialReference(4326)

# loop through states
states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
for st in states:

    # Set local variables
    inRaster = "Z:/Users/rilllydi/MidwestSALUS/Soils/States/" + st + "_soil_A"
    inMaskData = "Z:/Users/rilllydi/MidwestSALUS/CDL/Years2010-15/Region/" + Cr + "/cdl_corn1"

    # Execute ExtractByMask
    outExtractByMask = ExtractByMask(inRaster, inMaskData)
    print "now saving the extraction"
    # Save the output 
    output = "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/" + Cr + "/" + st + "_soil_" + Cr + ".tif"
    outExtractByMask.save(output)
    
    #################################################################################

    # Convert to WGS84
    wgs84out = output.replace(".tif","_WGS84.tif")
    arcpy.ProjectRaster_management(output,wgs84out,sr)