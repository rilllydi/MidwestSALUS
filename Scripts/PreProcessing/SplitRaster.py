# split raster

import arcpy
from arcpy import env
#from arcpy.sa import *
import os
import csv

Cr = "CornOnly"
basedir = "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/" + Cr + "/"

# Set environment settings
#env.workspace = basedir
#env.scratchWorkspace = r"Z:/Users/rilllydi/MidwestSALUS/scratch.gdb"
#env.overwriteOutput = True

states = ['mi','oh','in','il','ia','sd','mn','mo','wi']
#states = ['wi']

for st in states:
    rasterfile = basedir + st + "_soil_" + Cr + "_WGS84.tif"
    splitfolder = basedir + st + "_split_55/"
    name = st + "_soil_" + Cr + "_WGS84_"
    
    print rasterfile
    print splitfolder
    print name
    
    if not os.path.exists(splitfolder):
        os.makedirs(splitfolder)
    
    print "created folder"
    
    nodataval = '-32768'
    ### Split the raster (must use "NEAREST" option so the pixel values will not change!)
    #arcpy.SplitRaster_management(rasterfile, splitfolder, name, "SIZE_OF_TILE",\
    #                             "TIFF", "NEAREST",nodata_value=nodataval) # default tile size is 128 X 128
    arcpy.SplitRaster_management(rasterfile, splitfolder, name, "NUMBER_OF_TILES",\
                                 "TIFF", "NEAREST",num_rasters="5 5",nodata_value=nodataval) # 9 rasters (3 in the x direction, 3 in the y direction)