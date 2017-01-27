# Mosiac multiple rasters together (run within arcgis)

import arcpy
from arcpy import env
from arcpy.sa import *

# What is the crop type? (should be the name of the folder and the crop in the filenames also)
Cr = "corn"
Cr2 = "CornOnly"
Cr3 = "Corn"

env.workspace = r"Z:/Users/rilllydi/MidwestSALUS/Scratch"
#env.scratchWorkspace = r"Z:/Users/rilllydi/MidwestSALUS/scratch.gdb" # DO not use otherwise error
env.overwriteOutput = True

f0 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2010_" + Cr
f1 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2011_" + Cr
f2 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2012_" + Cr
f3 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2013_" + Cr
f4 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2014_" + Cr
f5 = "Z:/Users/rilllydi/MidwestSALUS/CDL/IndivYears/Region/" + Cr2 + "/cdl_2015_" + Cr

target = f0.replace("IndivYears","Years2010-15")
target = target.replace(Cr2 + "/cdl_2010_" + Cr,Cr3 + "/cdl_" + Cr)

# copy an existing raster (f0) to the output directory
arcpy.CopyRaster_management(f0,target)
target2 = target.replace("cdl_" + Cr,"cdl_" + Cr + "1")

##Background value: 0
##Nodata value: 0
##Merge to the exisiting raster (f0) which has been renamed
rasters = "\""+f1+";"+f2+";"+f3+";"+f4+";"+f5+"\" "
arcpy.Mosaic_management(rasters,target,"LAST","FIRST","0", "0", "", "", "")

###################
# Reclassify the new raster so that all the crops have a value of 1
outCon = Con(target, 1)
outCon.save(target2)