#!/usr/bin/Rscript
args = commandArgs(TRUE)

options(warn=1)

# Arguments: state, scenario, results folder, variable (summary or yearly variable), outdirectory for raster
# test if there is are 4 arguments: if not, return an error
if (length(args) != 1) {
  stop("1 argument must be supplied (state)", call.=FALSE)
} else {
  st <- args[1]
}

# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The SALUS results (per MUKEY and weather ID) are matched to each pixel grid point
# A new raster is created for each state and this raster should contain multiple attributes (all variables within the input file)
# Have the HPC do runs for each state (9), for each scenario (7), and for each yearly variable!

# Note using a statewide raster may crash R
# Note you need the associated .vat.dbf file along with the original raster for this to work

# Input Data:
#   - raster with MUKEY values
#   - polygon shapefile of the NLDAS weather grid (with the attribute table including the centroid coordinates)
#   - csv file of the SALUS results along with the associated MUKEY and Weather ID

# Process:
# 1. Read in the raster along with the associated dbf file (without the dbf file you get incorrect MUKEY values!)
# 2. Create a new raster of just the MUKEY values
# 3. Create a data table of the pixel coordinates and associated MUKEY value
# 4. Overlay the data table the the weather grid shapefile
# 5. Write out the unique weather grid coordinates and MUKEY values in a csv file (Output Data)
# (completed with all rasters in the list of files)

#VERY HELPFUL: http://ncss-tech.github.io/AQP/soilDB/gSSURGO-SDA.html

# My idea is to create a unique id for the WxID and Mukey combinations  and create an attribute table based on that! So it would be ID, wx_y, wx_x, MUKEY, GWADavg, CWADmin...
# Then create a raster based on ID
# (.tif and .tif.vat.dbf)

library(sp) # already installed
library(raster) # installed

library(rgeos) # installed
library(rgdal) # installed
library(maptools) # installed

library(data.table) # installed
library(foreign) # installed
library(plyr) # installed

############################################ The function!
myfunc <- function(inraster){
  
  print(inraster)
  # Read the raster
  r <- raster(inraster)
  #str(r) # To get more information about the raster
  r <- ratify(r)
  rat <- levels(r)[[1]] # get the raster attribute table from the raster 
  dbf_file <- gsub(".TIF",".TIF.vat.dbf",inraster) # read the .vat.dbf file which contains the raster attributes
  mu <- read.dbf(dbf_file, as.is=TRUE)
  names(mu)[toupper(names(mu)) == 'VALUE'] <- 'ID' # some rasters have Value, some VALUE so compare upper case
  mu$MUKEY <- as.integer(mu$MUKEY)
  #print("mu") # looks ok!
  #print(head(mu))
  #print(head(rat))
  
  rat.new <- join(rat, mu, by='ID', type='left') # join the actual attributes with the rat
  levels(r) <- rat.new # set the joined dataframe as the rat for the raster
  #print("RAT NEW")
  #print(head(rat.new)) # looks ok!
  
  r.mu <- deratify(r, att='MUKEY') # THIS IS CORRECT!
  #print("R.MU")
  #print(r.mu) # this looks ok!
  
  MUKEY<-extract(r.mu,1:ncell(r.mu))
  #print("MUKEY")
  #print(MUKEY)
  coord<-xyFromCell(r.mu,1:ncell(r.mu))
  pixels<-as.data.table(cbind(coord,MUKEY)) # define the pixels (coordinate and MUKEY value)
  #print("HSOULD HAVE MUKEY")
  #print(head(pixels))
  
  ##########################
  #b <- pixels
  
  coordinates(pixels) <- c("x", "y") # make this data table spatial 
  proj4string(pixels) <- proj4string(poly)
  df <- over(pixels, poly) # intersect the pixels with the Weather polygon to determine the closest weather grid to each pixel (very slow!)
  
  #bdf <- df
  
  pixels <- as.data.table(as.data.frame(pixels))
  #print(head(df))
  pixels[,wxID_x:=df$CENTROID_X] # new column with centroids (centroids are already columns in the df (they are not being calculated in this script))
  pixels[,wxID_y:=df$CENTROID_Y]
  #print(head(df$CENTROID_X))
  #print("Centroids")
  #print(head(pixels)) # looks ok
  #print(head(pixels))
  
  # Create unique ID for each pixel.
  rat <- pixels[complete.cases(pixels),] # return pixels with no missing values!
  #print("HERE BEFORE")
  #print(nrow(rat))
  #print(head(rat))
  rat[,c("x","y"):=NULL] # remove columns x and y
  #ratb <- rat
  setkeyv(rat,c("MUKEY","wxID_y","wxID_x")) 
  rat[ , Value := .GRP, by = key(rat)]
  setkeyv(rat,"Value") 
  rat <- unique(rat) # attribute table
  #print("first rat")
  #print(nrow(rat))
  #print(head(rat))
  
  #ratp <- rat
  
  # Format the RAT so that arcgis can read it (need unique ID, Value, and Count)
  # Join the unID with the raster (pixels)
  ratID <- rat[,c("Value","MUKEY","wxID_y","wxID_x"),with = FALSE] #create ratID
  setkeyv(ratID,c("MUKEY","wxID_y","wxID_x"))
  setkeyv(pixels,c("MUKEY","wxID_y","wxID_x"))
  pixels <- ratID[pixels]
  pixels[,Count:=.N, by=Value]
  
  setkeyv(rat,c("MUKEY","wxID_y","wxID_x"))
  pixels2 <- pixels[,c("MUKEY","wxID_y","wxID_x","Count"),with = FALSE]
  setkeyv(pixels2,c("MUKEY","wxID_y","wxID_x"))
  pixels2 <- unique(pixels2)
  rat <- pixels2[rat]
  
  #print("third rat")
  #print(head(rat))
  
  # need a Count and Value column for ArcGIS to read the raster correctly
  setcolorder(rat, c("Count",colnames(rat)[!(colnames(rat) %in% c("Count"))]))
  setcolorder(rat, c("Value",colnames(rat)[!(colnames(rat) %in% c("Value"))]))
  #print("REORDERED RAT")
  #print(head(rat))
  
  ras <- rasterFromXYZ(pixels[,c("x", "y", "Value"),with = FALSE])
  projection(ras) <- mycrs
  
  #print("rasterized!")
  #return (list(ras,rat,rat2))
  return (list(ras,rat))
  
}
#################################################################################################################

# Shapefile Data

#st <- 'ia'

rootdir <- "/mnt/home/rilllydi/Midwest/Soils_in_CDL/"
outdir <- "/mnt/home/rilllydi/Midwest/Results/Baseline_Rasters/"
wxpoly <- "/mnt/home/rilllydi/Midwest/WeatherPoly/"

wxname <- paste(wxpoly,st,"_WxPoly.shp",sep="")
poly <- readShapePoly(wxname)
mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(poly) <- mycrs

# Find all the raster data files 
folder <- paste(rootdir,st,"_split_55/",sep="")
#print(folder)
inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
#print(inraster)

####################
result <- list()
  
result <- lapply(inraster, myfunc)
  
rasters <- lapply(result,"[[", 1) 
ratlist <- lapply(result,"[[", 2) 
  
rout <- paste(outdir,st,"_base.TIF",sep="")
ratout <- gsub(".TIF",".TIF.vat.dbf",rout)

#print("HERE")
#print(rasters)
# mosaic rasters together and write
st.mosaicargs <- rasters
st.mosaicargs$fun <- mean
strast <- do.call(mosaic, st.mosaicargs)
#print("about to write raster")
writeRaster(strast, rout, datatype='INT2S', overwrite=TRUE)

#print("about to write rat")
rat <- rbindlist(ratlist)
write.dbf(as.data.frame(rat),ratout)
  
print("DONE")
