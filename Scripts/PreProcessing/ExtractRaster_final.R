# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The output is a data table of the unique combinations of weather grid point coordinates and MUKEY values of the SSURGO soils

# Note using a statewide raster may crash R
# Note this process was checked by creating rasters for multiple steps and comparing with the original raster (all.equal(testout,r.mu))
# Note you need the associated .vat.dbf file along with the original raster for this to work

# Input Data:
#   - raster with MUKEY values
#   - polygon shapefile of the NLDAS weather grid (with the attribute table including the centroid coordinates)

# Process:
# 1. Read in the raster along with the associated dbf file (without the dbf file you get incorrect MUKEY values!)
# 2. Create a new raster of just the MUKEY values
# 3. Create a data table of the pixel coordinates and associated MUKEY value
# 4. Overlay the data table the the weather grid shapefile
# 5. Write out the unique weather grid coordinates and MUKEY values in a csv file (Output Data)
# (completed with all rasters in the list of files)

#VERY HELPFUL: http://ncss-tech.github.io/AQP/soilDB/gSSURGO-SDA.html

library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(foreign)
library(plyr)

#install.packages('data.table')
setwd("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly")

############################################
# MY FUNCTION
############################################
myfunc <- function(inraster){

  print(inraster)
  # Read the raster
  r <- raster(inraster)
  #str(r) # To get more information about the raster
  r <- ratify(r)
  rat <- levels(r)[[1]]
  dbf_file <- gsub(".TIF",".TIF.vat.dbf",inraster)
  mu <- read.dbf(dbf_file, as.is=TRUE)
  names(mu)[1] <- 'ID'
  mu$MUKEY <- as.integer(mu$MUKEY)

  rat.new <- join(rat, mu, by='ID', type='left')
  levels(r) <- rat.new
  r.mu <- deratify(r, att='MUKEY') # THIS IS CORRECT!!!!!!!!!!!!!!!!!!! YOU CAN TRUST THESE VALUES
  
  MUKEY<-extract(r.mu,1:ncell(r.mu))
  coord<-xyFromCell(r.mu,1:ncell(r.mu))
  pixels<-as.data.table(cbind(coord,MUKEY))

  #gout <- paste('Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/Test/',state,'_wx_11_mukey_fromdbf.tif',sep="")
  #writeRaster(r.mu, gout, overwrite=TRUE) # CORRECT
  
  #output <- paste('Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/Test/',state,'_wx_11_pixels_fromdbf.csv',sep="")
  #write.table(pixels,file = output,sep = ",")
  
  #all.equal(testout,r.mu) # TRUE!!!!!!!!!!!!
  
  ##########################
  
  coordinates(pixels) <- c("x", "y")
  proj4string(pixels) <- proj4string(poly)
  df <- over(pixels, poly)
  pixels <- as.data.table(as.data.frame(pixels))
  pixels[,wx_x:=df$CENTROID_X]
  pixels[,wx_y:=df$CENTROID_Y]
  
  pixels <- na.omit(pixels)
  setkeyv(pixels,c("MUKEY","wx_y","wx_x")) 
  uniqRes <- unique(pixels)
  uniqRes <- uniqRes[,c("MUKEY","wx_y","wx_x"),with = FALSE]
  
  return(uniqRes)
}
###### END FUNCTION ######

########################################### 
# MAIN PROGRAM
########################################### 
#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('mi','wi','oh','sd','mn','mo')
#st <- 'il'

for (state in st) {
  print("Made it to a new state")
  
  # Shapefile Data
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(poly) <- mycrs
  
  # Find all the raster data files 
  folder <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_split_55/",sep="")
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  
  #inraster <- ("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/in_split_55/in_soil_CornOnly_WGS84_11.TIF")
  ####################
  # Run the function on all the raster files
  result <- list()
  result <- lapply(inraster, myfunc)
  print("made it through the function")

  # Once we have unique data frames for each split raster, we can combine the data frames and again take the unique values
  # This final data frame (one per state) will be written to a csv file.
  finalout <- rbindlist(result)
  setkeyv(finalout,c("MUKEY","wx_y","wx_x")) 
  finaluniq <- unique(finalout)
  output <- paste('Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/',state,'_wx_soil_unique.csv',sep="")
  header <- "ID,MUKEY,WX_X,WX_Y\n"
  cat(header, file=output)
  write.table(finaluniq,file = output,sep = ",", append=TRUE, col.names=FALSE)
  
}
