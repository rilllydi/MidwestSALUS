# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The SALUS results (per MUKEY and weather ID) are matched to each pixel grid point
# A new raster for average GWAD, CWAD, and max NLCC is create for each state.

# Note using a statewide raster may crash R???
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

library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(rgeos)
library(foreign)
library(plyr)
#install.packages('data.table')
setwd("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly")


############################################ The function!
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
  
  r.mu <- deratify(r, att='MUKEY') # THIS IS CORRECT!
  
  MUKEY<-extract(r.mu,1:ncell(r.mu))
  coord<-xyFromCell(r.mu,1:ncell(r.mu))
  pixels<-as.data.table(cbind(coord,MUKEY))
  
  ##########################
  
  coordinates(pixels) <- c("x", "y")
  proj4string(pixels) <- proj4string(poly)
  df <- over(pixels, poly)
  pixels <- as.data.table(as.data.frame(pixels))
  pixels[,wx_x:=df$CENTROID_X]
  pixels[,wx_y:=df$CENTROID_Y]
  
  ###backup <- pixels
  
  setkeyv(pixels,c("MUKEY","wx_y","wx_x")) 
  pixels <- SALUSres[pixels] # join the data with the SALUS results
  
  return (pixels)
  
}
#################################################################################################################

# Shapefile Data
SC = "SC1"
#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('ia')

for (state in st) {
  
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(poly) <- mycrs
  
  # Find all the raster data files 
  folder <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_split_55/",sep="")
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  ###inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/ia_split_55/ia_soil_CornOnly_WGS84_17.TIF"
  
  ############################ Read and add the SALUS results to the data frame
  
  csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
  chunkpatt = paste(state,"_._",SC,"_finalresults_avg.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  SALUStables <- list()
  
  # Read as data table
  rfun <- function(inSALUS){
    dat <- fread(inSALUS, header=TRUE)
    return(dat)
  }
  
  SALUStables <- lapply(inSALUS,rfun)
  SALUSres <- rbindlist(SALUStables)
  
  SALUSres <- fread("Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/ia_fakeSALUSdata.csv", header=TRUE)
  
  setkeyv(SALUSres,c("MUKEY","wxID_y","wxID_x")) 
  
  ####################
  result <- list()
  
  system.time(result <- lapply(inraster, myfunc)) #started 3:38, ended 

  result2 <- rbindlist(result)
  
  gras <- rasterFromXYZ(result2[,c("x", "y", "avgGWAD"),with = FALSE]) # CREATE RASTER STACK? AND THEN WRITE OUT RASTER STACK IN ONE LINE?
  cras <- rasterFromXYZ(result2[,c("x", "y", "avgCWAD"),with = FALSE])
  nras <- rasterFromXYZ(result2[,c("x", "y", "maxNLCC"),with = FALSE])
  
  projection(gras) <- mycrs
  projection(cras) <- mycrs
  projection(nras) <- mycrs
  
}