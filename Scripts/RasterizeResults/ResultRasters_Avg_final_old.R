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

# My idea is to create a unique id for the WxID and Mukey combinations  and create an attribute table based on that! So it would be ID, wx_y, wx_x, MUKEY, GWADavg, CWADmin...
# Then create a raster based on ID
# (.tif and .tif.vat.dbf)

library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(rgeos)
library(foreign)
library(plyr)

options(warn=1)
#install.packages('data.table')
#setwd("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly")

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
  print("mu") # looks ok!
  #print(head(mu))
  #print(head(rat))
  
  rat.new <- join(rat, mu, by='ID', type='left') # join the actual attributes with the rat
  levels(r) <- rat.new # set the joined dataframe as the rat for the raster
  print("RAT NEW")
  #print(head(rat.new)) # looks ok!
  
  r.mu <- deratify(r, att='MUKEY') # THIS IS CORRECT!
  print("R.MU")
  #print(r.mu) # this looks ok!
  
  MUKEY<-extract(r.mu,1:ncell(r.mu))
  #print("MUKEY")
  #print(MUKEY)
  coord<-xyFromCell(r.mu,1:ncell(r.mu))
  pixels<-as.data.table(cbind(coord,MUKEY)) # define the pixels (coordinate and MUKEY value)
  print("HSOULD HAVE MUKEY")
  print(head(pixels))
  
  ##########################
  
  coordinates(pixels) <- c("x", "y") # make this data table spatial 
  proj4string(pixels) <- proj4string(poly)
  plot(pixels)
  plot(poly)
  df <- over(pixels, poly) # intersect the pixels with the Weather polygon to determine the closest weather grid to each pixel (very slow!)
  pixels <- as.data.table(as.data.frame(pixels))
  pixels[,wx_x:=df$CENTROID_X] # new column with centroids (centroids are already columns in the df (they are not being calculated in this script))
  pixels[,wx_y:=df$CENTROID_Y]
  print(head(df$CENTROID_X))
  print("Centroids")
  print(head(pixels)) # looks ok
  pixels[,wx_y:=round(wx_y,4)]
  pixels[,wx_y:=round(wx_y,4)]
  print(head(pixels))
  
  setkeyv(pixels,c("MUKEY","wx_y","wx_x")) 
  print("ABOUT TO JOIN")
  print(head(SALUSres)) # this looks ok
  pixels <- SALUSres[pixels] # join the pixels with the SALUS results based on MUKEY and the weather grid location
  print("AFTER JOIN")
  print(head(pixels))
  print(head(SALUSres))
  
  # Get final version of RAT to write out with the raster:
  print("PIXELS")
  print(head(pixels))
  # Create unique ID for each pixel.
  rat <- pixels[complete.cases(pixels),] # return pixels with no missing values!
  print("HERE BEFORE")
  print(length(rat))
  print(head(rat))
  rat[,c("x","y"):=NULL] # remove columns x and y?
  print("HERE")
  print(length(rat))
  print(head(rat))
  setkeyv(rat,c("MUKEY","wxID_y","wxID_x")) 
  rat[ , Value := .GRP, by = key(rat)]
  setkeyv(rat,"Value") 
  rat <- unique(rat) # attribute table
  print("first rat")
  print(length(rat))
  print(head(rat))
  
  # Join the unID with the raster (pixels)
  ratID <- rat[,c("Value","MUKEY","wxID_y","wxID_x"),with = FALSE]
  setkeyv(ratID,c("MUKEY","wxID_y","wxID_x"))
  setkeyv(pixels,c("MUKEY","wxID_y","wxID_x"))
  pixels <- ratID[pixels]
  pixels[,Count:=.N, by=Value]
  
  #print("second rat")
  #print(rat)
  #setkeyv(rat,c("MUKEY","wxID_y","wxID_x"))
  #pixels2 <- pixels[,c("MUKEY","wxID_y","wxID_x","Count"),with = FALSE]
  #setkeyv(pixels2,c("MUKEY","wxID_y","wxID_x"))
  #pixels2 <- unique(pixels2)
  #rat <- pixels2[rat]
  
  print("third rat")
  print(head(rat))
  
  # join pixels with yearly data
  #pixels2 <- pixels[,c("Value","Count","MUKEY","wxID_y","wxID_x","x","y"),with = FALSE]
  #setkeyv(pixels2,c("MUKEY","wxID_y","wxID_x")) 
  #pixels2 <- pixels2[SALUSyear]
  #setkeyv(pixels2,c("MUKEY","wxID_y","wxID_x")) 
  #rat2 <- unique(pixels2[complete.cases(pixels2),]) # attribute table
  #rat2[,c("x","y"):=NULL]
  
  # need a Count and Value column for ArcGIS to read the raster correctly
  #setcolorder(rat, c("Count",colnames(rat)[!(colnames(rat) %in% c("Count"))]))
  #setcolorder(rat, c("Value",colnames(rat)[!(colnames(rat) %in% c("Value"))]))
  #print("REORDERED RAT")
  #print(head(rat))
  
  ras <- rasterFromXYZ(pixels[,c("x", "y", "Value"),with = FALSE])
  projection(ras) <- mycrs

  print("rasterized!")
  #return (list(ras,rat,rat2))
  return (list(ras,rat))
  
}
#################################################################################################################

# Shapefile Data
SC = "SC3"
#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('mi')
rootdir <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/"
outdir <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/State_Rasters/"
csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"

for (state in st) {
  
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(poly) <- mycrs
  
  # Find all the raster data files 
  folder <- paste(rootdir,state,"_split_55/",sep="")
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  inraster <- c("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/mi_split_55/mi_soil_CornOnly_WGS84_8.TIF","Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/mi_split_55/mi_soil_CornOnly_WGS84_9.TIF")
  print(inraster)
  
  ############################ Read and add the SALUS results to the data frame
  chunkpatt = paste(state,"_._",SC,"_yearly_GWAD.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  SALUStables <- list()

  # Read as data table
  rfun <- function(inSALUS){
    dat <- fread(inSALUS, header=TRUE)
    return(dat)
  }
  
  SALUStables <- lapply(inSALUS,rfun)
  print(inSALUS)
  SALUSres <- rbindlist(SALUStables)
  
  #SALUSres <- fread("Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/ia_fakeSALUSdata.csv", header=TRUE)
  #SALUSyear <- fread("Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/ia_fakeSALUSdata.csv", header=TRUE)
  
  setkeyv(SALUSres,c("MUKEY","wxID_y","wxID_x")) 
  #setkeyv(SALUSyear,c("MUKEY","wxID_y","wxID_x")) 
  
  ####################
  result <- list()
  
  system.time(result <- lapply(inraster, myfunc))
  
  rasters <- lapply(result,"[[", 1) 
  ratlist <- lapply(result,"[[", 2) 
  #rat2list <- lapply(result,"[[", 3) 
  
  rout <- paste(outdir,state,"_summary.TIF",sep="")
  #rout2 <- paste(outdir,state,"_yearly.TIF",sep="")
  ratout <- gsub(".TIF",".TIF.vat.dbf",rout)
  #rat2out <- gsub(".TIF",".TIF.vat.dbf",rout2)
  
  # mosaic rasters together and write
  st.mosaicargs <- rasters
  st.mosaicargs$fun <- mean
  strast <- do.call(mosaic, st.mosaicargs)
  writeRaster(strast, rout, datatype='INT2S', overwrite=TRUE)
  #writeRaster(strast, rout2, datatype='INT2S', overwrite=TRUE)
  
  rat <- rbindlist(ratlist)
  write.dbf(as.data.frame(rat),ratout)
  
  #rat2 <- rbindlist(rat2list)
  #write.dbf(as.data.frame(rat2),rat2out)

}