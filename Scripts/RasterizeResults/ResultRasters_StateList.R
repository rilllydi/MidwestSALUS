
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(rgeos)
#install.packages('data.table')
setwd("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly")

############################################
# combine wx_y and wx_x into 1 column
wxIDfunc <- function(V3,V2){
  WxID <- paste(V3,"N_",V2,"W",sep="")
  return(WxID)
}
############################################
# Read as data table
rfun <- function(inSALUS){
  dat <- fread(inSALUS, header=TRUE)
  return(dat)
}
###########################################
# mosaic rasters together and write
mosaicme <- function(result,gout){
  g.mosaicargs <- result
  g.mosaicargs$fun <- mean
  grast <- do.call(mosaic, g.mosaicargs)
  writeRaster(grast, gout, overwrite=TRUE)
}
############################################ The function!
myfunc <- function(inraster,num){
  r <- raster(inraster)
  dat <- as.data.table(as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=FALSE, row.names = NULL, optional = FALSE))
  colnames(dat)[num+2] <- 'MUKEY'
  dat <- dat[,c("MUKEY","y","x"),with = FALSE]
  coordinates(dat) <- c("x","y")
  proj4string(dat) <- proj4string(poly)
  
  df <- over(dat, poly)
  dat <- as.data.table(as.data.frame(dat))
  dat[,wx_x:=df$CENTROID_X]
  dat[,wx_y:=df$CENTROID_Y]
  dat[,WxID:=wxIDfunc(wx_y,abs(wx_x))]
  
  setkey(dat,WxID,MUKEY)
  dat <- SALUSres[dat] # join the data with the SALUS results
  
  #coordinates(dat) <- c("x","y")
  

  return (dat)
  
}
#################################################################################################################

# Shapefile Data
SC = "SC1"
#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('ia')

for (state in st) {
  
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  proj4string(poly) <- CRS("+init=epsg:4326")
  
  if (state == 'ia'){
    num <- 6
  }else{
    num <- 4
  }
  
  # Find all the raster data files 
  folder <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_split_55/",sep="")
  print(folder)
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  #inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/ia_split_55/ia_soil_CornOnly_WGS84_8.TIF"
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  #inraster <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_soil_CornOnly_WGS84.TIF",sep="")
  #inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/ia_split_55/ia_soil_CornOnly_WGS84_8.TIF"
  
  ############################################ Read and add the SALUS results to the data frame
  
  csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
  chunkpatt = paste(state,"_._",SC,"_finalresults_avg.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  print(inSALUS)
  SALUStables <- list()
  SALUStables <- lapply(inSALUS,rfun)
  SALUSres <- rbindlist(SALUStables)
  setkey(SALUSres,WxID,MUKEY)
  mylist <- list(SALUSres,csv2,csv3,csv4,csv5,csv6,csv7,csv8,csv9,csv10,csv11)
  SALUSbig <- rbindlist(mylist)
  SALUSres <- SALUSbig
  ####################
  result <- list()

  system.time(result <- lapply(inraster, myfunc, num))
  
  stateres <- rbindlist(result) #statewide data table
  
  gras <- rasterFromXYZ(stateres[,c("x", "y", "avgGWAD"),with = FALSE]) # CREATE RASTER STACK? AND THEN WRITE OUT RASTER STACK IN ONE LINE?
  cras <- rasterFromXYZ(stateres[,c("x", "y", "avgCWAD"),with = FALSE])
  nras <- rasterFromXYZ(stateres[,c("x", "y", "maxNLCC"),with = FALSE])
  
  projection(gras) <- mycrs
  projection(cras) <- mycrs
  projection(nras) <- mycrs
  
  gout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/AvgGWAD_",state,"list.tif",sep="")
  cout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_CWAD/AvgCWAD_",state,"list.tif",sep="")
  nout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_NLCC/MaxNLCC_",state,"list.tif",sep="")
  
  writeRaster(grast, gout, overwrite=TRUE)
  writeRaster(crast, cout, overwrite=TRUE)
  writeRaster(nrast, nout, overwrite=TRUE)
  
  
}