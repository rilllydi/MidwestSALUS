# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The output is a data frame of the pixel coordinates, associated weather grid point coordinates, and STMUKEY value (STMUKEY is the MUKEY value with the state abbreviation in front)

# Run with split rasters

# Then the results from SALUS are added to the data frame (any pixels without results are given the value -9999)
# the data frame is used to create a raster (1 raster for each variable)

# Then the new rasters are mosaicked together for the state

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
myfunc <- function(x,num){
  
  # Read the raster
  r <- raster(x)
  
  # Find the coordinates of each pixel (centroids)
  d <- data.frame(coordinates(r)[!is.na(values(r)),])
  
  # Find the corresponding weather grid point for each pixel
  coordinates(d) <- c("x", "y")
  
  proj4string(d) <- proj4string(poly)
  d$wx <- over(d, poly) # This part is slow 
  
  # Extract the MUKEY INTEGER attribute from the raster for each pixel 
  #   and create a data frame that includes the coordinates of the pixel, its associated weather grid point, and the STMUKEY.
  finalout <- data.table(coordinates(d), d$wx.CENTROID_X, d$wx.CENTROID_Y,
                         extract(r, d,df=TRUE, factor=TRUE)[num]) # This part is slow 
  # combine wx_y and wx_x into 1 column
  finalout[,WxID:=wxIDfunc(V3,abs(V2))]
  
  setkey(finalout,WxID,MUKEY)

  # Left join on finalout
  finalout <- SALUSres[finalout]
  finaloutdf <- as.data.frame(finalout)
  
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  finaloutdf[is.na(finaloutdf)] <- -9999
  x <- finaloutdf
  e <- extent(x)
  coordinates(x) <- ~x+y
  
  ######################
  #datafields <- list('avgCWAD','avgGWAD') #"GWAD_1980","GWAD_1981","GWAD_1982","GWAD_1983","GWAD_1984","GWAD_1985","GWAD_1986",
  #"GWAD_1987","GWAD_1988","GWAD_1989", "GWAD_1990","GWAD_1991","GWAD_1992","GWAD_1993","GWAD_1994","GWAD_1995","GWAD_1996",
  #"GWAD_1997","GWAD_1998","GWAD_1999", "GWAD_2000","GWAD_2001","GWAD_2002","GWAD_2003","GWAD_2004","GWAD_2005","GWAD_2006",
  #"GWAD_2007","GWAD_2008","GWAD_2009","GWAD_2010","GWAD_2011","GWAD_2012","GWAD_2013","GWAD_2014","GWAD_2015")
  
  #for (field in datafields) {
  rout <- raster(e, ncol=ncol(x), nrow=nrow(x),res=res(r), crs=mycrs)
  gout <- rasterize(x,rout,field='avgGWAD')
  cout <- rasterize(x,rout,field='avgCWAD')
  nout <- rasterize(x,rout,field='maxNLCC')
  print("rasterized!")
  return (list(gout,cout,nout))
  
    #if (grepl("GWAD",field)){ #field == "*GWAD*"){
    #  outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/",state,"_",field,".tif",sep="")
    #  gout <- rasterize(x,rout,field=field)
    #} else if (grepl("CWAD",field)){ #(field == "*CWAD*"){
    #  outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_CWAD/",state,"_",field,".tif",sep="")
    #  cout <- rasterize(x,rout,field=field)
    #} else {
    #  outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_NLCC/",state,"_",field,".tif",sep="")
    #  nout <- rasterize(x,rout,field=field)
    #}
    #writeRaster(rout, outraster, overwrite=TRUE)
  #}
}

#################################################################################################################

# Shapefile Data

#wxgrid <- read.table("Z:/Users/rilllydi/MidwestSALUS/Weather/WxGrid.txt",header=TRUE)
#coordinates(wxgrid) <- ~ lon + lat
#wxgridsp <- SpatialPoints(wxgrid)

SC = "SC1"

#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('ia')

for (state in st) {
  
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  proj4string(poly) <- CRS("+init=epsg:4326")
  
  # for MUKEY INTEGER NUMBER Instead of STMUKEY
  if (state == 'ia'){
    num <- 6
  }else{
    num <- 4
  }
  
  # Find all the raster data files 
  folder <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_split_55/",sep="")
  print(folder)
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/ia_split_55/ia_soil_CornOnly_WGS84_8.TIF"
  
  ############################################ Read and add the SALUS results to the data frame
  
  #SALUSoutlist <- list()
  #for (chunk in 0:11){
  #for (chunk in 0:0){
  #  inSALUS <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/",state,"_",chunk,"_",SC,"_finalresults.csv",sep="")
  #  #SALUSres <- lapply(inSALUS,fread(inSALUS, header=TRUE))
  #  SALUStables <- lapply(inSALUS,rfun) 
  #}
  csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
  chunkpatt = paste(state,"_._",SC,"_finalresults_avg.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  #print(inSALUS)
  SALUStables <- list()
  SALUStables <- lapply(inSALUS,rfun)
  SALUSres <- rbindlist(SALUStables)
  SALUSres <- SALUSres[, MUKEY:=as.character(MUKEY)] # THIS IS NECESSARY
  setkey(SALUSres,WxID,MUKEY)
  
  ####################
  # Run the function on all the raster files
  result <- list()
  gresult <- list()
  cresult <- list()
  nresult <- list()
  
  result <- lapply(inraster, myfunc, num)
  gresult <- lapply(result,"[[", 1) 
  cresult <- lapply(result,"[[", 2) 
  nresult <- lapply(result,"[[", 3) 

  gout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/AvgGWAD_",state,".tif",sep="")
  cout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_CWAD/AvgCWAD_",state,".tif",sep="")
  nout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_NLCC/MaxNLCC_",state,".tif",sep="")
  
  mosaicme(gresult,gout)
  mosaicme(cresult,cout)
  mosaicme(nresult,nout)
  
}

################### End of Loop ###################################
# 1 avgGWAD, 1 avgCWAD, 1 maxNLCC raster

# For all the output rasters
folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/"
GWADras <- lapply(list.files(folder, pattern="AvgGWAD_.*\\.tif$", full.names=TRUE),raster)
GWADout <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_GWAD/AvgGWAD_Midwest.tif"
mosaicme(GWADras, GWADout)

folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_CWAD/"
CWADras <- lapply(list.files(folder, pattern="AvgCWAD_.*\\.tif$", full.names=TRUE),raster)
CWADout <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_CWAD/AvgCWAD_Midwest.tif"
mosaicme(CWADras, CWADout)

folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_NLCC/"
NLCCras <- lapply(list.files(folder, pattern="MaxNLCC_.*\\.tif$", full.names=TRUE),raster)
NLCCout <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_NLCC/MaxNLCC_Midwest.tif"
mosaicme(NLCCras, NLCCout)

######################################################################## EXTRA

# Create raster layers (stores some or all pizel values) and pull them into a raster stack? or raster brick is more efficient

#table1 <- finalout[,c("x","y","MUKEY"),with = FALSE]
#table2 <- finalout[,c("x","y","avgGWAD"),with = FALSE]
#table3 <- finalout[,c("x","y","avgCWAD"),with = FALSE]

#projection(table1) <- CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#coordinates(table1) <- ~ x + y
#xy <- coordinates(table1)
#gridded(table1) <- TRUE
#raster(table1) # THIS MIGHT CRASH IF BY STATE
