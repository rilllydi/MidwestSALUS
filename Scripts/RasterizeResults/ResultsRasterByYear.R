# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The output is a data frame of the pixel coordinates, associated weather grid point coordinates, and STMUKEY value (STMUKEY is the MUKEY value with the state abbreviation in front)

# Try to run one state at a time instead of in split rasters

# Then the results from SALUS are added to the data frame (any pixels without results are given the value -9999)
# the data frame is used to create a raster

# Note using a statewide raster may crash R

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
  #print("reading raster")
  inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/mo_split_55/mo_soil_CornOnly_WGS84_0.TIF"
  x <- inraster
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
  
  mycrs="+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  finaloutdf[is.na(finaloutdf)] <- -9999
  x <- finaloutdf
  e <- extent(r)
  coordinates(x) <- ~x+y
  ######################
  datafields <- list("GWAD_1980") #,"GWAD_1981","GWAD_1982") #,"GWAD_1983","GWAD_1984","GWAD_1985","GWAD_1986",
  #"GWAD_1987","GWAD_1988","GWAD_1989", "GWAD_1990","GWAD_1991","GWAD_1992","GWAD_1993","GWAD_1994","GWAD_1995","GWAD_1996",
  #"GWAD_1997","GWAD_1998","GWAD_1999", "GWAD_2000","GWAD_2001","GWAD_2002","GWAD_2003","GWAD_2004","GWAD_2005","GWAD_2006",
  #"GWAD_2007","GWAD_2008","GWAD_2009","GWAD_2010","GWAD_2011","GWAD_2012","GWAD_2013","GWAD_2014","GWAD_2015")
  gout <- list()
  
  rastfunc <- function(field) {
    rout <- raster(extent(r), ncol=ncol(x), nrow=nrow(x),res=res(r))
    projection(rout) <- mycrs
    gout <- rasterize(x,rout,field=field)
  }
  rastfunc2 <- function(field) {
    gout2 <- rasterize(x,r,field=field)
  }
  
  
  system.time(rastfunc(field)) # user = 326.04
  plot(gout)
  
  system.time(rastfunc2(field)) # user = 324.55
  plot(gout2, add=T)
  
  
  gout <- lapply(datafields,rastfunc)
  return (gout)
  

}


# writeRaster(gout, "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/test_ia_rast5.tif", overwrite=TRUE)

#################################################################################################################

# Shapefile Data

SC = "SC1"
yearlist <- c(1980:2015)

#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('mo')

for (state in st) {
  wxname <- paste("Z:/Users/rilllydi/MidwestSALUS/Weather/",state,"_WxPoly.shp",sep="")
  poly <- readShapePoly(wxname)
  proj4string(poly) <- CRS("+init=epsg:4326")
  
  # for MUKEY INTEGER NUMBER Instead of STMUKEY
  if (state == 'ia'){
    num <- 5
  }else{
    num <- 4
  }
  
  print("Made it to a new state")
  
  # Find all the raster data files 
  folder <- paste("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/",state,"_split_55/",sep="")
  print(folder)
  inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)
  #inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/mo_split_55/mo_soil_CornOnly_WGS84_0.TIF"
  
  ############################################ Read and add the SALUS results to the data frame
  
  csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
  chunkpatt = paste(state,"_._",SC,"_finalresults_GWADyearly.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  SALUStables <- list()
  SALUStables <- lapply(inSALUS,rfun)
  SALUSres <- rbindlist(SALUStables)
  SALUSres <- SALUSres[, MUKEY:=as.character(MUKEY)] # THIS IS NECESSARY
  setkey(SALUSres,WxID,MUKEY)
  
  ####################
  # Run the function on all the raster files
  result <- list()
  result <- lapply(inraster, myfunc, num)
  # resultcopy <- result

  finalfunc <- function(year){
    num <- year - 1979
    resultyr <- lapply(result,"[[", num) 
    output <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/GWAD_",state,"_",year,".tif",sep="")
    mosaicme(resultyr,output)
  }
  
  yearlist <- list(1980)
  lapply(yearlist,finalfunc)

}

mosMidwest <- function(year) {
  folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/"
  yearpatt <- paste("GWAD_.*",year,".tif$",sep="")
  GWADras <- lapply(list.files(folder, pattern=yearpatt, full.names=TRUE),raster)
  GWADout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_GWAD/Midwest_GWAD_",year,".tif",sep="")
  mosaicme(GWADras, GWADout)
}

lapply(yearlist,mosMidwest)