# This script extracts data from the Soil/CDL raster. It also uses the Weather polygon shapefile to determine weather grid points for each raster pixel.
# The output is a data frame of the pixel coordinates, associated weather grid point coordinates, and STMUKEY value (STMUKEY is the MUKEY value with the state abbreviation in front)

# Try to run one state at a time instead of in split rasters

# Then the results from SALUS are added to the data frame (any pixels without results are given the value -9999)
# the data frame is used to create a state raster

# Note rasterizing the statewide data frame TAKES A LONG TIME

library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(rgeos)
#install.packages('data.table')
setwd("Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly")

############################################ The function!
myfunc <- function(x,num){
  
  # Read the raster
  r <- raster(x)
  #vv <- as.data.table(as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=TRUE))
  #vv <- vv[,c("x","y","mo_soil_CornOnly_WGS84_0_MUKEY"),with = FALSE]
  #vvsp <- vv[,c("x","y"),with = FALSE]
  #vvsp <- SpatialPoints(vvsp)
  
  #coordinates(vv) <- ~ x + y
  # for each x and y in vv, find the nearest weather grid points (wxgrid)
  #vv$nearwx <- apply(gDistance(vvsp, wxgridsp, byid=TRUE), 1, which.min)
  
  # Find the coordinates of each pixel (centroids)
  d <- data.frame(coordinates(r)[!is.na(values(r)),])

  # Find the corresponding weather grid point for each pixel
  coordinates(d) <- c("x", "y")

  proj4string(d) <- proj4string(poly)
  d$wx <- over(d, poly) # This part is slow 

  # Extract the MUKEY INTEGER attribute from the raster for each pixel 
  #   and create a data frame that includes the coordinates of the pixel, its associated weather grid point, and the STMUKEY.
  result <- data.table(coordinates(d), d$wx.CENTROID_X, d$wx.CENTROID_Y,
                       extract(r, d,df=TRUE, factor=TRUE)[num]) # This part is slow 
  return(result)
}
########################################### Run the function

# Shapefile Data
poly = readShapePoly("Z:/Users/rilllydi/MidwestSALUS/Weather/MidwestWxPoly.shp")
proj4string(poly)=CRS("+init=epsg:4326")

#wxgrid <- read.table("Z:/Users/rilllydi/MidwestSALUS/Weather/WxGrid.txt",header=TRUE)
#coordinates(wxgrid) <- ~ lon + lat
#wxgridsp <- SpatialPoints(wxgrid)

SC = "SC1"

#st <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
st <- c('mo')

for (state in st) {
  
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
  ####################
  # Run the function on all the raster files
  result <- list()
  result <- lapply(inraster, myfunc, num)
  
  print("made it through the function. here is the list of results:")
  print(result)
  
  # Once we have unique data frames for each split raster, we can combine the data frames and again take the unique values
  # This final data frame (one per state) will be written to a csv file.
  finalout <- rbindlist(result) #use.names=TRUE
  head(finalout)
  
  # combine wx_y and wx_x into 1 column
  wxIDfunc <- function(V3,V2){
    WxID <- paste(V3,"N_",V2,"W",sep="")
    return(WxID)
  }
  finalout[,WxID:=wxIDfunc(V3,abs(V2))]
  
  ############################################ Read and add the SALUS results to the data frame
  # Read as data table
  rfun <- function(inSALUS){
    dat <- fread(inSALUS, header=TRUE)
    return(dat)
  }
  #SALUSoutlist <- list()
  #for (chunk in 0:11){
  #for (chunk in 0:0){
  #  inSALUS <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/",state,"_",chunk,"_",SC,"_finalresults.csv",sep="")
  #  #SALUSres <- lapply(inSALUS,fread(inSALUS, header=TRUE))
  #  SALUStables <- lapply(inSALUS,rfun) 
  #}
  csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
  chunkpatt = paste(state,"_._",SC,"_finalresults.csv$",sep="")
  inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
  #print(inSALUS)
  SALUStables <- list()
  SALUStables <- lapply(inSALUS,rfun)
  SALUSres <- rbindlist(SALUStables)
  SALUSres <- SALUSres[, MUKEY:=as.character(MUKEY)] # THIS IS NECESSARY
  
  setkey(SALUSres,WxID,MUKEY)
  setkey(finalout,WxID,MUKEY)
  
  # Left join on finalout
  finalout <- SALUSres[finalout]
  finaloutdf <- as.data.frame(finalout)
  
  mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  inras <- raster(inraster[1])
  finaloutdf[is.na(finaloutdf)] <- -9999
  x <- finaloutdf
  e <- extent(x)
  coordinates(x) <- ~x+y
  
  ######################
  datafields <- list('avgCWAD','avgGWAD') #"GWAD_1980","GWAD_1981","GWAD_1982","GWAD_1983","GWAD_1984","GWAD_1985","GWAD_1986",
  #"GWAD_1987","GWAD_1988","GWAD_1989", "GWAD_1990","GWAD_1991","GWAD_1992","GWAD_1993","GWAD_1994","GWAD_1995","GWAD_1996",
  #"GWAD_1997","GWAD_1998","GWAD_1999", "GWAD_2000","GWAD_2001","GWAD_2002","GWAD_2003","GWAD_2004","GWAD_2005","GWAD_2006",
  #"GWAD_2007","GWAD_2008","GWAD_2009","GWAD_2010","GWAD_2011","GWAD_2012","GWAD_2013","GWAD_2014","GWAD_2015")

  for (field in datafields) {
    r <- raster(e, ncol=ncol(x), nrow=nrow(x),res=res(inras), crs=mycrs)
    r <- rasterize(x,r,field=field)
    if (grepl("GWAD",field)){ #field == "*GWAD*"){
      outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_GWAD/",state,"_",field,".tif",sep="")
    } else if (grepl("CWAD",field)){ #(field == "*CWAD*"){
      outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_CWAD/",state,"_",field,".tif",sep="")
    } else {
      outraster <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_State_NLCC/",state,"_",field,".tif",sep="")
    }
    writeRaster(r, outraster, overwrite=TRUE)
  }
  
}
################### End of Loop ###################################
# 1 avgGWAD, 1 avgCWAD, 1 maxNLCC raster
# 1 raster for GWAD for each year (36 years)
# 1 raster for CWAD for each year (36 years)

# For all the output rasters
folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster/raster_State_GWAD/"
GWADras <- list.files(folder, pattern="*avgGWAD*.TIF$", full.names=TRUE)
GWADmos <- do.call(mosaic,GWADras)
GWADout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_GWAD/AvgGWAD_Midwest.tif",sep="")
writeRaster(GWADmos, GWADout, overwrite=TRUE)

folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster/raster_State_CWAD/"
CWADras <- list.files(folder, pattern="*avgCWAD*.TIF$", full.names=TRUE)
CWADmos <- do.call(mosaic,GWADras)
CWADout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_CWAD/AvgCWAD_Midwest.tif",sep="")
writeRaster(CWADmos, CWADout, overwrite=TRUE)

folder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster/raster_State_NLCC/"
NLCCras <- list.files(folder, pattern="*maxNLCC*.TIF$", full.names=TRUE)
NLCCmos <- do.call(mosaic,GWADras)
NLCCout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_NLCC/MaxNLCC_Midwest.tif",sep="")
writeRaster(NLCCmos, NLCCout, overwrite=TRUE)

#x$filename <- 'test.tif'
#x$overwrite <- TRUE
#m <- do.call(merge, GWADras) # or reduce if merge doesn't work
#m <- do.call(merge, GWADras)
#m <- do.call(merge, GWADras)
#Reduce(function(x, y) merge(x, y), GWADras)

# For each year
for (year in 1980:2015){
  print(year)
  patt = paste("*GWAD_",year,"*.tif$",sep="")
  ras <- 
  paste("GWADras_",year,sep="") <- list.files(folder, pattern=patt, full.names=TRUE)
  GWADmos <- do.call(mosaic,GWADras)
  GWADout <- paste("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/raster_Midwest_GWAD/AvgGWAD_Midwest.tif",sep="")
  writeRaster(GWADmos, GWADout, overwrite=TRUE)
  
}

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
