##########################################################################################################################
##################################################### WINNER ############################################################
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
###############################################################

inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/ia_split_55/ia_soil_CornOnly_WGS84_8.TIF"
r <- raster(inraster)
poly <- readShapePoly("Z:/Users/rilllydi/MidwestSALUS/Weather/ia_WxPoly.shp")
proj4string(poly) <- CRS("+init=epsg:4326")
mycrs="+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
csvfolder <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/"
state<-'ia'
SC <- 'SC1'
field <- 'avgGWAD'
#chunkpatt = paste(state,"_._",SC,"_finalresults_GWADyearly.csv$",sep="")
#list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
inSALUS <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/ia_5_SC1_finalresults_avg.csv"
SALUSres <- fread(inSALUS, header=TRUE)
SALUSres <- SALUSres[, MUKEY:=as.character(MUKEY)] # THIS IS NECESSARY
setkey(SALUSres,WxID,MUKEY)

dat <- as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=FALSE)
names(dat)[names(dat) == 'ia_soil_CornOnly_WGS84_8_MUKEY'] <- 'MUKEY'
dat <- dat[c("MUKEY","x","y")]
coordinates(dat) <- ~ x + y
proj4string(dat) <- proj4string(poly)
datc <- dat
dat$wx <- over(dat, poly) # This part is slow 
dat <- as.data.table(as.data.frame(dat))
dat[,wx_x:=wx.CENTROID_X]
dat[,wx_y:=wx.CENTROID_Y]
dat <- dat[,c("MUKEY","y","x","wx_x","wx_y"),with = FALSE]
dat[,WxID:=wxIDfunc(wx_y,abs(wx_x))]

setkey(dat,WxID,MUKEY)
dat <- SALUSres[dat]
#dat <- as.data.frame(dat)
#dat[is.na(dat)] <- -9999
coordinates(dat) <- ~ x + y
######################

system.time(coordinates(dat) <- ~ x + y)
system.time(coordinates(dat) <- cbind(x,y))

rastfunc <- function(r) {
  #dat1 <- dat[c("MUKEY")]
  dat1$MUKEY <- as.integer(dat1$MUKEY)
  #datt <- as.data.table(dat1)
  #coordinates(dat) <- ~ x + y
  gout1 <- rasterFromXYZ(dat) #values=TRUE FIGURE OUT HOW TO GET THE VALUES YOU WANT! CREATE RASTER STACK? AND THEN WRITE OUT RASTER STACK IN ONE LINE?
  projection(gout1) <- mycrs
  plot(gout1)
  return(gout1)
}
# FAILS
rastfunc2 <- function(r) {
  r3 <- raster(extent(r), ncol=ncol(dat), nrow=nrow(dat),res=res(r))  
  r3 <- raster(nrow=nrow(dat), ncol=ncol(dat))
  r3 <- raster(extent(r), res=res(r))
  cells <- cellFromXY(r,coordinates(dat))
  r[cells] <- dat[1,]
  projection(r3) <- mycrs
  return(r3)
}

# CHECK WHICH RASTER CONVERSION IS FASTER
system.time(rastfunc(r)) #5.70 seconds
system.time(rastfunc2(r)) #

#Raster:::.quickStack

writeRaster(gout1, "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/test_ia_xyz.tif", overwrite=TRUE)
##########################################################################################################################
##########################################################################################################################
# CHECK FOR A FASTER OVERLAY
library(rgeos)

system.time(pp1 <- over(datc, poly)) # 137.58
system.time(pp2 <- datc[apply(gIntersects(datc, poly, byid=TRUE), 2, any),]) 
system.time(pp3 <- datc[apply(gIntersects(datc, poly, byid=TRUE), 1, sum)==1,])
system.time(pp4 <- gIntersects(datc, poly)) #
#system.time(pp5 <- point.in.polygon)

# THEN CHECK IF YOU CAN OVERLAY WITH DATA TABLE INSTEAD
myover2 <- function(r){
  dat <- as.data.table(as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=FALSE))
  dat[,MUKEY:=mo_soil_CornOnly_WGS84_16_MUKEY]
  dat <- dat[,c("MUKEY","y","x"),with = FALSE]
  coordinates(dat) <- ~ x + y
  proj4string(datt) <- proj4string(poly)
  #dat$wx <- over(dat[1:2], poly) # This part is slow
  ndatt <- over(datt, poly) # This part is slow
  datt[,new:=over(datt, poly)] # This part is slow
  datt[,new:=as.data.table(over(datt, poly))] # This part is slow
  
}


myover <- function(r){
  dat <- as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=FALSE)
  names(dat)[names(dat) == 'ia_soil_CornOnly_WGS84_8_MUKEY'] <- 'MUKEY'
  dat <- dat[c("MUKEY","x","y")]
  coordinates(dat) <- ~ x + y
  proj4string(dat) <- proj4string(poly)
  dat$wx <- over(dat, poly) # This part is slow
}

system.time(myover(r))
system.time(myover2(r))

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
###############################################################

inraster <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/States/CornOnly/mo_split_55/mo_soil_CornOnly_WGS84_16.TIF"
r <- raster(inraster)

# ORIGINAL # user: 87.02 seconds
orig <- function(r){
  d <- data.frame(coordinates(r)[!is.na(values(r)),])
  coordinates(d) <- c("x", "y")
  proj4string(d) <- proj4string(poly)
  d$wx <- over(d, poly) # This part is slow 
  df1 <- data.table(coordinates(d), d$wx.CENTROID_X, d$wx.CENTROID_Y,
                    extract(r, d,df=TRUE, factor=TRUE)[4]) # This part is slow 
}

poly <- readShapePoly("Z:/Users/rilllydi/MidwestSALUS/Weather/mo_WxPoly.shp")
proj4string(poly) <- CRS("+init=epsg:4326")
system.time(orig(r))

# 1.1 User = 86.71
func11 <- function(r){
  vv <- as.data.table(as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=TRUE))
  vv[,MUKEY:=mo_soil_CornOnly_WGS84_16_MUKEY]
  vv <- vv[,c("MUKEY","y","x"),with = FALSE]
  d <- as.data.frame(vv)
  coordinates(d) <- ~ x + y
  proj4string(d) <- proj4string(poly)
  d$wx <- over(d, poly) # This part is slow 
  d <- as.data.table(as.data.frame(d))
}
system.time(func1.1(r))

# 1.2 FAILS
func12 <- function(r){
  vv <- as.data.table(as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=TRUE))
  vv[,MUKEY:=mo_soil_CornOnly_WGS84_16_MUKEY]
  vv <- vv[,c("MUKEY","y","x"),with = FALSE]
  d <- as.data.frame(vv)
  coordinates(d) <- ~ x + y
  proj4string(d) <- proj4string(poly)
  d$wx <- over(d, poly) # This part is slow
  vv[,newx:=d$wx.CENTROID_X]
  vv[,newy:=d$wx.CENTROID_Y]
}
system.time(func12(r))

# 1.3 User = 80.92
func13 <- function(r){
  vv <- as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  names(vv)[names(vv) == 'mo_soil_CornOnly_WGS84_16_MUKEY'] <- 'MUKEY'
  vv <- vv[c("MUKEY","x","y")]
  coordinates(vv) <- ~ x + y
  proj4string(vv) <- proj4string(poly)
  vv$wx <- over(vv, poly) # This part is slow 
  d <- as.data.table(as.data.frame(vv))
  d[,wx_x:=wx.CENTROID_X]
  d[,wx_y:=wx.CENTROID_Y]
  df3 <- d[,c("MUKEY","y","x","wx_x","wx_y"),with = FALSE]
}
system.time(func13(r))

# MAKE SURE ORIGINAL AND 1.3 give the same results
df31 <- df3[,c("x","y","V2","V3","MUKEY"),with=FALSE]
identical(df1, df31) # YES THEY ARE THE SAME


#


############ ORIGINAL (REMOVE NA)
func14 <- function(r){
  dat <- as.data.frame(r, xy=TRUE, centroids=TRUE, na.rm=TRUE)
  names(dat)[names(dat) == 'ia_soil_CornOnly_WGS84_8_MUKEY'] <- 'MUKEY'
  dat <- dat[c("MUKEY","x","y")]
  coordinates(dat) <- ~ x + y
  proj4string(dat) <- proj4string(poly)
  dat$wx <- over(dat, poly) # This part is slow 
  dat <- as.data.table(as.data.frame(dat))
  dat[,wx_x:=wx.CENTROID_X]
  dat[,wx_y:=wx.CENTROID_Y]
  dat <- dat[,c("MUKEY","y","x","wx_x","wx_y"),with = FALSE]
  dat[,WxID:=wxIDfunc(wx_y,abs(wx_x))]
  
  setkey(dat,WxID,MUKEY)
  dat <- SALUSres[dat]
  copy <- as.data.frame(dat)
  dat <- as.data.frame(dat)
  #copy <- dat
  dat[is.na(dat)] <- -9999
  coordinates(dat) <- ~ x + y
  ######################
  
  #rastfunc <- function(fieldat) {
  projection(r) <- mycrs
  #rout <- raster(extent(r), ncol=ncol(dat), nrow=nrow(dat),res=res(r))
  #projection(rout) <- mycrs
  #print("rasterizing")
  #gout <- rasterize(dat,rout,field='avgGWAD')
  gout <- rasterize(dat,r,field=field)
  #  return(gout)
  #}
  #rastfunc('avgGWAD')
}
system.time(func14(r)) #434.05


##########################################################################################################################

# 3
s <- stack( sapply(1:5, function(i) setValues(r, rnorm(ncell(r), i, 3) )) )
s[1:3]<-NA
vals<-extract(s,1:ncell(s))
coord<-xyFromCell(s,1:ncell(s))
combine<-cbind(coord,vals)
write.table(combine,"xyvalues.txt") 

############################### PARALLEL PROCESSING ##########################################
# 4
library(snow)
z=vector('list',4)
z=1:4
system.time(lapply(z,function(x) Sys.sleep(1)))
cl<-makeCluster(###YOUR NUMBER OF CORES GOES HERE ###,type="SOCK")
  system.time(clusterApply(cl, z,function(x) Sys.sleep(1)))
  stopCluster(cl)
  