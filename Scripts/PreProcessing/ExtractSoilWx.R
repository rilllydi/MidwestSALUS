
library(raster)
library(rgdal)
library(sp)
library(maptools)
library(rgeos)
library(data.table)
library(rgeos)
library(foreign)
library(plyr)

rootdir <- "Z:/Users/rilllydi/MidwestSALUS/Soils_in_CDL/Baseline/"
outdir <- "Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/"

states <- c('mi','wi','oh','in','il','ia','sd','mn','mo')
#states <- ('mi')
for (st in states){
  print(st)
  invat <- paste(rootdir,st,"_base.TIF.vat.dbf",sep="")
  print(invat)
  out = paste(outdir,st,"_wx_soil_unique.csv",sep="")
  
  vat <- as.data.table(read.dbf(invat, as.is=TRUE))
  print(names(vat))
  
  vat <- vat[, c('MUKEY','wxID_y','wxID_x'), with=FALSE]
  vat <- unique(vat)
  
  
  write.csv(as.data.frame(vat),out,row.names=FALSE)
}