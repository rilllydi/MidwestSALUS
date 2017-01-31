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
  rat <- levels(r)[[1]]
  dbf_file <- gsub(".TIF",".TIF.vat.dbf",inraster)
  mu <- read.dbf(dbf_file, as.is=TRUE)
  names(mu)[names(mu) == 'VALUE'] <- 'ID'
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
  
  setkeyv(pixels,c("MUKEY","wx_y","wx_x")) 
  pixels <- SALUSres[pixels] # join the data with the SALUS results
  
  # Create unique ID for each pixel.
  rat <- pixels[complete.cases(pixels),]
  rat[,c("x","y"):=NULL]
  setkeyv(rat,c("MUKEY","wxID_y","wxID_x")) 
  rat[ , Value := .GRP, by = key(rat)]
  setkeyv(rat,"Value") 
  rat <- unique(rat) # attribute table
  
  # Join the unID with the raster (pixels)
  ratID <- rat[,c("Value","MUKEY","wxID_y","wxID_x"),with = FALSE]
  setkeyv(ratID,c("MUKEY","wxID_y","wxID_x"))
  setkeyv(pixels,c("MUKEY","wxID_y","wxID_x"))
  pixels <- ratID[pixels]
  pixels[,Count:=.N, by=Value]
  
  #setkeyv(rat,c("MUKEY","wxID_y","wxID_x"))
  #pixels2 <- pixels[,c("MUKEY","wxID_y","wxID_x","Count"),with = FALSE]
  #setkeyv(pixels2,c("MUKEY","wxID_y","wxID_x"))
  #pixels2 <- unique(pixels2)
  #rat <- pixels2[rat]
  
  setcolorder(rat, c("Count",colnames(rat)[!(colnames(rat) %in% c("Count"))]))
  setcolorder(rat, c("Value",colnames(rat)[!(colnames(rat) %in% c("Value"))]))

  ras <- rasterFromXYZ(pixels[,c("x", "y", "Value"),with = FALSE])
  projection(ras) <- mycrs

  return (list(ras,rat))
  
}
#################################################################################################################

# Shapefile Data
# Arguments: state, scenario, results folder, variable (summary or yearly variable)

# test if there is are 4 arguments: if not, return an error
if (length(args)==4) {
  stop("4 arguments must be supplied (input file).n", call.=FALSE)
} else {
  st <- args[1]
  SC <- args[2]
  resfolder <- args[3]
  resvar <- args[4]
  outdir <- args[5]
}

#st <- 'ia'
#SC <- 'SC1'
#resfolder <- "Summary_1979_2016"
#resvar <- "summary" # "yearly_GWAD"

rootdir <- "/mnt/home/rilllydi/Midwest/Soils_in_CDL/"
outdir <- paste("/mnt/home/rilllydi/Midwest/Results/",outdir,sep="")
csvfolder <- paste("/mnt/home/rilllydi/Midwest/Results/",resfolder,"/",sep="")
wxpoly <- "/mnt/home/rilllydi/Midwest/WeatherPoly/"

wxname <- paste(wxpoly,state,"_WxPoly.shp",sep="")
poly <- readShapePoly(wxname)
mycrs=CRS("+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(poly) <- mycrs

# Find all the raster data files 
folder <- paste(rootdir,state,"_split_55/",sep="")
inraster <- list.files(folder, pattern="*.TIF$", full.names=TRUE)

############################ Read and add the SALUS results to the data frame
chunkpatt = paste(state,"_._",SC,"_",resvar,".csv$",sep="") # get all chunks
inSALUS <- list.files(csvfolder, pattern=chunkpatt, full.names=TRUE)
SALUStables <- list()

# Read as data table
rfun <- function(inSALUS){
  dat <- fread(inSALUS, header=TRUE)
  return(dat)
}

SALUStables <- lapply(inSALUS,rfun)
SALUSres <- rbindlist(SALUStables)

setkeyv(SALUSres,c("MUKEY","wxID_y","wxID_x")) 

####################
result <- list()

result <- lapply(inraster, myfunc)

rasters <- lapply(result,"[[", 1) 
ratlist <- lapply(result,"[[", 2) 

rout <- paste(outdir,state,"_summary.TIF",sep="")
ratout <- gsub(".TIF",".TIF.vat.dbf",rout)

# mosaic rasters together and write
st.mosaicargs <- rasters
st.mosaicargs$fun <- mean
strast <- do.call(mosaic, st.mosaicargs)
writeRaster(strast, rout, datatype='INT2S', overwrite=TRUE)

rat <- rbindlist(ratlist)
write.dbf(as.data.frame(rat),ratout)