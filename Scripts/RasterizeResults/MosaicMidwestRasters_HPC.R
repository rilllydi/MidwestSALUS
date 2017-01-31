#!/usr/bin/Rscript
args = commandArgs(TRUE)

# This script gathers the state rasters and mosaics them into a Midwest raster!

# test if there is are 3 arguments: if not, return an error
if (length(args)!=3) {
  stop("3 arguments must be supplied", call.=FALSE)
} else {
  SC <- args[1]
  resvar <- args[2] # yearly_GWAD or endvalue 
  resdir <- args[3] # State_Rasters/Yearly/ or State_Rasters/Endof2016/
}

# Main Function
myfunc <- function(inraster){
  r <- raster(inraster) 
  r <- ratify(r)
  rat <- levels(r)[[1]]
  print(rat)
  return (list(ras,rat))
}

# Find all the state rasters for the variable we want
rootdir <- paste("/mnt/home/rilllydi/Midwest/Results/",resdir,sep="")
patt <- paste("*",SC,"_",resvar,".TIF$",sep="")
print(patt)
inraster <- list.files(rootdir, pattern=patt, full.names=TRUE)
print(inraster)

# Call the main function
result <- list()
result <- lapply(inraster, myfunc)
rasters <- lapply(result,"[[", 1) 
ratlist <- lapply(result,"[[", 2) 
print("succesfully ran function")

# Name the output raster
outdir <- "/mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/"
rout <- paste(outdir,SC,"_",resvar,".TIF",sep="")
ratout <- gsub(".TIF",".TIF.vat.dbf",rout)

# mosaic rasters together and write
st.mosaicargs <- rasters
st.mosaicargs$fun <- mean
strast <- do.call(mosaic, st.mosaicargs)
writeRaster(strast, rout, datatype='INT2S', overwrite=TRUE)

# write out raster attribute table
rat <- rbindlist(ratlist)
write.dbf(as.data.frame(rat),ratout)