#!/usr/bin/Rscript
args = commandArgs(TRUE)

options(warn=1)

# Arguments: state, scenario, results folder, variable (summary or yearly variable), outdirectory for raster
# test if there is are 4 arguments: if not, return an error
if (length(args) != 3) {
  stop("3 arguments must be supplied (SC, resdir, resvar)", call.=FALSE)
} else {
  SC <- args[1]
  resdir <- args[2]
  resvar <- args[3]
}

library(data.table)
library(foreign)

rootdir <- "/mnt/home/rilllydi/Midwest/Results/"

SALUSres = paste(rootdir, resdir, "/Midwest_", SC, "_", resvar, ".csv",sep="")

inraster <- paste("/mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_", SC, "_", resvar, ".tif",sep="")
dbf_file <- gsub(".tif", ".tif.vat.dbf", inraster) # read the .vat.dbf file which contains the raster attributes

vat <- as.data.table(read.dbf(dbf_file, as.is=TRUE))
print(names(vat))

setkeyv(vat,c("MUKEY", "wxID_y", "wxID_x")) 

print("about to read in results")
SALUSres <- fread(SALUSres, header=TRUE)

setkeyv(SALUSres,c("MUKEY", "wxID_y", "wxID_x"))

print("about to join")
vat <- SALUSres[vat]

# need a Count and Value column for ArcGIS to read the raster correctly
setcolorder(vat, c("Count",colnames(vat)[!(colnames(vat) %in% c("Count"))]))
setcolorder(vat, c("Value",colnames(vat)[!(colnames(vat) %in% c("Value"))]))

write.dbf(as.data.frame(vat),dbf_file)

print("done")
