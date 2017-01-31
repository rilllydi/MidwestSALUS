#!/usr/bin/Rscript
args = commandArgs(TRUE)

options(warn=1)

# Arguments: state, scenario, results folder, variable (summary or yearly variable), outdirectory for raster
# test if there is are 4 arguments: if not, return an error
if (length(args) != 2) {
  stop("2 argument must be supplied (state, SC)", call.=FALSE)
} else {
  st <- args[1]
  SC <- args[2]
}

library(data.table)
library(foreign)

inraster <- paste("/mnt/home/rilllydi/Midwest/Results/Baseline_Rasters/",st,"_",SC,"_endvalue.TIF",sep="")
rootdir <- "/mnt/home/rilllydi/Midwest/Results/"

#outavg = rootdir + "Summary_1979_2016/" + st + "_SC" + str(SC) + "_summary.csv" 
#outyr = rootdir + "Yearly_1979_2016/" + st + "_SC" + str(SC) + "_yearly_" + var + ".csv"  
endvalues = paste(rootdir,"End_of_2016/",st,"_",SC,"_endvalue.csv",sep="")

dbf_file <- gsub(".TIF",".TIF.vat.dbf",inraster) # read the .vat.dbf file which contains the raster attributes
vat <- as.data.table(read.dbf(dbf_file, as.is=TRUE))
print(names(vat))
#print(head(vat))

setkeyv(vat,c("MUKEY","wxID_y","wxID_x")) 

print("about to read in results")
SALUSres <- fread(endvalues, header=TRUE)

setkeyv(SALUSres,c("MUKEY","wxID_y","wxID_x"))

#print(head(SALUSres)) # this looks ok
print("about to join")
vat <- SALUSres[vat]

# need a Count and Value column for ArcGIS to read the raster correctly
setcolorder(vat, c("Count",colnames(vat)[!(colnames(vat) %in% c("Count"))]))
setcolorder(vat, c("Value",colnames(vat)[!(colnames(vat) %in% c("Value"))]))

write.dbf(as.data.frame(vat),dbf_file)

print("done")
