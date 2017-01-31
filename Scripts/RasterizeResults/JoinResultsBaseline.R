# THis is a test script to join the raster attribute table through the .vat.dbf file with SALUS results

library(data.table)

inraster <- "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/State_Rasters/ia_baseline.TIF"

dbf_file <- gsub(".TIF",".TIF.vat.dbf",inraster) # read the .vat.dbf file which contains the raster attributes
vat <- as.data.table(read.dbf(dbf_file, as.is=TRUE))
print(names(vat))
print(head(vat))

setkeyv(vat,c("MUKEY","wxID_y","wxID_x")) 
print("ABOUT TO JOIN")

SALUSres <- fread("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/results_csv/ia_SC3_endof2016.csv", header=TRUE)

setkeyv(SALUSres,c("MUKEY","wxID_y","wxID_x"))

print(head(SALUSres)) # this looks ok
vat <- SALUSres[vat]

# need a Count and Value column for ArcGIS to read the raster correctly
setcolorder(vat, c("Count",colnames(vat)[!(colnames(vat) %in% c("Count"))]))
setcolorder(vat, c("Value",colnames(vat)[!(colnames(vat) %in% c("Value"))]))

write.dbf(as.data.frame(vat),dbf_file)