# This script write the shell scripts to run the R script to create state rasters for each state, scenario, and variable (writeState)
# This script also write the shell scripts to run the R script to mosaic the state rasters into 1 Midwest Raster (writeMidwest)

from __future__ import print_function
import csv
import string
import sys
import time
import datetime
import os

def writeState (OutName, st, SC):

  ScratchDir = "/mnt/scratch/rilllydi/Midwest"
  OutFile = open(OutName, "w")
  OutFile.write("#!/bin/sh -login\n")
  OutFile.write("#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=20gb\n")
  OutFile.write("#PBS -N " + OutName + "\n")
  OutFile.write("#PBS -j oe\n")
  OutFile.write("\n")
  OutFile.write("module load GNU/4.4.5\n")
  OutFile.write("module load OpenMPI/1.4.3\n")
  OutFile.write("module load R/3.2.0\n")
  OutFile.write("module load GEOS\n")
  OutFile.write("module load GDAL\n")
  OutFile.write("export R_LIBS_USER=~/R/library\n")
  OutFile.write("\n")
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results/Baseline_Rasters\n")
  OutFile.write("cp "+st+"_base.tif /mnt/home/rilllydi/Midwest/Results/State_Rasters/"+st+"_"+SC+"_endvalue.tif\n")
  OutFile.write("cp "+st+"_base.TIF.vat.dbf /mnt/home/rilllydi/Midwest/Results/State_Rasters/"+st+"_"+SC+"_endvalue.TIF.vat.dbf\n")
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results\n")
  OutFile.write("Rscript --no-save JoinResultsBaseline.R " + st + " " + SC + "\n")
  OutFile.close()

####################################################

def writeMidwest (OutName, SC, resvar, resdir):

  ScratchDir = "/mnt/scratch/rilllydi/Midwest"
  OutFile = open(OutName, "w")
  OutFile.write("#!/bin/sh -login\n")
  OutFile.write("#PBS -l nodes=1:ppn=1,walltime=1:00:00,mem=20gb\n")
  OutFile.write("#PBS -N " + OutName + "\n")
  OutFile.write("#PBS -j oe\n")
  OutFile.write("\n")
  OutFile.write("module load GNU/4.4.5\n")
  OutFile.write("module load OpenMPI/1.4.3\n")
  OutFile.write("module load R/3.2.0\n")
  OutFile.write("module load GEOS\n")
  OutFile.write("module load GDAL\n")
  OutFile.write("export R_LIBS_USER=~/R/library\n")
  OutFile.write("\n")
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results/Baseline_Rasters\n")
  OutFile.write("cp Midwest_base.tif /mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_" + SC + "_" + resvar + ".tif\n")
  OutFile.write("cp Midwest_base.tif.vat.dbf /mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_" + SC + "_" + resvar + ".tif.vat.dbf\n")
  OutFile.write("cp Midwest_base.tif.aux.xml /mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_" + SC + "_" + resvar + ".tif.aux.xml\n")
  OutFile.write("cp Midwest_base.tif.ovr /mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_" + SC + "_" + resvar + ".tif.ovr\n")
  OutFile.write("cp Midwest_base.tfw /mnt/home/rilllydi/Midwest/Results/Midwest_Rasters/Midwest_" + SC + "_" + resvar + ".tfw\n")
  OutFile.write("\n")
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results\n")
  OutFile.write("Rscript --no-save JoinMidwestRaster.R " + SC + " " + resdir + " " + resvar + "\n")
  OutFile.close()

####################################################

#states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
scenarios = ["SC1","SC2","SC3","SC4","SC5","SC6","SC7"]

#BatchFileName = "/mnt/home/rilllydi/Midwest/Results/ShellScripts/allsubJoinRasterResults.bat"
#BatchFile = open(BatchFileName, "w")

#for st in states:
#  for SC in scenarios:
#    OutName = st + "_" + SC + "_JoinRasterResults.sh"
    #writeState(OutName, st, SC)
    #print(OutName)
    #BatchFile.write("qsub " + OutName + "\n")

###############################################################

BatchFileName = "/mnt/home/rilllydi/Midwest/Results/ShellScripts/allsubJoinMidwest.bat"
BatchFile = open(BatchFileName, "w")

yearvariables = ["GWAD"]
variables = ["GWAD","CWAD","NLCC","C_ResOrgBl","C_SloOrgBl","C_ActOrgBl","C_CO2","C_FertBl","NOAD","DRNC","PREC"] # NIAD

for SC in scenarios:
    for var in yearvariables:
	resvar = "yearly_" + var
	resdir = "Yearly_1979_2016"
	OutName = "Midwest_" + SC + "_" + resvar + ".sh"
        writeMidwest(OutName, SC, resvar, resdir)
	print(OutName)
	BatchFile.write("qsub " + OutName + "\n")

    resvar = "endvalue"
    resdir = "End_of_2016"
    OutName = "Midwest_" + SC + "_" + resvar + ".sh"
    writeMidwest(OutName, SC, resvar, resdir)
    print(OutName)
    BatchFile.write("qsub " + OutName + "\n")

    resvar = "summary"
    resdir = "Summary_1979_2016"
    OutName = "Midwest_" + SC + "_" + resvar + ".sh"
    writeMidwest(OutName, SC, resvar, resdir)
    print(OutName)
    BatchFile.write("qsub " + OutName + "\n")
