# This script write the shell scripts to run the R script to create state rasters for each state, scenario, and variable (writeState)
# This script also write the shell scripts to run the R script to mosaic the state rasters into 1 Midwest Raster (writeMidwest)

from __future__ import print_function
import csv
import string
import sys
import time
import datetime
import os

def writeState (OutName, st, SC, resfolder, resvar, outdir):

  ScratchDir = "/mnt/scratch/rilllydi/Midwest"
  OutFile = open(OutName, "w")
  OutFile.write("#!/bin/sh -login\n")
  OutFile.write("#PBS -l nodes=1:ppn=1,walltime=6:00:00,mem=20gb\n")
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
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results\n")
  OutFile.write("Rscript --no-save ResultRasters_Avg_final_HPC.R " + st + " " + SC + " " + resfolder + " " + resvar + " " + outdir + "\n")
  #OutFile.write("R < ResultRasters_Avg_final_HPC.R --no-save " + st + " " + SC + " " + resfolder + " " + resvar + "\n")
  OutFile.close()

####################################################

def writeMidwest (OutName, SC, resvar, resdir):

  ScratchDir = "/mnt/scratch/rilllydi/Midwest"
  OutFile = open(OutName, "w")
  OutFile.write("#!/bin/sh -login\n")
  OutFile.write("#PBS -l nodes=1:ppn=1,walltime=6:00:00,mem=20gb\n")
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
  OutFile.write("cd /mnt/home/rilllydi/Midwest/Results\n")
  OutFile.write("Rscript --no-save MosaicMidwestRasters_HPC.R " + SC + " " + resvar + " " + resdir + "\n")
  OutFile.close()

####################################################

states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
scenarios = ["SC1","SC2","SC3","SC4","SC5","SC6","SC7"]
yearvariables = ["GWAD"]
variables = ["GWAD","CWAD","NLCC","C_ResOrgBl","C_SloOrgBl","C_ActOrgBl","C_CO2","C_FertBl","NOAD","DRNC","PREC"] # NIAD

BatchFileName = "/mnt/home/rilllydi/Midwest/Results/ShellScripts/allsubR.bat"
BatchFile = open(BatchFileName, "w")

for st in states:
  for SC in scenarios:
    #OutName = st + "_" + SC + "_summary.sh"
    #resvar = "summary"
    #resfolder = "Summary_1979_2016"
    #writeOut(OutName, st, SC, resfolder, resvar)
    #print(OutName)
    #BatchFile.write("qsub " + OutName + "\n")

    for var in yearvariables:
	OutName = st + "_" + SC + "_" + var + "_yearly.sh"
	resvar = "yearly_" + var
	resfolder = "Yearly_1979_2016"
	outdir = "State_Rasters/Yearly/"
        writeState(OutName, st, SC, resfolder, resvar, outdir)

	print(OutName)
	BatchFile.write("qsub " + OutName + "\n")

    OutName = st + "_" + SC + "_endof2016.sh"
    resvar = "endvalue"
    resfolder = "End_of_2016"
    outdir = "State_Rasters/End_of_2016/"
    writeState(OutName, st, SC, resfolder, resvar, outdir)
    print(OutName)
    BatchFile.write("qsub " + OutName + "\n")

###############################################################

BatchFileName = "/mnt/home/rilllydi/Midwest/Results/ShellScripts/allsubMidwest.bat"
BatchFile = open(BatchFileName, "w")

for SC in scenarios:
    for var in yearvariables:
	OutName = "Midwest_" + SC + "_" + var + "_yearly.sh"
	resvar = "yearly_" + var
	resdir = "State_Rasters/Yearly/"
        writeMidwest(OutName, SC, resvar, resdir)

	print(OutName)
	BatchFile.write("qsub " + OutName + "\n")

    OutName = "Midwest_" + SC + "_endof2016.sh"
    resvar = "endvalue"
    resdir = "State_Rasters/End_of_2016/"
    writeMidwest(OutName, SC, resvar, resdir)
    print(OutName)
    BatchFile.write("qsub " + OutName + "\n")


