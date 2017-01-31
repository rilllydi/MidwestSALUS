# This script write the shell scripts to run the R script to create state rasters for each state, scenario, and variable (writeState)
# This script also write the shell scripts to run the R script to mosaic the state rasters into 1 Midwest Raster (writeMidwest)

from __future__ import print_function
import csv
import string
import sys
import time
import datetime
import os

####################################################

def writeBaseline (OutName, st):

  ScratchDir = "/mnt/scratch/rilllydi/Midwest"
  OutFile = open(OutName, "w")
  OutFile.write("#!/bin/sh -login\n")
  OutFile.write("#PBS -l nodes=1:ppn=1,walltime=10:00:00,mem=20gb\n")
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
  OutFile.write("Rscript --no-save MidwestRaster_base.R " + st + "\n")
  OutFile.close()

####################################################

states = ['wi','mi','oh','in','il','ia','sd','mn','mo']

BatchFileName = "/mnt/home/rilllydi/Midwest/Results/ShellScripts/allsubBase.bat"
BatchFile = open(BatchFileName, "w")

for st in states:
    OutName = st + "_Baseline.sh"
    writeBaseline(OutName, st)
    print(OutName)
    BatchFile.write("qsub " + OutName + "\n")
