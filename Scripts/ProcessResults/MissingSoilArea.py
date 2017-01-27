# This script determines the percent of area in the state that we do not have MUKEY characteristic data available

import os, sys
import numpy as np
import csv
import pandas as pd

# Remove the Soil Mukeys from the MukeyList that are not present in the XML file
#states = ['wi','mi','oh','sd','mn','mo','il','in','ia']
states = ['ia']
out = open("Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/MissingSoils/All_MissingSoilsArea.csv",'a')
out.write("State,Percent_Missing_rea\n")
Midtotal = 0.0
Mid_area = 0.0
for st in states:
	areafile = "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/MissingSoils/"+st+"_Mukey_Area.csv"
	soilfile = "Z:/Users/rilllydi/MidwestSALUS/SALUSresults/CornOnly/MissingSoils/"+st+"_missingSoils.csv"
	areas = pd.read_csv(areafile)
	areas = areas.dropna()
	soils = pd.read_csv(soilfile)
	soils = soils.dropna()
	st_area = areas["MU_Area"].sum()
	Mid_area = Mid_area + st_area
	
	#join the data tables and get the sum of the area of the MUKEY
	soilarea = pd.merge(soils, areas, left_on='MUKEY', right_on='MUKEYNUM')
	
	#soilarea["perarea"] = (soilarea["MU_Area"] / st_area) * 100
	#soilarea.to_csv(out2, sep=',', header=False)
	
	total = soilarea["MU_Area"].sum()
	totalper = float(total) / float(st_area)
	Midtotal = Midtotal + total
	out.write(st + "," + '{:.2%}'.format(totalper) + "\n")
Midtotalper = float(Midtotal) / Mid_area
out.write('Midwest,{:.2%}'.format(Midtotalper) + "\n")