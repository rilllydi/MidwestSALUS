# This script reads the output from the R script (MUKEY, Wx X coord, Wx Y coord)
# and then sorts by weather grid and writes the weather grid ( 1 per line) and its associated Mukeys as columns to a csv file

#/cygdrive/c/Python27/ArcGIS10.2/python

# Import system modules
import os, sys
import numpy as np
import csv
import pandas as pd

#############################
# This function finds the nearest wx grid point to given x,y
#states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
states = ['wi','mi','oh','sd','mn','mo']
#states = ['il']
for st in states:
    
    infile = "Z:/Users/rilllydi/MidwestSALUS/Wx_for_Soils/" + st + "_wx_soil_unique.csv"
    out_final = infile.replace(".csv","_SALUS.csv")
    df = pd.read_csv(infile) #,index_col=False)
    df = df.dropna()
    df["MUKEY"] = st.upper() + df["MUKEY"].astype(str)
    # Need to change MUKEY to STMUKEY (st.upper())
    
    df2 = (df.groupby(["WX_X","WX_Y"]).apply(lambda x: x["MUKEY"].tolist()).reset_index())
    dval = pd.DataFrame(df2[0].tolist(), )
    df2 = df2.drop(0, axis=1)
    data_final = pd.concat([df2, dval], axis=1)
    data_final.to_csv(out_final, sep=',', header=False)