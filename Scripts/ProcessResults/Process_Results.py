# This script reads the seasonal SALUS output and writes out the yearly and overall SALUS results to csv files
# This script uses pandas library
# Needs to be updated depending on the variables included in the runs!!! (update to varlist only)

import sys
import os
import pandas as pd
import numpy as np

################################################################################
#scrdir = "/mnt/scratch/rilllydi/Midwest/"
#rootdir = "/mnt/home/rilllydi/Midwest/Results/"
rootdir = "C:/Users/Lydia Rill/Documents/MidwestSALUStest/"
scrdir = "C:/Users/Lydia Rill/Documents/MidwestSALUStest/"

varlist = ["GWAD","CWAD","NLCC"]

#states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
states = ['ia']

# Split the title into wx_x, wx_y, and MUKEY (and then drop Title)
def splitTitle(df):
	splits = df["Title"].str.split("|")
	df['MUKEY'] = splits.str[1]
	df["MUKEY"] = df["MUKEY"].str.replace(st.upper(),"")
	
	wx = splits.str[2].str.split("N_")
	df['wxID_y'] = wx.str[0]
	df["wxID_x"] = wx.str[1]
	df["wxID_x"] = df["wxID_x"].str.replace("W","")
	df["wxID_x"] = "-" + df["wxID_x"]
	df = df.drop(['Title'], axis=1)
	return(df)

# Rename the data frame columns so when you merge them you know which is which			
def renamecol(df,new):
	for var in varlist:
		newvar = new + var
		df = df.rename(columns={var:newvar})
	return(df)

# loop through all the result files (9 states and 12 chunks)	
for st in states:
			SC = 1
			chunk = 2
	#for SC in range(1,8): 
		#for chunk in range(0,11):
			
			result = scrdir + st + "_" + str(chunk) + "_SC" + str(SC) + "_seasonal.csv"
			outavg = rootdir + "Summary_1979_2016/" + st + "_" + str(chunk) + "_SC" + str(SC) + "_summary.csv" 

			# Use pandas to read in the SALUS results
			df = pd.read_csv(result, index_col=None)
			
			# Cols: ExpID,RID,RcID,Year,DOY,Title,CWAD,GWAD,NLCC, etc...
			# NEED TO DOUBLE CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			if SC <= 5:
				df = df.loc[df.RcID == 1]
			elif SC == 6: # corn, soybean, and cover crop
				#dfs = df.loc[df.RcID == 3]
				df = df.loc[df.RcID == 1]
			elif SC == 7: # corn, soybean, winter wheat, and cover crop
				#dfs = df.loc[df.RcID == 3]
				#dfw = df.loc[df.RcID == 4]
				df = df.loc[df.RcID == 1]

			# Drop ExpID,RID,RcID,Year,DOY
			df2 = df.drop(['ExpID','RID','RcID','DOY'], axis=1)
			df = df.drop(['ExpID','RID','RcID','Year','DOY'], axis=1)
			
			# calculate the avg, min, and max for each experiment across the years
			meandf = df.groupby(['Title'], as_index=False).mean()
			meandf = renamecol(meandf, "avg") # call the renamecol function
			mindf = df.groupby(['Title'], as_index=False).min()
			mindf = renamecol(mindf, "min") # call the renamecol function
			maxdf = df.groupby(['Title'], as_index=False).max()
			maxdf = renamecol(maxdf, "max") # call the renamecol function
			mdf = pd.merge(meandf, mindf, on="Title")
			mdf = pd.merge(mdf, maxdf, on="Title")
			
			mdf = splitTitle(mdf) # call the splitTitle function									
			mdf.to_csv(outavg, index=False, float_format="%.2f")
			
			# For each variable write out the data with each year as a column
			for var in varlist:
				outyr = rootdir + "Yearly_1979_2016/" + st + "_" + str(chunk) + "_SC" + str(SC) + "_yearly_" + var + ".csv"  
				df3 = pd.pivot_table(df2, index="Title", columns="Year", values=var)
				df3 = df3.reset_index() # turn the index into a column #df.reset_index(level=0, inplace=True)
				df3 = splitTitle(df3) # call the splitTitle function
				df3.to_csv(outyr, index=False, float_format="%.2f")