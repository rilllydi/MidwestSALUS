# This script reads the seasonal SALUS output and writes out the yearly and overall SALUS results to csv files
# This script uses pandas library
# Needs to be updated depending on the variables included in the runs!!! (update to varlist only)
# One state at a time! No chunks in the output!

import sys
import os
import pandas as pd
import numpy as np

################################################################################
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

def main(st,SC,chunk,vartype,var):
	result = scrdir + st + "_" + str(chunk) + "_SC" + str(SC) + "_seasonal.csv"

	# Use pandas to read in the SALUS results
	df = pd.read_csv(result, index_col=None)
		
	# Cols: ExpID,RID,RcID,Year,DOY,Title,CWAD,GWAD,NLCC, etc...
	# NEED TO DOUBLE CHECK THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if SC <= 5:
		df = df.loc[df.RcID == 1]
	elif SC == 6: # corn, soybean, and cover crop
		#dfs = df.loc[df.RcID == 3] # soybean
		df = df.loc[df.RcID == 1]
	elif SC == 7: # corn, soybean, winter wheat, and cover crop
		#dfs = df.loc[df.RcID == 3] # soybean
		#dfw = df.loc[df.RcID == 4] # wheat
		df = df.loc[df.RcID == 1]

	if vartype == "summary":
		# Drop ExpID,RID,RcID,Year,DOY
		df2 = df.drop(['ExpID','RID','RcID','Year','DOY'], axis=1)
		
		# calculate the avg, min, and max for each experiment across the years
		meandf = df2.groupby(['Title'], as_index=False).mean()
		meandf = renamecol(meandf, "avg") # call the renamecol function
		mindf = df2.groupby(['Title'], as_index=False).min()
		mindf = renamecol(mindf, "min") # call the renamecol function
		maxdf = df2.groupby(['Title'], as_index=False).max()
		maxdf = renamecol(maxdf, "max") # call the renamecol function
		mdf = pd.merge(meandf, mindf, on="Title")
		mdf = pd.merge(mdf, maxdf, on="Title")
			
		mdf = splitTitle(mdf) # call the splitTitle function									
		return(mdf)
		#mdf.to_csv(outavg, index=False, float_format="%.2f")
			
	elif vartype == "yearly":
		# For each variable write out the data with each year as a column
		df3 = pd.pivot_table(df, index="Title", columns="Year", values=var)
		df3 = df3.reset_index() # turn the index into a column #df.reset_index(level=0, inplace=True)
		df3 = splitTitle(df3) # call the splitTitle function
		return(df3)
		#df3.to_csv(outyr, index=False, float_format="%.2f")

	elif vartype == "endof2016":
		# For each variable write out the data for the last day in the seasonal file (2016 at harvest or covercrop planting) to get the cumulative or final values
		# NEED TO CHECK FOR OTHER SCENARIOS BESIDES 1 - 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DONT WANT TO SUBSET RCID
		# find the last days
		#df4 = df.loc[:,(df.Year == 2016) & (df.DOY == 300)] # the last line in the data frame
		df4 = df[(df["Year"] == 2016) & (df["DOY"] == 300)] # the last line in the data frame
		df4 = df4.drop(['ExpID','RID','RcID','Year','DOY'], axis=1)
		#print(df4)
		df4 = splitTitle(df4) # call the splitTitle function
		return(df4)
		#df4.to_csv(outcumul, index=False, float_format="%.2f")

###################################################################################
scrdir = "/mnt/scratch/rilllydi/Midwest/"
rootdir = "/mnt/home/rilllydi/Midwest/Results/"

#allvarlist = ["GWAD","CWAD","NLCC","C_ResOrgBl","C_SloOrgBl","C_ActOrgBl","C_CO2","C_FertBl","NAID","NOAD","DRNC","PREC"] # NIAD
varlist = ["GWAD"] 

states = ['wi','mi','oh','in','il','ia','sd','mn','mo']
#states = ['mo']
#states = ['mi']


# loop through all the result files (9 states and 12 chunks)	
for SC in range(3,5): # 1 - 8
	dfyearlylist = []
	dfendlist = []
	dfsummlist = []
	for st in states:
		for chunk in range(0,11): # 0 - 10
			print st, ":", "SC", SC, "chunk:", chunk
			dfyearlylist.append(main(st,SC,chunk,"yearly","GWAD"))
			dfendlist.append(main(st,SC,chunk,"endof2016",""))
			dfsummlist.append(main(st,SC,chunk,"summary",""))

		# now join the data frames together for the whole state
		dfyearly = pd.concat(dfyearlylist) # join='outer', ignore_index=True)
		dfend = pd.concat(dfendlist)
		dfsumm = pd.concat(dfsummlist)

		dfyearly = dfyearly.drop_duplicates()
		dfend = dfend.drop_duplicates()
		dfsumm = dfsumm.drop_duplicates()

		var = "GWAD"
		outavg = rootdir + "Summary_1979_2016/" + st + "_SC" + str(SC) + "_summary.csv" 
		outyr = rootdir + "Yearly_1979_2016/" + st + "_SC" + str(SC) + "_yearly_" + var + ".csv"  
		outcumul = rootdir + "End_of_2016/" + st + "_SC" + str(SC) + "_endvalue.csv"  

		dfyearly.to_csv(outyr, index=False, float_format="%.2f")
		dfend.to_csv(outcumul, index=False, float_format="%.2f")
		dfsumm.to_csv(outavg, index=False, float_format="%.2f")
