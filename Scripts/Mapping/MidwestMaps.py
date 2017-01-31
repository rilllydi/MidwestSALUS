# This script is originally Ruben's and has been modified by Lydia
# This script should create maps from the Midwest rasters for each SALUS scenario and for each variable within the raster

import arcpy, numpy, os, shutil, time, sys
from arcpy import env
import arcpy.mapping as mapping

############################################
def writeLyr(rast, name):
	fullrast = os.path.join(rootdir,rast)
	if os.path.isfile(fullrast):
		arcpy.env.workspace = rootdir
		outlyr = rast.replace(".tif",".tif.lyr")
		outlyr = os.path.join(lyr_folder,outlyr)
		lyrname = SC + "_" + name
		
		# Execute MakeFeatureLayer (save raster to lyr symbology file)
		arcpy.MakeRasterLayer_management(rast, lyrname, rootdir) 
		# do i want a layer per variable?? Or should I loop through the variables when making the maps? (can use where clause here)
		
		# can i adjust the layer symbology now?
		
		
		# Execute SaveToLayerFile (write the lyr file to the computer )
		arcpy.SaveToLayerFile_management(lyrname, outlyr, "ABSOLUTE")
	else:
		print "This raster does not exist:", fullrast

############################################
def makeMaps(rast,name):
	fullrast = os.path.join(rootdir,rast)
	if os.path.isfile(fullrast):


		#cropAccList = [("alf", "Alfalfa"), ("dry", "Beans"), ("cot", "Cotton"), ("cor", "Corn"),
		#				   ("pea", "Peas"), ("sil", "Corn Silage"), ("soy", "Soybeans"), ("wht", "Wheat")]
		#yieldStabList = [("level", "Yield Level"), ("_stab", "Yield Stability"), ("t_var", "Temporal Variance")]

		mxd = mapping.MapDocument(MXD) # get empty arcgis file
		df = mapping.ListDataFrames(mxd, "*")[0]

		for elm in mapping.ListLayoutElements(mxd, "TEXT_ELEMENT"): # loop through the elements in the mxd and set them
			print "element: ", elm.name
			if elm.name == "Title":
				elm.text = "SALUS simulated VARIABLE, " + SC
			#if elm.name == "Unstable":
			#	elm.elementPositionX = -5
			#if elm.name == "total_area":
			#	elm.elementPositionX = -5

		layer = mapping.Layer(fullrast) # add the raster layer to the mxd
		mapping.AddLayer(df, layer, "TOP")

		arcpy.RefreshTOC() # refresh
		arcpy.RefreshActiveView()

		# Load the symbology layer (from a file or create one from the raster on the fly)
		#lyrname = SC + "_" + name
		## Execute MakeFeatureLayer (save raster to lyr symbology file)
		##arcpy.MakeRasterLayer_management(rast, lyrname, rootdir)
		#outlyr = rast.replace(".tif",".tif.lyr")
		outlyr = "mi_SC3_endvalue.tif.lyr"
		outlyr = os.path.join(lyr_folder,outlyr)
		sym_lay = mapping.Layer(outlyr) # lyrname (using lines above) # or I can skip the line above and call the full path name (if lyr file is already made)
		
		lay_name = "Variable (units)"
		#mapping.AddLayer(df, layer, "TOP")
		updateLayer = mapping.ListLayers(mxd, "", df)[0]
		update_layer(df, updateLayer, sym_lay, True)
		apply_sym(updateLayer, sym_lay) # added symbology_only = True so that the layer becomes classified (symbology)??
		updateLayer.name = lay_name
				
		print "Added layer."
		
		#style_stab = mapping.ListStyleItems("USER_STYLE", "Legend Items")[0]
		
		legend = mapping.ListLayoutElements(mxd, "LEGEND_ELEMENT")[0] # updates all legend items to the chosen style
		for lyr in legend.listLegendItemLayers(): #[:3]:
		#	legend.updateItem(lyr, style_stab)
			print "in legend: ", lyr.name
			if "States" in lyr.name or "Counties" in lyr.name:
				print "states is in legend"
				legend.removeItem(lyr)
		
		##############
		endvarlist = ["avgGWAD","avgCWAD"]
		for var in endvarlist:
			name = var + "_" + name
			
			for lyr in mapping.ListLayers(mxd, "", df)[:-1]:  	# -1 to avoid boundary layer              						
				lyr.visible = True
				print "layer name: ", lyr.name
				
				# change the symbology value field
				if lyr.symbologyType == "RASTER_CLASSIFIED": #, then create a variable reference to the RasterClassifiedSymbology class for that layer (lyrSymbolClass = lyr.symbology)
					print "RASTER IS CLASSIFIED"
					lyr.symbology.valueField = var
					# first value is the minimum
					#lyr.symbology.classBreakValues = [1, 60, 118, 165, 255]
					# one less label than the classBreakValues
					#lyr.symbology.classBreakLabels = ["1 to 60", "61 to 118", 
					 #                   "119 to 165", "166 to 255"]
					#lyr.symbology.classBreakDescriptions = ["Class A", "Class B",
					 #                         "Class C", "Class D"]
					# lyr.symbology.excludedValues = '0'
				
				
				#name_temp = lyr.name
				#name =  "{}.jpg".format(name_temp.replace(" ", "_"))
				arcpy.RefreshTOC()                                          # refresh the TOC and active view
				arcpy.RefreshActiveView()
				outjpeg = os.path.join(outdir, name)
				mapping.ExportToJPEG(mxd, outjpeg, "PAGE_LAYOUT", resolution = 300, jpeg_quality = 100) #export to jpeg
				lyr.visible = False                                         # switch off layer
				arcpy.RefreshTOC()                                          # refresh again
				arcpy.RefreshActiveView()
			print "Created jpeg's of layers."
			
			# Save a copy of the mxd with the layer in it to the MXDs directory
			new_mxd_name = name + ".mxd"
			new_mxd = os.path.join(newMXDdir, new_mxd_name)
			arcpy.RefreshTOC()                  # refresh again
			arcpy.RefreshActiveView()
			mxd.saveACopy(new_mxd)
		##################
		
		# Save the mapp and its data to a single compressed .mpkx file
		#MPK_name = name + ".mpk"
		#MPK_name = os.path.join(MPKdir, MPK_name)
		#print new_mxd
		#print MPK_name
		#arcpy.PackageMap_management(new_mxd, MPK_name, "PRESERVE", "CONVERT_ARCSDE", "#", "ALL") # crashes the script :(

		del mxd
		del df
		print "Done with ", rast
		print
		arcpy.Delete_management("in_memory")
	
	else:
		print "This raster does not exist:", fullrast
###########################################################
rootdir = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\Midwest_Rasters" # location of rasters
lyr_folder = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\DataforMapping\LyrFiles" # where your symbology (lyr) files are/will be located
MXD = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\DataforMapping\CRMapper.mxd" # The MXD (empty map formatted to what you want)
outdir = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\Maps" # output directory for the maps
newMXDdir = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\Maps\MXDs" # output directory for the maps
MPKdir = r"Z:\Users\rilllydi\MidwestSALUS\SALUSresults\Maps\MPKs" # output directory for the maps

env.workspace = rootdir
env.overwriteOutput = True	

#sym_folder = r"C:\Users\ulbrichr\Desktop\python_data\02_arcgis\yield_symbology" # where your symbology (lyr) files are located
#MXD = r"C:\Users\ulbrichr\Desktop\basemap_only_ppt_promissing_3.mxd" # THIS IS IMPORTANT. look at benajamins: CRMapper.mxd

# shorthand for commands
apply_sym = arcpy.ApplySymbologyFromLayer_management
update_layer = mapping.UpdateLayer

start_time = time.time()
print "Start script..."
print

# For each SALUS scenario load the Midwest rasters (end of year and GWAD yearly)
#scenarios = ["SC1","SC2","SC3","SC4","SC5","SC6","SC7"]
scenarios = ["SC1"]
for SC in scenarios:
	# get the end of year raster:
	infile1 = SC + "_endof2016.TIF"
	infile2 = SC + "_yearly_GWAD.TIF"
	
	################ TESTING
	infile1 = "test11_1ras.tif"
	################
	
	#writeLyr(infile1,"endof2016")
	#writeLyr(infile2,"GWAD")

	makeMaps(infile1,"endof2016")

m, s = divmod((time.time() - start_time), 60)
h, m = divmod(m, 60)
print "Done! Time needed: %02d:%02d:%02d" % (h, m, s)
print

