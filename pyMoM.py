import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import shutil
import subprocess

def extractAndTidyCurated(experimentDate,experimentFolderPath,preprocFolder,strainToPosDic,rx0,rxPos):

    posToStrainDic={'Pos'+str(e):key for key in strainToPosDic.keys() for e in strainToPosDic[key]}

    for root, dirs, files in os.walk("%s"%(experimentFolderPath+"/"+preprocFolder)):
        for fileName in files:
            #print(fileName)
            if re.search(rx0,fileName): # if a relevant CSV file is found
                outputFolder=re.search(rx0,fileName).group(1) #name of the folder
                fullPath=os.path.join(root,fileName) #full path to the csv file, should NOT contain `curated`
                outputFolderPath=fullPath.split("/")
                outputFolderPath='/'.join(outputFolderPath[:-1])
                #if re.search(rxOrientation,fullPath): #in the full path there should be the orientation ./date_bottom|up
                #orientation=re.search(rxOrientation,fullPath).group(1) #save orientation
                if re.search(rxPos,outputFolder): # if the file was found in a desired folder
                    pos=re.search(rxPos,outputFolder).group(1)
                    if pos in list(posToStrainDic.keys()): # just another way to control for the list of positions that are considered
                        strain=posToStrainDic[re.search(rxPos,outputFolder).group(1)] #saves strain based on the position
                        outputFolder=outputFolder[:-4]
                        outputPath="%s/curated_data/%s_%s_curated/%s_curated"%(experimentFolderPath,experimentDate,strain,outputFolder) #date_bottom|top_strain_curated
                    #print(outputFolder)
                        if "curated" not in fullPath: #avoid copying already curated files
                            try: 
                                shutil.copytree(outputFolderPath,outputPath)
                                print("Copying %s to %s."%(outputFolderPath,outputPath))
                            except:
                                print("Error for %s: cannot copy curated folder. Curated folder certainly already exists."%(outputFolderPath))
                        
def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]
    
    
                   
  
