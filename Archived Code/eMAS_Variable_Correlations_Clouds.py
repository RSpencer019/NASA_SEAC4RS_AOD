#------------------------------------------------------------------------------
# Name:        eMAS (SEAC4RS) Variable Correlations
# Description: Generates statistics on each granule variable in SEAC4RS
#               campaign to be plotted for correlations
#
# Author:      Robert S. Spencer
#
# Created:     6/15/2016
# Python:      2.7
#------------------------------------------------------------------------------

import os
import numpy as np
from pyhdf.SD import SD, SDC
import pandas as pd
import csv
from datetime import datetime
import time as tm
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

start_time = tm.time()
cwd = os.getcwd()

eMAS_file, AOD_mean, AOD_var, AOD_range, AOD_elevated, solar_zenith, cloud_percent, altitude = ([] for i in range(8))

def file_search(ext):        
    complete = False
    subdir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data'     ###### SPECIFY DIRECTORY!
    currentfolderdict = {'CWD':subdir}
    while complete is False:      
        newfolderdict = {}      
        for folder in currentfolderdict:
            os.chdir(currentfolderdict[folder])                        
            for item in os.listdir(currentfolderdict[folder]):
                if os.path.isdir(item) == True:
                    newfolderdict[item] = os.path.abspath(item)
                else:
                    if item[-len(ext):] == ext:
                        statistics(os.path.abspath(item))
                        cloudstats(os.path.abspath(item))
        if len(newfolderdict) == 0:
            complete = True
        else:
            currentfolderdict = newfolderdict
    os.chdir(cwd)



def cloudstats(eMAS_file):
    cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'
    eMAS_ID = eMAS_file[-45:-37]
    for cloud_file in os.listdir(cloud_dir):
        if cloud_file[-45:-37] == eMAS_ID:
            cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)

            ### Height ###
            DATAFIELD_NAME = 'Cloud_Top_Height'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Height = data3D[:,:].astype('float')
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Height = np.ma.masked_array(Height, Height == fv)
            scalefactor=attrs['scale_factor']
            sf = scalefactor[0]
            Height *= sf
            Height /= 1000 # to km

            Height_mean.append(np.mean(Height))
            




    
def statistics(FILE_NAME):
    print '--Running--'
    # Read HDF file
    hdf = SD(FILE_NAME, SDC.READ)
    
    # AOD
    dataset = hdf.select('Optical_Depth_Land_And_Ocean')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    scalefactor=attrs['Scale_factor']
    fv = fillvalue[0]
    sf = scalefactor[0]
    AOD = dataset[:,:].astype('float')
    AOD[AOD == fv] = np.nan
    AOD *= sf 
    
    # Solar Zenith Angle
    dataset = hdf.select('Solar_Zenith')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    scalefactor=attrs['Scale_factor']
    fv = fillvalue[0]
    sf = scalefactor[0]
    sol_zen = dataset[:,:].astype('float')
    sol_zen[sol_zen == fv] = np.nan
    sol_zen *= sf
    
    # Aerosol_Cldmask_Land_Ocean
    dataset = hdf.select('Aerosol_Cldmask_Land_Ocean')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    fv = fillvalue[0]
    cld_msk = dataset[:,:].astype('float')
    cld_msk[cld_msk == fv] = np.nan
    cld_perc = float(len(cld_msk[cld_msk == 0])) / float( len(cld_msk[cld_msk == 0]) + len(cld_msk[cld_msk == 1]))
    
    # Aircraft_Altitude
    dataset = hdf.select('Aircraft_Altitude')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    scalefactor=attrs['Scale_factor']
    fv = fillvalue[0]
    sf = scalefactor[0]
    alt = dataset[:].astype('float')
    alt[alt == fv] = np.nan
    alt *= sf
    
    
    ### APPEND ###
    eMAS_file.append(FILE_NAME)
    AOD_mean.append(np.nanmean(AOD))
    AOD_var.append(np.nanvar(AOD))
    AOD_range.append(np.nanmax(AOD)-np.nanmin(AOD))
    AOD_elevated.append(float(len(AOD[AOD>0.5])) / (len(AOD[AOD>-2]) + 1))
    solar_zenith.append(np.nanmean(sol_zen))
    cloud_percent.append(cld_perc)
    altitude.append(np.nanmean(alt))





    
file_search('hdf')


# Join the attributes together
compiled = [eMAS_file,AOD_mean,AOD_var,AOD_range,AOD_elevated,solar_zenith,cloud_percent,altitude]

# Generate csv labels
labels = 'eMAS_file,AOD_mean,AOD_var,AOD_range,AOD_elevated,solar_zenith,cloud_percent,altitude'

# Write to file
if os.path.exists(cwd+'/eMAS_Statistics.csv'):
    os.remove(cwd+'/eMAS_Statistics.csv')
with open(cwd+'/eMAS_Statistics.csv', "w") as output:  
    output.write(labels+'\n')  
    writer = csv.writer(output, lineterminator='\n')  
    writer.writerows(np.transpose(compiled))


print("--- %s seconds ---" % (tm.time() - start_time))



eMAS_Stats = pd.read_csv(cwd+'/eMAS_Statistics.csv')

fig = plt.figure()
plt.figtext(0.5,0.025,'Effects of Solar Zenith Angle',ha='center')
pl1 = fig.add_subplot(2,2,1)
pl1.scatter(eMAS_Stats.solar_zenith, eMAS_Stats.AOD_mean)
plt.title('AOD Mean')
pl2 = fig.add_subplot(2,2,2)
pl2.scatter(eMAS_Stats.solar_zenith, eMAS_Stats.AOD_var)
plt.title('AOD Variance')
pl3 = fig.add_subplot(2,2,3)
pl3.scatter(eMAS_Stats.solar_zenith, eMAS_Stats.AOD_range)
plt.title('AOD Range')
pl4 = fig.add_subplot(2,2,4)
pl4.scatter(eMAS_Stats.solar_zenith, eMAS_Stats.AOD_elevated)
plt.title('AOD Elevated')

fig2 = plt.figure()
plt.figtext(0.5,0.025,'Effects of Clouds Present',ha='center')
pl1 = fig2.add_subplot(2,2,1)
pl1.scatter(eMAS_Stats.cloud_percent, eMAS_Stats.AOD_mean)
plt.title('AOD Mean')
pl2 = fig2.add_subplot(2,2,2)
pl2.scatter(eMAS_Stats.cloud_percent, eMAS_Stats.AOD_var)
plt.title('AOD Variance')
pl3 = fig2.add_subplot(2,2,3)
pl3.scatter(eMAS_Stats.cloud_percent, eMAS_Stats.AOD_range)
plt.title('AOD Range')
pl4 = fig2.add_subplot(2,2,4)
pl4.scatter(eMAS_Stats.cloud_percent, eMAS_Stats.AOD_elevated)
plt.title('AOD Elevated')

fig3 = plt.figure()
plt.figtext(0.5,0.025,'Effects of Altitude',ha='center')
pl1 = fig3.add_subplot(2,2,1)
pl1.scatter(eMAS_Stats.altitude, eMAS_Stats.AOD_mean)
plt.title('AOD Mean')
pl2 = fig3.add_subplot(2,2,2)
pl2.scatter(eMAS_Stats.altitude, eMAS_Stats.AOD_var)
plt.title('AOD Variance')
pl3 = fig3.add_subplot(2,2,3)
pl3.scatter(eMAS_Stats.altitude, eMAS_Stats.AOD_range)
plt.title('AOD Range')
pl4 = fig3.add_subplot(2,2,4)
pl4.scatter(eMAS_Stats.altitude, eMAS_Stats.AOD_elevated)
plt.title('AOD Elevated')