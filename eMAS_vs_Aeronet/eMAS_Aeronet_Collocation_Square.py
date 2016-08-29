#------------------------------------------------------------------------------
# Name:        eMAS vs Aeronet Collocation (SEAC4RS)
# Description: Generates collocation compiled.csv file for eMAS and Aeronet
#               observations
#
# Author:      Robert S. Spencer
#
# Created:     6/15/2016
# Python:      2.7
#------------------------------------------------------------------------------

"""
Generates: Compiled.csv, Aeronet_Sites.txt
"""

import os
import numpy as np
from pyhdf.SD import SD, SDC
import pandas as pd
import csv
from datetime import datetime
import time as tm
from sklearn.metrics import r2_score


start_time = tm.time()
cwd = os.getcwd()

if os.path.exists(cwd+'/Aeronet_Sites.txt'):
    os.remove(cwd+'/Aeronet_Sites.txt')


# INPUT: specify collocation criteria!
sample_radius = 3 # km
time_delta = pd.Timedelta('30 min')

# initiate final attributes
timestamp, location, longitude, latitude, elevation = ([] for i in range(5))
HDFfile, longitude_index, latitude_index, cval_eMAS, nval_eMAS, meanval_eMAS = ([] for i in range(6))
lev20file, nval_aeronet, meanval_aeronet = ([] for i in range(3))
meanval_aeronet_550_intrp, meanval_aeronet_550_intrp_r2 = ([] for i in range(2))



def generate_aeronet_meta(item2): 
    csvfile = open(item2, 'rb')
    csvreader = csv.reader(csvfile, delimiter=',')
    n = 0
    t = []
    for i in csvreader:
        t.append(i)
        if n == 3:
            break
        n += 1
    sites = open(cwd+'/Aeronet_Sites.txt','a')
    if os.path.getsize(cwd+'/Aeronet_Sites.txt') == 0:
        sites.write('Site_Name,Longitude(decimal_degrees),Latitude(decimal_degrees),Elevation(meters)\n')
    sites.write(t[2][0][9:]+','+t[2][1][5:]+','+t[2][2][4:]+','+t[2][3][5:]+'\n')
    sites.close()
    


# Searches through specified files in directory
def file_search(ext, current_location):      # lev20 or hdf input
    complete = False
    if ext == 'hdf':
        subdir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data'     ###### SPECIFY DIRECTORY!
    else:
        subdir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/Aeronet_Data'  ###### SPECIFY DIRECTORY!
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
                        if ext == 'hdf':
                            collocation(item)
                        if ext == 'lev20':
                            #print os.path.abspath(item)
                            generate_aeronet_meta(os.path.abspath(item))
                    else:
                        if current_location in item:
                            print item, '  true  ', current_location
                            #temporal_collocation(os.path.abspath(item))

        if len(newfolderdict) == 0:
            complete = True
        else:
            currentfolderdict = newfolderdict
    os.chdir(cwd)



# Main collocation function
def collocation(FILE_NAME):
    
    ## Read Aeronet    
    # Read aeronet site list
    global aeronet_sites
    aeronet_sites = pd.read_csv(cwd+'/Aeronet_Sites.txt', index_col=0)
    
    ## Read eMAS 
    hdf = SD(FILE_NAME, SDC.READ)
    
    # Read geolocation dataset.
    lat = hdf.select('Latitude')
    lon = hdf.select('Longitude')
    global eMAS_lat, eMAS_lon
    eMAS_lat = lat[:,:].astype(float) * 0.0001
    eMAS_lon = lon[:,:].astype(float) * 0.0001
    global eMAS_lat_max, eMAS_lat_min, eMAS_lon_max, eMAS_lon_min
    eMAS_lat_max = eMAS_lat.max()
    eMAS_lat_min = eMAS_lat.min()
    eMAS_lon_max = eMAS_lon.max()
    eMAS_lon_min = eMAS_lon.min()
    
    # Read dataset.
    global eMAS_data
    DATAFIELD_NAME='Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS_data = data3D[:,:].astype(float) * 0.001

    # Handle fill value.
    global fv
    attrs = data3D.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    fv = fillvalue[0] * 0.001
    
    #data[data == fv] = np.nan
    eMAS_data = np.ma.masked_array(eMAS_data, eMAS_data == fv)    
    
    spatial_collocation(FILE_NAME)
      


def spatial_collocation(FILE_NAME):        
    
    for i in aeronet_sites.index:      
        if (
        (aeronet_sites['Latitude(decimal_degrees)'][i] >= eMAS_lat_min ) &
        (aeronet_sites['Latitude(decimal_degrees)'][i] <= eMAS_lat_max ) &
        (aeronet_sites['Longitude(decimal_degrees)'][i] >= eMAS_lon_min ) &
        (aeronet_sites['Longitude(decimal_degrees)'][i] <= eMAS_lon_max )
        ):
            aeronet_lat = aeronet_sites.ix[i][1]
            aeronet_lon = aeronet_sites.ix[i][0]
    
            distance = ( np.abs(eMAS_lon - aeronet_lon) + np.abs(eMAS_lat - aeronet_lat) )
            t = []
            for z in range(len(eMAS_lat[0])-1):
                t.append(eMAS_lat[0][z] - eMAS_lat[0][z+1])
            for z in range(len(eMAS_lat)-1):
                t.append(eMAS_lat[z][0] - eMAS_lat[z+1][0])
            for z in range(len(eMAS_lon[0])-1):
                t.append(eMAS_lon[0][z] - eMAS_lon[0][z+1])
            for z in range(len(eMAS_lon)-1):
                t.append(eMAS_lon[z][0] - eMAS_lon[z+1][0])
            global dist_mean
            dist_mean = np.mean(np.abs(t))
            dist_max = np.sqrt(2 * max(t) ** 2)
            dist_min = np.min(distance)            
            if dist_min <= dist_max:
            
                ind = np.where(distance == dist_min)
                global cval_ind_lon, cval_ind_lat                  
                cval_ind_lon = ind[0][0]
                cval_ind_lat = ind[1][0]
                
                """# APPEND #""" 
                cval_eMAS.append(eMAS_data[cval_ind_lon][cval_ind_lat])
                longitude_index.append(cval_ind_lon)
                latitude_index.append(cval_ind_lat)
                location.append(i)
                latitude.append(aeronet_lat)                      
                longitude.append(aeronet_lon)
                elevation.append(aeronet_sites.ix[i][2])
                HDFfile.append(os.path.abspath(FILE_NAME))                         
                eMAS_meanval()
                
                print FILE_NAME , "    " , i
                
                ###################### Initiates the Aeronet query
                working_directory = os.getcwd()
                global current_timestamp
                current_timestamp = pd.Timestamp(datetime(int(FILE_NAME[-36:-32]), int(FILE_NAME[-32:-30]), int(FILE_NAME[-30:-28]),int(FILE_NAME[-27:-25]),int(FILE_NAME[-25:-23])))
                timestamp.append(current_timestamp)                
                current_location = i+'.lev20'
                file_search('temporal', current_location)
                os.chdir(working_directory)
                ######################


# spatial averaging of eMAS
def eMAS_meanval():
    sample_radius_pixels = int(round(sample_radius / (dist_mean * 111.32)))
    samples = []
    for i in range(-sample_radius_pixels,sample_radius_pixels+1):
        for j in range(-sample_radius_pixels,sample_radius_pixels+1):
            if 0 <= (cval_ind_lon+i) & (cval_ind_lon+i) < len(eMAS_data):
                if 0 <= (cval_ind_lat+j) & (cval_ind_lat+j) < len(eMAS_data[0]):
                    if eMAS_data[cval_ind_lon+i][cval_ind_lat+j] != fv:
                        samples.append(eMAS_data[cval_ind_lon+i][cval_ind_lat+j])
    """# APPEND #""" 
    nval_eMAS.append(len(samples))         
    meanval_eMAS.append(np.mean(samples))                    


# compiles the aeronet data
def temporal_collocation(item):   
    aeronet_data = pd.read_csv(item,header=4)
    aeronet_data['Date(dd-mm-yy)'] = aeronet_data['Date(dd-mm-yy)'].str.replace(':','-')
    aeronet_data['Timestamp'] = pd.to_datetime(aeronet_data['Date(dd-mm-yy)']+' '+aeronet_data['Time(hh:mm:ss)'],dayfirst=True)   
    aeronet_samples = aeronet_data[(aeronet_data['Timestamp'] > current_timestamp - time_delta) & (aeronet_data['Timestamp'] < current_timestamp + time_delta)]
    aeronet_samples = aeronet_samples[aeronet_samples.columns[aeronet_samples.columns.str.startswith('AOT')]]  
    """# APPEND #""" 
    nval_aeronet.append(len(aeronet_samples))
    meanval_aeronet.append(np.mean(aeronet_samples))
    lev20file.append(item)
    meanval_aeronet_550_intrp.append(aeronet_550_interpolation(np.mean(aeronet_samples))[0])
    meanval_aeronet_550_intrp_r2.append(aeronet_550_interpolation(np.mean(aeronet_samples))[1])



def aeronet_550_interpolation(sample):

    aeronet_col = list(sample.index)
    aeronet_wl = []
    for item in aeronet_col:
        aeronet_wl.append(int(item[4:]))
    aeronet_aod = []    
    for item in aeronet_col:
        aeronet_aod.append(sample[item])
    if np.any(np.isfinite(aeronet_aod)):
        
        aeronet_aod = np.array(np.log10(aeronet_aod))
        aeronet_wl = np.array(np.log10(aeronet_wl))
        idx = np.isfinite(aeronet_aod)
        z = np.polyfit(aeronet_wl[idx],aeronet_aod[idx],3)
        p = np.poly1d(z)
        r2 = r2_score(aeronet_aod[idx], p(aeronet_wl[idx]))
        return 10**(p(np.log10(550))), r2
    else:
        return np.nan, np.nan


# Initiation of collocation functions
file_search('lev20', 'null')
file_search('hdf','null')

# Generation of compiled.csv file
meanval_aeronet_transposed = np.transpose(meanval_aeronet)

compiled = [timestamp, location, longitude, latitude, elevation, longitude_index, latitude_index, HDFfile, lev20file, cval_eMAS, nval_eMAS, meanval_eMAS, nval_aeronet, meanval_aeronet_550_intrp, meanval_aeronet_550_intrp_r2]
for i in meanval_aeronet_transposed:
    compiled.append(i)

labels = 'timestamp,location,longitude,latitude,elevation,longitude(index),latitude(index),eMAS_file,Aeronet_file,cval_eMAS,nval_eMAS,meanval_eMAS_550,nval_aeronet,meanval_aeronet_550_intrp,meanval_aeronet_550_intrp_r2'
for i in meanval_aeronet[0].index:  
    labels = labels + ',meanval_aeronet_' + i[4:]

if os.path.exists(cwd+'/Compiled.csv'):
    os.remove(cwd+'/Compiled.csv')
with open(cwd+'/Compiled.csv', "w") as output:  
    output.write(labels+'\n')  
    writer = csv.writer(output, lineterminator='\n')  
    writer.writerows(np.transpose(compiled))


print("--- %s seconds ---" % (tm.time() - start_time))

