#------------------------------------------------------------------------------
# Name:        MODIS Location Search
# Description: Search MODIS directory for modis granules that fall within
#               a eMAS granule
#
# Author:      Robert S. Spencer
#
# Created:     6/20/2016
# Python:      2.7
#------------------------------------------------------------------------------


### INPUT: specify collocation criteria and data directories! ###
time_delta = 30 # min
eMAS_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data'
MODIS_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data'
ext = 'hdf'
#------------------------------------------------------------------------------


import os
from pyhdf.SD import SD, SDC
from datetime import datetime, timedelta

start_time = tm.time()
cwd = os.getcwd()


# Generates list of desired files from directory tree
def file_search(data_dir): 
    complete = False    
    currentfolderdict = {'CWD':data_dir}
    data_files = []    
    while complete is False:      
        newfolderdict = {}      
        for folder in currentfolderdict:
            os.chdir(currentfolderdict[folder])                        
            for item in os.listdir(currentfolderdict[folder]):
                if os.path.isdir(item) == True:
                    newfolderdict[item] = os.path.abspath(item)
                elif item[-len(ext):] == ext:
                    data_files.append(os.path.abspath(item))
        if len(newfolderdict) == 0:
            complete = True
        else:
            currentfolderdict = newfolderdict
    os.chdir(cwd)
    return data_files


# Main collocation function
def collocation():
    
    # Extract data from eMAS files
    for eMAS_item in eMAS_files:

        # Extract timestamp from filename
        eMAS_year = int(eMAS_item[-36:-32])
        eMAS_month = int(eMAS_item[-32:-30])
        eMAS_day = int(eMAS_item[-30:-28])
        eMAS_hour = int(eMAS_item[-27:-25])
        eMAS_min = int(eMAS_item[-25:-23])
        eMAS_date = datetime(eMAS_year,eMAS_month,eMAS_day,eMAS_hour,eMAS_min)        
        
        # Temporal collocation date range
        date_min = eMAS_date - timedelta(minutes=time_delta)
        date_max = eMAS_date + timedelta(minutes=time_delta)

        # Read geolocation dataset.
        hdf = SD(eMAS_item, SDC.READ)
        
        # Extract spatial information from dataset
        eMAS_lat = hdf.select('Latitude')
        eMAS_lon = hdf.select('Longitude')
        eMAS_lat = eMAS_lat[:,:].astype(float) * 0.0001
        eMAS_lon = eMAS_lon[:,:].astype(float) * 0.0001
        eMAS_lat_max = eMAS_lat.max()
        eMAS_lat_min = eMAS_lat.min()
        eMAS_lon_max = eMAS_lon.max()
        eMAS_lon_min = eMAS_lon.min()
                
        # Extract data from MODIS files
        for MODIS_item in MODIS_files:

            # Extract timestamp from filename
            MODIS_year = int(MODIS_item[-34:-30])
            MODIS_jday = int(MODIS_item[-30:-27])
            MODIS_hour = int(MODIS_item[-26:-24])
            MODIS_min = int(MODIS_item[-24:-22])
            MODIS_date = datetime(MODIS_year,1,1,MODIS_hour,MODIS_min) + timedelta(MODIS_jday-1)    
        
            if (MODIS_date > date_min) & (MODIS_date < date_max):
                
                hdf = SD(MODIS_item, SDC.READ)
                # Extract spatial information from dataset
                MODIS_lat = hdf.select('Latitude')
                MODIS_lon = hdf.select('Longitude')
                MODIS_lat = MODIS_lat[:,:]
                MODIS_lon = MODIS_lon[:,:]
                MODIS_lat_max = MODIS_lat.max()
                MODIS_lat_min = MODIS_lat.min()
                MODIS_lon_max = MODIS_lon.max()
                MODIS_lon_min = MODIS_lon.min() 
                                
                # Handles the PROBLEM with granules crossing the 180th meridian
                if (MODIS_lat_max < 85) & (MODIS_lat_min > -85) & (MODIS_lon_max < 175) & (MODIS_lon_min > -175):
                    # Spatial Collocation
                    if (MODIS_lat_max > eMAS_lat_min) & (MODIS_lat_min < eMAS_lat_max) & (MODIS_lon_max > eMAS_lon_min) & (MODIS_lon_min < eMAS_lon_max):
                        
                        ####### APPEND ######                 
                        compiled.append([str(eMAS_month)+'/'+str(eMAS_day), eMAS_item, MODIS_item])
                        print ('---MATCH---')
                        
    
compiled = []    

# Initiate functions
eMAS_files = file_search(eMAS_dir)
MODIS_files = file_search(MODIS_dir)
collocation()

# Output
labels = 'day,eMAS_file,MODIS_file'
if os.path.exists('eMAS_vs_MODIS_Compiled.csv'):
    os.remove('eMAS_vs_MODIS_Compiled.csv')
with open('eMAS_vs_MODIS_Compiled.csv', "w") as output:  
    output.write(labels+'\n')  
    writer = csv.writer(output, lineterminator='\n')  
    writer.writerows(compiled)


print("--- %s seconds ---" % (tm.time() - start_time))




''' SCRATCH CODE

if (MODIS_lat_max != -999.0) & (MODIS_lat_min != -999.0) & (MODIS_lon_max != -999.0) & (MODIS_lon_min != -999.0):

'''