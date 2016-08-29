#------------------------------------------------------------------------------
# Name:        MODIS Location Search
# Description: Search MODIS directory for modis granules that fall within
#               a specified spatial range (currently set for continental US)
#
# Author:      Robert S. Spencer
#
# Created:     6/20/2016
# Python:      2.7
#------------------------------------------------------------------------------

import os
from pyhdf.SD import SD, SDC

# INPUT: Specify spatial area!
latmin = 20
latmax = 50
lonmin = -130
lonmax = -70

# INPUT: Specify search directory!
directory = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/'            
            
for item in os.listdir(directory):
    FILE_NAME = directory + item
    hdf = SD(FILE_NAME, SDC.READ)
    lat = hdf.select('Latitude')
    latitude = lat[:]
    lon = hdf.select('Longitude')
    longitude = lon[:]
    latmax_M = latitude.max()
    latmin_M = latitude.min()
    lonmax_M = longitude.max()
    lonmin_M = longitude.min()
    
    if (latmax_M < 85) & (latmin_M > -85) & (lonmax_M < 175) & (lonmin_M > -175):
        
        if (latmax_M > latmin) & (latmin_M < latmax) & (lonmax_M > lonmin) & (lonmin_M < lonmax):
            
            print '\'' + directory + item + '\','

    
