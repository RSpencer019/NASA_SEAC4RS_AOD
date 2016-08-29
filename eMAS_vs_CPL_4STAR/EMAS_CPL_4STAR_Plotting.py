#------------------------------------------------------------------------------
# Name:        HDF4 Data Extraction
# Description: NASA Aerosol Project
#
# Author:      Robert S. Spencer
#
# Created:     6/13/2016
# Python:      3.5
#------------------------------------------------------------------------------


import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
import h5py


fig.clear()

lowerval = 0.0
upperval = 0.75

# Draw an equidistant cylindrical projection using the low resolution
# coastline database.
m = Basemap(projection='cyl', resolution='l',
#            llcrnrlat=35.5, urcrnrlat = 38,
#            llcrnrlon=-87, urcrnrlon = -85)
            llcrnrlat=20, urcrnrlat = 55,     # EXTENT: USA
            llcrnrlon=-125, urcrnrlon = -65)
m.bluemarble(alpha=0.75)
m.drawcoastlines(linewidth=0.5)
m.drawstates(linewidth=0.75)
m.drawparallels(np.arange(-90., 90., 10.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 180., 10.), labels=[0, 0, 0, 1])
fig = plt.gcf()
#plt.title('{0}\n {1} at H20PrsLvls=11'.format(FILE_NAME, DATAFIELD_NAME))




#### MODIS ####
'''
""" Retrieve list from MODIS_Search.py """
FILE_NAMES_M = [
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.1710.006.2014115095519.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.1715.006.2014115100401.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.1845.006.2014115101029.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.1850.006.2014115101234.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.1855.006.2014115100008.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2020.006.2014115095848.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2025.006.2014115100414.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2030.006.2014115100329.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2200.006.2014115095831.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2205.006.2014115101501.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/MODIS_Data/Aqua/228/MYD04_L2.A2013228.2210.006.2014115095806.hdf'
]


for item in FILE_NAMES_M:    

    # Open MODIS hdf file.
    FILE_NAME = item
    hdf = SD(FILE_NAME, SDC.READ)
    
    # Read dataset.
    DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    MODIS = data3D[:,:]       # 1: 550 nm
    
    # Read geolocation dataset.
    lat = hdf.select('Latitude')
    latitude = lat[:]
    lon = hdf.select('Longitude')
    longitude = lon[:]
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["_FillValue"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0]
    
    #data[data == fv] = np.nan
    MODIS = np.ma.masked_array(MODIS, MODIS == fv)
    
    # Plot MODIS
    m.pcolormesh(longitude, latitude, MODIS*0.001,alpha=1,vmin=lowerval,vmax=upperval)
'''



#### eMAS ####

FILE_NAMES_eMAS = []     ## Retrieve flight
eMAS_directory = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Sept/Emas_18Sept/'
for item in os.listdir(eMAS_directory):
    FILE_NAMES_eMAS.append(eMAS_directory+item)
    
''' ALTERNATIVE
""" Retrieve list from Compiled_Cleaned.txt """
FILE_NAMES = [
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_03_20130830_1759_1827_20160627_1603.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_07_20130830_1917_1930_20160627_1604.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_09_20130830_1939_1952_20160627_1605.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_10_20130830_2000_2009_20160627_1605.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_11_20130830_2016_2023_20160627_1606.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_19_20130830_2209_2223_20160627_1608.hdf',
'/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_21_20130830_2236_2250_20160627_1608.hdf'
]'''

for item in FILE_NAMES_eMAS:    

    # Open eMAS hdf file.
    FILE_NAME = item
    hdf = SD(FILE_NAME, SDC.READ)
  
    # Read dataset.
    DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float) * 0.001
    
    # Read geolocation dataset.
    lat = hdf.select('Latitude')
    latitude = lat[:,:].astype(float) * 0.0001
    lon = hdf.select('Longitude')
    longitude = lon[:,:].astype(float) * 0.0001
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0] * 0.001
    
    #data[data == fv] = np.nan
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)
    
    # Plot eMAS
    res = 1
    m.pcolormesh(longitude[::res][::res], latitude[::res][::res], eMAS[::res][::res],vmin=lowerval,vmax=upperval)
    

    

'''
#### Cloud CPL ####

# Open CPL file
FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/CPL_ATB_13954_16aug13.hdf5'
hdf = h5py.File(FILE_NAME,'r')

# Plot Clouds
x = hdf['Longitude'][hdf['Layer_Type'][:,0]==3]
y = hdf['Latitude'][hdf['Layer_Type'][:,0]==3]
m.scatter(x,y,c='white',lw=0, s=30)
'''


#### 4STAR ####
# Read 4STAR data
fourstar = pd.DataFrame.from_csv('/Users/rsspenc3/Desktop/SEAC4RS/DATA/4STAR_Data/SEAC4RS-4STAR-AOD-CWV_DC8_20130918_R2.ict',header=78)

# Plot 4STAR Data
x = fourstar['Longitude'][(fourstar['   AOD0550']!=-9999.0)]
y = fourstar['Latitude'][(fourstar['   AOD0550']!=-9999.0)]
z = fourstar['   AOD0550'][(fourstar['   AOD0550']!=-9999.0)]
c1 = m.scatter(x,y,c=z,lw=0, s=20, vmin=lowerval,vmax=upperval)


fig.colorbar(c1,ticks=[0,0.25,0.5,0.75])

"""
# Save plot.
pngfile = "{0}.py.png".format(FILE_NAME)
fig.savefig(pngfile,dpi=200)
"""



