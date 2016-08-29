#------------------------------------------------------------------------------
# Name:        eMAS (SEAC4RS) vs. MODIS & Aeronet Plotting
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


#fig.clear()

DATAFIELD_NAME='Optical_Depth_Land_And_Ocean'

lowerval = 0.0
upperval = 0.75

# Draw an equidistant cylindrical projection using the low resolution
# coastline database.
m = Basemap(projection='cyl', resolution='l',
#            llcrnrlat=35.5, urcrnrlat = 38,
#            llcrnrlon=-87, urcrnrlon = -85)
            llcrnrlat=26, urcrnrlat = 38.5,     # EXTENT: USA
            llcrnrlon=-99, urcrnrlon = -82)
m.bluemarble(alpha=0.6)
m.drawcoastlines(linewidth=0.5)
m.drawstates(linewidth=0.75)
m.drawparallels(np.arange(-90., 90., 4), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 180., 4), labels=[0, 0, 0, 1])
fig = plt.gcf()
#plt.title('SEAC4RS AOD ("{0}") vs. Aeronet'.format(DATAFIELD_NAME))



'''
#### MODIS ####

""" Retrieve list from MODIS_Search.py """
FILE_NAMES = [
'Aug_30/MOD04_L2.A2013242.1540.006.2015069175744.hdf',
'Aug_30/MOD04_L2.A2013242.1545.006.2015069180306.hdf',
'Aug_30/MOD04_L2.A2013242.1550.006.2015069175440.hdf',
'Aug_30/MOD04_L2.A2013242.1720.006.2015069180633.hdf',
'Aug_30/MOD04_L2.A2013242.1725.006.2015069180700.hdf',
'Aug_30/MOD04_L2.A2013242.1730.006.2015069175548.hdf',
'Aug_30/MOD04_L2.A2013242.1855.006.2015069180114.hdf',
'Aug_30/MOD04_L2.A2013242.1900.006.2015069175932.hdf',
'Aug_30/MOD04_L2.A2013242.1905.006.2015069175730.hdf',
'Aug_30/MOD04_L2.A2013242.2035.006.2015069175614.hdf',
'Aug_30/MOD04_L2.A2013242.2040.006.2015069175634.hdf',
'Aug_30/MOD04_L2.A2013242.2215.006.2015069175852.hdf',
'Aug_30/MOD04_L2.A2013242.2220.006.2015069175547.hdf',
'Aug_30/MOD04_L2.A2013242.2225.006.2015069180459.hdf',
'Aug_30/MOD04_L2.A2013242.2355.006.2015069180236.hdf'
]

for item in FILE_NAMES:    

    # Open MODIS hdf file.
    FILE_NAME = 'MODIS_Data/' + item
    hdf = SD(FILE_NAME, SDC.READ)
    
    # Read dataset.
    DATAFIELD_NAME='Corrected_Optical_Depth_Land'
    data3D = hdf.select(DATAFIELD_NAME)
    MODIS = data3D[1,:,:]       # 1: 550 nm
    
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

""" Retrieve list from Compiled_Cleaned.csv """

#FILE_NAMES = [
#'Aug/Emas_30Aug/eMASL2AER_13959_05_20130830_1854_1909_20160627_1604.hdf'
#]

data = pd.DataFrame.from_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
FILE_NAMES = data['eMAS_file']


for item in FILE_NAMES:    

    # Open eMAS hdf file.
#    FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/eMAS_Data/' + item
    FILE_NAME = item
    
    hdf = SD(FILE_NAME, SDC.READ)
  
    # Read dataset.
    DATAFIELD_NAME='Optical_Depth_Land_And_Ocean'
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

    


#### All Aeronet Sites ####

# Read Aeronet sites
aeronet = pd.DataFrame.from_csv('Aeronet_Sites_All.txt',header=1)
x = aeronet['Longitude(decimal_degrees)']
y = aeronet['Latitude(decimal_degrees)']
# Plot Aeronet Sites
plt.scatter(x,y,s=30,c='black')




#### Collocated Aeronet Sites ####

# Read Aeronet sites
aeronet = pd.DataFrame.from_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
x = aeronet['longitude']
y = aeronet['latitude']
n = aeronet['location']
t = aeronet.index
z1 = aeronet['meanval_eMAS_550']
z2 = aeronet['meanval_aeronet_550_intrp']

# Plot Aeronet Sites
plt.scatter(x,y,c=z1,s=600,vmin=lowerval,vmax=upperval)
plt.scatter(x,y,c=z2,s=200,vmin=lowerval,vmax=upperval)
'''
for i, txt in enumerate(n):
    plt.annotate(txt, (x[i]+.1,y[i]+.1))    
for i, txt in enumerate(t):
    plt.annotate(txt, (x[i]+.1,y[i]))
'''

m.colorbar(location='right',ticks=[0,0.25,0.5,0.75])





"""
# Save plot.
pngfile = "{0}.py.png".format(FILE_NAME)
fig.savefig(pngfile,dpi=200)



##### Aeronet Collocation LEIGH ####
#
## Read Aeronet sites
#aeronet = pd.DataFrame.from_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_AERONET/validation_with_aeronet/mapss_emas_seac4rs_alldata_2p5km_radius.csv',header=0)
#x = aeronet['Longitude']
#y = aeronet['Latitude']
#z1 = aeronet['mean_AOD0550']
#z2 = aeronet['mean_AOD0550intrp']
#
## Plot Aeronet Sites
#plt.scatter(x,y,c=z1,s=600,vmin=lowerval,vmax=upperval)
#plt.scatter(x,y,c=z2,s=200,vmin=lowerval,vmax=upperval)

"""



