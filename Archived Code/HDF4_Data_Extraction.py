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
from mpl_toolkits.basemap import Basemap
import numpy as np
from pyhdf.SD import SD, SDC

# Open file.
FILE_NAME = 'eMAS_vs_Aeronet/eMAS_Data/eMASL2CLD_13948_02_20130801_2212_2225_20150825_1343.hdf'
hdf = SD(FILE_NAME, SDC.READ)

# List available SDS datasets.
#print (hdf.datasets())


# Read dataset.
DATAFIELD_NAME='Surface_Temperature'
data3D = hdf.select(DATAFIELD_NAME)
data = data3D[:,:]

# Read geolocation dataset.
lat = hdf.select('PixelLatitude')
latitude = lat[:,:]
lon = hdf.select('PixelLongitude')
longitude = lon[:,:]

lat_min = 0.99 * latitude.min()
lat_max = 1.01 * latitude.max()
lon_min = 1.01 * longitude.min()
lon_max = 0.99 * longitude.max()

# Handle fill value.
attrs = data3D.attributes(full=1)
fillvalue=attrs["_FillValue"]

# fillvalue[0] is the attribute value.
fv = fillvalue[0]
#data[data == fv] = np.nan
data = np.ma.masked_array(data, data == fv)

# Draw an equidistant cylindrical projection using the low resolution
# coastline database.
m = Basemap(projection='cyl', resolution='l',
            llcrnrlat=lat_min, urcrnrlat = lat_max,
            llcrnrlon=lon_min, urcrnrlon = lon_max)
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(-90., 90., 0.5), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 180., 0.5), labels=[0, 0, 0, 1])
m.pcolormesh(longitude, latitude, data)
cb = m.colorbar()
cb.set_label('Unit:%')

plt.title('{0}\n {1} at H20PrsLvls=11'.format(FILE_NAME, DATAFIELD_NAME))
fig = plt.gcf()
# Show the plot window.
# plt.show()
"""
# Save plot.
pngfile = "{0}.py.png".format(FILE_NAME)
fig.savefig(pngfile,dpi=200)
"""