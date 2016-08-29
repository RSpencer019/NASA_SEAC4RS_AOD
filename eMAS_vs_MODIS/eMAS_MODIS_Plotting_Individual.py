#------------------------------------------------------------------------------
# Name:        eMAS (SEAC4RS) vs. MODIS - Individual Module Plotting
# Description: NASA Aerosol Project
#
# Author:      Robert S. Spencer
#
# Created:     07/01/2016
# Python:      3.5
#------------------------------------------------------------------------------

plot_buffer = 1.0


import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_MODIS/eMAS_vs_MODIS_Compiled.csv',header=0)
grouped = data['MODIS_file'].groupby(data['eMAS_file'])

for item in grouped:
    print '------'
    eMAS_file = item[0]
    MODIS_files = item[1]
    
    # Read dataset
    hdf = SD(eMAS_file, SDC.READ)
    DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float) * 0.001
    
    # Read geolocation dataset.
    eMAS_lat = hdf.select('Latitude')
    eMAS_lon = hdf.select('Longitude')
    eMAS_lat = eMAS_lat[:,:].astype(float) * 0.0001
    eMAS_lon = eMAS_lon[:,:].astype(float) * 0.0001
    eMAS_lat_max = eMAS_lat.max()
    eMAS_lat_min = eMAS_lat.min()
    eMAS_lon_max = eMAS_lon.max()
    eMAS_lon_min = eMAS_lon.min()
    
    # Define eMAS border outline            
    lon_corners = [eMAS_lon[0,0],eMAS_lon[-1,0],eMAS_lon[-1,-1],eMAS_lon[0,-1],eMAS_lon[0,0]]
    lat_corners = [eMAS_lat[0,0],eMAS_lat[-1,0],eMAS_lat[-1,-1],eMAS_lat[0,-1],eMAS_lat[0,0]]
            
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0] * 0.001
    
    #data[data == fv] = np.nan
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)
        
        
    
    fig = plt.figure()
    fig.clear()
    lowerval = 0.0
    upperval = 0.75
    
    # Draw an equidistant cylindrical projection using the low resolution
    # coastline database.
    m = Basemap(projection='cyl', resolution='l',
                llcrnrlat = eMAS_lat_min-plot_buffer, urcrnrlat = eMAS_lat_max+plot_buffer,   # Extent: eMAS
                llcrnrlon = eMAS_lon_min-plot_buffer, urcrnrlon = eMAS_lon_max+plot_buffer)
#                llcrnrlat = 20, urcrnrlat = 50,     # EXTENT: USA
#                llcrnrlon = -125, urcrnrlon = -70)
#                llcrnrlat = -90, urcrnrlat = 90,     # EXTENT: WORLD
#                llcrnrlon = -180, urcrnrlon = 180)
    m.bluemarble(alpha=0.75)
    m.drawcoastlines(linewidth=0.75)
    m.drawparallels(np.arange(-90., 90., 1), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 180., 1), labels=[0, 0, 0, 1])
    fig = plt.gcf()
    plt.title('eMAS vs MODIS (%s)'%(eMAS_file[-55:]))
    
    
    # Map out background MODIS files    
    for MODIS_file in MODIS_files:
        
        # Read dataset.
        hdf = SD(MODIS_file, SDC.READ)
        DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
        data3D = hdf.select(DATAFIELD_NAME)
        MODIS = data3D[1,:,:]       # 1: 550 nm
        
        # Read geolocation dataset.
        MODIS_lat = hdf.select('Latitude')
        MODIS_lon = hdf.select('Longitude')
        MODIS_lat = MODIS_lat[:,:]
        MODIS_lon = MODIS_lon[:,:]
        
        # Handle fill value.
        attrs = data3D.attributes(full=1)
        fillvalue=attrs["_FillValue"]
        
        # fillvalue[0] is the attribute value.
        fv = fillvalue[0]
        
        #data[data == fv] = np.nan
        MODIS = np.ma.masked_array(MODIS, MODIS == fv)
        
        # Plot MODIS
        m.pcolormesh(MODIS_lon, MODIS_lat, MODIS*0.001,alpha=1,vmin=lowerval,vmax=upperval)
        
    # Overlay eMAS
    res = 1
    m.pcolormesh(eMAS_lon[::res][::res], eMAS_lat[::res][::res], eMAS[::res][::res],vmin=lowerval,vmax=upperval)
    m.plot(lon_corners,lat_corners,'k-',linewidth=0.75)            

    # Save figure
    pngfile = "Plots/{0}.png".format(eMAS_file[-55:])
    fig.savefig(pngfile,dpi=200)
    plt.close()
        