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
cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'

import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC


data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_MODIS/eMAS_vs_MODIS_3KTerra_Compiled.csv',header=0)
grouped = data['MODIS_file'].groupby(data['eMAS_file'])

for item in grouped:
    print '------'
    eMAS_file = item[0]
    MODIS_files = item[1]
    
    ###--------------- eMAS AOD ------------------###

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
    cent = [eMAS.shape[0]/2,eMAS.shape[1]/2]         
    lon_corners = [eMAS_lon[0,0],eMAS_lon[cent[0],0],eMAS_lon[-1,0],eMAS_lon[-1,-1],eMAS_lon[cent[0],-1],eMAS_lon[0,-1],eMAS_lon[0,0]]
    lat_corners = [eMAS_lat[0,0],eMAS_lat[cent[0],0],eMAS_lat[-1,0],eMAS_lat[-1,-1],eMAS_lat[cent[0],-1],eMAS_lat[0,-1],eMAS_lat[0,0]]
            
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0] * 0.001
    
    #data[data == fv] = np.nan
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)
        
    
    
    ###--------------- RGB Projections ------------------###

        
    def scaled_refl(refl):
        result = 0. + 3.8393 * refl - 2.1626e-2 * refl**2 + 4.1278e-5 * refl**3
        return result
    
    MAXVALUES = 65535
    color_scale = np.zeros(MAXVALUES)
    band_scale = np.zeros(MAXVALUES)
    
    
    for k in range(0, MAXVALUES):
    
        val1 = k - int(MAXVALUES/2) + 0.5 -1

        if val1 <= 0. :
            band_scale[k] = 0.
        else :
            if val1 >= 11000. :  
                band_scale[k] = 255.
            else :
                band_scale[k] = 0.0 + (val1 - 0.) * (255. - 0.) / ( 11000. - 0.) 


        color_scale[k] = scaled_refl(k*255. / MAXVALUES*1.0)	
        if (color_scale[k] > 255.) :
            color_scale[k] = 255.
    
    for k in range(0, MAXVALUES):
    	
        intensity = band_scale[k] / 255.
        
        new_intensity = color_scale[int(intensity*(MAXVALUES-1)+0.5)] / 255.
    
        if intensity > 0. and intensity < 1. : 
            band_scale[k] = band_scale[k] * new_intensity / intensity
        if band_scale[k] > 255. :
            band_scale[k] = 255.
    
        		
    #---------- Read MODIS RGB Data ----------#
        
        
    for MODIS_file in MODIS_files:
        if 'Terra' in MODIS_file:
            fold1 = 'Terra/'
            fold2 = 'TerraGeo/'
        elif 'Aqua' in MODIS_file:
            fold1 = 'Aqua/'
            fold2 = 'AquaGeo/'
        timestamp = MODIS_file[-30:-22]
        
        
        for rgb_item in os.listdir('/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/'+fold1):
            if rgb_item[-30:-22] == MODIS_file[-30:-22]:
                rgb_file = rgb_item
        for rgb_item in os.listdir('/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/'+fold2):
            if rgb_item[-30:-22] == MODIS_file[-30:-22]:
                rgbgeo_file = rgb_item
                
        myd021km_name = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/' + fold1 + rgb_file
        myd03_name = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/' + fold2 + rgbgeo_file
                
        myd021km = SD(myd021km_name, SDC.READ)
        myd03 = SD(myd03_name, SDC.READ)
        
        ht = 2030 #1015
        wid = 1354 #677
        
        
        #---------- Reflectances (Bands 1, 4 and 3) from MODIS ----------#
        
        sds = myd021km.select('EV_250_Aggr1km_RefSB')
        Band_01 = sds.get(start=(0,0,0), count=(1, ht, wid), stride =(1,1,1)) *5.31009e-5
        
        sds = myd021km.select('EV_500_Aggr1km_RefSB')
        Band_04 = sds.get(start=(1,0,0), count=(1, ht, wid), stride =(1,1,1)) *3.42507e-5
        
        sds = myd021km.select('EV_500_Aggr1km_RefSB')
        Band_03 = sds.get(start=(0,0,0), count=(1, ht, wid), stride =(1,1,1)) *3.67092e-5
        
        
        sds = myd03.select('SolarZenith')
        miu0_read = sds.get(start=(0,0), count=(ht, wid), stride =(1,1)) *0.01
        miu0 = np.cos(miu0_read/180.*np.pi)
        
        sds = myd03.select('Latitude')
        latitude = sds.get(start=(0,0), count=(ht, wid), stride =(1,1))
        
        sds = myd03.select('Longitude')
        longitude = sds.get(start=(0,0), count=(ht, wid), stride =(1,1))
        
        
        #---------- RGB Matrix ----------#
        rgb_ModisScale = np.zeros((ht, wid,3), dtype=np.uint8)
        
        scale = np.zeros((ht,wid), dtype=np.int)
        scale[:,:] = Band_01[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
        scale[scale > MAXVALUES] = MAXVALUES-1
        scale[scale < 0] = 0
        rgb_ModisScale[:,:,0] = band_scale[scale[:,:]]
        
        scale = np.zeros((ht,wid), dtype=np.int)
        scale[:,:] = Band_04[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
        scale[scale > MAXVALUES] = MAXVALUES-1
        scale[scale < 0] = 0
        rgb_ModisScale[:,:,1] = band_scale[scale[:,:]]
        
        scale = np.zeros((ht,wid), dtype=np.int)
        scale[:,:] = Band_03[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
        scale[scale > MAXVALUES] = MAXVALUES-1
        scale[scale < 0] = 0
        rgb_ModisScale[:,:,2] = band_scale[scale[:,:]]
        
        rgb_ModisScale[rgb_ModisScale > 255] = 255
        rgb_ModisScale[rgb_ModisScale < 0] = 0
                
        r = rgb_ModisScale[:,:, 0]
        g = rgb_ModisScale[:,:, 1]
        b = rgb_ModisScale[:,:, 2]
        rgb = np.array([r,g,b]).T
        color_tuple = rgb.transpose((1,0,2)).reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0
        
        
        
    #---------- Read eMAS RGB Data ----------#
        
    eMAS_RGB_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/eMAS/'
    for item in os.listdir(eMAS_RGB_dir):
        if item[-35:-27] == eMAS_file[-45:-37]:
            eMASL1B_name = eMAS_RGB_dir + item
            eMASL1B = SD(eMASL1B_name, SDC.READ)
            
    print eMASL1B_name
    #---------- Reflectances (Bands 3, 2 and 1) from eMAS ----------#
    
    sds = eMASL1B.select('CalibratedData')
#    ht = sds.dimensions()['NumberOfScanlines']
#    wid = sds.dimensions()['NumberOfPixels']
    Band_03 = sds[::4,2,::4] * 0.1 * np.pi / 1532.2
    Band_02 = sds[::4,1,::4] * 0.1 * np.pi / 1856.6
    Band_01 = sds[::4,0,::4] * 0.1 * np.pi / 2004.7
    
    ht = Band_03.shape[0]
    wid = Band_03.shape[1]
    
    sds2 = eMASL1B.select('SolarZenithAngle')
    miu0_read = sds2.get(start=(0,0), count=(ht, wid), stride =(4,4))
    miu0 = np.cos(miu0_read/180.*np.pi)
    
    sds = eMASL1B.select('PixelLatitude')
    latitude_eMAS = sds.get(start=(0,0), count=(ht, wid), stride =(4,4))
    
    sds = eMASL1B.select('PixelLongitude')
    longitude_eMAS = sds.get(start=(0,0), count=(ht, wid), stride =(4,4))
    
    
    
    #---------- RGB Matrix ----------#
    rgb_ModisScale = np.zeros((ht, wid,3), dtype=np.uint8)
    
    scale = np.zeros((ht,wid), dtype=np.int)
    scale[:,:] = Band_03[0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
    scale[scale > MAXVALUES] = MAXVALUES-1
    scale[scale < 0] = 0
    rgb_ModisScale[:,:,0] = band_scale[scale[:,:]]
    
    scale = np.zeros((ht,wid), dtype=np.int)
    scale[:,:] = Band_02[0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
    scale[scale > MAXVALUES] = MAXVALUES-1
    scale[scale < 0] = 0
    rgb_ModisScale[:,:,1] = band_scale[scale[:,:]]
    
    scale = np.zeros((ht,wid), dtype=np.int)
    scale[:,:] = Band_01[0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
    scale[scale > MAXVALUES] = MAXVALUES-1
    scale[scale < 0] = 0
    rgb_ModisScale[:,:,2] = band_scale[scale[:,:]]
    
    rgb_ModisScale[rgb_ModisScale > 255] = 255
    rgb_ModisScale[rgb_ModisScale < 0] = 0
    
    r = rgb_ModisScale[:,:, 0]
    g = rgb_ModisScale[:,:, 1]
    b = rgb_ModisScale[:,:, 2]
    rgb = np.array([r,g,b]).T
    color_tuple_eMAS = rgb.transpose((1,0,2)).reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0
    
    
        

    ###----------- PLOTTING -----------###
        
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
#    m.bluemarble(alpha=0.75)
    m.drawcoastlines(linewidth=0.75)
    m.drawparallels(np.arange(-90., 90., 1), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(-180., 180., 1), labels=[0, 0, 0, 1])
    fig = plt.gcf()
    plt.figtext(0.1,0.025,'eMAS vs MODIS (%s)'%(eMAS_file[-55:]),size='small')
    

    # RGB MODIS
    x, y = m(longitude, latitude) # compute map projection coordinates
    indexarray = np.zeros(shape = x.shape,dtype=np.int)
    for i in range(len(indexarray)):
        for j in range(len(indexarray[i])):
            indexarray[i][j] = abs(( len(indexarray[0]) / 2 ) - j)
    m.scatter(x,y, c=color_tuple, s=(indexarray / 25 + 20), edgecolor='none')


    # RGB eMAS
    x, y = m(longitude_eMAS, latitude_eMAS) # compute map projection coordinates       ######## !!!!!!!!!!
    m.scatter(x,y, c=color_tuple_eMAS, s=20, edgecolor='none')



    # Map out background MODIS files    
    for MODIS_file in MODIS_files:
        
        # Read dataset.
        hdf = SD(MODIS_file, SDC.READ)
        DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
        data3D = hdf.select(DATAFIELD_NAME)
        MODIS = data3D[:,:]       # 1: 550 nm
        
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
    m.plot(lon_corners,lat_corners,'k-',linewidth=1.0)            
    
    # Save figure
    pngfile = "Plots_3K_RGB/Terra/{0}.png".format(eMAS_file[-55:-4])
    fig.savefig(pngfile,dpi=200)
    plt.close()
        