#------------------------------------------------------------------------------
# Name:        Distance to Cloud Generator
# Description: Generates the distance to cloud from cloud mask
#
# Author:      Robert S. Spencer
#
# Created:     7/11/2016
# Python:      2.7
#------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import time as tm

start_time = tm.time()

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'

#y = data['longitude(index)']
#x = data['latitude(index)']
n = data['location']
#t = data.index
HDFfile = data['eMAS_file']
indices = len(n)
#z1 = data['meanval_eMAS_550']
#z2 = data['meanval_aeronet_550_intrp']

for loc in range(indices):           # 28,29    # Modify for parallel computing !!!!!!!!
    print 'Computing file...'    
    
    eMAS_file = HDFfile[loc]
    eMAS_ID = eMAS_file[-45:-37]
    cloud_hdf = ''
    for cloud_file in os.listdir(cloud_dir):
        if cloud_file[-45:-37] == eMAS_ID:
            print cloud_file
            cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)    
    
    print eMAS_file
    print cloud_hdf

    # Solar Azimuth Angle
    dataset = cloud_hdf.select('SolarAzimuthAngle')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['_FillValue']
    fv = float(fillvalue[0])
    scale_factor=attrs['scale_factor']
    sf = float(scale_factor[0])
    add_offset=attrs['add_offset']
    aoff = float(add_offset[0])
    sol_azim = dataset[:,:].astype(float)
    sol_azim[sol_azim == fv] = np.nan
    sol_azim = ( np.nanmean(sol_azim) + aoff ) * sf
    print sol_azim

    # Sensor Azimuth Angle
    dataset = cloud_hdf.select('SensorAzimuthAngle')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['_FillValue']
    fv = float(fillvalue[0])
    scale_factor=attrs['scale_factor']
    sf = float(scale_factor[0])
    add_offset=attrs['add_offset']
    aoff = float(add_offset[0])
    sen_azim = dataset[:,:].astype(float)
    sen_azim[sen_azim == fv] = np.nan                   # Aircraft's right side
    sen_azim = ( sen_azim[0,1] + aoff ) * sf           # Top of Granule direction
    print sen_azim

    # Aerosol_Cldmask_Land_Ocean
    dataset = cloud_hdf.select('Cloud_Top_Height')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['_FillValue']
    fv = float(fillvalue[0])
    cld_msk = dataset[:,:].astype(float)
    # handle the values along the boundaries (not sure why they exist...)
    cld_msk[1] = fv
    cld_msk[-2] = fv
    cld_msk[:,1] = fv
    cld_msk[:,-2] = fv
    # convert to mask from cloud height dataset
    cld_msk[cld_msk > -1] = 0
    cld_msk[cld_msk == fv] = 1
    
    # Initiate distance / direction rasters
    cld_dist = np.empty([cld_msk.shape[0],cld_msk.shape[1]])
    cld_dist.fill(np.nan)
    cld_direction = np.empty([cld_msk.shape[0],cld_msk.shape[1]])
    cld_direction.fill(np.nan)

    rows = cld_msk.shape[0]
    cols = cld_msk.shape[1]
    print "total rows: ", rows

    for i in range(rows):
        n = 0
        if i % 100 == 0:
            print "row: ", i
        for j in range(cols):
            if cld_msk[i,j] == 1: # if clear
                while True:
    
                    # Determines the next step size, s
                    # Step size gets added to search radias, n
                    # Optimized: Once cloud is found, the next pixel starts off with n - s instead of 0
                    if n < 10:
                        s = 1
                    elif n < 50:
                        s = 5
                    elif n < 100:
                        s = 10
                    elif n < 358: # half the swath width
                        s = 50
                    if n >= 358:
                        cld_dist[i,j] = n
                        n -= s
                        break

                    if n>i:
                        rowl = 0
                    else:
                        rowl = i - n
                    if n>j:
                        coll = 0
                    else:
                        coll = j - n
                    if n>rows-i:
                        rowu = rows
                    else:
                        rowu = i + n + 1
                    if n>cols-j:
                        colu = cols
                    else:
                        colu = j + n + 1
                                        

                    if 0 in cld_msk[rowl:rowl+s+1,coll:colu]: # TOP

                        iraw, jraw = np.where(cld_msk[rowl:rowl+s+1,coll:colu]==0)
                        i2 = rowl + iraw
                        j2 = coll + jraw
                        di = i2 - i
                        dj = j2 - j
                        dist = np.sqrt( di ** 2 + dj ** 2 )
                        di = di[np.argmin(dist)]
                        dj = dj[np.argmin(dist)]
                        cld_dist[i,j] = dist.min()

                        if dj == 0:
                            dj += 0.01
                        cld_azim = np.arctan(float(di)/float(dj)) / (2 * np.pi) * 360
                        if (dj < 0):
                            cld_azim += 180
                        if (di < 0) & (dj > 0):
                            cld_azim += 360
                        cld_direction[i,j] = ( sen_azim + cld_azim ) - sol_azim                                            
                        if cld_direction[i,j] > 360:
                            cld_direction[i,j] -= 360  
                        if cld_direction[i,j] < 0:
                            cld_direction[i,j] += 360  
                        n -= s
                        break

                    elif 0 in cld_msk[rowl:rowu,colu-s-1:colu]: # RIGHT
                        
                        iraw, jraw = np.where(cld_msk[rowl:rowu,colu-s-1:colu]==0)
                        i2 = rowl + iraw
                        j2 = colu-s-1 + jraw
                        di = i2 - i
                        dj = j2 - j
                        dist = np.sqrt( di ** 2 + dj ** 2 )
                        di = di[np.argmin(dist)]
                        dj = dj[np.argmin(dist)]
                        cld_dist[i,j] = dist.min()
                        
                        if dj == 0:
                            dj += 0.01
                        cld_azim = np.arctan(float(di)/float(dj)) / (2 * np.pi) * 360
                        if (dj < 0):
                            cld_azim += 180
                        if (di < 0) & (dj > 0):
                            cld_azim += 360
                        cld_direction[i,j] = ( sen_azim + cld_azim ) - sol_azim
                        if cld_direction[i,j] > 360:
                            cld_direction[i,j] -= 360  
                        if cld_direction[i,j] < 0:
                            cld_direction[i,j] += 360  
                        n -= s
                        break

                    elif 0 in cld_msk[rowu-s-1:rowu,coll:colu]: # BOTTOM
                        
                        iraw, jraw = np.where(cld_msk[rowu-s-1:rowu,coll:colu]==0)
                        i2 = rowu-s-1 + iraw
                        j2 = coll + jraw
                        di = i2 - i
                        dj = j2 - j
                        dist = np.sqrt( di ** 2 + dj ** 2 )
                        di = di[np.argmin(dist)]
                        dj = dj[np.argmin(dist)]
                        cld_dist[i,j] = dist.min()
                        
                        if dj == 0:
                            dj += 0.01
                        cld_azim = np.arctan(float(di)/float(dj)) / (2 * np.pi) * 360
                        if (dj < 0):
                            cld_azim += 180
                        if (di < 0) & (dj > 0):
                            cld_azim += 360
                        cld_direction[i,j] = ( sen_azim + cld_azim ) - sol_azim
                        if cld_direction[i,j] > 360:
                            cld_direction[i,j] -= 360  
                        if cld_direction[i,j] < 0:
                            cld_direction[i,j] += 360  
                        n -= s
                        break

                    elif 0 in cld_msk[rowl:rowu,coll:coll+s+1]: # LEFT
                        
                        iraw, jraw = np.where(cld_msk[rowl:rowu,coll:coll+s+1]==0)
                        i2 = rowl + iraw
                        j2 = coll + jraw
                        di = i2 - i
                        dj = j2 - j
                        dist = np.sqrt( di ** 2 + dj ** 2 )
                        di = di[np.argmin(dist)]
                        dj = dj[np.argmin(dist)]
                        cld_dist[i,j] = dist.min()
                        
                        if dj == 0:
                            dj += 0.01
                        cld_azim = np.arctan(float(di)/float(dj)) / (2 * np.pi) * 360
                        if (dj < 0):
                            cld_azim += 180
                        if (di < 0) & (dj > 0):
                            cld_azim += 360
                        cld_direction[i,j] = ( sen_azim + cld_azim ) - sol_azim
                        if cld_direction[i,j] > 360:
                            cld_direction[i,j] -= 360  
                        if cld_direction[i,j] < 0:
                            cld_direction[i,j] += 360  
                        n -= s
                        break

                    n += s

            if cld_msk[i,j] == 0:
                cld_dist[i,j] = 0

    cld_direction = np.abs( cld_direction - 180 )
    
    def rebin(a, shape):
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).mean(-1).mean(1)
    
    rowsnip = cld_dist.shape[0]%10
    colsnip = cld_dist.shape[1]%10

    cld_dist_lowres = rebin(cld_dist[:len(cld_dist)-rowsnip,:len(cld_dist[0])-colsnip],[cld_dist.shape[0]/10,cld_dist.shape[1]/10])    
    cld_direction_lowres = rebin(cld_direction[:len(cld_direction)-rowsnip,:len(cld_direction[0])-colsnip],[cld_direction.shape[0]/10,cld_direction.shape[1]/10])
    print cld_dist.shape
    print cld_dist_lowres.shape    
    
    np.savetxt("DATA/Dist_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), cld_dist_lowres, delimiter=",")
    np.savetxt("DATA/Dir_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), cld_direction_lowres, delimiter=",")
    
print("--- %s seconds ---" % (tm.time() - start_time))