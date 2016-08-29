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

for loc in range(2,indices,3):                # Modify for parallel computing !!!!!!!!
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
    
                
    cld_dist = np.empty([cld_msk.shape[0],cld_msk.shape[1]])
    cld_dist.fill(np.nan)
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
    
                    if 0 in cld_msk[rowl:rowu,coll:coll+s+1]: # LEFT
                        cld_dist[i,j] = n
                        n -= s
                        break
                    if 0 in cld_msk[rowl:rowu,colu-s-1:colu]: # RIGHT
                        cld_dist[i,j] = n
                        n -= s
                        break
                    if 0 in cld_msk[rowl:rowl+s+1,coll:colu]: # TOP
                        cld_dist[i,j] = n
                        n -= s
                        break
                    if 0 in cld_msk[rowu-s-1:rowu,coll:colu]: # BOTTOM
                        cld_dist[i,j] = n
                        n -= s
                        break

                    n += s

            if cld_msk[i,j] == 0:
                cld_dist[i,j] = 0
    
    
    def rebin(a, shape):
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).mean(-1).mean(1)
    
    rowsnip = cld_dist.shape[0]%10
    colsnip = cld_dist.shape[1]%10

    cld_dist_lowres = rebin(cld_dist[:len(cld_dist)-rowsnip,:len(cld_dist[0])-colsnip],[cld_dist.shape[0]/10,cld_dist.shape[1]/10])    
    
    print cld_dist.shape
    print cld_dist_lowres.shape    
    
    np.savetxt("Dist_Cloud_Rasters_L2CLD2/{0}.csv".format(eMAS_file[-55:-4]), cld_dist_lowres, delimiter=",")
    
print("--- %s seconds ---" % (tm.time() - start_time))