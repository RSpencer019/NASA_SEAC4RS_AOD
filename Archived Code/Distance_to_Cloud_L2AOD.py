#------------------------------------------------------------------------------
# Name:        Distance to Cloud Generator
# Description: Generates the distance to cloud from cloud mask
#
# Author:      Robert S. Spencer
#
# Created:     7/11/2016
# Python:      2.7
#------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import time as tm

start_time = tm.time()

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)

#y = data['longitude(index)']
#x = data['latitude(index)']
n = data['location']
#t = data.index
HDFfile = data['eMAS_file']
indices = len(n)
#z1 = data['meanval_eMAS_550']
#z2 = data['meanval_aeronet_550_intrp']

for loc in range(0,indices):
    print 'Computing file...'    
    eMAS_file = HDFfile[loc]
    
    hdf = SD(eMAS_file, SDC.READ)
    print eMAS_file
    # Aerosol_Cldmask_Land_Ocean
    dataset = hdf.select('Aerosol_Cldmask_Land_Ocean')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    fv = fillvalue[0]
    cld_msk = dataset[:,:].astype('float')
    cld_msk[cld_msk == fv] = 1
        
                
    cld_dist = np.empty([cld_msk.shape[0],cld_msk.shape[1]])
    cld_dist.fill(np.nan)
    rows = cld_msk.shape[0]
    cols = cld_msk.shape[1]
    
    for i in range(rows):
        if i % 100 == 0:
            print i
        for j in range(cols):
            n = 0
            if cld_msk[i,j] == 1: # if clear
                while True:
    
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
    
                    if 0 in cld_msk[rowl:rowu,coll:colu]:
                        cld_dist[i,j] = n
                        break
    
                    n += 1
                    if n > rows:
                        cld_dist[i,j] = n
                        break

            if cld_msk[i,j] == 0:
                cld_dist[i,j] = 0


    def rebin(a, shape):
        sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
        return a.reshape(sh).mean(-1).mean(1)

    cld_dist_lowres = rebin(cld_dist,[cld_dist.shape[0]/10,cld_dist.shape[1]/10])    
    
    print cld_dist.shape
    print cld_dist_lowres.shape       
    
    np.savetxt("Dist_Cloud_Rasters/{0}.csv".format(eMAS_file[-55:-4]), cld_dist_lowres, delimiter=",")
    
print("--- %s seconds ---" % (tm.time() - start_time))