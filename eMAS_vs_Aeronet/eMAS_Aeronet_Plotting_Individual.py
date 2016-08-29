#------------------------------------------------------------------------------
# Name:        eMAS (SEAC4RS) vs. MODIS & Aeronet Plotting
# Description: NASA Aerosol Project
#
# Author:      Robert S. Spencer
#
# Created:     6/28/2016
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

lowerval = 0.0
upperval = 0.75

data = pd.DataFrame.from_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)

y = data['longitude(index)']
x = data['latitude(index)']
n = data['location']
t = data.index
HDFfile = data['eMAS_file']
indices = len(x)
z1 = data['meanval_eMAS_550']
z2 = data['meanval_aeronet_550_intrp']

for loc in range(indices):
    
    fig = plt.figure()
    fig.clear()
    plt.figtext(0.4,0.075,n[loc])
    pl1 = fig.add_subplot(1,2,1)
    pl1.axes.get_xaxis().set_visible(False)
    pl1.axes.get_yaxis().set_visible(False)
    plt.title('Cloud Mask')
    
    # Read Cloud Mask dataset.
    hdf = SD(HDFfile[loc], SDC.READ) 
    DATAFIELD_NAME='Aerosol_Cldmask_Land_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float)
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0]
    
    #data[data == fv] = np.nan
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)
    
    # Plot Mask
    res = 10
    c2 = pl1.imshow(eMAS[::res,::res],vmin=lowerval,vmax=upperval,aspect='equal',cmap=plt.get_cmap('Greens'))

    
    #######    
    
    
    pl2 = fig.add_subplot(1,2,2)
    pl2.axes.get_xaxis().set_visible(False)
    pl2.axes.get_yaxis().set_visible(False)
    
    # Read eMAS dataset.  
    hdf = SD(HDFfile[loc], SDC.READ)
    DATAFIELD_NAME='Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float) * 0.001
    
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    
    # fillvalue[0] is the attribute value.
    fv = fillvalue[0] * 0.001
    
    #data[data == fv] = np.nan
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)
        
    # Plot eMAS
    res = 1
    c1 = pl2.imshow(eMAS[::res,::res],vmin=lowerval,vmax=upperval,aspect='equal')
    
  
    fig.colorbar(c1,ticks=[0,0.25,0.5,0.75])
    
    
    # Plot Aeronet Sites
    plt.scatter(x[loc],y[loc],c=z1[loc],s=120,vmin=lowerval,vmax=upperval)
    plt.scatter(x[loc],y[loc],c=z2[loc],s=40,vmin=lowerval,vmax=upperval)
    
    
    fig.subplots_adjust(left=0.3,right=0.7)

    '''
    for i, txt in enumerate(n):
        plt.annotate(txt, (x[i]+.1,y[i]+.1))    
    for i, txt in enumerate(t):
        plt.annotate(txt, (x[i]+.1,y[i]))
    
    '''

    
    plt.xlim([0,eMAS.shape[1]])
    plt.ylim([eMAS.shape[0],0])  
    
    plt.title('AOD')
    # Save plot.
    pngfile = "Plots_Individual/{0}_{1}.png".format(n[loc],t[loc])
    fig.savefig(pngfile,dpi=200)
    plt.close()
    
