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

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
grouped_data = data.groupby(data['eMAS_file'])

for item in grouped_data:
    eMAS_file = item[0]
    Aeronet_files = item[1]['Aeronet_file']
    n = item[1]['location']
    z1 = item[1]['meanval_eMAS_550']
    z2 = item[1]['meanval_aeronet_550_intrp']
    y = item[1]['longitude(index)']
    x = item[1]['latitude(index)']

    
    
    
    fig = plt.figure()
    fig.clear()
    pl1 = fig.add_subplot(1,2,1)
    pl1.axes.get_xaxis().set_visible(False)
    pl1.axes.get_yaxis().set_visible(False)
    plt.title('Cloud Fraction')
        
    # Read Cloud Mask dataset.
    hdf = SD(eMAS_file, SDC.READ) 
    DATAFIELD_NAME='Aerosol_Cloud_Fraction_Land'
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
    res = 1
    c2 = pl1.imshow(eMAS[::res,::res],aspect='equal',cmap=plt.get_cmap('Greens_r'))


    
    #######    
    
    
    pl2 = fig.add_subplot(1,2,2)
    pl2.axes.get_xaxis().set_visible(False)
    pl2.axes.get_yaxis().set_visible(False)
    
    # Read eMAS dataset.  
    hdf = SD(eMAS_file, SDC.READ)
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
    for loc in Aeronet_files.index:
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
    pngfile = "Plots_Individual_Merged/{0}.png".format(eMAS_file[-55:])
    fig.savefig(pngfile,dpi=200)
    plt.close()
    
