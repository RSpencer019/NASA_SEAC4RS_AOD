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

breaker = 0
for item in grouped_data:
    if breaker == 1000:
        break
    breaker += 1

    eMAS_file = item[0]
    Aeronet_files = item[1]['Aeronet_file']
    n = item[1]['location']
    z1 = item[1]['meanval_eMAS_550']
    z2 = item[1]['meanval_aeronet_550_intrp']
    y = item[1]['longitude(index)']
    x = item[1]['latitude(index)']

       
    
    fig = plt.figure()
    fig.clear()
    pl1 = fig.add_subplot(1,3,1)
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
    
    
    pl2 = fig.add_subplot(1,3,2)
    pl2.axes.get_xaxis().set_visible(False)
    pl2.axes.get_yaxis().set_visible(False)
    
    # Read eMAS dataset.  
    hdf = SD(eMAS_file, SDC.READ)
    DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float) * 0.001
    
    
    hour_start = eMAS_file[-27:-25]
    minute_start = eMAS_file[-25:-23]
    hour_end = eMAS_file[-22:-20]
    minute_end = eMAS_file[-20:-18]
    time_start = float(hour_start) + (float(minute_start) / 60)
    time_end = float(hour_end) + (float(minute_end) / 60)


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
    
    

    
    #### Cloud CPL ####
    
    # Find CPL file
    CPL_file = 'N/A'
    CPL_directory = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/'
    for item in os.listdir(CPL_directory):
        if 'OP' in item:
            if eMAS_file[-59:-56].lower() in item: # Filter Month Aug
                if eMAS_file[-61:-59] in item[-12:-10]: # Filter Day
                    CPL_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/'+item
            elif eMAS_file[-60:-57].lower() in item: # Filter Month Sept  
                if eMAS_file[-62:-60] in item[-12:-10]: # Filter Day
                    CPL_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/'+item
            
    if CPL_file != 'N/A':
        # Open CPL file
    #    CPL_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/CPL_OP_13955_19aug13.hdf5'
        hdf5 = h5py.File(CPL_file,'r')
        
        hour = hdf5['Hour'][:].astype(float)
        minute = hdf5['Minute'][:].astype(float)
        time = hour + (minute / 60)
        clouds = hdf5['Layer_Type'][:,0]
        CPL = np.array([time, clouds])
        CPL2 = CPL[:,(CPL[0] > time_start) & (CPL[0] < time_end)]
        
        # Plot Clouds    
        cpl_y = np.linspace(1,len(eMAS),len(CPL2[0]))
        cpl_x = np.array([len(eMAS[0]) / 2] * len(cpl_y))
        pl2.scatter(cpl_x[CPL2[1]==3]+2, cpl_y[CPL2[1]==3],linewidths=0,s=0.75,c='0.0',marker='s')
        pl2.scatter(cpl_x[CPL2[1]==3]-2, cpl_y[CPL2[1]==3],linewidths=0,s=0.75,c='0.0',marker='s')    
        pl2.set_ylim(0,2)    
            
    
    
    
    # Plot Aeronet Sites
    for loc in Aeronet_files.index:
        pl2.scatter(x[loc],y[loc],c=z1[loc],s=120,vmin=lowerval,vmax=upperval)
        pl2.scatter(x[loc],y[loc],c=z2[loc],s=40,vmin=lowerval,vmax=upperval)
    
    
#    fig.subplots_adjust(left=0.3,right=0.7)

    '''
    for i, txt in enumerate(n):
        plt.annotate(txt, (x[i]+.1,y[i]+.1))    
    for i, txt in enumerate(t):
        plt.annotate(txt, (x[i]+.1,y[i]))
    
    '''

    
    pl2.set_xlim([0,eMAS.shape[1]])
    pl2.set_ylim([eMAS.shape[0],0])  
    
    pl2.set_title('AOD')
    fig.text(0.1,0.02,CPL_file,size='x-small')

    ### Cloud Optical Depth ###
    cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_RGB_Clouds/'
    eMAS_ID = eMAS_file[-45:-37]
    for cloud_file in os.listdir(cloud_dir):
        if cloud_file[-45:-37] == eMAS_ID:
            cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)
            DATAFIELD_NAME = 'Cloud_Optical_Thickness'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Clouds = data3D[:,:]
            Clouds_lat = cloud_hdf.select('PixelLatitude')
            Clouds_lon = cloud_hdf.select('PixelLongitude')
            Clouds_lat = Clouds_lat[:,:]
            Clouds_lon = Clouds_lon[:,:]
            
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Clouds = np.ma.masked_array(Clouds, Clouds == fv)
    
    pl3 = fig.add_subplot(1,3,3)
    pl3.axes.get_xaxis().set_visible(False)
    pl3.axes.get_yaxis().set_visible(False)  
    pl3.set_axis_bgcolor('black')
    res = 1   
    c2 = pl3.imshow(Clouds.astype(float)/Clouds.max(),cmap=plt.get_cmap('pink'))

    fig.colorbar(c2,ticks=[0,0.25,0.5,0.75,1.0])
    fig.tight_layout()
    # Save plot.
    pngfile = "Plots_Individual_Merged_CPL_Clouds/{0}.png".format(eMAS_file[-55:])
    fig.savefig(pngfile,dpi=200)
    plt.close()
    
