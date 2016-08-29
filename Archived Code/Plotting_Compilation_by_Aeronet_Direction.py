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
from scipy.stats import gaussian_kde
from mpl_toolkits.basemap import Basemap
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
import h5py

# AOD scale
lowerval = 0.0
upperval = 0.75

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
HDFfile = data['eMAS_file']

for loc in range(len(HDFfile)): #28, 29):   #

    print 'Plotting...'
       
    eMAS_file = HDFfile[loc]
    eMAS_ID = eMAS_file[-45:-37]
           
    fig = plt.figure(figsize=(16,8))
    fig.clear()



    ###--------------- RGB ---------------###

        
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
            

    #---------- Read eMAS RGB Data ----------#
        
    eMAS_RGB_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/eMAS/'
    for item in os.listdir(eMAS_RGB_dir):
        if item[-35:-27] == eMAS_file[-45:-37]:
            eMASL1B_name = eMAS_RGB_dir + item
            
    eMASL1B = SD(eMASL1B_name, SDC.READ)
    print eMASL1B_name
    #---------- Reflectances (Bands 3, 2 and 1) from eMAS ----------#
    
    sds = eMASL1B.select('CalibratedData')
    Band_03 = sds[:,2,:] * 0.1 * np.pi / 1532.2
    Band_02 = sds[:,1,:] * 0.1 * np.pi / 1856.6
    Band_01 = sds[:,0,:] * 0.1 * np.pi / 2004.7
    
    ht = Band_03.shape[0]
    wid = Band_03.shape[1]
    
    sds2 = eMASL1B.select('SolarZenithAngle')
    miu0_read = sds2.get(start=(0,0), count=(ht, wid), stride =(1,1))
    miu0 = np.cos(miu0_read/180.*np.pi)
    
    
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
    
    r = rgb_ModisScale[:,:, 0].T
    g = rgb_ModisScale[:,:, 1].T
    b = rgb_ModisScale[:,:, 2].T
    rgb = np.array([r,g,b]).T
    
    # Plot RGB eMAS
    pl1 = fig.add_subplot(1,11,1)    
    pl1.axes.get_xaxis().set_visible(False)
    pl1.axes.get_yaxis().set_visible(False)  
    pl1.imshow(rgb)
    pl1.set_title('RGB',fontsize='medium')



    ###--------------- Cloud Properties ---------------###

    cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'
    for cloud_file in os.listdir(cloud_dir):
        if cloud_file[-45:-37] == eMAS_ID:
            cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)

            ### COD ###
            DATAFIELD_NAME = 'Cloud_Optical_Thickness'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Clouds = data3D[:,:].astype('float')
            '''
            Clouds_lat = cloud_hdf.select('PixelLatitude')                         #### !!!!
            Clouds_lon = cloud_hdf.select('PixelLongitude')
            Clouds_lat = Clouds_lat[:,:]
            Clouds_lon = Clouds_lon[:,:]
            '''
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Clouds = np.ma.masked_array(Clouds, Clouds == fv)
            scalefactor=attrs['scale_factor']
            sf = scalefactor[0]
            Clouds *= sf

            ### Height ###
            DATAFIELD_NAME = 'Cloud_Top_Height'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Height = data3D[:,:].astype('float')
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Height = np.ma.masked_array(Height, Height == fv)
            scalefactor=attrs['scale_factor']
            sf = scalefactor[0]
            Height *= sf
            Height /= 1000 # to km

            ### Phase ###
            DATAFIELD_NAME = 'Cloud_Phase_Optical_Properties'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Phase = data3D[:,:].astype('float')
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Phase = np.ma.masked_array(Phase, (Phase == fv) | (Phase == 0) | (Phase == 1))
            scalefactor=attrs['scale_factor']
            sf = scalefactor[0]
            Phase *= sf

            ### Temp ###
            DATAFIELD_NAME = 'Cloud_Top_Temperature'
            data3D = cloud_hdf.select(DATAFIELD_NAME)
            Temp = data3D[:,:].astype('float')
            # Handle fill value.
            attrs = data3D.attributes(full=1)
            fillvalue=attrs["_FillValue"]
            fv = fillvalue[0]
            Temp = np.ma.masked_array(Temp, Temp == fv)
            scalefactor=attrs['scale_factor']
            sf = scalefactor[0]
            offset=attrs['add_offset']
            ofs = offset[0]            
            Temp += ofs
            Temp *= sf

    
    # Plot Cloud Depth
    pl2 = fig.add_subplot(1,11,2)
    pl2.axes.get_xaxis().set_visible(False)
    pl2.axes.get_yaxis().set_visible(False)  
    pl2.set_axis_bgcolor('darkolivegreen')
    res = 1   
    c2 = pl2.imshow(Clouds,cmap=plt.get_cmap('pink_r'))
#    fig.colorbar(c2)
    pl2.set_title('COD',fontsize='medium')

    # Plot Cloud Height
    pl3 = fig.add_subplot(1,11,3)
    pl3.axes.get_xaxis().set_visible(False)
    pl3.axes.get_yaxis().set_visible(False)  
    pl3.set_axis_bgcolor('darkolivegreen')
    res = 1   
    c3 = pl3.imshow(Height,cmap=plt.get_cmap('pink_r'))
#    fig.colorbar(c3)
    pl3.set_title('Cld Height (km)',fontsize='medium')

    # Plot Cloud Phase
    pl4 = fig.add_subplot(1,11,4)
    pl4.axes.get_xaxis().set_visible(False)
    pl4.axes.get_yaxis().set_visible(False)  
    pl4.set_axis_bgcolor('darkolivegreen')
    res = 1   
    c4 = pl4.imshow(Phase,vmin=1.5,vmax=4.5,cmap=plt.get_cmap('Blues_r',3))
    formatter = plt.FuncFormatter(lambda val, loc: ['n/a','clear','Liquid','Ice','Unknown'][val])
#    plt.colorbar(c4,ticks=[2,3,4],format = formatter)    
    pl4.set_title('Cld Phase',fontsize='medium')
    

    # Plot Cloud Temp
    pl5 = fig.add_subplot(1,11,5)
    pl5.axes.get_xaxis().set_visible(False)
    pl5.axes.get_yaxis().set_visible(False)  
    pl5.set_axis_bgcolor('darkolivegreen')
    res = 1   
    c5 = pl5.imshow(Temp,cmap=plt.get_cmap('hot_r'))
#    fig.colorbar(c5)
    pl5.set_title('Cld Temp (K)',fontsize='medium')



    ###--------------- Cloud Mask ---------------###  

    # Read Cloud Mask dataset.
    hdf = SD(eMAS_file, SDC.READ) 
    DATAFIELD_NAME='Aerosol_Cldmask_Land_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    eMAS = data3D[:,:].astype(float)
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    fv = fillvalue[0]
    eMAS = np.ma.masked_array(eMAS, eMAS == fv)

    # Plot Mask
    pl6 = fig.add_subplot(1,11,6)
    pl6.axes.get_xaxis().set_visible(False)
    pl6.axes.get_yaxis().set_visible(False)
    pl6.set_title('Cld_Msk',fontsize='medium')
    pl6.imshow(eMAS,aspect='equal',cmap=plt.get_cmap('Greens'))
    
    
    
    ###--------------- eMAS AOD ---------------### 

    # Read eMAS dataset. 
    hdf = SD(eMAS_file, SDC.READ)
    dataset = hdf.select('Optical_Depth_Land_And_Ocean')
    attrs = dataset.attributes(full=1)
    fillvalue=attrs['Fill_Val']
    scalefactor=attrs['Scale_factor']
    fv = fillvalue[0]
    sf = scalefactor[0]
    aod = dataset[:,:].astype('float')
    aod[aod == fv] = np.nan
    aod *= sf 
    '''
    DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
    data3D = hdf.select(DATAFIELD_NAME)
    aod = data3D[:,:].astype(float) * 0.001
    '''
    hour_start = eMAS_file[-27:-25]
    minute_start = eMAS_file[-25:-23]
    hour_end = eMAS_file[-22:-20]
    minute_end = eMAS_file[-20:-18]
    time_start = float(hour_start) + (float(minute_start) / 60)
    time_end = float(hour_end) + (float(minute_end) / 60)
    '''
    # Handle fill value.
    attrs = data3D.attributes(full=1)
    fillvalue=attrs["Fill_Val"]
    fv = fillvalue[0] * 0.001
    aod = np.ma.masked_array(aod, aod == fv)
    '''
    # Plot eMAS
    pl8 = fig.add_subplot(1,11,8)
    pl8.axes.get_xaxis().set_visible(False)
    pl8.axes.get_yaxis().set_visible(False)
    pl8.set_xlim([0,aod.shape[1]])
    pl8.set_ylim([aod.shape[0],0]) 
    res = 1
    c6 = pl8.imshow(aod[::res,::res],vmin=lowerval,vmax=upperval,aspect='equal')
#    fig.colorbar(c6,ticks=[0,0.25,0.5,0.75])
    pl8.set_title('AOD',fontsize='medium')
    
    
    
    ###--------------- Cloud CPL ---------------### 
    
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
        # Read CPL file
        hdf5 = h5py.File(CPL_file,'r')
        hour = hdf5['Hour'][:].astype(float)
        minute = hdf5['Minute'][:].astype(float)
        second = hdf5['Second'][:].astype(float)
        time = hour + (minute / 60) + (second / 3600)
        clouds = hdf5['Layer_Type'][:,0]
        CPL = np.array([time, clouds])
        CPL2 = CPL[:,(CPL[0] > time_start) & (CPL[0] < time_end)]
        
        # Plot CPL Clouds    
        cpl_y = np.linspace(1,len(aod),len(CPL2[0]))
        cpl_x = np.array([len(aod[0]) / 2] * len(cpl_y))
        pl8.scatter(cpl_x[CPL2[1]==3]+2, cpl_y[CPL2[1]==3],linewidths=0,s=1.0,c='w',marker='s')
        pl8.scatter(cpl_x[CPL2[1]==3]-2, cpl_y[CPL2[1]==3],linewidths=0,s=1.0,c='w',marker='s')  
        pl8.scatter(cpl_x[CPL2[1]==3]+3, cpl_y[CPL2[1]==3],linewidths=0,s=2.0,c='k',marker='s')
        pl8.scatter(cpl_x[CPL2[1]==3]-3, cpl_y[CPL2[1]==3],linewidths=0,s=2.0,c='k',marker='s') 
        fig.text(0.75,0.2,CPL_file[-25:-5],size='x-small')    
    
    


    ###--------------- Direction to Cloud Scatter ---------------###
    
    cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
    cld_dir = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dir_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')

    print eMAS_file
    # Generate scatter data
    xx = np.ndarray.flatten(cld_dist)
    yy = np.ndarray.flatten(aod)
    zz = np.ndarray.flatten(cld_dir)

    x = xx[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)] / 2 #km
    y = yy[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]
    z = zz[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]
    
    # Plot
    pl7 = fig.add_subplot(1,11,7)
    pl7.axes.get_xaxis().set_visible(False)
    pl7.axes.get_yaxis().set_visible(False)
    cld_dir_plt = cld_dir
    cld_dir_plt[cld_dir == 0] = np.nan
    c7 = pl7.imshow(cld_dir_plt,cmap = 'jet_r')
    pl7.set_title('Sol_Dir',fontsize='medium')
#    plt.colorbar(c7) 

    # Direction to Cloud
    dist = 10  # Use all points within this distance (km)
    pl9 = fig.add_subplot(2,4,4)
    pl9.scatter(z[x<dist],y[x<dist],s=x[x<dist])
    plt.xlim(0,180)
    plt.xticks([0,90,180])
    plt.ylim(-0.05,1)
    plt.yticks([0.0,0.25,0.5,0.75,1.0])
    plt.xlabel('Solar Direction (deg)')
    plt.ylabel('AOD')

    
    
    ###--------------- Direction to Cloud Scatter (CHARACTERIZATION) ---------------###
  
    t = np.polyfit(z[x<dist],y[x<dist],1)
    p = np.poly1d(t)
    r2 = r2_score(y[x<dist], p(z[x<dist]))
    #plt.plot([0,180],p([0,180]),'r-',lw=2)
    fig.text(0.75,0.4,'Total Solar Bias: %.3f at %i km'%(-t[0]*180,dist),size='small')
    t = np.polyfit(z[(x<dist) & (z<90)],y[(x<dist) & (z<90)],1)
    p = np.poly1d(t)
    r2 = r2_score(y[(x<dist) & (z<90)], p(z[(x<dist) & (z<90)]))
    plt.plot([0,90],p([0,90]),'r-',lw=2)

    t = np.polyfit(z[(x<dist) & (z>90)],y[(x<dist) & (z>90)],1)
    p = np.poly1d(t)
    r2 = r2_score(y[(x<dist) & (z>90)], p(z[(x<dist) & (z>90)]))
    plt.plot([90,180],p([90,180]),'r-',lw=2)

    plt.axvline(x=90,c='k')

    ###--------------- Save Plot ---------------###

    # Save plot.
#    fig.tight_layout()
    pngfile = "Plots_Compilation_by_Aeronet_Direction/{0}_L2CLD.png".format(eMAS_file[-55:-4])
    fig.savefig(pngfile,dpi=200)
    plt.close()


