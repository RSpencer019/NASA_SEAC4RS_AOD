#------------------------------------------------------------------------------
# Name:        eMAS (SEAC4RS) Plotting Compilation by Aeronet Collocations
# Description: NASA Aerosol Project
#
# Author:      Robert S. Spencer
#
# Created:     7/27/2016
# Python:      2.7
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


if os.path.exists('Background_AOD.csv'):
    os.remove('Background_AOD.csv')

# AOD scale
lowerval = 0.0
upperval = 0.75

data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Compiled_Cleaned.csv',header=0)
grouped_data = data.groupby(data['eMAS_file'])

breaker = 0
for item in grouped_data:
    if breaker == 1000:
        break
    breaker += 1
    print 'Plotting...'
    
    eMAS_file = item[0]
    eMAS_ID = eMAS_file[-45:-37]
    Aeronet_files = item[1]['Aeronet_file']
    n = item[1]['location']
    z1 = item[1]['meanval_eMAS_550']
    z2 = item[1]['meanval_aeronet_550_intrp']
    n_aero = item[1]['nval_aeronet']
    var_aero = item[1]['variance_aeronet']
    y = item[1]['longitude(index)']
    x = item[1]['latitude(index)']
           
    fig = plt.figure(figsize=(18,8))
    fig.clear()



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

    hour_start = eMAS_file[-27:-25]
    minute_start = eMAS_file[-25:-23]
    hour_end = eMAS_file[-22:-20]
    minute_end = eMAS_file[-20:-18]
    time_start = float(hour_start) + (float(minute_start) / 60)
    time_end = float(hour_end) + (float(minute_end) / 60)

    
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
        


    ###--------------- Distance to Cloud Scatter ---------------###
    
    cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
    print eMAS_file
    # Generate scatter data
    xx = np.ndarray.flatten(cld_dist)
    yy = np.ndarray.flatten(aod)
    x = xx[np.isfinite(xx) & np.isfinite(yy)] / 20 #km
    y = yy[np.isfinite(xx) & np.isfinite(yy)]
    
    # Calculate the point density
    xy = np.vstack([x,y])
    z = np.zeros(shape=x.shape)
    if xy[0].mean() != xy[0,1]:
        z = gaussian_kde(xy)(xy)
    
        # Sort the points by density, so that the densest points are plotted last
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]
    

    # Plot
    pl1 = fig.add_subplot(2,3,1)
    c1 = pl1.scatter(x,y,c=z,s=5,edgecolor='',cmap=plt.get_cmap('hot'))
    plt.xlim(xmin=0)
    plt.ylim(-0.05,1)
    plt.xlabel('Distance to Cloud (km)')
    plt.ylabel('AOD')
    plt.axvline(x=10,lw=0.75,c='k',ls='--')
#    plt.colorbar(c1) 

    aod2 = y
    dist2 = x
    dens2 = z

    
    ###--------------- Distance to Cloud Scatter (CHARACTERIZATION) ---------------###
    
    cld_dist = np.round(cld_dist)
    cld_dist = np.ndarray.flatten(cld_dist)
    aod = np.ndarray.flatten(aod)
    df = pd.DataFrame(data={'aod':aod,'cld_dist':cld_dist})
    df_g = (pd.groupby(df['aod'],df['cld_dist']))
    df_gm = df_g.mean()
    df_gm = np.array(df_gm)
    
    background_aod = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)])
    near_cloud_aod = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)][::-1])
    range_aod = near_cloud_aod - background_aod
    
    aeronet_aod = np.mean(z2)
    eMAS_aod = np.mean(z1)
    aeronet_n = np.sum(n_aero)
    aeronet_var = np.mean(var_aero)    
    
    near_cloud_perc = 100 * ( near_cloud_aod - aeronet_aod ) / aeronet_aod
    background_perc = 100 * ( background_aod - aeronet_aod ) / aeronet_aod
    eMAS_sampel_perc = 100 * ( eMAS_aod - aeronet_aod ) / aeronet_aod

    plt.axhline(y=aeronet_aod,lw=2.0,c='g')
    plt.axhline(y=background_aod, lw=1.75, c='b')
    plt.axhline(y=near_cloud_aod, lw=1.5, c='r')
    
    plt.figtext(0.15,0.4,'AOD Range: %.2f'%(range_aod))
    plt.figtext(0.15,0.38,r'Aeronet ($n: %i, \sigma: %.2f$)'%(aeronet_n, np.sqrt(aeronet_var)),color='g')
    plt.figtext(0.15,0.36,r'Near Cloud AOD ($\epsilon: %.1f$%%)'%(near_cloud_perc),color='r')
    plt.figtext(0.15,0.34,r'Background AOD ($\epsilon: %.1f$%%)'%(background_perc),color='b')
    plt.figtext(0.15,0.32,r'Collocation AOD ($\epsilon: %.1f$%%)'%(eMAS_sampel_perc),color='k')
 


    ### Compile Background AOD csv ###
    sites = open('Background_AOD.csv','a')
    if os.path.getsize('Background_AOD.csv') == 0:
        sites.write('eMAS_file,Background_AOD,Ave_Aeronet_AOD\n')
    sites.write(eMAS_file+','+str(background_aod)+','+str(aeronet_aod)+'\n')
    sites.close()




    ###--------------- Direction to Cloud Scatter ---------------###
    
    cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
    cld_dir = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dir_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')

    print eMAS_file
    # Generate scatter data
    xx = np.ndarray.flatten(cld_dist)
    yy = np.ndarray.flatten(aod)
    zz = np.ndarray.flatten(cld_dir)

    x = xx[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)] / 20 #km
    y = yy[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]
    z = zz[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]
    
    # Direction to Cloud
    dist = 1  # Use all points within this distance (km)
    pl2 = fig.add_subplot(2,3,2)
    pl2.scatter(z[x<dist],y[x<dist],s=x[x<dist])
    plt.xlim(0,180)
    plt.xticks([0,90,180])
    plt.ylim(-0.05,1)
    plt.yticks([0.0,0.25,0.5,0.75,1.0])
    plt.xlabel('Solar Direction (deg)')
    plt.ylabel('AOD')

    
    
    ###--------------- Direction to Cloud Scatter (CHARACTERIZATION) ---------------###
    if len(y[x<dist]) != 0:
        t = np.polyfit(z[x<dist],y[x<dist],1)
        p = np.poly1d(t)
        r2 = r2_score(y[x<dist], p(z[x<dist]))
        #plt.plot([0,180],p([0,180]),'r-',lw=2)
        fig.text(0.4,0.4,'Total Solar Bias: %.3f at %i km'%(-t[0]*180,dist))

    if len(y[(x<dist) & (z<90)]) != 0:
        t = np.polyfit(z[(x<dist) & (z<90)],y[(x<dist) & (z<90)],1)
        p = np.poly1d(t)
        r2 = r2_score(y[(x<dist) & (z<90)], p(z[(x<dist) & (z<90)]))
        plt.plot([0,90],p([0,90]),'r-',lw=2)

    if len(y[(x<dist) & (z>90)]) != 0:
        t = np.polyfit(z[(x<dist) & (z>90)],y[(x<dist) & (z>90)],1)
        p = np.poly1d(t)
        r2 = r2_score(y[(x<dist) & (z>90)], p(z[(x<dist) & (z>90)]))
        plt.plot([90,180],p([90,180]),'r-',lw=2)

    plt.axvline(x=90,c='k')






    ###--------------- Direction and Distance to Cloud Graphing ---------------###
    
    cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
    cld_dir = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dir_Cloud_Rasters_L2CLD/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')

    # Generate layers data
    xx = np.ndarray.flatten(cld_dist)
    yy = np.ndarray.flatten(aod)
    zz = np.ndarray.flatten(cld_dir)

    x = xx[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)] / 20 #km
    y = yy[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]
    z = zz[np.isfinite(xx) & np.isfinite(yy) & np.isfinite(zz)]

    aod_front = y[z<45]
    dist_front = x[z<45]
    aod_back = y[z>135]
    dist_back = x[z>135]
    aod_side = y[(45<z) & (z<135)]
    dist_side = x[(45<z) & (z<135)]

    dist_front = np.ceil(dist_front * 4) / 4  # bins into 0-0.25, 0.25-0.5, etc
    data_front = pd.DataFrame(np.transpose([dist_front,aod_front]), columns=['dist','aod'])
    grouped_front = data_front.groupby(by='dist')
    mean_front = np.array(grouped_front.mean()).flatten()
    std_front = np.array(grouped_front.std()).flatten()
    bins_front = np.array(grouped_front.std().index)
    

    dist_back = np.ceil(dist_back * 4) / 4  # bins into 0-0.25, 0.25-0.5, etc
    data_back = pd.DataFrame(np.transpose([dist_back,aod_back]), columns=['dist','aod'])
    grouped_back = data_back.groupby(by='dist')
    mean_back = np.array(grouped_back.mean()).flatten()
    std_back = np.array(grouped_back.std()).flatten()
    bins_back = np.array(grouped_back.std().index)
    

    dist_side = np.ceil(dist_side * 4) / 4  # bins into 0-0.25, 0.25-0.5, etc
    data_side = pd.DataFrame(np.transpose([dist_side,aod_side]), columns=['dist','aod'])
    grouped_side = data_side.groupby(by='dist')
    mean_side = np.array(grouped_side.mean()).flatten()
    std_side = np.array(grouped_side.std()).flatten()
    bins_side = np.array(grouped_side.std().index)

    # Plot
    pl3 = fig.add_subplot(2,3,3)
    
    plt.fill_between(bins_front, mean_front-std_front, mean_front+std_front, alpha = 0.35, color='r')
    plt.fill_between(bins_back, mean_back-std_back, mean_back+std_back, alpha = 0.35, color='b')
    
    plt.plot(bins_front, mean_front, color='r')
    plt.plot(bins_back, mean_back, color='b')
    plt.plot(bins_side, mean_side, color='g')

    plt.axvline(x=10,lw=0.75,c='k',ls='--')
    plt.xlim(0.5,15)
    plt.ylim(-0.05,1)
    plt.yticks([0.0,0.25,0.5,0.75,1.0])
    plt.xlabel('Distance to Cloud (km)')
    plt.ylabel('AOD')

    # Aeronet
    plt.axhline(y=aeronet_aod,lw=1.2,c='k')

    ### Solar Zenith Angle ###
    stats = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_Statistics.csv')
    sol_zen = float(stats.solar_zenith.loc[stats['eMAS_file'] == eMAS_file])
    plt.figtext(0.4,0.3,'Solar Zenith Angle: %.2f degrees'%(sol_zen))
    plt.figtext(0.4,0.28,'Solar Side',color='r')
    plt.figtext(0.4,0.26,'Shadow Side',color='b')
    plt.figtext(0.4,0.24,'Neutral Sides',color='g')



    ###--------------- Distance Distribution ---------------###
    pl4 = fig.add_subplot(2,3,6)
    x_cum = np.linspace(1,len(x[x<20]),len(x[x<20])) / len(x[x<20])
    plt.plot(sorted(x[x<20]),x_cum,color='k')
    plt.fill_between(sorted(x[x<20]),0,x_cum,alpha = 0.2)
    plt.xlim(0,15)
    plt.ylim(0,1)
    plt.yticks([0.0,0.25,0.5,0.75,1.0])
    plt.xlabel('Distance to Cloud (km)')
    plt.ylabel('Cummulative')
    plt.grid()

    x_cum = np.linspace(1,len(x[(z<45)&(x<20)]),len(x[(z<45)&(x<20)])) / len(x[x<20])
    plt.plot(sorted(x[(z<45)&(x<20)]),x_cum,color='r',lw=2)

    x_cum = np.linspace(1,len(x[(z>135)&(x<20)]),len(x[(z>135)&(x<20)])) / len(x[x<20])
    plt.plot(sorted(x[(z>135)&(x<20)]),x_cum,color='b',lw=2)

    x_cum = np.linspace(1,len(x[(45<z)&(z<135)&(x<20)]),len(x[(45<z)&(z<135)&(x<20)])) / len(x[x<20])
    plt.plot(sorted(x[(45<z)&(z<135)&(x<20)]),x_cum,color='g',lw=2)

    ###--------------- Save Plot ---------------###

    # Save plot.
    pngfile = "Plots_Compilation_by_Aeronet/{0}_Graphs.png".format(eMAS_file[-55:-4])
    fig.savefig(pngfile,dpi=200)
    plt.close()


    ####################### NEW SCATTER PLOT w/ Density Contours and Histograms ############

    fig2 = plt.figure(1, figsize=(8, 8))

    ###### Scatter plot with histogram axes ######

    from matplotlib.ticker import NullFormatter

    nullfmt = NullFormatter()         # no labels

    x = dist2
    y = aod2
    z = dens2

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.004

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.1]
    rect_histy = [left_h, bottom, 0.1, height]

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # no histogram axes
    axHistx.axes.get_yaxis().set_visible(False)
    axHisty.axes.get_xaxis().set_visible(False) 
    axHistx.axis('off')
    axHisty.axis('off')

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.scatter(x, y, c=z, s=3, edgecolor='', vmin = z.mean(), vmax = z.max()*3, cmap=plt.get_cmap('Blues_r'))

    # now determine nice limits by hand:
    xbinwidth = 0.1
    ybinwidth = 0.005
    xlim = 15
    ylim = 0.75

    axScatter.set_xlim((0, xlim))
    axScatter.set_ylim((-0.05, ylim))

    xbins = np.arange(0, xlim + xbinwidth, xbinwidth)
    ybins = np.arange(-0.05, ylim + ybinwidth, ybinwidth)

    axHistx.hist(x, bins=xbins, color='navy')
    axHisty.hist(y, bins=ybins, orientation='horizontal', color='navy')

    axHistx.set_xlim(axScatter.get_xlim())
    axHisty.set_ylim(axScatter.get_ylim())


    ###### Scatter plot with density contours ######

    H, xedges, yedges = np.histogram2d(x,y,bins=50)
    extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
    cset1 = axScatter.contour(H.T,extent=extent,colors='chartreuse',linewidths=0.75)


    ########################################################################################



    ###--------------- Save Plot ---------------###

    # Save plot.
    pngfile = "Plots_Compilation_by_Aeronet/{0}_Graphs2.png".format(eMAS_file[-55:-4])
    fig2.savefig(pngfile,dpi=200)
    plt.close()
    

