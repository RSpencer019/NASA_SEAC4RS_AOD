#------------------------------------------------------------------------------
# Name:        SEAC4RS Campaign Aggregation
# Description: Aggregates AOD and Cloud Property statistics for the entire campaign
#
# Author:      Robert S. Spencer
#
# Created:     8/9/2016
# Python:      2.7
#------------------------------------------------------------------------------

import os
import numpy as np
import pandas as pd
from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.stats.mstats import mode
import time as tm
import csv

start_time = tm.time()

directory = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data'
path_names = []
ext = 'hdf'
cwd = os.getcwd()

# Searches for specified files in data directory
def file_search():
    complete = False
    currentfolderdict = {'CWD':directory}
    while complete is False:      
        newfolderdict = {}      
        for folder in currentfolderdict:
            os.chdir(currentfolderdict[folder])                        
            for item in os.listdir(currentfolderdict[folder]):
                if os.path.isdir(item) == True:
                    newfolderdict[item] = os.path.abspath(item)
                elif item[-len(ext):] == ext:
                    path_names.append(os.path.abspath(item))
        if len(newfolderdict) == 0:
            complete = True
        else:
            currentfolderdict = newfolderdict

file_search()
os.chdir(cwd)

cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'

# Initialize the compiled file
if os.path.exists(cwd+'/SEAC4RS_Aggregation_Compiled.csv'):
    os.remove(cwd+'/SEAC4RS_Aggregation_Compiled.csv')
# Generate csv labels
labels = 'eMAS_ID,COD_sect_ave,COD_sect_std,COD_sect_n,Height_sect_ave,Height_sect_std,Height_sect_n,Phase_sect_ave,Phase_sect_std,Phase_sect_n,Phase_sect_mode,Phase_sect_mode_perc,Temp_sect_ave,Temp_sect_std,Temp_sect_n,background_aod,near_cloud_aod,background_aod_sect,near_cloud_aod_sect,solar_zen'
output = open(cwd+'/SEAC4RS_Aggregation_Compiled.csv', "w")
writer = csv.writer(output, lineterminator='\n')
output.write(labels+'\n')

for i in range(0, len(path_names), 3):
    print "================="
    print 'Computing file...' 

    # Initialize
    background_aod, near_cloud_aod, background_aod_sect, near_cloud_aod_sect, solar_zen = (np.nan for i in range(5))
    COD_sect_ave, COD_sect_std, COD_sect_n, Height_sect_ave, Height_sect_std, Height_sect_n, Phase_sect_ave, Phase_sect_std, Phase_sect_n, Phase_sect_mode, Phase_sect_mode_perc, Temp_sect_ave, Temp_sect_std, Temp_sect_n = (np.nan for i in range(14))

    #####------------------------------------------------#####
    #####--------------- Load Datasets ------------------#####
    #####------------------------------------------------#####


    # Parent filename (all other datasets follow this)
    eMAS_file = path_names[i]
    eMAS_ID = eMAS_file[-45:-37]

    if os.path.exists("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD_ALL/{0}.csv".format(eMAS_file[-55:-4])):

        ######------ Read eMAS dataset ------######
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

        dataset = hdf.select('Solar_Zenith')
        attrs = dataset.attributes(full=1)
        fillvalue=attrs['Fill_Val']
        scalefactor=attrs['Scale_factor']
        sf = scalefactor[0]
        solar_zen = dataset[:,:].mean()
        solar_zen *= sf

        length = aod.shape[0]
        width = aod.shape[1]
        pixels = length * width


        ######------ Read Direction & Distance Data ------######
        cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dist_Cloud_Rasters_L2CLD_ALL/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
        cld_dir = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/DATA/Dir_Cloud_Rasters_L2CLD_ALL/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')
        # Throw away possible shadow contaminated pixels
        cld_dist[(cld_dir > 135) & (cld_dist < 40)] = np.nan # 40 pixels is about 2 km





        ######------ Read Cloud Properties Data ------######
        for cloud_file in os.listdir(cloud_dir):
            if cloud_file[-45:-37] == eMAS_ID:
                cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)

        ### COD ###
        DATAFIELD_NAME = 'Cloud_Optical_Thickness'
        data3D = cloud_hdf.select(DATAFIELD_NAME)
        COD = data3D[:,:].astype('float')
        # Handle fill value.
        attrs = data3D.attributes(full=1)
        fillvalue=attrs["_FillValue"]
        fv = fillvalue[0]
        COD[COD==fv] = np.nan
        scalefactor=attrs['scale_factor']
        sf = scalefactor[0]
        COD *= sf

        ### Height ###
        DATAFIELD_NAME = 'Cloud_Top_Height'
        data3D = cloud_hdf.select(DATAFIELD_NAME)
        Height = data3D[:,:].astype('float')
        # Handle fill value.
        attrs = data3D.attributes(full=1)
        fillvalue=attrs["_FillValue"]
        fv = fillvalue[0]
        Height[Height==fv] = np.nan
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
        Phase[(Phase == fv) | (Phase == 0) | (Phase == 1)] = np.nan
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
        Temp[Temp==fv] = np.nan
        scalefactor=attrs['scale_factor']
        sf = scalefactor[0]
        offset=attrs['add_offset']
        ofs = offset[0]            
        Temp -= ofs
        Temp *= sf

        print float(len(aod[np.isfinite(aod)])) / pixels

        #####------------------------------------------------#####
        #####--------------- Whole Strip Stats --------------#####
        #####------------------------------------------------#####

        if float(len(aod[np.isfinite(aod)])) / pixels > 0.1:
            print 'Strip: TRUE'
            cld_dist2 = np.round(cld_dist)
            cld_dist2 = np.ndarray.flatten(cld_dist2)
            aod_flat = np.ndarray.flatten(aod)
            df = pd.DataFrame(data={'aod':aod_flat,'cld_dist':cld_dist2})
            df_g = (pd.groupby(df['aod'],df['cld_dist']))
            df_gm = df_g.mean()
            df_gm = np.array(df_gm)
            
            background_aod = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)])
            near_cloud_aod = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)][::-1]+1)


            #####------------------------------------------------#####
            #####--------------- Fragmented Strips --------------#####
            #####------------------------------------------------#####

            # Fragment strips into square boxes along track (excess track length is not used)
            sections = length/width
            pixels = width * width
            print "Sections:",sections
            print "-----------------"

            for section in range(sections):

                # Fragment individual datasets
                aod_sect = aod[section*width:(section+1)*width,:]

                COD_sect = COD[section*width*10:(section+1)*width*10,:]
                Height_sect = Height[section*width*10:(section+1)*width*10,:]
                Phase_sect = Phase[section*width*10:(section+1)*width*10,:]
                Temp_sect = Temp[section*width*10:(section+1)*width*10,:]

                cld_dist_sect = cld_dist[section*width:(section+1)*width,:]
                '''cld_dir_sect = cld_dir[section*width:(section+1)*width,:]'''

                
                print float(len(aod_sect[np.isfinite(aod_sect)])) / pixels
                if float(len(aod_sect[np.isfinite(aod_sect)])) / pixels > 0.1:
                    print 'Granule: TRUE'

                    # Calculate AOD statistics
                    cld_dist2 = np.round(cld_dist_sect)
                    cld_dist2 = np.ndarray.flatten(cld_dist2)
                    aod_flat = np.ndarray.flatten(aod_sect)
                    df = pd.DataFrame(data={'aod':aod_flat,'cld_dist':cld_dist2})
                    df_g = (pd.groupby(df['aod'],df['cld_dist']))
                    df_gm = df_g.mean()
                    df_gm = np.array(df_gm)

                    background_aod_sect = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)])
                    near_cloud_aod_sect = np.average(df_gm[np.isfinite(df_gm)],weights=df_g.max().index[np.isfinite(df_gm)][::-1]+1)


                    if float(len(COD_sect[np.isfinite(COD_sect)])) / pixels > 0.1:
                        # Calculate Cloud Property Statistics
                        COD_sect_ave = np.nanmean(COD_sect)
                        COD_sect_std = np.nanstd(COD_sect)
                        COD_sect_n = len(COD_sect[np.isfinite(COD_sect)])

                        Height_sect_ave = np.nanmean(Height_sect)
                        Height_sect_std = np.nanstd(Height_sect)
                        Height_sect_n = len(Height_sect[np.isfinite(Height_sect)])

                        Phase_sect_ave = np.nanmean(Phase_sect)
                        Phase_sect_std = np.nanstd(Phase_sect)
                        Phase_sect_n = len(Phase_sect[np.isfinite(Phase_sect)])
                        Phase_sect_mode = mode(Phase_sect.flatten())[0][0]
                        Phase_sect_mode_perc = mode(Phase_sect.flatten())[1][0] / Phase_sect_n

                        Temp_sect_ave = np.nanmean(Temp_sect)
                        Temp_sect_std = np.nanstd(Temp_sect)
                        Temp_sect_n = len(Temp_sect[np.isfinite(Temp_sect)])

                        # Join the attributes together
                        compiled = [eMAS_ID, COD_sect_ave, COD_sect_std, COD_sect_n, Height_sect_ave, Height_sect_std, Height_sect_n, Phase_sect_ave, Phase_sect_std, Phase_sect_n, Phase_sect_mode, Phase_sect_mode_perc, Temp_sect_ave, Temp_sect_std, Temp_sect_n, background_aod, near_cloud_aod, background_aod_sect, near_cloud_aod_sect, solar_zen]
                        # Write to file
                        writer.writerow(compiled)

output.close()
