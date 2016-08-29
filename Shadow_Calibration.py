#------------------------------------------------------------------------------
# Name:        Shadow Calibration
# Description: NASA Aerosol Project.
#
# Author:      Robert S. Spencer
#
# Created:     08/08/2016
# Python:      2.7
#------------------------------------------------------------------------------


import os
import h5py
from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#### eMAS AOD ####

eMAS_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_10_20130830_2000_2009_20160627_1605.hdf'
eMAS_ID = eMAS_file[-45:-37]
hdf = SD(eMAS_file, SDC.READ)
  
# Read dataset.
DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
data3D = hdf.select(DATAFIELD_NAME)
aod = data3D[:,:].astype(float) * 0.001

# Handle fill value.
attrs = data3D.attributes(full=1)
fillvalue=attrs["Fill_Val"]
fv = fillvalue[0] * 0.001
aod[aod == fv] = np.nan

# Read Datetime
hour_start = eMAS_file[-27:-25]
minute_start = eMAS_file[-25:-23]
hour_end = eMAS_file[-22:-20]
minute_end = eMAS_file[-20:-18]
time_start = float(hour_start) + (float(minute_start) / 60)
time_end = float(hour_end) + (float(minute_end) / 60)

# Cast into original pixel resolution
aod_highres = np.empty([aod.shape[0]*10,aod.shape[1]*10])
for i in range(len(aod_highres)):
    for j in range(len(aod_highres[i])):
        aod_highres[i,j] = aod[i/10,j/10]

# Solar Zenith
stats = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_Statistics.csv')
sol_zen = float(stats.solar_zenith.loc[stats['eMAS_file'] == eMAS_file])



#### eMAS RGB ####



#####-----------------------------------#####
#####--------------- RGB ---------------#####
#####-----------------------------------#####
    
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




#####------------------------------------------------#####
#####--------------- Cloud Properties ---------------#####
#####------------------------------------------------#####

cloud_dir = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Clouds/'
for cloud_file in os.listdir(cloud_dir):
    if cloud_file[-45:-37] == eMAS_ID:
        cloud_hdf = SD(cloud_dir+cloud_file, SDC.READ)

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

        ### COD ###
        DATAFIELD_NAME = 'Cloud_Optical_Thickness'
        data3D = cloud_hdf.select(DATAFIELD_NAME)
        Clouds = data3D[:,:].astype('float')
        # Handle fill value.
        attrs = data3D.attributes(full=1)
        fillvalue=attrs["_FillValue"]
        fv = fillvalue[0]
        Clouds = np.ma.masked_array(Clouds, Clouds == fv)
        scalefactor=attrs['scale_factor']
        sf = scalefactor[0]
        Clouds *= sf




#### Plotting ####

fig = plt.figure(figsize=(10,16))

# RGB
pl1 = fig.add_subplot(4,1,1)
pl1.axes.get_xaxis().set_visible(False)
pl1.axes.get_yaxis().set_visible(False)  
pl1.imshow(np.transpose(rgb,(1,0,2)))
plt.xlim(xmin=0,xmax=aod_highres.shape[0])
plt.ylim(ymin=aod_highres.shape[1],ymax=0)

# AOD
pl2 = fig.add_subplot(4,1,2)
pl2.axes.get_xaxis().set_visible(False)
pl2.axes.get_yaxis().set_visible(False)  
pl2.imshow(np.transpose(rgb,(1,0,2)))
c2 = plt.imshow(aod_highres.T,vmin=0,vmax=0.75)
#plt.colorbar(c2)
plt.xlim(xmin=0,xmax=aod_highres.shape[0])
plt.ylim(ymin=aod_highres.shape[1],ymax=0)

# Plot Cloud Height
pl3 = fig.add_subplot(4,1,3)
pl3.axes.get_xaxis().set_visible(False)
pl3.axes.get_yaxis().set_visible(False)  
pl3.set_axis_bgcolor('darkolivegreen')
res = 1   
c3 = pl3.imshow(Height.T,cmap=plt.get_cmap('pink_r'))
#fig.colorbar(c3)

# Plot Cloud Depth
pl4 = fig.add_subplot(4,1,4)
pl4.axes.get_xaxis().set_visible(False)
pl4.axes.get_yaxis().set_visible(False)  
pl4.set_axis_bgcolor('darkolivegreen')
res = 1   
c4 = pl4.imshow(Clouds.T,cmap=plt.get_cmap('pink_r'))
fig.colorbar(c4)

plt.show()