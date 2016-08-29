#------------------------------------------------------------------------------
# Name:        eMAS vs CPL Collocations
# Description: NASA Aerosol Project. Uses indices of premade matchfiles (Paolo Veglio)
#
# Author:      Robert S. Spencer
#
# Created:     08/05/2016
# Python:      2.7
#------------------------------------------------------------------------------


import os
import h5py
from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt


'''
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
'''


aeronet = 0.3

# eMASL2AER_13959_03_20130830_1759_1827_20160627_1603.hdf

#### MATCHFILE ####

FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Match_Examples/Wisconsin_Match_Example/Match_CPL.EMAS.COLLOCATION.20130830.195947.hdf'
hdf = SD(FILE_NAME, SDC.READ)
  
# Read CPL Index dataset.
DATAFIELD_NAME='master_idx'
data = hdf.select(DATAFIELD_NAME)
CPL_idx = data[0,:]

# Read eMAS Along Index dataset.
DATAFIELD_NAME='follower_idx_alongtrack'
data = hdf.select(DATAFIELD_NAME)
eMAS_along_idx = data[:,:]

# Read eMAS Across Index dataset.
DATAFIELD_NAME='follower_idx_acrosstrack'
data = hdf.select(DATAFIELD_NAME)
eMAS_across_idx = data[:,:]

# Read eMAS Across Number dataset.
DATAFIELD_NAME='eMAS_Num_FOV'
data = hdf.select(DATAFIELD_NAME)
eMAS_across_number = data[0,:]

# Read Datetime
hour_start = FILE_NAME[-10:-8]
minute_start = FILE_NAME[-8:-6]
second_start = FILE_NAME[-6:-4]
time_start = float(hour_start) + (   (float(minute_start) + round(float(second_start) / 60) ) / 60   )



#### eMAS AOD ####

eMAS_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_10_20130830_2000_2009_20160627_1605.hdf'
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

# Initiate the new AOD collocation DS
new_aod_length = len(eMAS_across_number)
new_aod_width = int(max(eMAS_across_number))
new_aod = np.empty([new_aod_length, new_aod_width])
new_aod.fill(np.nan)

# Cast values into new AOD DS
for n in range(new_aod_length):
	for j in range(int(eMAS_across_number[n])):
		if eMAS_along_idx[n,j] > aod_highres.shape[0]:
			break
		new_aod[n,j] = aod_highres[int(eMAS_along_idx[n,j])-1, int(eMAS_across_idx[n,j])-1]


# Reduce casted array into single dimension
AOD_centerline = np.nanmean(new_aod, axis=1)





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








#### Cloud CPL ####

# Open CPL file
FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/CPL_ATB_13959_30aug13.hdf5'
hdf5 = h5py.File(FILE_NAME,'r')

# Read Ground Height dataset.
Gnd_Hgt = hdf5['Gnd_Hgt'][:]
Gnd_Hgt = Gnd_Hgt[CPL_idx.astype(int)]
Bin_Alt = hdf5['Bin_Alt']
Gnd_Ind = np.empty(shape=Gnd_Hgt.shape)
for vert_column in range(len(Gnd_Hgt)):
	Gnd_Ind[vert_column] = np.argmin(np.abs(Bin_Alt - Gnd_Hgt[vert_column]))



# CPL Data
clouds = hdf5['Layer_Type'][:,:]
cloud_layers = clouds[CPL_idx.astype(int)]
clouds = cloud_layers.max(axis=1)
ATB = hdf5['ATB_532'][:,:]
ATB_profile = ATB[CPL_idx.astype(int)]

ATB_profile_norm = np.copy(ATB_profile)
ATB_profile_norm[ATB_profile_norm > 0.012] = np.nan
ATB_sum_norm = np.nansum(ATB_profile_norm, axis=1)

ATB_sum = np.nansum(ATB_profile[:,:730], axis=1)
ATB_sum_cleaned = np.copy(ATB_sum)
ATB_sum_cleaned[clouds==3] = np.nan


ATB_depth = np.empty(shape=ATB_sum.shape)
for vertical_column in range(len(ATB_profile)):
    total = 0
    for del_z in range(len(ATB_profile[vertical_column])):
        total += ATB_profile[vertical_column,del_z]
        if total > ATB_sum[vertical_column] * 0.9:
            ATB_depth[vertical_column] = del_z
            break



#### Plotting ####

fig = plt.figure(figsize=(10,16))

# RGB
pl1 = fig.add_subplot(5,1,1)
pl1.axes.get_xaxis().set_visible(False)
pl1.axes.get_yaxis().set_visible(False)  
pl1.imshow(np.transpose(rgb,(1,0,2)))
cpl_x = np.linspace(1,len(aod_highres),len(clouds))
cpl_y = np.array([len(aod_highres[0]) / 2] * len(cpl_x))
pl1.scatter(cpl_x[clouds==3], cpl_y[clouds==3]+20,linewidths=0,s=10.0,c='k',marker='s')
pl1.scatter(cpl_x[clouds==3], cpl_y[clouds==3]-20,linewidths=0,s=10.0,c='k',marker='s')
plt.xlim(xmin=0,xmax=aod_highres.shape[0])
plt.ylim(ymin=aod_highres.shape[1],ymax=0)

# AOD
pl2 = fig.add_subplot(5,1,2)
pl2.axes.get_xaxis().set_visible(False)
pl2.axes.get_yaxis().set_visible(False)  
pl2.imshow(np.transpose(rgb,(1,0,2)))
c2 = plt.imshow(aod_highres.T,vmin=0,vmax=0.75)
#plt.colorbar(c2)
cpl_x = np.linspace(1,len(aod_highres),len(clouds))
cpl_y = np.array([len(aod_highres[0]) / 2] * len(cpl_x))
pl2.scatter(cpl_x[clouds==3], cpl_y[clouds==3]+20,linewidths=0,s=10.0,c='k',marker='s')
pl2.scatter(cpl_x[clouds==3], cpl_y[clouds==3]-20,linewidths=0,s=10.0,c='k',marker='s')
plt.xlim(xmin=0,xmax=aod_highres.shape[0])
plt.ylim(ymin=aod_highres.shape[1],ymax=0)

# CPL Layers
pl3 = fig.add_subplot(5,1,3)
c3 = pl3.imshow(cloud_layers.T,aspect='auto',interpolation='none')
pl3.axes.get_xaxis().set_visible(False)
pl3.axes.get_yaxis().set_visible(False)  
#plt.colorbar(c3)
'''
# Graphs
pl3 = fig.add_subplot(5,2,5)
pl3.scatter(ATB_sum_cleaned,AOD_centerline)
'''
# CPL ATB Profile
pl4 = fig.add_subplot(5,1,4)
pl4.axes.get_xaxis().set_visible(False)
pl4.axes.get_yaxis().set_visible(False)  
c4 = plt.imshow(ATB_profile.T,vmin=0,vmax=0.014)
plt.xlim(xmin=0,xmax=ATB_profile.shape[0])
plt.ylim(ymin=750,ymax=600)
#plt.colorbar(c4)
plt.plot(Gnd_Ind,color='w')

# AOD CPL Comparison
pl5 = fig.add_subplot(5,1,5)
pl5.axes.get_xaxis().set_visible(False)
pl5.axes.get_yaxis().set_visible(False)  
pl5.fill_between(x=range(len(ATB_sum)),y1=0,y2=ATB_sum,alpha=0.1,color='r')
pl5.fill_between(x=range(len(ATB_sum)),y1=0,y2=ATB_sum_cleaned,alpha=0.6,color='r')
#pl5.fill_between(x=range(len(ATB_depth)),y1=0,y2=(900-ATB_depth)/100 - 1,alpha=0.6,color='r')
pl5.fill_between(x=range(len(AOD_centerline)),y1=0,y2=AOD_centerline,alpha=0.6,color='b')
pl5.bar(left=range(len(clouds)),height=clouds.astype(float)/5-0.4,color='k')
plt.axhline(y=aeronet,ls='--',lw=2,c='k')
plt.xlim(xmin=0,xmax=len(clouds))

fig.tight_layout()
plt.show()




