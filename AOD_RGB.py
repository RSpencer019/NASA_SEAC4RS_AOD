import numpy as np
import os
from pyhdf.SD import SD, SDC
import matplotlib.pyplot as plt

eMAS_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Sept/Emas_11Sept/eMASL2AER_13964_13_20130911_1949_2002_20160628_0952.hdf'


# AOD scale
lowerval = 0.0
upperval = 0.75


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

# Plot RGB eMAS
fig = plt.figure()    
pl1 = fig.add_subplot(111)
pl1.axes.get_xaxis().set_visible(False)
pl1.axes.get_yaxis().set_visible(False)  
pl1.imshow(rgb)












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
aod_highres = np.empty([aod.shape[0]*10,aod.shape[1]*10])
for i in range(len(aod_highres)):
    for j in range(len(aod_highres[i])):
        aod_highres[i,j] = aod[i/10,j/10]

c7 = pl1.imshow(aod_highres,vmin=lowerval,vmax=upperval,aspect='equal')
#fig.colorbar(c7,ticks=[0,0.25,0.5,0.75])
pl1.set_xlim([0,aod_highres.shape[1]])
pl1.set_ylim([aod_highres.shape[0],0])
pl1.set_title('AOD',fontsize='medium')



pngfile = "{0}_AOD_RGB.png".format(eMAS_file[-55:-4])
fig.savefig(pngfile,dpi=1000)
plt.close()
