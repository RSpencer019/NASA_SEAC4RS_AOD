#!/usr/bin/env python

'''
Copyright (c) 2013 Marchant Benjamin
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap


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
		

# scale_image function References:
#    - http://www.idlcoyote.com/ip_tips/brightmodis.html) 
#    - http://code.google.com/p/ccnworks/source/browse/trunk/modis/true.py) 

#---------- Read Data ----------#

myd021km_name = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/Aqua/MYD021KM.A2013218.1635.006.2013219181718.hdf'
myd03_name = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/RGB/AquaGeo/MYD03.A2013218.1635.006.2013219155551.hdf'
output_name = '/Users/rsspenc3/Desktop/rgb.png'

print myd021km_name
print myd03_name
print output_name

myd021km = SD(myd021km_name, SDC.READ)
myd03 = SD(myd03_name, SDC.READ)

ht = 1015
wid = 677


#---------- Reflectances (Bands 1, 4 and 3) from MODIS ----------#

sds = myd021km.select('EV_250_Aggr1km_RefSB')
Band_01 = sds.get(start=(0,0,0), count=(1, ht, wid), stride =(1,2,2)) *5.31009e-5

sds = myd021km.select('EV_500_Aggr1km_RefSB')
Band_04 = sds.get(start=(1,0,0), count=(1, ht, wid), stride =(1,2,2)) *3.42507e-5

sds = myd021km.select('EV_500_Aggr1km_RefSB')
Band_03 = sds.get(start=(0,0,0), count=(1, ht, wid), stride =(1,2,2)) *3.67092e-5


sds2 = myd03.select('SolarZenith')
miu0_read = sds2.get(start=(0,0), count=(ht, wid), stride =(2,2)) *0.01
miu0 = np.cos(miu0_read/180.*np.pi)

sds = myd03.select('Latitude')
latitude = sds.get(start=(0,0), count=(ht, wid), stride =(2,2))

sds = myd03.select('Longitude')
longitude = sds.get(start=(0,0), count=(ht, wid), stride =(2,2))



#---------- RGB Matrix ----------#
rgb_ModisScale = np.zeros((ht, wid,3), dtype=np.uint8)

scale = np.zeros((ht,wid), dtype=np.int)
scale[:,:] = Band_01[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
scale[scale > MAXVALUES] = MAXVALUES-1
scale[scale < 0] = 0
rgb_ModisScale[:,:,0] = band_scale[scale[:,:]]

scale = np.zeros((ht,wid), dtype=np.int)
scale[:,:] = Band_04[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
scale[scale > MAXVALUES] = MAXVALUES-1
scale[scale < 0] = 0
rgb_ModisScale[:,:,1] = band_scale[scale[:,:]]

scale = np.zeros((ht,wid), dtype=np.int)
scale[:,:] = Band_03[0, 0:ht, 0:wid]/miu0[0:ht,0:wid]*10000. + (MAXVALUES/2) + 0.5
scale[scale > MAXVALUES] = MAXVALUES-1
scale[scale < 0] = 0
rgb_ModisScale[:,:,2] = band_scale[scale[:,:]]

rgb_ModisScale[rgb_ModisScale > 255] = 255
rgb_ModisScale[rgb_ModisScale < 0] = 0


Xcoord = longitude
Ycoord = latitude

ysize = latitude.shape[0]
xsize = latitude.shape[1]

# set these to 0 for SEVIRI
lon0 = longitude[ysize/2, xsize/2]
lat0 = latitude[ysize/2, xsize/2]

print lon0, lat0
print ysize, xsize

# 8 for polar + np.trunc
# 0.5 or 1 for airborne

parallels = np.trunc(np.arange(Ycoord.min(), Ycoord.max(), 8)) # was 8 for both
meridians = np.trunc(np.arange(Xcoord.min(), Xcoord.max(), 8))

print parallels
print meridians

mymap = Basemap(projection='ortho', lon_0 = lon0, lat_0 = lat0, #) 
# 	llcrnrx=-1.6e5,llcrnry=-1.6e5,urcrnrx=1.6e5, urcrnry=1.6e5, resolution='i') # MASTER, EMAS
 	llcrnrx=-1.6e6,llcrnry=-1.6e6,urcrnrx=1.6e6, urcrnry=1.6e6) # for polar
# nothing for geostationary, do not set vars

x, y = mymap(Xcoord, Ycoord) # compute map projection coordinates

r = rgb_ModisScale[:,:, 0]
g = rgb_ModisScale[:,:, 1]
b = rgb_ModisScale[:,:, 2]

rgb = np.array([r,g,b]).T
color_tuple = rgb.transpose((1,0,2)).reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0

#---------- Plot Image ----------#
fig = plt.figure(figsize=(7, 7))
ax1 = plt.subplot(111)

m = mymap.scatter(x,y, c=color_tuple, s=3.0, edgecolor='none')

mymap.drawcoastlines()
mymap.drawmapboundary(fill_color='silver')
mymap.fillcontinents(color='gainsboro', lake_color='silver', zorder=0)
mymap.drawparallels(parallels, labels=[1,0,0,1])
mymap.drawmeridians(meridians, labels=[1,0,0,1])


plt.savefig(output_name, dpi=200 )


#print "display 3"

plt.show()


