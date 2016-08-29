#!/usr/bin/env python


from pyhdf.SD import SD, SDC
from scipy.misc import bytescale

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.colors import ListedColormap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import BoundaryNorm
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.basemap import Basemap

import math


def setup_other_plots2(plot_data, minval, maxval, use_log) :

	Xcoord_use = Xcoord[plot_data > 0.]
	Ycoord_use = Ycoord[plot_data > 0.]
	Zcoord = plot_data[plot_data > 0.]

	x, y = mymap(Xcoord_use, Ycoord_use) # compute map projection coordinates
	z = Zcoord
	
	mybins = 15
	cmap = plt.get_cmap('hot')


	if (use_log == 0) : 
		levels = MaxNLocator(nbins=mybins).tick_values(minval, maxval)
		norm = BoundaryNorm(levels, ncolors=cmap.N)
		cs1 = mymap.scatter(x, y, c=z, s=0.05, cmap = cmap, norm=norm, edgecolor='none') # was 0.05, use 1 for aircraft

	mymap.drawcoastlines()
	mymap.drawmapboundary(fill_color='silver')
	mymap.fillcontinents(color='gainsboro', lake_color='silver', zorder=0)
	mymap.drawparallels(parallels, labels=[1,0,0,1], fontsize = 8)
	mymap.drawmeridians(meridians, labels=[1,0,0,1], fontsize = 8)


		
	return cs1



def setup_other_plots(plot_data, minval, maxval, use_log) :

	Xcoord_use = Xcoord[plot_data == 2]
	Ycoord_use = Ycoord[plot_data == 2]
	x, y = mymap(Xcoord_use, Ycoord_use) # compute map projection coordinates
	cs1 = mymap.scatter(x, y, c='dodgerblue', s=0.075, edgecolor='none') # use 1 for aircraft, 0.05 for polar

	Xcoord_use = Xcoord[plot_data == 3]
	Ycoord_use = Ycoord[plot_data == 3]
	x, y = mymap(Xcoord_use, Ycoord_use) # compute map projection coordinates
	cs2 = mymap.scatter(x, y, c='lightcyan', s=0.075, edgecolor='none') # use 1 for aircraft, 0.05 for polar

	Xcoord_use = Xcoord[plot_data == 4]
	Ycoord_use = Ycoord[plot_data == 4]
	x, y = mymap(Xcoord_use, Ycoord_use) # compute map projection coordinates
	cs3 = mymap.scatter(x, y, c='limegreen', s=0.075, edgecolor='none') # use 1 for aircraft, 0.05 for polar

	return cs1
	
	
def align_plot_other(ax1, text_title):

	mymap.drawcoastlines()
	mymap.drawmapboundary(fill_color='silver')
	mymap.fillcontinents(color='gainsboro', lake_color='silver', zorder=0)
	mymap.drawparallels(parallels, labels=[1,0,0,1], fontsize = 8)
	mymap.drawmeridians(meridians, labels=[1,0,0,1], fontsize = 8)
	plt.title(text_title, fontsize=10)

	box = ax1.get_position()
	ax1.set_position([box.x0*1.05, box.y0, box.width, box.height])

	axColor1 = plt.axes([box.x0*1.05 + box.width*0.9, box.y0+0.2, 0.01, box.height-0.20])
	cmap = ListedColormap(['dodgerblue', 'lightcyan', 'limegreen'])
	bounds = [2,3,4,5]
	norm = BoundaryNorm(bounds, cmap.N)
	cb1 = ColorbarBase(axColor1, cmap=cmap, norm=norm, orientation='vertical')
	cax = cb1.ax
	cax.set_yticklabels(['Liquid', 'Ice', 'Undetermined'])

#---------- Read Data ----------#

import sys

num_args = len(sys.argv)
if (num_args < 4) : 
	print "Insufficient arguments"
	exit()

mod03_name = sys.argv[1]
mod06_name_1263 = sys.argv[2]
output_name = sys.argv[3]
granule_title = sys.argv[4]

print mod03_name
print mod06_name_1263
print output_name


#---------- lat/lon data ----------#
mod03 = SD(mod03_name, SDC.READ)
sds = mod03.select('Latitude')
#sds = mod03.select('PixelLatitude') #EMAS
latitude = sds[:,:]

sds = mod03.select('Longitude')
#sds = mod03.select('PixelLongitude') #EMAS
longitude = sds[:,:]

#---------- science data -----------#
mod06_1263 = SD(mod06_name_1263, SDC.READ)

sds = mod06_1263.select('Cloud_Phase_Optical_Properties')
phase = sds[:,:]

sds = mod06_1263.select('cloud_top_temperature_1km')
#sds = mod06_1263.select('Cloud_Top_Temperature') #EMAS
ctt = (sds[:,:]+15000)*0.01 
ctt[ctt < 100.] = -999.

sds = mod06_1263.select('cloud_top_height_1km')
#sds = mod06_1263.select('Cloud_Top_Height') #EMAS
cth = sds[:,:]


Xcoord = longitude
Ycoord = latitude

ysize = latitude.shape[0]
xsize = latitude.shape[1]

lon0 = longitude[ysize/2, xsize/2]
lat0 = latitude[ysize/2, xsize/2]

# 8 for polar + np.trunc
# 0.5 or 1 for airborne
parallels = np.trunc(np.arange(Ycoord.min(), Ycoord.max(), 8.))
meridians = np.trunc(np.arange(Xcoord.min(), Xcoord.max(), 8.))

mymap = Basemap(projection='ortho', lon_0 = lon0, lat_0 = lat0, 
# 	llcrnrx=-1.6e5,llcrnry=-1.6e5,urcrnrx=1.6e5, urcrnry=1.6e5, resolution='i') # MASTER, EMAS
	llcrnrx=-1.6e6,llcrnry=-1.6e6,urcrnrx=1.6e6, urcrnry=1.6e6)



#---------- Plot Images ----------#
fig = plt.figure(figsize=(14,7))
ax1 = plt.subplot(131)

phase_use = phase

cs1 = setup_other_plots( phase_use, 2, 4, 0) # phase plot
align_plot_other(ax1, 'COP Phase ' + granule_title)

ax1 = plt.subplot(132)
cs1 = setup_other_plots2( ctt, 100, 325, 0) # CTT plot

ax1 = plt.subplot(133)
cs1 = setup_other_plots2( cth, 100, 20000, 0) # CTH plot


plt.savefig(output_name, dpi=300 )


plt.show()


exit()





