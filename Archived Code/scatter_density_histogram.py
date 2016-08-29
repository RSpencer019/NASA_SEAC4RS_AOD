###### Scatter plot with histogram axes ######

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

'''
# the random data
x = np.random.randn(10000)
y = np.random.randn(10000)
'''

# import data

data = pd.read_csv('aggregated_AOD_dist.csv')
x = np.array(data['agg_dist'])
y = np.array(data['agg_aod'])


nullfmt = NullFormatter()         # no labels

# definitions for the axes
left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
bottom_h = left_h = left + width + 0.004

rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.1]
rect_histy = [left_h, bottom, 0.1, height]

# start with a rectangular Figure
plt.figure(1, figsize=(8, 8))

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
axScatter.scatter(x, y, color='navy',s=0.5)

# now determine nice limits by hand:
xbinwidth = 0.1
ybinwidth = 0.1
xlim = 20
ylim = 5

axScatter.set_xlim((0, xlim))
axScatter.set_ylim((-0.05, ylim))

xbins = np.arange(0, xlim + xbinwidth, xbinwidth)
ybins = np.arange(-0.05, ylim + ybinwidth, ybinwidth)

axHistx.hist(x, bins=xbins, color='navy')
axHisty.hist(y, bins=ybins, orientation='horizontal', color='navy')

axHistx.set_xlim(axScatter.get_xlim())
axHisty.set_ylim(axScatter.get_ylim())


###### Scatter plot with density contours ######

H, xedges, yedges = np.histogram2d(x,y)
extent = [xedges[0],xedges[-1],yedges[0],yedges[-1]]
cset1 = axScatter.contour(H,extent=extent,colors='chartreuse',linewidths=2)


plt.show()