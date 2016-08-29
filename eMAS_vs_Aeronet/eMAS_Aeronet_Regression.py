## Linear Trendline Plotting ##

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap


plt.clf()

compiled = pd.read_csv('/Volumes/NASA_Spencer/SEAC4RS/eMAS_vs_Aeronet/Compiled_3K.csv',header=0)

'''
plt.subplot(2,1,1)
m = Basemap(projection='cyl', resolution='l',
            llcrnrlat=10, urcrnrlat = 55,
            llcrnrlon=-130, urcrnrlon = -60)
m.drawcoastlines(linewidth=0.75)
m.drawcountries(linewidth=0.25)
m.drawparallels(np.arange(-90., 90., 15.), labels=[1, 0, 0, 0])
m.drawmeridians(np.arange(-180., 180., 15.), labels=[0, 0, 0, 1])
plt.scatter(compiled['longitude'],compiled['latitude'],s=80)
plt.title('Aeronet Sites')
'''


#### by Collocation ####

plt.subplot(2,1,1)
x = compiled['meanval_aeronet_550_intrp'].replace(-9999.0,np.nan)
y = compiled['meanval_eMAS_550'].replace(-9999.0,np.nan)
z = compiled['nval_eMAS'].replace(-9999.0,np.nan) #* compiled['nval_aeronet'].replace(-9999.0,np.nan)
xerr = compiled['std_aeronet']
yerr = compiled['std_eMAS']
print ( max(z), np.mean(z), min(z) )

#plt.scatter(x,y,s=4*z)
plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o', elinewidth = 0.25)
plt.title('eMAS vs. Aeronet (AOD @ 550nm)')
plt.xlabel('Aeronet')
plt.ylabel('eMAS')
plt.plot([0,1],[0,1],'k-')
plt.ylim([0,1])
plt.xlim((0,1))

print (min(z), max(z))
#plt.plot(sorted(z))
z_lower_limit = 10
idx = np.isfinite(x) & np.isfinite(y) & (x<1.0) & (z>z_lower_limit)
pf = np.polyfit(x[idx], y[idx], 1)  # w = z[idx]
p = np.poly1d(pf)
r2 = r2_score(y[idx], p(x[idx]))

plt.plot([0,1],p([0,1]),'k--',linewidth=3)
plt.figtext(0.55,0.6,r'$n > %i$'%(z_lower_limit)+'\n'+r'$y = %.2fx+%.2f$'%(pf[0],pf[1])+'\n'+r'$r^2=%.3f$'%(r2),size='x-large')
#plt.figtext(0.2,0.45,r'$z > %i$'%(200), size='large')
fig = plt.gcf()

print (len(z), len(z[np.isfinite(x) & np.isfinite(y)]),len(z[np.isfinite(x) & np.isfinite(y) & (z>z_lower_limit)]))

# highlight the regression points
index = (z > z_lower_limit)
#plt.scatter(x[index],y[index],s=60,edgecolor='r',linewidths=3)  #s=z[index]*4




#### by Background ####

compiled = pd.read_csv('/Volumes/NASA_Spencer/SEAC4RS/Background_AOD.csv',header=0)

plt.subplot(2,1,2)
x = compiled['Ave_Aeronet_AOD'].replace(-9999.0,np.nan)
y = compiled['Background_AOD'].replace(-9999.0,np.nan)

plt.scatter(x,y,s=50,c='b')
plt.title('"Background" eMAS vs. Aeronet (AOD @ 550nm)')
plt.xlabel('Aeronet')
plt.ylabel('eMAS')
plt.plot([0,1],[0,1],'k-')
plt.ylim([0,1])
plt.xlim((0,1))

pf = np.polyfit(x, y, 1) 
p = np.poly1d(pf)
r2 = r2_score(y, p(x))

plt.plot([0,1],p([0,1]),'r--',linewidth=4)
plt.figtext(0.55,0.15,r'$y = %.2fx+%.2f$'%(pf[0],pf[1])+'\n'+r'$r^2=%.3f$'%(r2),size='x-large')
fig = plt.gcf()

# highlight the regression points

'''
pngfile = "Plots_Overall/Regression_Error.png"
fig.savefig(pngfile,dpi=300)
'''



