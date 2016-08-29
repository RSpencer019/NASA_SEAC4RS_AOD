# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 16:36:24 2016

@author: rsspenc3
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

fig.clear()

aeronet_data = pd.read_csv('/Users/rsspenc3/Desktop/SEAC4RS/DATA/Aeronet_Data/920801_160611_IMPROVE-MammothCave.lev20',header=4)
aeronet_col = list(aeronet_data.columns[3:19])
aeronet_wl = []
for item in aeronet_col:
    aeronet_wl.append(int(item[4:]))
aeronet_aod = []    
for item in aeronet_col:
    aeronet_aod.append(aeronet_data[item][0])


aeronet_aod = np.array(np.log10(aeronet_aod))
aeronet_wl = np.array(np.log10(aeronet_wl))
idx = np.isfinite(aeronet_aod)

z = np.polyfit(aeronet_wl[idx],aeronet_aod[idx],2)
p = np.poly1d(z)
r2 = r2_score(aeronet_aod[idx], p(aeronet_wl[idx]))

wl = np.linspace(np.log10(300),np.log10(1700),100)

plt.scatter(aeronet_wl[idx],aeronet_aod[idx])
plt.title('log(AOD) vs. log(lambda)')
plt.grid()
plt.axis('tight')
plt.plot(wl,p(wl),'r--')
plt.figtext(0.4,0.75,'Polynomial Coefficients: %.2f, %.2f, %.2f'%(z[0],z[1],z[2])+'\n'+r'$r^2=%.3f$'%(r2),size='large')


print 10**(p(np.log10(550)))
plt.scatter(np.log10(550),p(np.log10(550)),c='green',s=200)
