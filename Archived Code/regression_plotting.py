# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:06:19 2016

@author: rsspenc3
"""

import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import numpy as np
import pandas as pd

aeronet = pd.DataFrame.from_csv('Aeronet_Sites_OLD.txt',header=1)

z = np.polyfit(aeronet['Latitude(decimal_degrees)'],aeronet['Longitude(decimal_degrees)'],1)
p = np.poly1d(z)
r2 = r2_score(aeronet['Longitude(decimal_degrees)'], p(aeronet['Latitude(decimal_degrees)']))

plt.scatter(aeronet['Latitude(decimal_degrees)'],aeronet['Longitude(decimal_degrees)'])
plt.title('Title')
plt.grid()
plt.axis('tight')
plt.plot([25,55],p([25,55]),'r--')
plt.ylim(ymin=75,ymax=100)
plt.xlim(xmin=25,xmax=55)
plt.figtext(0.2,0.75,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')

