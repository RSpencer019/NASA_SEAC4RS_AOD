from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import pandas as pd
import math
import numpy as np
import heapq


#cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/Dist_Cloud_Rasters/{0}.csv".format(eMAS_file[-55:-4]), delimiter=',')

cld_dist = np.genfromtxt("/Users/rsspenc3/Desktop/SEAC4RS/Dist_Cloud_Rasters/eMASL2AER_13958_19_20130827_2112_2133_20160627_1558.csv", delimiter=',')

# AOD
eMAS_file = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_27Aug/eMASL2AER_13958_19_20130827_2112_2133_20160627_1558.hdf'
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



def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

   
cld_dist = np.round(rebin(cld_dist,[cld_dist.shape[0]/10,cld_dist.shape[1]/10]))

cld_dist = np.ndarray.flatten(cld_dist)
aod = np.ndarray.flatten(aod)

plt.scatter(cld_dist,aod,c='b')

df = pd.DataFrame(data={'aod':aod,'cld_dist':cld_dist})

df_g = df.groupby('cld_dist')
df_i = []
for i in df_g.indices:
    df_i.append(i)

df_g = (pd.groupby(df['aod'],df['cld_dist']))

df_gl = df_g.nlargest(10)
df_v = []
for i in df_i:
    df_v.append(df_gl[i].mean())

x = np.array(df_i)
y = np.array(df_v)

background_aod = np.average(df_g.mean(),weights=df_g.mean().index)
z = np.polyfit(x,1/y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))
peak_aod = 1/p(0)
range_aod = peak_aod - background_aod

plt.scatter(x,y,c='r')
plt.plot([0,x.max()],[background_aod,background_aod],c='b')
plt.plot([0,x.max()],[peak_aod,peak_aod],c='b')
plt.title('Title')
plt.grid()
plt.axis('tight')
plt.plot(x,1/p(x),'r--')

