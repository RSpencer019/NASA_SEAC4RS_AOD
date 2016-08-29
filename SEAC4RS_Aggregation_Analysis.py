import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import numpy as np

data = pd.read_csv('SEAC4RS_Aggregation_Compiled.csv', header=0)

fig = plt.figure()

'''
### Enhancement vs. Background AOD ###

x = data['background_aod_sect']
y = data['near_cloud_aod_sect']-data['background_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl1 = fig.add_subplot(1,2,1)
plt.scatter(x, y)
plt.title('Enhancement vs. Background AOD')
plt.xlabel('Background AOD')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,1.2],p([0,1.2]),'r--')
plt.figtext(0.25,0.95,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')


### Enhancement vs. Solar Zenith ###

x = data['solar_zen']
y = data['near_cloud_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl2 = fig.add_subplot(1,2,2)
plt.scatter(x, y)
plt.title('Enhancement vs. Solar Zenith')
plt.xlabel('Solar Zenith')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,90],p([0,90]),'r--')
plt.figtext(0.7,0.95,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')



'''


### Enhancement vs. Background AOD ###

x = data['COD_sect_ave']
y = data['near_cloud_aod_sect']-data['background_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl1 = fig.add_subplot(2,2,1)
plt.scatter(x, y)
plt.title('Enhancement vs. Cloud Properties')
plt.xlabel('Cloud Properties')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,100],p([0,100]),'r--')
plt.figtext(0.25,0.95,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')


### Enhancement vs. Background AOD ###

x = data['Height_sect_ave']
y = data['near_cloud_aod_sect']-data['background_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl2 = fig.add_subplot(2,2,2)
plt.scatter(x, y)
plt.title('Enhancement vs. Cloud Properties')
plt.xlabel('Cloud Properties')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,100],p([0,100]),'r--')
plt.figtext(0.7,0.95,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')


### Enhancement vs. Background AOD ###

x = data['Phase_sect_ave']
y = data['near_cloud_aod_sect']-data['background_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl3 = fig.add_subplot(2,2,3)
plt.scatter(x, y)
plt.title('Enhancement vs. Cloud Properties')
plt.xlabel('Cloud Properties')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,100],p([0,100]),'r--')
plt.figtext(0.28,0.03,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')


### Enhancement vs. Background AOD ###

x = data['Temp_sect_ave']
y = data['near_cloud_aod_sect']-data['background_aod_sect']

z = np.polyfit(x,y,1)
p = np.poly1d(z)
r2 = r2_score(y, p(x))

pl4 = fig.add_subplot(2,2,4)
plt.scatter(x, y)
plt.title('Enhancement vs. Cloud Properties')
plt.xlabel('Cloud Properties')
plt.ylabel('Enhancement Near Clouds')
plt.grid()
plt.axis('tight')
plt.plot([0,100],p([0,100]),'r--')
plt.figtext(0.7,0.03,r'$y = %.2fx + %.2f$'%(z[0],z[1])+'\n'+r'$r^2=%.3f$'%(r2),size='large')




plt.show()