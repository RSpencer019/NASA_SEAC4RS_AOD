import h5py
from pyhdf.SD import SD, SDC
import numpy as np
import matplotlib.pyplot as plt



# eMASL2AER_13959_03_20130830_1759_1827_20160627_1603.hdf

#### eMAS ####

FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/eMAS_Data/Aug/Emas_30Aug/eMASL2AER_13959_03_20130830_1759_1827_20160627_1603.hdf'

hdf = SD(FILE_NAME, SDC.READ)
  
# Read dataset.
DATAFIELD_NAME='Image_Optical_Depth_Land_And_Ocean'
data3D = hdf.select(DATAFIELD_NAME)
eMAS = data3D[:,35:36].astype(float) * 0.001

# Handle fill value.
attrs = data3D.attributes(full=1)
fillvalue=attrs["Fill_Val"]

# fillvalue[0] is the attribute value.
fv = fillvalue[0] * 0.001

#data[data == fv] = np.nan
eMAS = np.ma.masked_array(eMAS, eMAS == fv)

# Plot
x = np.linspace(1,len(eMAS),len(eMAS))

x3 = np.linspace(1,ATB.shape[1],len(eMAS))

#fig.clear()
plt.plot(eMAS)


hour_start = FILE_NAME[-27:-25]
minute_start = FILE_NAME[-25:-23]
hour_end = FILE_NAME[-22:-20]
minute_end = FILE_NAME[-20:-18]
time_start = float(hour_start) + (float(minute_start) / 60)
time_end = float(hour_end) + (float(minute_end) / 60)





#### Cloud CPL ####

# Open CPL file
FILE_NAME = '/Users/rsspenc3/Desktop/SEAC4RS/DATA/CPL_Data/CPL_ATB_13959_30aug13.hdf5'
hdf5 = h5py.File(FILE_NAME,'r')

hour = hdf5['Hour'][:].astype(float)
minute = hdf5['Minute'][:].astype(float)
time = hour + (minute / 60)
clouds = hdf5['Layer_Type'][:,0]
CPL = np.array([time, clouds])
CPL2 = CPL[:,(CPL[0] > time_start) & (CPL[0] < time_end)]

ATB = hdf5['ATB_532'][:,:]
ATB = ATB[(time > time_start) & (time < time_end),:].T

c = plt.imshow(ATB,vmin=0,vmax=0.014)
plt.colorbar(c)
plt.axhline(y=200,lw=4,color='r')

# Plot Clouds
x2 = np.linspace(1,len(eMAS),len(CPL2[0]))
plt.fill(x2, CPL2[1]-2.75)
plt.ylim(0,2)

