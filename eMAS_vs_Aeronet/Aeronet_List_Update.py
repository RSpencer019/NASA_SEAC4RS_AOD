

import os

for i in os.listdir('/Users/rsspenc3/Downloads/AOT/LEV20/ALL_POINTS/'):
	n = 0
	for j in os.listdir('/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Aeronet_Data/'):
		if j[14:] in i:
			print ('exists', i, j)
			n += 1
	if n == 0:
		os.rename('/Users/rsspenc3/Downloads/AOT/LEV20/ALL_POINTS/'+i,'/Users/rsspenc3/Desktop/SEAC4RS/eMAS_vs_Aeronet/Aeronet_Data/'+i)