# SEAC4RS Python Scripts #

SEAC4RS_Aggregation.py
	Input: eMAS AOD, Distance, Direction, eMAS Cloud Product
	Output: Fragments each strip in the SEAC4RS campaign and summarizes the AOD and cloud properties in each

Distance_Direction_to_Cloud_L2CLD_Scaled.py
	Input: Uses the cloud mask from the eMAS cloud data product
	Output: creates the distance and direction rasters for aod retrievals in entire campaign

Distance_Direction_to_Cloud_L2CLD.py
	Output: Uses the cloud mask from the eMAS cloud data product
	Input: creates the distance and direction rasters for aod retrievals in aeronet collocated files

Plotting_Compilation_by_Aeronet_Graphs.py
	Input: Distance & Direction Rasters, eMAS AOD, aeronet collocation file (compiled_cleaned.csv)
	Output: Graph plots of eMAS vs Distance and Direction of Clouds

Plotting_Compilation_by_Aeronet_Strips.py
	Input: Distance & Direction Rasters, eMAS AOD, eMAS Clouds, CPL, aeronet collocation file (compiled_cleaned.csv)
	Output: Plots image strips of eMAS, Distance and Direction of Clouds, Cloud Properties, CPL

CPL_Collocation.py
	Input: eMAS AOD, CPL Matchfiles, CPL
	Output: Provides plots and strips of eMAS vs the CPL backscatter profile.

MODIS_Search.py
	Input: MODIS (hdf)
	Output: List of MODIS files by spatial search only

eMas_Aeronet_Collocation.py
	Input: eMAS (hdf), Aeronet (lev20)
	Output: Compiled.csv, Aeronet_Sites.txt

eMAS_Aeronet_Plotting_Whole.py
	Input: Compiled.csv, eMAS files
	Output: Mapping Plot (PNG) of the US

eMAS_Aeronet_Plotting_Individual.py
	Input: Compiled.csv, eMAS files
	Output: Mapping Plot (PNG) of individual Aeronet sites

eMAS_Aeronet_Plotting_Individual_Merged.py
	Input: Compiled.csv, eMAS files
	Output: Mapping Plot (PNG) of individual eMAS files w/ all Aeronet sites

eMAS_Aeronet_Plotting_Individual_Merged_CPL.py
	Input: Compiled.csv, eMAS files
	Output: Mapping Plot (PNG) of individual eMAS files w/ all Aeronet sites and a CPL Centerline

Aeronet_Regression.py
	Input: Compiled.csv
	Output: Regression Plot (PNG)

EMAS_CPL_4STAR_Plotting.py
	Input: eMAS, CPL, 4STAR, MODIS data. 
	Output: Mapping Plot (PNG) for entire day

eMAS_vs_MODIS_Aeronet.py
	Output: Map - Individual eMAS w/ MODIS and Aeronet collocation

eMAS_vs_MODIS_Collocation.py
	Input: eMAS, MODIS
	Output: eMAS_vs_MODIS_Compiled.csv

eMAS_CPL_Profile.py


Misc:
display_*.py
	reference files from Gala to do plot eMAS rgb and cloud properties
Aeronet_AOD550_Interpolation.py
	Method used to interpolate a 550 nm from Aeronet
File_Reduce.py
HDF4_Data_Extraction.py