# -*- coding: utf-8 -*-
"""


Data - model comparison:
    
Comparing modeled with measured d18Oc values, using species specific settings

@author: Anne
"""

import netCDF4
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os
import math

# for statistics
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import ttest_ind

# for seasons
from datetime import timedelta
from datetime import date, datetime
from operator import attrgetter
import pandas as pd
from netCDF4 import num2date, date2num

# for plotting
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% Functions 

class ParticleData(object):
    def __init__(self):
        self.id = []

def load_particles_file(fname, varnames):
        T = ParticleData()
        pfile = Dataset(fname, 'r')
        T.id = pfile.variables['trajectory'][:]
        for v in varnames:
            setattr(T, v, pfile.variables[v][:])
        return T

## ---------Temperature analysis
        
def find_dwell(DwellTemp, Depth):
    obsmax = int(len(DwellTemp[0])) #(number of observations)
    trajmax = int(len(DwellTemp[:,0]))    # (number of particles and their trajectories)

    for i in range(trajmax):    
        for j in range(obsmax):
            if Depth[i,j] > 501: #301: # 51: # 101: # 
                DwellTemp[i,j]= np.nan
                
    AvTemp = np.nanmean(DwellTemp, axis= 1)
    Stdev = np.nanstd(DwellTemp, axis=1)
    return DwellTemp, AvTemp, obsmax, trajmax

## --------Final locations and value of variable at final location
    
def find_final(trajmax, variable, trajlon, trajlat):
    Final_variable= np.zeros(trajmax)
    FinalLon  = np.zeros(trajmax)
    FinalLat  = np.zeros(trajmax)

    for i in range(trajmax):
        new_mask = (np.isnan(variable[i,:]))           # making vector with masked elements for each non-nan value
        if len(variable[i,:][~new_mask]) == 0:
            Final_variable[i] == 0
        else:
            Final_variable[i] = variable[i,:][~new_mask][-1]    # applying the masked vector to the row, making each nan false
            FinalLon[i]  = trajlon[i,:][~new_mask][-1]
            FinalLat[i]  = trajlat[i,:][~new_mask][-1]
    
    Final_variable[Final_variable == 0]= 'nan'                # getting rid of the 0 values (from the initial filling of 'empty 'array')
    FinalLon[FinalLon == 0]= 'nan'                 
    FinalLat[FinalLat == 0]= 'nan'    
    return Final_variable, FinalLon, FinalLat



def makeKDE_array(values):
    def create_filterdata(data):
        data= ma.getdata(data)
        filterdata = data[~np.isnan(data)]
        return filterdata
    values= create_filterdata(values)
    kernel = stats.gaussian_kde(values)
    xx = np.linspace(-2, 3, 10000)
    x = kernel(xx)
    return x

def temp_grouping(Av_variable, trajmax):
    Temp_group= np.zeros(trajmax)
    Warm_species= []
    Trans_species= []
    Cold_species = []
    Warmquantile = []
    Coldquantile = []
    Tempcount=np.zeros(3)
    quantiles = np.nanquantile(Av_variable, [0.10, 0.25, 0.5, 0.75, 0.90])
    for i in range(trajmax):
        if Av_variable[i] >= 22:
            Temp_group[i] = 3
            Warm_species.append(Av_variable[i])
            #Tempcount[0]=Tempcount[0]+1
        if Av_variable[i] >= quantiles[4]:
            Warmquantile.append(Av_variable[i])
            Tempcount[0]=Tempcount[0]+1
        if Av_variable[i] <22 and Av_variable[i] >= 12:
            Temp_group[i] = 2
            Trans_species.append(Av_variable[i])
        if Av_variable[i] < quantiles[4] and Av_variable[i] >= quantiles[0]:
            Tempcount[1]=Tempcount[1]+1
        if Av_variable[i] <12:
            Temp_group[i] = 1
            Cold_species.append(Av_variable[i])
        if Av_variable[i] < quantiles[0]:
            Coldquantile.append(Av_variable[i])
            Tempcount[2]=Tempcount[2]+1
        if Av_variable[i] == 'nan':
            Temp_group[i] = 0
    return Temp_group, Warm_species, Trans_species, Cold_species, Warmquantile, Coldquantile, Tempcount




#%%  Imporing data
file_d18OR = "C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/d18O/Ruber_d18O_T_traj_and_local_dd50_sp200_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc"
file_d18OP = "C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/d18O/Pachy_d18O_T_traj_and_local_dd100_sp250_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc"
file_d18OI = "C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/d18O/Inflata_d18O_T_traj_and_local_dd300_sp250_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc"
file_d18OIdeep = "C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/d18O/InflataDeep_d18O_T_traj_and_local_dd500_sp250_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc"

#file_d18OG = 'C:/Users/Anne/Desktop/Datasets_thesis/d18o/d18O_T_traj_and_local_dd50_sp200_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc'


Tr = load_particles_file(file_d18OR,['time','lon', 'lat', 'age', 'temp', 'z', 'd18oc', 'd18ow', 'loctemp'])    
Tp = load_particles_file(file_d18OP,['time','lon', 'lat', 'age', 'temp', 'z', 'd18oc', 'd18ow', 'loctemp'])
Ti = load_particles_file(file_d18OI,['time','lon', 'lat', 'age', 'temp', 'z', 'd18oc', 'd18ow', 'loctemp'])
Tid = load_particles_file(file_d18OIdeep,['time','lon', 'lat', 'age', 'temp', 'z', 'd18oc', 'd18ow', 'loctemp'])

#%% Using functions
# make sure to set the right depth in the find_dwell function

# For Ruber: (Depth =51)
Dwelld18O_arrR, Avd18O_arrR, obsmax, trajmax, = find_dwell(Tr.d18ow, Tr.z)
Dwelld18Oc_arrR, Avd18Oc_arrR, obsmax, trajmax, = find_dwell(Tr.d18oc, Tr.z)
indi = np.where(Tr.time.mask == True)
Dwelld18O_arrR[indi]= "nan"

Final_d18OR, FinalLonR, FinalLatR = find_final(trajmax, Dwelld18O_arrR, Tr.lon, Tr.lat)

# For Pachyderma: (Depth = 101)
Dwelld18O_arrP, Avd18O_arrP, obsmax, trajmax, = find_dwell(Tp.d18ow, Tp.z)
Dwelld18Oc_arrP, Avd18Oc_arrP, obsmax, trajmax, = find_dwell(Tp.d18oc, Tp.z)
indi = np.where(Tp.time.mask == True)
Dwelld18O_arrP[indi]= "nan"

Final_d18OP, FinalLonP, FinalLatP = find_final(trajmax, Dwelld18O_arrP, Tp.lon, Tp.lat)


# For Inflata: (Depth = 301)
Dwelld18O_arrI, Avd18O_arrI, obsmax, trajmax, = find_dwell(Ti.d18ow, Ti.z)
Dwelld18Oc_arrI, Avd18Oc_arrI, obsmax, trajmax, = find_dwell(Ti.d18oc, Ti.z)
indi = np.where(Ti.time.mask == True)
Dwelld18O_arrI[indi]= "nan"

Final_d18OI, FinalLonI, FinalLatI = find_final(trajmax, Dwelld18O_arrI, Ti.lon, Ti.lat)

# For Inflata deep : (Depth = 501)
Dwelld18O_arrId, Avd18O_arrId, obsmax, trajmax, = find_dwell(Tid.d18ow, Tid.z)
Dwelld18Oc_arrId, Avd18Oc_arrId, obsmax, trajmax, = find_dwell(Tid.d18oc, Tid.z)
indi = np.where(Tid.time.mask == True)
Dwelld18O_arrId[indi]= "nan"

Final_d18OId, FinalLonId, FinalLatId = find_final(trajmax, Dwelld18O_arrId, Tid.lon, Tid.lat)



#%% 
path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/"  
a = [-0.59, -0.39, -0.57, -0.38, 0.09] # measurements on G. Ruber

x = [2.39, 2.28, 1.8, 1.26, 1.23] # measurements on G. Inflata


fig = plt.figure()
kde_data1 = makeKDE_array(Avd18Oc_arrR)
kde_data2 = makeKDE_array(Avd18Oc_arrP)
kde_data3 = makeKDE_array(Avd18Oc_arrI)
kde_data4 = makeKDE_array(Avd18Oc_arrId)
xx = np.linspace(-2, 3, 10000)
plt.plot(xx, kde_data1, linewidth=2)
plt.plot(xx, kde_data2, linewidth=2)
plt.plot(xx, kde_data3, linewidth=2)
plt.plot(xx, kde_data4, linewidth=2)
for i in range(len(x)):
    plt.axvline(x[i], ymax=0.3, c= 'r', linestyle='dashed', linewidth=1 )
for i in range(len(a)):
    plt.axvline(a[i], ymax=0.3, c= 'b', linestyle='dashed', linewidth=1)

plt.xlabel('$\delta^{18}O$c', fontsize=12)
plt.ylabel('Probability density function', fontsize=10)
plt.tick_params(labelsize=10)
plt.xlim(-1.5, 2.5)

plt.legend(("G. Ruber", "N. Pachyderma", 'G. Inflata', 'G. Inflata deep'))
plt.tight_layout()
fig.savefig(path+'modeldatacomp2.png',dpi= 600, bbox_inches='tight')
plt.show()




