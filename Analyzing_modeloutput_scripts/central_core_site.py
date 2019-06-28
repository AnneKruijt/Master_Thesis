# -*- coding: utf-8 -*-
"""
Created on Mon May 20 23:34:47 2019

Investigation of central core site (-36.5N, -52.5E )

- average temperature signal at origin locations in map
- PDF's of the temperature range and the warm and cold tails

@author: Anne
"""

import netCDF4
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import os

# for statistics
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import ttest_ind

# for plotting
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% 
#----------------------FUNCTIONS--------------------

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
            if Depth[i,j] > 51:
                DwellTemp[i,j]= np.nan
                
    AvTemp = np.nanmean(DwellTemp, axis= 1)
    Stdev = np.nanstd(DwellTemp, axis=1)
    return DwellTemp, AvTemp, Stdev, obsmax, trajmax

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
    xx = np.linspace(-10, 30, 10000)
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



# ----------Plotting Functions-------------------
    
def Uruguayan_margin():
    mp = plt.axes(projection=ccrs.PlateCarree())
    mp.set_xticks([-65, -55, -45, -35], crs=ccrs.PlateCarree())
    mp.set_yticks([-50, -40, -30, -20, -10], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    mp.xaxis.set_major_formatter(lon_formatter)
    mp.yaxis.set_major_formatter(lat_formatter)
    mp.grid(linewidth=1, color='black', alpha=0.3, linestyle='--')
    mp.set_extent([-60, -30, -50, -10], ccrs.PlateCarree()) #this sets the dimensions of your map (west-east-norht-south extent)
    mp.coastlines()
    mp.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='orangered'))

#%%
#----------Import data for temperature investigation------------------------------------------------#

file = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/for_countsSink_Dwell_pluslocal_dd50_sp200_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc'


nc = Dataset(file, 'r')
nc.variables
T = load_particles_file(file,['time','lon', 'lat', 'age', 'temp', 'z', 'loctemp'])

#%% Using functions to compute relevant variables

DwellTemp_arr, AvTemp_arr, Stdev_arr, obsmax, trajmax, = find_dwell(T.temp, T.z)
DwellTemploc_arr, AvTemploc_arr, Stdevloc_arr, obsmax, trajmax, = find_dwell(T.loctemp, T.z)

indi = np.where(T.time.mask == True)
DwellTemp_arr[indi]= "nan"
Final_temp, FinalLon, FinalLat = find_final(trajmax, DwellTemp_arr, T.lon, T.lat)

#%% 
# -------------Plotting average temperature in map --------------------
    
 
name = str('Average dwelling temperatures')
path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/"  


sinks = np.where(T.z > 50.) #100.)
fig = plt.figure()
Uruguayan_margin()

plt.scatter(FinalLon[:], FinalLat[:], c=AvTemp_arr[:], s=1)
plt.set_cmap('viridis')
cbar = plt.colorbar()
plt.clim(10,25)
cbar.ax.xaxis.set_label_position('top')
cbar.ax.set_xlabel('[$^\circ$C]')
plt.scatter(T.lon[sinks], T.lat[sinks], c= 'k', s=1)

plt.title(name + '\n', fontweight='bold',fontsize= 10)    
fig.savefig(path+name + 'lat_'+ str(T.lat[0,0]) + 'lon_'+str(T.lon[0,0]) + 'depth_'+ str(T.z[0,0]) +'.png',dpi= 600, bbox_inches='tight')
plt.show()

#%%  
# -------------Plotting PDF of the temperature range
quantiles = np.nanquantile(AvTemp_arr, [0.10, 0.25, 0.5, 0.75, 0.90])

fig = plt.figure()
kde_data1 = makeKDE_array(AvTemp_arr)
kde_data2 = makeKDE_array(AvTemploc_arr)
xx = np.linspace(-10, 30, 10000)
plt.plot(xx, kde_data1, linewidth=1)
plt.plot(xx, kde_data2, linewidth=1)

plt.axvline(x=quantiles[0], c= 'r', linestyle='dashed', linewidth=1 )
plt.axvline(x=quantiles[4], c= 'r', linestyle='dashed', linewidth=1)

plt.xlabel('Temperatures [$^\circ$C]', fontsize=10)
plt.ylabel('Propability density function', fontsize=10)
plt.tick_params(labelsize=10)
plt.xlim(10, 30)
plt.legend(("With advection", "Without advection"))
plt.title(name, fontweight='bold',fontsize= 10)
plt.tight_layout()
fig.savefig(path+ 'KDE'+ name + 'lat_'+ str(T.lat[0,0]) + 'lon_'+str(T.lon[0,0]) + 'depth_'+ str(T.z[0,0]) +'.png',dpi= 600, bbox_inches='tight')
plt.show()

#%%
# ----------- Plotting PDF's of tails of distribution --------------------

Temp_group, Warm_species, Trans_species, Cold_species, Warmquantile, Coldquantile, Tempcount = temp_grouping(AvTemp_arr, trajmax)
Temp_grouploc, Warm_speciesloc, Trans_speciesloc, Cold_speciesloc, Warmquantileloc, Coldquantileloc, Tempcountloc = temp_grouping(AvTemploc_arr, trajmax)


# KDE of the distribution of the tails
fig = plt.figure()
kde_data1 = makeKDE_array(Warmquantile)
kde_data2 = makeKDE_array(Warmquantileloc)
kde_data3 = makeKDE_array(Coldquantile)
kde_data4 = makeKDE_array(Coldquantileloc)
xx = np.linspace(10, 30, 10000)
plt.plot(xx, kde_data1, linewidth=1)
plt.plot(xx, kde_data2, linewidth=1)
plt.plot(xx, kde_data3, linewidth=1)
plt.plot(xx, kde_data4, linewidth=1)
plt.xlabel('Temperatures [$^\circ$C]', fontsize=10)
plt.ylabel('Propability density function', fontsize=10)
plt.title('Warm and cold tails of the distribution', fontweight='bold',fontsize= 10)
plt.tick_params(labelsize=10)
plt.xlim(18, 30)
plt.legend(("With advection", "Without advection", 'With advection', 'Without advection'))
plt.tight_layout()
fig.savefig(path+ 'TAILS'+ name + 'lat_'+ str(T.lat[0,0]) + 'lon_'+str(T.lon[0,0]) + 'depth_'+ str(T.z[0,0]) +'.png',dpi= 600, bbox_inches='tight')
plt.show()

