# -*- coding: utf-8 -*-
"""
Created on Tue May 21 11:13:01 2019

@author: Anne

Investigation of central core site (-36.5N, -52.5E )
Script that visualizes the total traveled distances of particles at Urugayan margin

- map showing final positions of particles and their total travelled distance
- histograms of spread in total distances

"""

#matplotlib inline
import netCDF4
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import math

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from matplotlib import cm

from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import ttest_ind

#%% 
# --------------------FUNCTIONS------------------------

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


#%% Importing data to be investigated

file = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/Distance_Sink_Dwell_dd50_sp200_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.0.nc'
nc = netCDF4.Dataset(file)
nc.variables

T = load_particles_file(file, ['lon', 'lat', 'age', 'distance', 'z'])

#%% Computing travelled distances

trajdist = nc.variables["distance"][:,:]
trajlon = nc.variables["lon"][:,:]
trajlat = nc.variables[ "lat"][:,:]
trajdepth = nc.variables["z"][:,:]

obsmax = int(len(trajdist[0])) # = (number of observations)
trajmax = int(len(trajdist[:,0]))    # = (number of particles and their trajectories)

#______Save final travelled distance of each particle and their final location ________________________________#

Finaldist = np.zeros(trajmax)
FinalLon  = np.zeros(trajmax)
FinalLat  = np.zeros(trajmax)
directdist= np.zeros(trajmax)           # for the direct distance from core site to final location
sinkdist = np.zeros((trajmax,obsmax))   # for the distance travelled while sinking
dist_direct = np.zeros((trajmax,obsmax))

# Obtaining total distance along trajectory and direct distance from core to origin location
for i in range(trajmax):
    new_mask = (np.isnan(trajdist[i,:]))           # making vector with masked elements for each non-nan value
    Finaldist[i] = trajdist[i,:][~new_mask][-1]    # applying the masked vector to the row, making each nan false
    FinalLon[i]  = trajlon[i,:][~new_mask][-1]
    FinalLat[i]  = trajlat[i,:][~new_mask][-1]
   
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    lat_dist = (FinalLat[i] - trajlat[0,0]) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    lon_dist = (FinalLon[i] - trajlon[0,0]) * 1.11e2 * math.cos(FinalLon[i] * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    directdist[i] += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))
    
# Computing total distance travelled while sinking        
for i in range(trajmax):
        for j in range(obsmax):
            # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
            dist_directlat = (trajlat[i,j] - trajlat[0,0]) * 1.11e2
            # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
            dist_directlon = (trajlon[i,j] - trajlon[0,0]) * 1.11e2 * math.cos(FinalLon[i] * math.pi / 180)
            # Calculate the total Euclidean distance travelled by the particle
            dist_direct[i,j] += math.sqrt(math.pow(dist_directlon, 2) + math.pow(dist_directlat, 2))
            if trajdepth[i,j]<50:
                sinkdist[i,:]=dist_direct[i,j]
                break

Finaldist[Finaldist == 0]= 'nan'  # getting rid of the 0 values (from the initial filling of 'empty 'array')
FinalLon[FinalLon == 0]= 'nan'                 
FinalLat[FinalLat == 0]= 'nan'                

#%%
# ---------------Plotting along path travelled distance in map -----------------------------------

path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/"  

fig = plt.figure()
Uruguayan_margin()


plt.scatter(FinalLon[:], FinalLat[:], c=Finaldist[:], s=1)
plt.set_cmap('spring')
cbar = plt.colorbar()
cbar.ax.xaxis.set_label_position('top')
cbar.ax.set_xlabel('km', fontsize=10 )
plt.scatter(-52.5, -36.5, c='k', s=5)
name = str('Total travelled distance \nfrom core to origin' + '\n')
plt.title(name, fontweight='bold',fontsize= 10)

fig.savefig(path  + 'Distance' +'lat_'+ str(T.lat[0,0]) + 'lon_'+str(T.lon[0,0]) + 'depth_'+ str(T.z[0,0])+'.png',dpi= 600, bbox_inches='tight')
plt.show()



#%%
#----------------Plotting the histogram of different travelled distances------------------------

def histogram_plot(variable, name):
    # masking the nan-objects, for creation of histogram
    fig = plt.figure()
    x = np.array(variable)                  # changing from list to array
    x = x[np.logical_not(np.isnan(x))]
    n, bins, patches = plt.hist(x, 30, density=1, facecolor='g')
    plt.title(name, fontweight='bold',fontsize= 13)
    plt.xlabel('Travelled distance (km)', fontsize=11)
    plt.ylabel('Propability density', fontsize=10)
    fig.savefig(path + 'Histogram_' + name +'.png',dpi= 600)
    plt.show()

sink= np.array(sinkdist[:,0]) 
name1 =str('Distance along trajectory')
name2= str('Direct distance')
name3 =str('Distance while sinking') 
histogram_plot(Finaldist, name1)
histogram_plot(directdist, name2)
histogram_plot(sink, name3)


