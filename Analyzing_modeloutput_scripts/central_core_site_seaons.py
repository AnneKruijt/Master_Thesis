# -*- coding: utf-8 -*-
"""

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

## ---------Seasonal analysis-----------------------

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

def second_to_date(trajmax, obsmax, time):
    
    Unix_starttime = '1970-01-01T00:00:00Z'
    formatted_timeUTZ = datetime.strptime(Unix_starttime, '%Y-%m-%dT%H:%M:%SZ')
    experiment_start = '2000-01-09T00:00:00Z'
    formatted_exptime = datetime.strptime(experiment_start, '%Y-%m-%dT%H:%M:%SZ')
    time_delta = formatted_exptime - formatted_timeUTZ
    delta_in_seconds = time_delta.total_seconds()
    
    dates=np.zeros((trajmax, obsmax), dtype='datetime64[s]') # dtype= object
    totalruntime = time[0,0]
    timenan = np.ma.filled(time, np.nan)   # some items of the time array are masked, here they are filled with nans
    for i in range(trajmax):
        for j in range(obsmax):
            if np.isnan(timenan[i,j]):   # replacing the nan with a  strange default value! 
                dates[i,j]= datetime(2100,12,12,0,0,0)  #the strange default value
            else:
                dates[i,j]=datetime.fromtimestamp(delta_in_seconds + totalruntime - timenan[i,j])
    return dates    # be aware that dates contains unrealistic default values, with the same index as the masked elements in the time array
    

def get_season(seasons, day):
    for season,(season_start, season_end) in seasons.items():
        if day>=season_start and day<= season_end:
            return season
    else:
        return '3'


def seasonvalue(trajmax, obsmax, dates):
    seasonval=np.zeros((trajmax, obsmax))
    for i in range(trajmax):
        for j in range(obsmax):
            day_array = dates[i,j]
            Y = day_array.astype(object).year
            seasons = {'1':(datetime(Y,6,21,0,0,0), datetime(Y,9,22,23,59,59)),
                       '2':(datetime(Y,9,23,0,0,0), datetime(Y,12,20,23,59,59)),
                       '4':(datetime(Y,3,21,0,0,0), datetime(Y,6,20,23,59,59))}
            seasonval[i,j]= get_season(seasons, day_array)
    return seasonval

 

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
Timecentral = T.time

#%% Using functions on data
DwellTemp_arr, AvTemp_arr, Stdev_arr, obsmax, trajmax, = find_dwell(T.temp, T.z)
DwellTemploc_arr, AvTemploc_arr, Stdevloc_arr, obsmax, trajmax, = find_dwell(T.loctemp, T.z)

dates = second_to_date(trajmax, obsmax, T.time)
seasonval = seasonvalue(trajmax, obsmax, dates)

indi = np.where(T.time.mask == True)
seasonval[indi]= "nan"

#%% plotting the ranges of temperatures for different seasons
path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/" 


def create_filterdata(data):
    data= ma.getdata(data)
    filterdata = data[~np.isnan(data)]
    return filterdata

winter = np.where(seasonval == 1)
spring = np.where(seasonval == 2)
summer = np.where(seasonval == 3)
autumn = np.where(seasonval == 4)


data1= DwellTemp_arr[winter]
data2= DwellTemp_arr[spring]
data3= DwellTemp_arr[summer]
data4= DwellTemp_arr[autumn]

data1b= DwellTemploc_arr[winter]
data2b= DwellTemploc_arr[spring]
data3b= DwellTemploc_arr[summer]
data4b= DwellTemploc_arr[autumn]

filterdata1= create_filterdata(data1)
filterdata2= create_filterdata(data2)
filterdata3= create_filterdata(data3)
filterdata4= create_filterdata(data4)

filterdata1b= create_filterdata(data1b)
filterdata2b= create_filterdata(data2b)
filterdata3b= create_filterdata(data3b)
filterdata4b= create_filterdata(data4b)
filterdata5= create_filterdata(DwellTemp_arr)
filterdata5b= create_filterdata(DwellTemploc_arr)
filterdata6 = create_filterdata(AvTemp_arr)

#%% Fancy boxplotting
# (source: https://matplotlib.org/examples/pylab_examples/boxplot_demo2.html)

from matplotlib.patches import Polygon

data = [filterdata1, filterdata1b, filterdata2, filterdata2b, filterdata3, filterdata3b, filterdata4, filterdata4b, filterdata5, filterdata5b]

fig, ax1 = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey',
               alpha=0.5)

# Hide these grid behind plot objects
ax1.set_axisbelow(True)
ax1.set_title('Range of dwelling temperatures for each season', fontweight='bold', fontsize=13 )
#ax1.set_xlabel('Seasons', fontsize=13)
ax1.set_ylabel('Dwelling Temperatures [$^\circ$C]', fontsize=13)


# Now fill the boxes with desired colors
boxColors = ['darkkhaki', 'royalblue']
numBoxes = len(data)
medians = list(range(numBoxes))
for i in range(numBoxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    boxCoords = list(zip(boxX, boxY))
    # Alternate between Dark Khaki and Royal Blue
    k = i % 2
    boxPolygon = Polygon(boxCoords, facecolor=boxColors[k])
    ax1.add_patch(boxPolygon)
    # Now draw the median lines back over what we just filled in
    med = bp['medians'][i]
    medianX = []
    medianY = []
    for j in range(2):
        medianX.append(med.get_xdata()[j])
        medianY.append(med.get_ydata()[j])
        plt.plot(medianX, medianY, 'k')
        medians[i] = medianY[0]
    # Finally, overplot the sample averages, with horizontal alignment
    # in the center of each box
    plt.plot([np.average(med.get_xdata())], [np.average(data[i])],
             color='w', marker='*', markeredgecolor='k')

# Set axes labels

TickNames = ['Winter', '', 'Spring', '', 'Summer', '', 'Autumn', '', 'All', ''] 
xtick =plt.setp(ax1, xticklabels=np.repeat(TickNames, 1))
plt.setp(xtick, rotation=0, fontsize=10)
fig.savefig(path + 'Av_T_seasons.png',dpi= 600)
plt.show()

