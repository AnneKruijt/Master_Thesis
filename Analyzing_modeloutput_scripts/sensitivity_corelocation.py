# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 13:10:47 2019

Script to investigate the spatial sensitivity on small scale (0.5degree resolution),
in the study area.

Statistical difference between points is computed.

@author: Anne
"""

import netCDF4
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma 

#for plotting
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# for statistical analysis:
from scipy import stats
from scipy.stats import ks_2samp
from scipy.stats import pearsonr
from scipy.stats import ks_2samp, ttest_ind
from scipy.stats import entropy




directory= 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/sensitivity/'


#________________LOADING ALL data files

file1 = directory + 'Sink_Dwell_dd50_sp200_lon[-53]_lat[-36]_depth[679]_nr10.nc'
file2= directory + 'Sink_Dwell_dd50_sp200_lon[-53]_lat[-36.5]_depth[1597.7]_nr10.nc'
file3= directory + 'Sink_Dwell_dd50_sp200_lon[-53]_lat[-37]_depth[2945.5]_nr10.nc'
file4= directory + 'Sink_Dwell_dd50_sp200_lon[-52.5]_lat[-36]_depth[2003.75]_nr10.nc'
file5= directory + 'Sink_Dwell_dd50_sp200_lon[-52.5]_lat[-36.5]_depth[2578.8]_nr10.nc'
file6= directory + 'Sink_Dwell_dd50_sp200_lon[-52.5]_lat[-37]_depth[3301.8]_nr10.nc'
file7= directory + 'Sink_Dwell_dd50_sp200_lon[-52]_lat[-36]_depth[2800.5]_nr10.nc'
file8= directory + 'Sink_Dwell_dd50_sp200_lon[-52]_lat[-36.5]_depth[3347.4]_nr10.nc'
file9= directory + 'Sink_Dwell_dd50_sp200_lon[-52]_lat[-37]_depth[3609.8]_nr10.nc'



#____________Functions_________________________________________________#


def find_dwell(fname):
    nc = netCDF4.Dataset(fname)
    DwellTemp = nc.variables["temp"][:,:]
    Depth= nc.variables["z"][:,:]
    obsmax = len(DwellTemp[0]) #(number of observations)
    trajmax = len(DwellTemp[:,0])    # (number of particles and their trajectories)
    
    obsmax = int(obsmax)
    trajmax = int(trajmax)

    for i in range(trajmax):
    
        for j in range(obsmax):
            if Depth[i,j] > 51:
                DwellTemp[i,j]= np.nan
                    
    nc.close()
    return DwellTemp


def chanceT(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a > 22).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance


def makeKDE_array(values):
    def create_filterdata(data):
        data= ma.getdata(data)
        filterdata = data[~np.isnan(data)]
        return filterdata
    values= create_filterdata(values)
    kernel = stats.gaussian_kde(values)
    xx = np.linspace(0, 30, 1000)
    x = kernel(xx)
    return x



#%%    
#________________COMPUTING FOR ALL CORE LOCATIONS_________________________________________________-

## Making a dictionary
    
files = file1, file2, file3, file4, file5, file6 , file7, file8, file9
all_files = {}
for file in files:
    all_files[file] = {}

Medians = np.zeros(len(files))
Mean = np.zeros(len(files))
count = 0
''' This doesnt work yet'''
for file in files:
    all_files[file]['DwellTemp'] = find_dwell(file)
    all_files[file]['AvTemp']=np.nanmean(all_files[file]['DwellTemp'], axis= 1)
    all_files[file]['Tchance']=chanceT(all_files[file]['AvTemp'])
    all_files[file]['KDE']=makeKDE_array(all_files[file]['AvTemp'])
    Medians[count]=np.nanmedian(all_files[file]['AvTemp'])
    Mean[count]=np.nanmean(all_files[file]['AvTemp'])
    count = count +1

    

#%%
#_______________Plotting KDE's to visualise differences between locations_______________________________________
path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/" 
    
kde1 = all_files[file4]['KDE']
kde2 = all_files[file6]['KDE']

xx = np.linspace(0, 30, 1000)
#n, bins, patches = plt.hist([filterdata2, filterdata2], 50, normed=1, facecolor='g')
fig = plt.figure()
plt.plot(xx, kde1)
plt.plot(xx, kde2)
plt.xlabel('Temperature' '[$^\circ$C]', fontsize=12)
plt.xlim(5, 30)
#plt.ylabel('Probability')
plt.legend(("Location 4", "Location 6"))
plt.text(7, 0.100, "ks-statistic = 0.27")
plt.text(7, 0.08800, "p-value = $1.44*10^{-39}$")
fig.savefig(path  + 'Loc4-6.png',dpi= 600, bbox_inches='tight')

#%%
#________________Performing statistical tests on KDE's of all sites... does this make sense>???
i=0
KSval=np.zeros((len(files), len(files)))
Pksval=np.zeros((len(files), len(files)))

for file in files:
    x = all_files[file]['KDE']
    i=i+1
    j=0
    for file in files:
        y= all_files[file]['KDE']
        KSval[i-1,j], Pksval[i-1,j]= stats.ks_2samp(x,y)
        j=j+1

#%% 
#________Testing the statistical difference between all point in the grid

i=0
Tval=np.zeros((len(files), len(files)))
Pval=np.zeros((len(files), len(files)))
KSval2=np.zeros((len(files), len(files)))
Pksval2=np.zeros((len(files), len(files)))
entro = np.zeros((len(files), len(files)))

for file in files:
    x = all_files[file]['AvTemp']
    new_mask = (np.isnan(x))
    x_nonan = x[~new_mask]
    i=i+1
    j=0
    for file in files:
        y= all_files[file]['AvTemp']
        new_mask = (np.isnan(x))
        y_nonan = y[~new_mask]
        Tval[i-1,j], Pval[i-1,j]= stats.ttest_ind(x_nonan,y_nonan)
        KSval2[i-1,j], Pksval2[i-1,j]= stats.ks_2samp(x_nonan,y_nonan)
        entro[i-1,j] = entropy(x_nonan, qk=y_nonan)
        j=j+1


#%%
#______Plotting________________________________#

path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/" 


def create_filterdata(file):
    #data= all_files[file]['Tchance']
    data= all_files[file]['AvTemp']
    data= ma.getdata(data)
    filterdata = data[~np.isnan(data)]
    return filterdata



filterdata1 = create_filterdata(file1)
filterdata2 = create_filterdata(file2)
filterdata3 = create_filterdata(file3)
filterdata4 = create_filterdata(file4)
filterdata5 = create_filterdata(file5)
filterdata6 = create_filterdata(file6)
filterdata7 = create_filterdata(file7)
filterdata8 = create_filterdata(file8)
filterdata9 = create_filterdata(file9)



data = [filterdata1, filterdata2, filterdata3, filterdata4, filterdata5, filterdata6, filterdata7, filterdata8, filterdata9]
fig, ax1 = plt.subplots(figsize=(10, 6))
plt.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = plt.boxplot(data, notch=0, sym='+', vert=1, whis=1.5)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['fliers'], color='red', marker='+')

# Add a horizontal grid to the plot, but make it very light in color
# so we can use it for reading data values but not be distracting
ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax1.set_axisbelow(True)
ax1.set_title("Range of foraminiferal temperature signals", fontweight='bold', fontsize=13 )
ax1.set_xlabel('Core sites', fontsize =13)
ax1.set_ylabel('Average Temperatures [$^\circ$C]', fontsize =13)

TickNames = ['1', '2','3', '4', '5', '6', '7', '8', '9'] 
xtick =plt.setp(ax1, xticklabels=np.repeat(TickNames, 1))
plt.setp(xtick, rotation=0, fontsize=10)
fig.savefig(path +'Av_T_allcores.png', dpi= 600)
plt.show()

