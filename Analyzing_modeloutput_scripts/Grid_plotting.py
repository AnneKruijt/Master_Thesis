# -*- coding: utf-8 -*-
"""

Script which computes:
    
- the chance of particles having an average temperature above or below a certain threshold, for each location in the grid
- the offset between the tails of the local temperature distribution and the distribution of advected particles, for each location in the grid

@author: Anne
"""

import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import os
from scipy import stats
from scipy.stats import pearsonr
from scipy.stats import ttest_ind
import numpy.ma as ma


#_______importing functions to use on the datasets_____________________________________

os.chdir("C:/Users/Anne/Documents/Studie Master/Thesis/python_oefenen")

from functions_for_seasons import second_to_date, get_season, seasonvalue
from functions_for_traj import load_particles_file, find_dwell, find_final, chanceT, value_ranges, temp_grouping
from functions_for_plotting import Uruguayan_margin, finalloc_plot, alldwell_plot, histogram_plot
os.chdir("C:/Users/Anne/Documents/Studie Master/Thesis/python_oefenen")

#%% Load data file

#For grid
#file = 'C:/Users/Anne/Desktop/Datasets_thesis/for_countsTest_100days_dd50_sp200_res1.nc'
file1 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid2_10Years_part1_dd50_sp200_res1.nc'
file2 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid2_10Years_part2_dd50_sp200_res1.nc'
file3 ='C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid_10Years_part3_dd50_sp200_res1.nc'
file4 ='C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid2_10Years_part4_dd50_sp200_res1.nc'
file5 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid_10Years_part5_dd50_sp200_res1.nc'
file6 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid_10Years_part6_dd50_sp200_res1.nc'
file7 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid2_10Years_part7_dd50_sp200_res1.nc'
file8 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid_10Years_part8_dd50_sp200_res1.nc'
file9 = 'C:/Users/Anne/Desktop/Datasets_thesis/For_FINAL_figures/grid/new/Grid_10Years_part9_dd50_sp200_res1.nc'

T = load_particles_file(file1,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T2 = load_particles_file(file2,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T3 = load_particles_file(file3,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T4 = load_particles_file(file4,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T5 = load_particles_file(file5,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T6 = load_particles_file(file6,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T7 = load_particles_file(file7,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T8 = load_particles_file(file8,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
T9 = load_particles_file(file9,['time','lon', 'lat', 'age', 'shellage', 'temp', 'z', 'loctemp', 'lon0', 'lat0'])
#


#%%

def find_dwell(DwellTemp, Depth):
    obsmax = int(len(DwellTemp[0])) #(number of observations)
    trajmax = int(len(DwellTemp[:,0]))    # (number of particles and their trajectories)

    for i in range(trajmax):    
        for j in range(obsmax):
            if Depth[i,j] > 51:
                DwellTemp[i,j]= np.nan
            
                
    AvTemp = np.nanmean(DwellTemp, axis= 1)
    DwellTempnan = np.ma.filled(DwellTemp, np.nan)
    AvTempnan= np.ma.filled(AvTemp, np.nan)
    return DwellTempnan, AvTempnan, obsmax, trajmax

print ('have you set the dwelldepht correctly???')


#%%

def find_sharedGridPoint(ds, lat1, lon1):
    '''
    ds = dataset
    lon1 and lon2 between -180, +180
    '''
    lat = ds.lat0[:,1]
    lon = ds.lon0[:,1]
    indices     =  np.where((lat == lat1 ) & (lon == lon1 ))[0]
    return indices


#%%     

def temp_grouping(Av_variable, trajmax):
    Warmquantile = []
    Coldquantile = []
    quantiles = np.nanquantile(Av_variable, [0.10, 0.25, 0.5, 0.75, 0.90])
    for i in range(trajmax):
        if Av_variable[i] >= quantiles[4]:
            Warmquantile.append(Av_variable[i])
        if Av_variable[i] < quantiles[0]:
            Coldquantile.append(Av_variable[i])
    return  Warmquantile, Coldquantile


def WarmQavdifffunc(Av_variable, Avloc_variable):
    trajmax = int(len(Av_variable))
    Warm, Cold =temp_grouping(Av_variable, trajmax)
    Warmloc, Coldloc = temp_grouping(Avloc_variable, trajmax)
    a =np.mean(Warm)
    b =np.mean(Warmloc)
    diff= a-b
    return diff

def ColdQavdifffunc(Av_variable, Avloc_variable):
    trajmax = int(len(Av_variable))
    Warm, Cold =temp_grouping(Av_variable, trajmax)
    Warmloc, Coldloc = temp_grouping(Avloc_variable, trajmax)
    a =np.mean(Cold)
    b =np.mean(Coldloc)
    diff= a-b
    return diff

def Gridcalc (T, LatStart, LonStart):
    
    DwellTemp_arr, AvTemp, obsmax, trajmax, = find_dwell(T.temp, T.z)
    DwellTemploc_arr, AvTemploc_arr, obsmax, trajmax, = find_dwell(T.loctemp, T.z)
    
    Data_Loc = {}
    loc=0
    WarmDifGrid= np.zeros((len(LatStart),len(LonStart)))
    ColdDifGrid= np.zeros((len(LatStart),len(LonStart)))
    WarmChanceGrid= np.zeros((len(LatStart),len(LonStart)))
    ColdChanceGrid= np.zeros((len(LatStart),len(LonStart)))  
    WarmChanceGrid22= np.zeros((len(LatStart),len(LonStart)))
    ColdChanceGrid10= np.zeros((len(LatStart),len(LonStart)))    
    WarmChanceGrid18= np.zeros((len(LatStart),len(LonStart)))
    ColdChanceGrid8= np.zeros((len(LatStart),len(LonStart)))
    WarmChanceGrid16= np.zeros((len(LatStart),len(LonStart)))

    for i in range(int(len(LatStart))):
        for j in range(int(len(LonStart))):
            la = LatStart[i]
            lo = LonStart[j]
            index = find_sharedGridPoint(T, la, lo)
            
            Data_Loc[loc]={}
            
            Data_Loc[loc]['AvT']= AvTemp[index]
            Data_Loc[loc]['AvTna']= AvTemploc_arr[index]
            WarmDifGrid[i,j]=WarmQavdifffunc((Data_Loc[loc]['AvT']), (Data_Loc[loc]['AvTna']))
            ColdDifGrid[i,j]=ColdQavdifffunc((Data_Loc[loc]['AvT']), (Data_Loc[loc]['AvTna']))
            Data_Loc[loc]['T_warmDiff']= WarmQavdifffunc((Data_Loc[loc]['AvT']), (Data_Loc[loc]['AvTna']))
            Data_Loc[loc]['T_coldDiff']= ColdQavdifffunc((Data_Loc[loc]['AvT']), (Data_Loc[loc]['AvTna']))
            Data_Loc[loc]['TchanceW']=chanceT(Data_Loc[loc]['AvT'])
            Data_Loc[loc]['TchanceWloc']=chanceT(Data_Loc[loc]['AvTna'])
            WarmChanceGrid[i,j] = chanceWT(Data_Loc[loc]['AvT'])
            ColdChanceGrid[i,j] = chanceCT(Data_Loc[loc]['AvT'])
            WarmChanceGrid22[i,j] = chanceWT22(Data_Loc[loc]['AvT'])
            ColdChanceGrid10[i,j] = chanceCT10(Data_Loc[loc]['AvT'])
            WarmChanceGrid18[i,j] = chanceWT18(Data_Loc[loc]['AvT'])
            ColdChanceGrid8[i,j] = chanceCT8(Data_Loc[loc]['AvT'])
            WarmChanceGrid16[i,j] = chanceWT16(Data_Loc[loc]['AvT'])            
            loc = loc+1
        loc = loc +1
        
    return WarmDifGrid, ColdDifGrid, WarmChanceGrid, ColdChanceGrid, WarmChanceGrid22, ColdChanceGrid10, WarmChanceGrid18, ColdChanceGrid8, WarmChanceGrid16

''' Below are the functions used to compute the chance of particles having an average temperature above or below a certain threshold'''

def chanceWT(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a > 20).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance


def chanceWT18(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a > 18).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance

def chanceWT22(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a > 22).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance

def chanceWT16(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a > 16).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance


def chanceCT8(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a < 8).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance

def chanceCT10(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a < 10).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance

def chanceCT(Temperatures):
    a= np.sort(Temperatures)
    T_interest = (a < 12).sum()
    total = len(a)
    Tchance = T_interest/total *100
    return Tchance
 



#%% Performing computations for all files, together representing the whole grid under investigation
    
LatStart1 = np.arange(-54.5, -44.5, 1)
LonStart1 = np.arange(300.5, 310.5, 1)
WarmDif1, ColdDif1, Warmc1, Coldc1, Warm22c1, Cold10c1, Warm18c1, Cold8c1, Warm16c1 = Gridcalc(T, LatStart1, LonStart1)

LatStart2 = np.arange(-44.5, -34.5, 1)
LonStart2 = np.arange(300.5, 310.5, 1)
WarmDif2, ColdDif2, Warmc2, Coldc2, Warm22c2, Cold10c2, Warm18c2, Cold8c2, Warm16c2  = Gridcalc(T2, LatStart2, LonStart2)

LatStart3 = np.arange(-34.5, -24.5, 1)
LonStart3 = np.arange(300.5, 310.5, 1)
WarmDif3, ColdDif3, Warmc3, Coldc3, Warm22c3, Cold10c3, Warm18c3, Cold8c3, Warm16c3 = Gridcalc(T3, LatStart3, LonStart3)

LatStart4 = np.arange(-54.5, -44.5, 1)
LonStart4 = np.arange(310.5, 320.5, 1)
WarmDif4, ColdDif4, Warmc4, Coldc4, Warm22c4, Cold10c4, Warm18c4, Cold8c4, Warm16c4 = Gridcalc(T4, LatStart4, LonStart4)


LatStart5 = np.arange(-44.5, -34.5, 1)
LonStart5 = np.arange(310.5, 320.5, 1)
WarmDif5, ColdDif5, Warmc5, Coldc5, Warm22c5, Cold10c5, Warm18c5, Cold8c5, Warm16c5 = Gridcalc(T5, LatStart5, LonStart5)

LatStart6 = np.arange(-34.5, -24.5, 1)
LonStart6 = np.arange(310.5, 320.5, 1)
WarmDif6, ColdDif6, Warmc6, Coldc6, Warm22c6, Cold10c6, Warm18c6, Cold8c6, Warm16c6= Gridcalc(T6, LatStart6, LonStart6)


LatStart7 = np.arange(-54.5, -44.5, 1)
LonStart7 = np.arange(320.5, 330.5, 1)
WarmDif7, ColdDif7, Warmc7, Coldc7, Warm22c7, Cold10c7, Warm18c7, Cold8c7, Warm16c7 = Gridcalc(T7, LatStart7, LonStart7)

LatStart8 = np.arange(-44.5, -34.5, 1)
LonStart8 = np.arange(320.5, 330.5, 1)
WarmDif8, ColdDif8, Warmc8, Coldc8, Warm22c8, Cold10c8, Warm18c8, Cold8c8, Warm16c8 = Gridcalc(T8, LatStart8, LonStart8)


LatStart9 = np.arange(-34.5, -24.5, 1)
LonStart9 = np.arange(320.5, 330.5, 1)
WarmDif9, ColdDif9, Warmc9, Coldc9, Warm22c9, Cold10c9, Warm18c9, Cold8c9, Warm16c9 = Gridcalc(T9, LatStart9, LonStart9)


#%%

'''Loop to add all together '''

LonTot = np.concatenate((LonStart1, LonStart4, LonStart7))
LatTot = np.concatenate((LatStart1, LatStart2, LatStart3))
WarmDifTot = np.zeros((len(LonTot),(len(LonTot))))
ColdDifTot = np.zeros((len(LonTot),(len(LonTot))))
WarmcTot= np.zeros((len(LonTot),(len(LonTot))))
ColdcTot = np.zeros((len(LonTot),(len(LonTot))))
WarmcTot22 = np.zeros((len(LonTot),(len(LonTot))))
WarmcTot18 = np.zeros((len(LonTot),(len(LonTot))))
WarmcTot16 = np.zeros((len(LonTot),(len(LonTot))))
ColdcTot10 = np.zeros((len(LonTot),(len(LonTot))))
ColdcTot8 = np.zeros((len(LonTot),(len(LonTot))))

for i in range(int(len(LonTot))):
    if LonTot[i] < 310.:
        WarmDifTot[:,i] = np.concatenate((WarmDif1[:,i],WarmDif2[:,i], WarmDif3[:,i]))
        ColdDifTot[:,i] = np.concatenate((ColdDif1[:,i],ColdDif2[:,i], ColdDif3[:,i]))
        WarmcTot[:,i] = np.concatenate((Warmc1[:,i],Warmc2[:,i], Warmc3[:,i]))
        ColdcTot[:,i] = np.concatenate((Coldc1[:,i],Coldc2[:,i], Coldc3[:,i]))
        WarmcTot22[:,i] = np.concatenate((Warm22c1[:,i],Warm22c2[:,i], Warm22c3[:,i]))
        ColdcTot10[:,i] = np.concatenate((Cold10c1[:,i],Cold10c2[:,i], Cold10c3[:,i]))
        WarmcTot18[:,i] = np.concatenate((Warm18c1[:,i],Warm18c2[:,i], Warm18c3[:,i]))
        ColdcTot8[:,i] = np.concatenate((Cold8c1[:,i],Cold8c2[:,i], Cold8c3[:,i]))
        WarmcTot16[:,i] = np.concatenate((Warm16c1[:,i],Warm16c2[:,i], Warm16c3[:,i]))        
    if LonTot[i] < 319 and LonTot[i] > 310:
        WarmDifTot[:,i] = np.concatenate((WarmDif4[:,i-10],WarmDif5[:,i-10], WarmDif6[:,i-10]))
        ColdDifTot[:,i] = np.concatenate((ColdDif4[:,i-10],ColdDif5[:,i-10], ColdDif6[:,i-10]))
        WarmcTot[:,i] = np.concatenate((Warmc4[:,i-10],Warmc5[:,i-10], Warmc6[:,i-10]))
        ColdcTot[:,i] = np.concatenate((Coldc4[:,i-10],Coldc5[:,i-10], Coldc6[:,i-10]))
        WarmcTot22[:,i] = np.concatenate((Warm22c4[:,i-10],Warm22c5[:,i-10], Warm22c6[:,i-10]))
        ColdcTot10[:,i] = np.concatenate((Cold10c4[:,i-10],Cold10c5[:,i-10], Cold10c6[:,i-10]))
        WarmcTot18[:,i] = np.concatenate((Warm18c4[:,i-10],Warm18c5[:,i-10], Warm18c6[:,i-10]))
        ColdcTot8[:,i] = np.concatenate((Cold8c4[:,i-10],Cold8c5[:,i-10], Cold8c6[:,i-10]))
        WarmcTot16[:,i] = np.concatenate((Warm16c4[:,i-10],Warm16c5[:,i-10], Warm16c6[:,i-10]))
    if LonTot[i] > 319:
        WarmDifTot[:,i] = np.concatenate((WarmDif7[:,i-20],WarmDif8[:,i-20], WarmDif9[:,i-20]))
        ColdDifTot[:,i] = np.concatenate((ColdDif7[:,i-20],ColdDif8[:,i-20], ColdDif9[:,i-20]))
        WarmcTot[:,i] = np.concatenate((Warmc7[:,i-20],Warmc8[:,i-20], Warmc9[:,i-20]))
        ColdcTot[:,i] = np.concatenate((Coldc7[:,i-20],Coldc8[:,i-20], Coldc9[:,i-20]))
        WarmcTot22[:,i] = np.concatenate((Warm22c7[:,i-20],Warm22c8[:,i-20], Warm22c9[:,i-20]))
        ColdcTot10[:,i] = np.concatenate((Cold10c7[:,i-20],Cold10c8[:,i-20], Cold10c9[:,i-20]))
        WarmcTot18[:,i] = np.concatenate((Warm18c7[:,i-20],Warm18c8[:,i-20], Warm18c9[:,i-20]))
        ColdcTot8[:,i] = np.concatenate((Cold8c7[:,i-20],Cold8c8[:,i-20], Cold8c9[:,i-20]))
        WarmcTot16[:,i] = np.concatenate((Warm16c7[:,i-20],Warm16c8[:,i-20], Warm16c9[:,i-20]))

#%% Plotting
        
import netCDF4
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

def Uruguayan_margin():
    mp = plt.axes(projection=ccrs.PlateCarree())
    mp.set_xticks([-65, -55, -45, -35], crs=ccrs.PlateCarree())
    mp.set_yticks([-70, -60, -50, -40, -30, -20, -10], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter()
    lat_formatter = LatitudeFormatter()
    mp.xaxis.set_major_formatter(lon_formatter)
    mp.yaxis.set_major_formatter(lat_formatter)
    mp.grid(linewidth=1, color='black', alpha=0.3, linestyle='--')
    mp.set_extent([-60, -29, -56, -24], ccrs.PlateCarree()) #this sets the dimensions of your map (west-east-norht-south extent)
    mp.coastlines()
    mp.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='orangered'))


path = "C:/Users/Anne/Desktop/Datasets_thesis/Final_figures/" 

fig = plt.figure()
Uruguayan_margin()
levels	 = np.arange(0, 105, 1)

plt.contourf(LonTot, LatTot, WarmcTot16[:,:], levels, cmap = 'viridis', transform=ccrs.PlateCarree())
cbar = plt.colorbar()
cbar.ax.xaxis.set_label_position('top')
cbar.ax.set_xlabel('[$^\circ$C]') #('% chance')   #('[$^\circ$C]')
plt.scatter(-52.5, -36.5, c= 'k')
#plt.title( 'Offset in the warm tail', fontweight='bold',fontsize= 13)
#fig.savefig(path  + 'Offsetwarm' +'Grid.png',dpi= 600, bbox_inches='tight')
#plt.title( 'Offset in the cold tail', fontweight='bold',fontsize= 13)
#fig.savefig(path  + 'Offsetcold' +'Grid.png',dpi= 600, bbox_inches='tight')
#plt.title( 'Cold dwellers, < 8 $\degree$C', fontweight='bold',fontsize= 13)
#fig.savefig(path  + 'Chancecold8' +'Grid.png',dpi= 600, bbox_inches='tight')
plt.title( 'Warm dwellers, > 16 $\degree$C', fontweight='bold',fontsize= 13)
fig.savefig(path  + 'Chancewarm16' +'Grid.png',dpi= 600, bbox_inches='tight')
plt.show()

#%% For obtaining location within the grid where offsets are largest

print (np.nanmax(WarmDifTot))
print ( np.nanmin(WarmDifTot))

print (np.nanmax(ColdDifTot))
print ( np.nanmin(ColdDifTot))

indexwarmmax =(np.where(WarmDifTot==np.nanmax(WarmDifTot)))
indexcoldmin =(np.where(ColdDifTot[0:15,:]==np.nanmin(ColdDifTot[0:15,:])))


Warmtestlat = LatTot[indexwarmmax[0]]
Warmtestlon = LonTot[indexwarmmax[1]]
print (Warmtestlat)
print (Warmtestlon-360)

Coldtestlat=LatTot[indexcoldmin[0]]
Coldtestlon=LonTot[indexcoldmin[1]]

print (Coldtestlat)
print (Coldtestlon-360)

