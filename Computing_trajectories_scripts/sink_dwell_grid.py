# -*- coding: utf-8 -*-
"""
Modification of script: Sink_dwell.py
- used to release particles at multiple locations in a grid simultaneously
- loading topography file to determine if points in a grid are on land or in the ocean

@author: Anne Kruijt
"""

''' Importing relevant modules'''

from parcels import (FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
from netCDF4 import Dataset
import sys

def GetOFESLandArray(filename, fieldname):
    """
    Function to return a Field with 1's at land points and 0's at ocean points, based on python basemap
    :param f: a field of the .nc file, but not a parcels field object! For OFES, land points are masked. This is used here!
    """
    pfile = Dataset(filename, 'r')
    Lon = pfile.variables['LONN1799_1800'][:]
    Lat = pfile.variables['LAT'][:]
    f = pfile.variables[fieldname][:]
    f = f[0,0,:,:] 
    Land=Field('Land',f,transpose=False, lon=Lon,lat=Lat)
    return Land
    
def snapshot_function(start, end, delta):
    """
    
    """
    curr = start
    result = []
    result.append('{:%Y%m%d}'.format(curr))
    while curr <= end:
     #   yield curr
        curr += delta
        result.append('{:%Y%m%d}'.format(curr))
    del result[-1]
    result = np.array(result)
    return result

## Characteristics of foraminifera
sp = 200. #The sinkspeed m/day
dd = 50. #The dwelling depth
ls = 30.0 # life span while dwelling, in days


##Data are located on the imau folder on gemini and results written to thesis_folder
directory = '/data2/imau/oceanparcels/hydrodynamic_data/OFESdata/OFES_0.1_HIND/DATA_3D/snap_3day'
directoryb = '/home/students/4082842/projects/Master_thesis/'
dirwrite ='/scratch/AnneK/'

## full grid spanning lontitudes from -60 to -30 and latitudes from -55 to -25
## split into 9 grids, for faster run-time of each separare script

#creating part of the grid : part 1
# minlon = -60
# maxlon = -50
# minlat = -55
# maxlat = -45
# res = 1    #resolution in degrees

# # # creating part of the grid : part 2
# minlon = -60
# maxlon = -50
# minlat = -45
# maxlat = -35
# res = 1    #resolution in degrees

# # creating part of the grid : part 3
# minlon = -60
# maxlon = -50
# minlat = -35
# maxlat = -25
# res = 1    #resolution in degrees

# # creating part of the grid : part 4
# minlon = -50
# maxlon = -40
# minlat = -55
# maxlat = -45
# res = 1    #resolution in degrees

# # creating part of the grid : part 5
# minlon = -50
# maxlon = -40
# minlat = -45
# maxlat = -35
# res = 1    #resolution in degrees

# # creating part of the grid : part 6
# minlon = -50
# maxlon = -40
# minlat = -35
# maxlat = -25
# res = 1    #resolution in degrees

# # creating part of the grid : part 7
# minlon = -40
# maxlon = -30
# minlat = -55
# maxlat = -45
# res = 1    #resolution in degrees

# # creating part of the grid : part 8
# minlon = -40
# maxlon = -30
# minlat = -45
# maxlat = -35
# res = 1    #resolution in degrees

# # creating part of the grid : part 9
minlon = -40
maxlon = -30
minlat = -35
maxlat = -25
res = 1    #resolution in degrees

latminind = (max(-75,minlat-30)+75)*10
latmaxind = (min(75,maxlat+30)+75)*10

grid = np.mgrid[int(minlon):int(maxlon):res,int(minlat):int(maxlat):res]+0.5
n=grid[0].size
lons = np.reshape(grid[0],n)
lats = np.reshape(grid[1],n)

#Delete the particles on the land
landfilename = directoryb + "topography_OFES.nc"
Land = GetOFESLandArray(landfilename, 'HT')
[lonsz,latsz]=[np.array([lo for lo, la in zip(lons,lats) if Land[0,0,la,lo]>dd  ]),np.array([la for lo, la in zip(lons,lats) if Land[0,0,la,lo]>dd ])]

lonsz[lonsz<0] = lonsz[lonsz<0]+360

if(not lonsz.size):
    sys.exit("Only land in the run with this idx")

times = np.array([datetime(2010, 6, 15) - delta(days=x) for x in range(0,int(365*10+1),3)]) 
time = np.empty(shape=(0));lons = np.empty(shape=(0));lats = np.empty(shape=(0)); dep =np.empty(shape=(0))
for i in range(len(times)):
    lons = np.append(lons,lonsz)
    lats = np.append(lats, latsz)
    time = np.append(time, np.full(len(lonsz),times[i]))


sys.stdout.flush()  
#%%
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

    sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    
    bfile = directoryb + 'topography_OFES.nc'

    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles, 'salin': sfiles, 'B':bfile}
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature', 'salin':'salinity', 'B': 'HT'}
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'B':{'lat': 'LAT', 'lon': 'LONN1799_1800', 'time': 'TIME', 'depth': 'LEV'}}
        
    #indices = {'lat': range(0,900), 'lon': range(2700,3600)}                 
    indices = {'lat': range(0,950)}
    #indices = {'lat': range(latminind, latmaxind)}
    
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=True)

    return fieldset
        
#%% Functions to be used
def SampleTemp(particle, fieldset, time):
    particle.temp = fieldset.temp[time, particle.depth, particle.lat, particle.lon]
        
def LocalConditions(particle, fieldset, time):
    lonbottom = particle.lon0
    latbottom = particle.lat0
    particle.loctemp = fieldset.temp[time, particle.depth, latbottom, lonbottom]    #does this work?
    
    
def Sink(particle, fieldset, time):
    if particle.depth > fieldset.dwellingdepth:
        particle.depth = particle.depth + fieldset.sinkspeed * dt
    else:
        particle.depth = fieldset.dwellingdepth
    


def Age(particle, fieldset, time):
    if particle.depth <= fieldset.dwellingdepth:
        particle.age = particle.age + math.fabs(dt)
        
    if particle.age > fieldset.maxage:
        particle.delete()

def ShellAge(particle, fieldset, time):
    particle.shellage = particle.shellage + math.fabs(dt)
    

def initials(particle, fieldset, time):
    if particle.shellage==0.:
        particle.depth = fieldset.B[time+dt, particle.depth, particle.lat, particle.lon]
        particle.lon0 = particle.lon
        particle.lat0 = particle.lat
        particle.depth0 = particle.depth
    
        
def DeleteParticle(particle, fieldset, time):
    particle.delete()

 
def run_corefootprintparticles(dirwrite,outfile,lonss,lonsz,latss,latsz,startdep,time):
    snapshots = snapshot_function(date(2000,1,9), date(2010, 6, 15),delta(days=3))
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.B.allow_time_extrapolation = True
    fieldset.add_constant('dwellingdepth', np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400.)
    fieldset.add_constant('maxage', 30.*86400)


    class ForamParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=np.nan)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin = Variable('salin', dtype=np.float32, initial=np.nan)
        lon0 = Variable('lon0', dtype=np.float32, initial=0.)
        lat0 = Variable('lat0', dtype=np.float32, initial=0.)
        depth0 = Variable('depth0',dtype=np.float32, initial=0., to_write=False)
        time0 = Variable('time0',dtype=np.float32, initial=0., to_write=False)
        shellage = Variable('shellage', dtype=np.float32, initial=0.)
        loctemp =Variable('loctemp', dtype=np.float32, initial=np.nan) 

    pset = ParticleSet.from_list(fieldset=fieldset, pclass=ForamParticle, lon=lonss.tolist(), lat=latss.tolist(), time = time)
    
    pfile = ParticleFile(dirwrite + outfile, pset, outputdt=delta(days=1))

    kernels = pset.Kernel(initials) + ShellAge + Age + Sink +AdvectionRK4_3D  + SampleTemp  + LocalConditions

    pset.execute(kernels, runtime=delta(days=3650), dt=delta(minutes=-5), output_file=pfile, recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile = "Grid_part_number" +'_dd'+str(int(dd)) +'_sp'+str(int(sp))+"_res"+str(res)

run_corefootprintparticles(dirwrite,outfile,lons,lonsz,lats,latsz,dep,time)

print 'Exection finished'