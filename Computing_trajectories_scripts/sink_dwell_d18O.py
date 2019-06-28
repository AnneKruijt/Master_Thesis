# -*- coding: utf-8 -*-
"""
Modification of script: Sink_dwell.py
- addition of kernels for computing d18Ow and d18Oc
- addition of salinity field necessary for the computation of d18O

Use to compute particle trajectory and store local water temperature, 
d18O of the water and d18O of foraminifera's calcite  along the path,
and at the core site location.

@author: Anne
"""

''' Importing relevant modules'''
from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
import sys

''' function used for creating an array of the 3-day snapshots of OFES output, for the required time frame '''
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

##Data are located on the IMAU folder on gemini and results written to the authors scratch directory
directory = '/data2/imau/oceanparcels/hydrodynamic_data/OFESdata/OFES_0.1_HIND/DATA_3D/snap_3day'
dirwrite ='/scratch/AnneK/'


## Characteristics of foraminifera
## Different settings are used for different species, for accurate prediction of d18Oc values and enable comparison with values obtained from geochemical analysis of core top sediments.

## Characters for Pachyderma
# sp = 250. #The sinkspeed m/day
# dd = 100. #The dwelling depth
# ls = 30.0 # life span while dwelling, in days

# ##Characters for Glutinata
# sp = 200. #The sinkspeed m/day
# dd = 80. #The dwelling depth
# ls = 30.0 # life span while dwelling, in days

# ##Characters for Ruber
# sp = 200. #The sinkspeed m/day
# dd = 50. #The dwelling depth
# ls = 15.0 # life span while dwelling, in days

##Characters for Inflata
sp = 250. #The sinkspeed m/day
dd = 300. #The dwelling depth
ls = 30.0 # life span while dwelling, in days

# location and depth
###Lons and Lats of interest: Uruguayan margin: domain=[-20.0, -60.0,-40.0,-80.0]
corelon = [-52.5] #52.91     #53.2  	   #-52.41		#53 max value 
corelat = [-36.5] #36.91     #36.5 	   #-36.23		#37 max value
coredepth = [2578.8] #3026.75  #1154.19   #2535		#5580 
bottomlon = -52.5
bottomlat = -36.5

# amount of years of snapshots to be loaded
snap = 10
print 'did you set the amount of snapshots correctly?'

#%%
''' function that loads the relevant velocity and temperature fields and creates a fieldset'''
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    
    sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    
    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles, 'salin': sfiles} 
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature', 'salin':'salinity'} 
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'}}
                       
    ## Different indices can be picked for prefered size of loaded field                 
    indices = {'lat': range(100,700), 'lon': range(2900,3600)}
    #indices = {'lat': range(0,950)} 
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=True)

    return fieldset

def LocalConditions(particle, fieldset, time):
    lonbottom = fieldset.bottomlon
    latbottom = fieldset.bottomlat
    particle.loctemp = fieldset.temp[time, particle.depth, latbottom, lonbottom]    #does this work?
       
def SampleTemp(particle, fieldset, time):
    particle.temp = fieldset.temp[time, particle.depth, particle.lat, particle.lon]

# computing the local d18Ow and d18Oc values above the core site
def Sampled18OLoc(particle, fieldset, time):
    lonbottom = fieldset.bottomlon
    latbottom = fieldset.bottomlat
    
    # constants for linear relations between d18Owater and salinity, depending on watermass
    # obtained from Le Grande 2006 (noaa website)
    St= 0.15
    It= -4.61

    Sss= 0.51
    Iss= -17.40

    Snad= 0.51
    Inad= -17.75

    Saab= 0.23
    Iaab= -8.11

    Smd= 0.42
    Imd= -14.38

    Tl = fieldset.temp[time, particle.depth, latbottom, lonbottom]
    Sl = fieldset.salin[time, particle.depth, latbottom, lonbottom]

    # identifying watermass (T and S ranges derrived from Emery, 2003)
    if Tl > 18.:
        particle.d18owL = St*Sl + It
    if Tl <18. and Tl > 4.:
        particle.d18owL = Sss*Sl + Iss
    if Tl < 4. and Tl > 1. and Sl < 35. and Sl > 34.8:
        particle.d18owL = Snad*Sl + Inad
    if Tl < 2. and Sl < 34.8 and Sl > 34.5:
        particle.d18owL = Saab*Sl + Iaab
    else:
        particle.d18owL = Smd*Sl + Imd

    #particle.d18ocL  = particle.d18owL - 0.27 - (Tl-14.9)/4.8 #Ruber equation
    #particle.d18ocL  = particle.d18owL - 0.27 - (Tl-16.5)/4.8 #Pachyderma equation
    particle.d18ocL =  particle.d18owL - 0.27 - (Tl-15.2)/4.6 #Inflata and Glutinata equation

# computing the d18Ow and d18Oc values along the trajectory  
def Sampled18O(particle, fieldset, time):
    St= 0.15
    It= -4.61

    Sss= 0.51
    Iss= -17.40

    Snad= 0.51
    Inad= -17.75

    Saab= 0.23
    Iaab= -8.11

    Smd= 0.42
    Imd= -14.38

    T = fieldset.temp[time, particle.depth, particle.lat, particle.lon]
    S = fieldset.salin[time, particle.depth, particle.lat, particle.lon]

    if T > 18.:
        particle.d18ow = St*S + It
    if T <18. and T > 4.:
        particle.d18ow = Sss*S + Iss
    if T < 4. and T > 1. and S < 35. and S > 34.8:
        particle.d18ow = Snad*S + Inad
    if T < 2. and S < 34.8 and S > 34.5:
        particle.d18ow = Saab*S + Iaab
    else:
        particle.d18ow = Smd*S + Imd
    
    particle.d18oc  = particle.d18ow - 0.27 - (T-14.9)/4.8 #Ruber equation
    particle.d18oc  = particle.d18ow - 0.27 - (T-16.5)/4.8 #Pachyderma equation
    particle.d18oc =  particle.d18ow - 0.27 - (T-15.2)/4.6 #Inflata and Glutinata equation
    
      
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


def DeleteParticle(particle, fieldset, time):
    particle.delete()

def run_corefootprintparticles(dirwrite,outfile):
    snapshots = snapshot_function(date(2000,1,9), date(2010, 6, 15),delta(days=3))
    # one needs to uses the precise dates of the snapshots in the ofes-data
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.add_constant('dwellingdepth',np.float(dd))
    fieldset.add_constant('sinkspeed', sp/86400)
    fieldset.add_constant('maxage', ls*86400)
    fieldset.add_constant('bottomlon', bottomlon)
    fieldset.add_constant('bottomlat', bottomlat)

    class ForamParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=fieldset.temp)
        age = Variable('age', dtype=np.float32, initial=0.)
        salin= Variable('salin', dtype=np.float32, initial=np.nan)
        d18ow = Variable('d18ow', dtype=np.float32, initial=0.)
        d18oc = Variable('d18oc', dtype=np.float32, initial=0.)
        loctemp = Variable('loctemp', dtype=np.float32, initial=0)
        d18owL = Variable('d18owL', dtype=np.float32, initial=0.)
        d18ocL = Variable('d18ocL', dtype=np.float32, initial=0.)
       

    pset = ParticleSet(fieldset=fieldset, pclass=ForamParticle, lon=corelon, lat=corelat,
                       depth=coredepth, time=fieldset.U.grid.time[-1],
                       repeatdt=delta(days=3 ))  
    pfile = ParticleFile(dirwrite + outfile, pset, outputdt=delta(days=1))  

    kernels = pset.Kernel(AdvectionRK4_3D) + Sink + SampleTemp + Age + Sampled18OLoc + LocalConditions

    pset.execute(kernels, dt=delta(minutes=-5), output_file=pfile,
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})    
    
#here one can set the name of the file to be produced
outfile = "Species_name_d18O_T_traj_and_local" + '_dd'+str(int(dd)) +'_sp'+str(int(sp)) + '_lon'+str(corelon) + '_lat' + str(corelat) + '_depth' + str(coredepth) + '_nr'+ str(snap)

run_corefootprintparticles(dirwrite,outfile)
print 'Exection finished'