# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 12:19:55 2018

@author: Anne
"""

from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D,
                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
import sys


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

##Data are located on the imau folder on gemini and results written to thesis_folder
directory = '/data2/imau/oceanparcels/hydrodynamic_data/OFESdata/OFES_0.1_HIND/DATA_3D/snap_3day'
dirwrite = '/scratch/AnneK/temp/'

#time = np.array([datetime(2010, 12, 25) - delta(days=x) for x in range(0,int(365*4+1),5)])  #what happens here: it makes amount of point in time at which particle will be released

##print maxlat, minlat, maxlon, minlon
###Lons and Lats of interest: Uruguayan margin: domain=[-20.0, -60.0,-40.0,-80.0]
#lons = [-53]   #lon=west_east
#lats = [-37]   #lat=north_south

#lons = np.array(lons*time.shape[0])
#lats = np.array(lats*time.shape[0])
#
#print 'lons: ',np.unique(lons)
#print 'lats shape: ',lats.shape
#
#sys.stdout.flush()


#%%
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

  # sfiles = [path.join(directory, 'y'+s[:4], 'salinity', "nest_1_"+s+"000000s.nc".format(s)) for s in snapshots]
    
    bfile = '/home/students/4082842/projects/Master_thesis/' + 'topography_OFES.nc'

    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles,'B':bfile} #'salin': sfiles, 'B':bfile}
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature', 'B': 'HT'} #'salin':'salinity', 'B': 'HT'}
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        #'salin':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'B':{'lat': 'LAT', 'lon': 'LONN1799_1800', 'time': 'TIME', 'depth': 'LEV'}} 
        
    indices = {'lat': range(0,950)}                 

    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=True)

    return fieldset



## In Sebille et al they scale the velocities, converting them from cm/s to m/s ...
#    fieldset.U.set_scaling_factor(0.01)  # convert from cm/s to m/s
#    fieldset.V.set_scaling_factor(0.01)  # convert from cm/s to m/s
#    fieldset.W.set_scaling_factor(0.01)  # convert from cm/s to m/s

#%%
    
def SampleTemp(particle, fieldset, time, dt):
    particle.temp = fieldset.temp[time, particle.lon, particle.lat, particle.depth]


def Sink(particle, fieldset, time, dt):
    if particle.depth > fieldset.dwellingdepth:
        particle.depth = particle.depth + fieldset.sinkspeed * dt
    else:
        particle.depth = fieldset.dwellingdepth


def Age(particle, fieldset, time, dt):
    if particle.depth <= fieldset.dwellingdepth:
        particle.age = particle.age + math.fabs(dt)
    if particle.age > fieldset.maxage:
        particle.delete()


def DeleteParticle(particle, fieldset, time, dt):
    particle.delete()

def run_corefootprintparticles(dirwrite,outfile):
    snapshots = snapshot_function(date(2009,6,14), date(2010, 6, 15),delta(days=3))
    fieldset = set_ofes_fieldset(snapshots)
    fieldset.add_periodic_halo(zonal=True)
    fieldset.B.allow_time_extrapolation = True
    fieldset.add_constant('dwellingdepth', 50.)
    fieldset.add_constant('sinkspeed', 200./86400)
    fieldset.add_constant('maxage', 30.*86400)

    corelon = [-53]
    corelat = [-37]
    coredepth = [5580]

    class ForamParticle(JITParticle):
        temp = Variable('temp', dtype=np.float32, initial=fieldset.temp)
        age = Variable('age', dtype=np.float32, initial=0.)

    pset = ParticleSet(fieldset=fieldset, pclass=ForamParticle, lon=corelon, lat=corelat,
                       depth=coredepth, time=fieldset.U.grid.time[-1],
                       repeatdt=delta(days=3 ))  # the new argument 'repeatdt' means no need to call pset.add() anymore in for-loop
    pfile = ParticleFile(outfile, pset, outputdt=delta(days=1))  # `interval` argument has changed to `outputdt`

    kernels = pset.Kernel(AdvectionRK4_3D) + Sink + SampleTemp + Age

    pset.execute(kernels, dt=delta(minutes=-5), output_file=pfile,
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile = "Newtry_loc_dwelling"
run_corefootprintparticles(dirwrite,outfile)

print 'Exection finished'
