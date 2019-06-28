# -*- coding: utf-8 -*-
"""

Script combing and adding to code from Dr. E. v Sebille and P. Nooteboom, to advect particles both while sinking and while dwelling, at a specific core location at the Uruguayan margin. It uses the python tool Parcels, available on: github.com/OceanParcels/parcels

The hydrodynamic ﬁelds that carry the foraminifera come from the OFES model (Masumoto et al., 2004) and can be accessed from http://apdrc.soest.hawaii.edu/datadoc/ofes/ncep_0.1_global_3day.php 

- the indexed area of loaded field can be varied (smaller makes script run faster)
- the amount of snapshots loaded can be varied (less makes script run faster)

@author: Anne Kruijt
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
dirwrite ='/scratch/AnneK/temp/'


## Characteristics of foraminifera
sp = 200. #The sinkspeed m/day
dd = 50. #The dwelling depth
ls = 30.0 # life span while dwelling, in days

# location and depth of the particle release site
###Lons and Lats of interest: Uruguayan margin: domain=[-20.0, -60.0,-40.0,-80.0]
corelon = [-52.5.] 
corelat = [-36.5.] 
coredepth = [2578.8]
bottomlon = -52.5
bottomlat = -36.5

# amount of years of snapshots to be loaded
snap = 10.
print 'did you set the amount of years correctly?'

#%%
''' function that loads the relevant velocity and temperature fields and creates a fieldset'''
def set_ofes_fieldset(snapshots):
    ufiles = [path.join(directory, 'y'+s[:4], 'u_vel', "nest_1_"+s+"000000u.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    vfiles = [path.join(directory, 'y'+s[:4], 'v_vel', "nest_1_"+s+"000000v.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103
    wfiles = [path.join(directory, 'y'+s[:4], 'w_vel', "nest_1_"+s+"000000w.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103   
    tfiles = [path.join(directory, 'y'+s[:4], 'temp', "nest_1_"+s+"000000t.nc".format(s)) for s in snapshots]#path.dirname(__file__)#0103    

    filenames = {'U': ufiles, 'V': vfiles, 'W': wfiles, 'temp': tfiles} 
    variables = {'U': 'zu', 'V': 'zv', 'W': 'zw', 'temp': 'temperature'} 
    
    dimensions = {      'U':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'V':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'W':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'},
                        'temp':{'lat': 'Latitude', 'lon': 'Longitude', 'time': 'Time', 'depth': 'Depth'}}
                       
    ## Different indices can be picked for prefered size of loaded field                 
    indices = {'lat': range(100,700), 'lon': range(2900,3600)} #latitude ranges from 0 to 1500 and longitude from 0 to 3600
    #indices = {'lat': range(0,950)} 
    fieldset = FieldSet.from_netcdf(filenames, variables, dimensions, indices = indices, allow_time_extrapolation=True)

    return fieldset


#%%
#for storing the temperature at the particle's location  
def SampleTemp(particle, fieldset, time):
    particle.temp = fieldset.temp[time, particle.depth, particle.lat, particle.lon]

#for storing the temperature right above the core site but at the particle's depth  
def LocalConditions(particle, fieldset, time):
    lonbottom = fieldset.bottomlon
    latbottom = fieldset.bottomlat
    particle.loctemp = fieldset.temp[time, particle.depth, latbottom, lonbottom]

#for letting the particle sink as long as it is located between the release site and its dwelling depth near the surface
def Sink(particle, fieldset, time):
    if particle.depth > fieldset.dwellingdepth:
        particle.depth = particle.depth + fieldset.sinkspeed * dt
    else:
        particle.depth = fieldset.dwellingdepth

#for computing the particles from when it has reached its dwelling depth until it has reached its maximum age
def Age(particle, fieldset, time):
    if particle.depth <= fieldset.dwellingdepth:
        particle.age = particle.age + math.fabs(dt)
    if particle.age > fieldset.maxage:
        particle.delete()


def DeleteParticle(particle, fieldset, time):
    particle.delete()

#for computing the distance traveled along the particle trajectory
def TotalDistance(particle, fieldset, time):
    # Calculate the distance in latitudinal direction (using 1.11e2 kilometer per degree latitude)
    lat_dist = (particle.lat - particle.prev_lat) * 1.11e2
    # Calculate the distance in longitudinal direction, using cosine(latitude) - spherical earth
    lon_dist = (particle.lon - particle.prev_lon) * 1.11e2 * math.cos(particle.lat * math.pi / 180)
    # Calculate the total Euclidean distance travelled by the particle
    particle.distance += math.sqrt(math.pow(lon_dist, 2) + math.pow(lat_dist, 2))

    particle.prev_lon = particle.lon  # Set the stored values for next iteration.
    particle.prev_lat = particle.lat
    
#for executing the script   
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
        loctemp = Variable('loctemp', dtype=np.float32, initial=0)
        distance = Variable('distance', initial=0., dtype=np.float32)  # the distance travelled
        prev_lon = Variable('prev_lon', dtype=np.float32, to_write=False,
                            initial=attrgetter('lon'))  # the previous longitude
        prev_lat = Variable('prev_lat', dtype=np.float32, to_write=False,
                            initial=attrgetter('lat'))  # the previous latitude.     

    pset = ParticleSet(fieldset=fieldset, pclass=ForamParticle, lon=corelon, lat=corelat,
                       depth=coredepth, time=fieldset.U.grid.time[-1],
                       repeatdt=delta(days=3 ))  
    pfile = ParticleFile(dirwrite + outfile, pset, outputdt=delta(days=1))  

    kernels = pset.Kernel(AdvectionRK4_3D) + Sink + SampleTemp + Age + LocalConditions + TotalDistance

    pset.execute(kernels, dt=delta(minutes=-5), output_file=pfile,
                 recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})


#here one can set the name of the file to be produced
outfile = "Sink_Dwell_pluslocal" + '_dd'+str(int(dd)) +'_sp'+str(int(sp)) + '_lon'+str(corelon) + '_lat' + str(corelat) + '_depth' + str(coredepth) + '_nr'+ str(snap)

run_corefootprintparticles(dirwrite,outfile)
print 'Exection finished'
