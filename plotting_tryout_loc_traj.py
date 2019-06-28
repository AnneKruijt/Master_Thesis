# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 17:38:18 2018

@author: Anne
"""
#from parcels import (FieldSet, ParticleSet, JITParticle, AdvectionRK4_3D,
#                     ErrorCode, ParticleFile, Variable)
from datetime import timedelta as delta
from datetime import date, datetime
import numpy as np
from os import path
import math
import sys
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
    
def make_plot(trajfile):
    from netCDF4 import Dataset
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap

    class ParticleData(object):
        def __init__(self):
            self.id = []

    def load_particles_file(fname, varnames):
        T = ParticleData()
        pfile = Dataset(fname, 'r')
        T.id = pfile.variables['trajectory'][:]   #not sure what is happening here
        for v in varnames:
            setattr(T, v, pfile.variables[v][:])
        return T

    T = load_particles_file(trajfile, ['lon', 'lat', 'temp', 'z'])
    m = Basemap(projection='merc', llcrnrlat=-65, urcrnrlat=-10.5, llcrnrlon=-70, urcrnrlon=-20, resolution='h')
    m.drawcoastlines()
    m.fillcontinents(color='burlywood')
    m.drawparallels(np.arange(-50, -20, 10), labels=[True, False, False, False])
    m.drawmeridians(np.arange(0, 40, 10), labels=[False, False, False, True])

    sinks = np.where(T.z > 50.)
    dwell = np.where(T.z == 50.)
    xs, ys = m(T.lon[dwell], T.lat[dwell])
    m.scatter(xs, ys, c=T.temp[dwell], s=5)
    cbar = plt.colorbar()
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.set_xlabel('[$^\circ$C]')

    xs, ys = m(T.lon[sinks], T.lat[sinks])
    m.scatter(xs, ys, c='k', s=5)
    xs, ys = m(T.lon[0, 0], T.lat[0, 0])
    m.plot(xs, ys, 'om')
    plt.show()    
    
    # I would only need  latitudes from -65 to -5 N and longitudes from -70 to 0 E
    # search the indices belonging  to these coordinates to load the appropriate field


##---------Importing data needed for plotting-------------##
  
#outfile = '/home/students/4082842/projects/Master_thesis/results/loc_dd50_sp200_lon-56_lat-40.nc' #(zonder nc !!?)
outfile = 'C:/Users/Anne/Desktop/loc_dd50_sp200_lon-56_lat-40.nc' #(zonder nc !!?)
#convert_IndexedOutputToArray(file_in=outfile+".nc", file_out=outfile+"_array.nc")
make_plot(outfile)