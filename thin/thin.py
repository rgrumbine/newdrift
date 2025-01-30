import sys
from math import *

import numpy as np
import netCDF4 as nc

def harcdis(lat1, lon1, lat2, lon2):
  return 25000.

class buoy:
  def __init__(self, lat = 0., lon = 0.):
    self.latitude = lat
    self.longitude = lon

  def distance(self, buoy2):
    return harcdis(self.latitude, self.longitude, buoy2.latitude, buoy2.longitude)

  def mindist(self, buoylist, toler):
    #Return True if mindist > toler
    #RG: can do a thinning test in latitude -- if greater than (toler/111.1e3) not near
    #RG: only need compute distances up to point of finding one closer than toler
    #RG: 
    nb = len(buoylist)
    delta_lat = toler/111.1e3
    retval = True
    for i in range(-1,-nb,-1):
      if (abs(buoylist[i].latitude - self.latitude) < delta_lat):
          return False
    return retval 
#    #dists = np.zeros((nb-max(0,nb-25000) ))
#    dists = np.zeros((min(nb,25000) ))
#    dists[0] = 1.e8
#    #for i in range(max(0,nb-25000), nb):
#    for i in range(1, min(nb,25000)):
#      dists[i] = harcdis(self.latitude, self.longitude, 
#                 buoylist[-i].latitude, buoylist[-i].longitude)
#    dmin = dists.min()
#    del dists
#    return dmin

  def may_be_ice(self, posteriori):
    if (self.latitude > -40. and self.latitude < 30.): return False
    delta = 1./12.
    #debug: print(self.longitude, self.latitude, end="")
    tlon = self.longitude
    if    tlon < 0:     tlon += 360.
    #while tlon >= 360.: tlon -= 360.
    if (tlon >= 360.):
        if (tlon >= 720):
            tlon -= 720.
        else:
            tlon -= 360.
    i = int( (tlon - delta/2.)/delta + 0.5)
    j = int( (90. - delta/2. - self.latitude)/delta + 0.5)
    #debug: if (i < 0 or j < 0 or i >= 4320 or j >= 2160 ):
    #debug:   print(" oops",tlon, self.latitude, i, j)
    #debug:   return False
    #debug: print(' i',i, 'j',j, " p ",posteriori[j,i], flush=True)
    if (posteriori[j,i] == 165):
      return True
    else:
      return False

#------------------------------------------------------------
collection = []

#read in skip grid (fixed fields.nc)
skip = nc.Dataset(sys.argv[1], 'r')
posteriori = skip.variables['posteriori'][:,:]
print(type(posteriori))
print("posteriori",posteriori.max(), posteriori.min(), posteriori.mean() )

#with skiles: ----------------------------------
skiles = nc.Dataset(sys.argv[2], 'r')
nbuoy = skiles.dimensions['nbouy'].size
print("nbuoy = ",nbuoy, flush=True)
lat   = skiles.variables['initial_latitude'][:]
lon   = skiles.variables['initial_longitude'][:]
for i in range(0,nbuoy):
  #debug: print("in skiles ",i,flush=True)
  tmp = buoy(lat[i], lon[i])
  collection.append(tmp)
  del tmp
print("after skiles",flush=True)

#with fullgrid: ----------------------------------
fullgrid = nc.Dataset(sys.argv[3], 'r')
nbuoy = fullgrid.dimensions['nbuoy'].size
lat   = fullgrid.variables['Initial_Latitude'][:]
lon   = fullgrid.variables['Initial_Longitude'][:]
toler = 20e3
for i in range(0, nbuoy):
  #debug:
  if ((i % 10000) == 0): print('in fullgrid',i,'nbuoy_tmp',len(collection),flush=True)
  tmp = buoy(lat[i], lon[i])
  if (tmp.may_be_ice(posteriori)):
    #if (tmp.mindist(collection) > toler):
    if (tmp.mindist(collection, toler)): #true if mindist > toler
      collection.append(tmp)
      #debug: print(len(collection),'buoys')
  del tmp
  #if (i == 750000): break

  
#------------------------------------------------------------
# Write out rtofs full buoy list
#write collection(lat,lon,k)
print("found",len(collection),"buoy points")
for i in range(0, len(collection)):
    print(i,collection[i].longitude, collection[i].latitude)

#------------------------------------------------------------
