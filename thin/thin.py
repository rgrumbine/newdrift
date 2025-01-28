import sys

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

  def mindist(self, buoylist):
    dists = np.zeros((len(buoylist)))
    for i in range(0, len(buoylist)):
      dists[i] = harcdis(self.latitude, self.longitude, 
                 buoylist[i].latitude, buoylist[i].longitude)
    dmin = dists.min()
    del dists
    return dmin

  def may_be_ice(self, posteriori):
    return True
#    i = fn(self.longitude, vs. 1/12th grid)
#    j = fn(self.latitude   vs  1/12th grid)
#    if (posteriori[j,i] == 165):
#      return True
#    else:

#------------------------------------------------------------
collection = []

#read in skip grid (fixed fields.nc)
skip = nc.Dataset(sys.argv[1], 'r')
posteriori = skip.variables['posteriori'][:,:]
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
  if (i%2000 == 0): print('in fullgrid',i,'nbuoy_tmp',len(collection),flush=True)
  tmp = buoy(lat[i], lon[i])
  if (tmp.may_be_ice(posteriori)):
    if (tmp.mindist(collection) > toler):
      collection.append(tmp)
  del tmp

print("found",len(collection),"buoy points")
  
#------------------------------------------------------------
# Write out rtofs full buoy list
#nbuoy
#write collection(lat,lon,k)

#------------------------------------------------------------
