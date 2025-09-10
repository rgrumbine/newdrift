import sys
from math import *

import numpy as np
import netCDF4 as nc

'''
Take an .nc buoy file and rewrite with new variable labels
Args: skiles_pts.nc outname.nc
Robert Grumbine 8 September 2025

'''

class buoy:
  def __init__(self, lat = 0., lon = 0.):
    self.latitude = lat
    self.longitude = lon

  def mindist(self, buoylist, toler):
    '''
      Return True if no buoy in list is within toler of the trial buoy (self) 
    '''
    #RG: can do a thinning test in latitude -- if greater than (toler/111.2e3) not near
    #RG: only need compute distances up to point of finding one closer than toler
  
    nb = len(buoylist)
    delta_lat = toler/111.2e3

    for i in range(-1,-nb,-1):
      if (abs(buoylist[i].latitude - self.latitude) < delta_lat):
        d = harcdis(self.latitude, self.longitude, buoylist[i].latitude, buoylist[i].longitude)
        if (d < toler):
          return False
    return True 

  def may_be_ice(self, posteriori):
    if (self.latitude > -40. and self.latitude < 30.): return False
    delta = 1./12.
    #debug: print(self.longitude, self.latitude, end="")
    tlon = self.longitude
    #if    tlon < 0:     tlon += 360.
    ##while tlon >= 360.: tlon -= 360.
    #if (tlon >= 360.):
    #    if (tlon >= 720):
    #        tlon -= 720.
    #    else:
    #        tlon -= 360.
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


#with skiles: ----------------------------------
skiles = nc.Dataset(sys.argv[1], 'r')
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
tlats=np.zeros((nbuoy))
tlons=np.zeros((nbuoy))
print("found",nbuoy,"buoy points")
for i in range(0, nbuoy):
    tlats[i] = collection[i].latitude
    tlons[i] = collection[i].longitude
    print(collection[i].longitude, collection[i].latitude, i)

# Open the file for output and establish its size:
ncfile = nc.Dataset(sys.argv[2], mode='w', format='NETCDF4')
nbuoy  = ncfile.createDimension('nbuoy', size=nbuoy)

#Generic global header info:   
ncfile.title = sys.argv[2]
ncfile.setncattr("Institution","NOAA/NWS/NCEP/MDC")

import datetime
tmp = datetime.datetime(2025,9,8)
ncfile.setncattr("date_created",tmp.strftime("%Y-%m-%d") )

#More specialized header:
ncfile.setncattr("contributor_name","Robert Grumbine")
ncfile.setncattr("contributor_email","Robert.Grumbine@noaa.gov")
ncfile.setncattr("creator_name","Robert Grumbine")
ncfile.setncattr("creator_email","Robert.Grumbine@noaa.gov")

# Buoy information --------------------------------------------
dtype = np.dtype('float32')

#For python 3.10 / netcdf 1.6.4
ncfile.createVariable('Initial_Longitude', dtype, dimensions=( nbuoy ) )
ncfile.createVariable('Initial_Latitude', dtype, dimensions=( nbuoy )  )

# At last, give lat-lons of points  ----------------------------------

ncfile.variables['Initial_Longitude'][:] = tlons
ncfile.variables['Initial_Latitude'][:] = tlats

# Save file --------------------------------------------
ncfile.close()

#------------------------------------------------------------
