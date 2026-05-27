'''
Take a file with full listing of lat-lon of points on model grid and thin to ~25 km
Args:
    seaice_fixed_fields.nc skiles.nc fullgrid.nc outfile.nc
Robert Grumbine 8 April 2025

'''
import sys
import datetime

import numpy as np
import netCDF4 as nc

from buoy import buoy

#------------------------------------------------------------
collection = []

#read in posteriori grid (seaice_fixed_fields.nc)
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
lon[lon < 0 ] += 360.
lon[lon >= 720 ] -= 720.
lon[lon >= 360 ] -= 360.

toler = 25e3
for i in range(0, nbuoy):
  #debug: if ((i % 10000) == 0): print('in fullgrid',i/1000,'k nbuoy_tmp',len(collection),flush=True)
  tmp = buoy(lat[i], lon[i])
  if (tmp.may_be_ice(posteriori)):
    #if (tmp.mindist(collection) > toler):
    if (tmp.mindist(collection, toler)): #true if mindist > toler
      collection.append(tmp)
      #debug: print(len(collection),'buoys')
  del tmp
  #if (i == 920000): break

#------------------------------------------------------------
nb = len(collection)
tlats = np.zeros((nb))
tlons = np.zeros((nb))
tid   = np.linspace(0,nb,1)

# Write out rtofs full buoy list
#write collection(lat,lon,k)
print("found",nb,"buoy points")
for i in range(0, nb):
    tlats[i] = collection[i].latitude
    tlons[i] = collection[i].longitude
    print(collection[i].longitude, collection[i].latitude, i)

# Open the file for output and establish its size:
ncfile = nc.Dataset(sys.argv[4], mode='w', format='NETCDF4')
nbuoy  = ncfile.createDimension('nbuoy', size=nb)

#Generic global header info:
ncfile.title = sys.argv[4]
ncfile.setncattr("institution","NOAA/NWS/NCEP/EMC")

tmp = datetime.datetime(2025,1,31)
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
