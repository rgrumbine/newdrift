'''
Take an .nc buoy file and rewrite with new variable labels
Args: skiles_pts.nc outname.nc
Robert Grumbine 8 September 2025
'''

import sys
import datetime

import numpy as np
import netCDF4 as nc

from buoy import buoy

#------------------------------------------------------------
collection = []


#with skiles: ----------------------------------
skiles = nc.Dataset(sys.argv[1], 'r')
nbuoy = skiles.dimensions['nbouy'].size
#debug: print("nbuoy = ",nbuoy, flush=True)
lat   = skiles.variables['initial_latitude'][:]
lon   = skiles.variables['initial_longitude'][:]
for i in range(0,nbuoy):
  #debug: print("in skiles ",i,flush=True)
  tmp = buoy(lat[i], lon[i])
  collection.append(tmp)
  del tmp
#debug: print("after skiles",flush=True)

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

tmp = datetime.datetime(2025,9,8)
ncfile.setncattr("date_created",tmp.strftime("%Y-%m-%d") )

#More specialized header:
ncfile.setncattr("contributor_name","Robert Grumbine")
ncfile.setncattr("contributor_email","Robert.Grumbine@noaa.gov")
ncfile.setncattr("creator_name","Robert Grumbine")
ncfile.setncattr("creator_email","Robert.Grumbine@noaa.gov")

# Buoy information --------------------------------------------
dtype = np.dtype('float32')

#For python 3.10 / netcdf 1.6.4 or later
ncfile.createVariable('Initial_Longitude', dtype, dimensions=( nbuoy ) )
ncfile.createVariable('Initial_Latitude', dtype, dimensions=( nbuoy )  )

# At last, give lat-lons of points  ----------------------------------

ncfile.variables['Initial_Longitude'][:] = tlons
ncfile.variables['Initial_Latitude'][:] = tlats

# Save file --------------------------------------------
ncfile.close()
#------------------------------------------------------------
