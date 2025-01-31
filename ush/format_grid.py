import sys
from datetime import date

import numpy as np
import netCDF4 
#-------------------------------------------------
"""
Read in the latitude and longitude from a 2ds_ice file from RTOFS
and write out a netcdf file with a 25:1 thinned set of points

Robert Grumbine
5 Nov 2024

"""

#-------------------------------------------------
fin = open(sys.argv[1],"r")

model = netCDF4.Dataset(sys.argv[1], "r")
nx = len(model.dimensions['X'])
ny = len(model.dimensions['Y'])
print (nx, ny)

tlons = model.variables['Longitude'][:,:]
tlats = model.variables['Latitude'][:,:]

print("lons ", tlons.max(), tlons.min() )
print("lats ", tlats.max(), tlats.min() )

ratio = 5
lats = np.zeros((int(nx/ratio)*int(ny/ratio)) ) 
lons = np.zeros((int(nx/ratio)*int(ny/ratio)) ) 
k = 0
for i in range(0, int(nx/ratio) ):
  for j in range(0, int(ny/ratio) ):
    #if (tlats[j*ratio, i*ratio] > 30. or tlats[j*ratio, i*ratio] < -45.):
      lats[k] = tlats[j*ratio, i*ratio]
      lons[k] = tlons[j*ratio, i*ratio]
      if (lons[k] > 720): lons[k] -= 720.
      if (lons[k] > 360): lons[k] -= 360.
      k += 1
print(k, (nx*ny)/(ratio*ratio))
npts = k

#debug: exit(0)

#-------------------------------------------------
# Open the file for output and establish its size:
ncfile = netCDF4.Dataset(sys.argv[2], mode='w', format='NETCDF4')
nbuoy  = ncfile.createDimension('nbuoy', size=npts)
#debug: 
print(nbuoy)
#debug: exit(0)

#Generic global header info:   --------------------------------------------
ncfile.title = sys.argv[2]
ncfile.setncattr("institution","NOAA/NWS/NCEP")

tmp = date.today()
ncfile.setncattr("date_created",tmp.strftime("%Y-%m-%d") )

#More specialized header:
ncfile.setncattr("contributor_name","Robert Grumbine")
ncfile.setncattr("contributor_email","Robert.Grumbine@noaa.gov")
ncfile.setncattr("creator_name","Robert Grumbine")
ncfile.setncattr("creator_email","Robert.Grumbine@noaa.gov")

# Buoy information --------------------------------------------
dtype = np.dtype('float32')

#For python 3.10 / netcdf 1.6.4
ncfile.createVariable('initial_longitude', dtype, dimensions=( nbuoy ) )
ncfile.createVariable('initial_latitude', dtype, dimensions=( nbuoy )  )

##For python 3.8 / netcdf 1.5.6 ?
#ncfile.createVariable('initial_longitude', dtype, ('nbuoy')  )
#ncfile.createVariable('initial_latitude', dtype, ('nbuoy')  )

#debug: exit(0)

# At last, give lat-lons of points  ----------------------------------

ncfile.variables['initial_longitude'][:] = lons
ncfile.variables['initial_latitude'][:] = lats

# Save file --------------------------------------------
ncfile.close()

print("output lons",lons.max(), lons.min() )
