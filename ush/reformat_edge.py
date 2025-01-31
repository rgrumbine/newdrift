import sys
from datetime import date

import numpy as np
import netCDF4 as nc
#-------------------------------------------------
"""
Read in a text file of lat-lon points, optionally with names (for future)
and write out a netcdf file with those as initial latitude-longitudes for
lagrangean drift guidance

Robert Grumbine
26 July 2023

"""

#-------------------------------------------------
fin = open(sys.argv[1],"r")

lats = []
lons = []
names = []
for line in fin:
  words = line.split()
  lats.append(float(words[0]) )
  lons.append(float(words[1]) )
  if (len(words) > 2):
    names.append(words[2])
  else:
    names.append("null")

npts = len(lats)
#debug: print("npts = ",npts)
#debug: for i in range(0,npts):
  #debug: print(i,lats[i],lons[i],names[i])

#-------------------------------------------------
# Open the file for output and establish its side:
ncfile = nc.Dataset(sys.argv[2], mode='w', format='NETCDF4')
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

ncfile.variables['initial_longitude'][:] = lons
ncfile.variables['initial_latitude'][:] = lats

# Save file --------------------------------------------
ncfile.close()
