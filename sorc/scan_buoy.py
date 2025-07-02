import sys

import netCDF4

buoys = netCDF4.Dataset(sys.argv[1], 'r')
nbuoy = len(buoys.dimensions['nbuoy'])
dist  = buoys.variables['Drift_Distance'][:]
bear  = buoys.variables['Drift_Bearing'][:]
ilat  = buoys.variables['Initial_Latitude'][:]
ilon  = buoys.variables['Initial_Longitude'][:]
flat  = buoys.variables['Final_Latitude'][:]
flon  = buoys.variables['Final_Longitude'][:] 

count = 0
for i in range(0,nbuoy):
    if (dist[i] != 0):
      count += 1
      print(i,ilat[i], ilon[i], dist[i], bear[i], flat[i], flon[i])

print(count,' buoys moved')
