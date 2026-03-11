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
    #if (dist[i] != 0 and ilat[i] < 1.e9 and dist[i] < 1.e9):
    if (ilat[i] < 1.e9 and dist[i] < 1.e9):
      count += 1
      print(i,f"{ilat[i]:.4f}", f"{ilon[i]:.4f}", f"{dist[i]:.4f}", f"{bear[i]:.4f}", f"{flat[i]:.4f}", f"{flon[i]:.4f}")

print(count,' buoys fcst')
