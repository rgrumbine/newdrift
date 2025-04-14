import sys
from math import *
import netCDF4

buoys = netCDF4.Dataset(sys.argv[1], 'r')
nx = len(buoys.dimensions['X'])
ny = len(buoys.dimensions['Y'])

aice  = buoys.variables['ice_coverage'][0,:,:]
uvel  = buoys.variables['ice_uvelocity'][0,:,:]
vvel  = buoys.variables['ice_vvelocity'][0,:,:]

count = 0
for i in range(0,nx):
  for j in range(0,ny):
    if (uvel[j,i] < 1.e30 and vvel[j,i] < 1.e30):
      tu = uvel[j,i]
      tv = vvel[j,i]
      sp = sqrt(tu*tu + tv*tv)
      if (sp > 0):
        count += 1
        print(i,j,aice[j,i], tu, tv, sp, flush=True)

print(count,' points with valid velocities')
