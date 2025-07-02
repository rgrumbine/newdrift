import sys
from math import *
import numpy as np

import netCDF4

buoys = netCDF4.Dataset(sys.argv[1], 'r')
nbuoy = len(buoys.dimensions['nbuoy'])
dist  = buoys.variables['Drift_Distance'][:]
bear  = buoys.variables['Drift_Bearing'][:]
ilat  = buoys.variables['Initial_Latitude'][:]
ilon  = buoys.variables['Initial_Longitude'][:]
flat  = buoys.variables['Final_Latitude'][:]
flon  = buoys.variables['Final_Longitude'][:] 

init_lat = np.zeros((nbuoy))
init_lon = np.zeros((nbuoy))
plot_bear  = np.zeros((nbuoy))
plot_dist  = np.zeros((nbuoy))

count = 0
for i in range(0,nbuoy):
    if (dist[i] != 0 and flat[i] < 1.e30 and flon[i] < 1.e30):
      init_lat[count] = ilat[i]
      init_lon[count] = ilon[i]
      plot_bear[count] = bear[i]
      plot_dist[count] = dist[i]
      count += 1
      #print(i,ilat[i], ilon[i], dist[i], bear[i], flat[i], flon[i])
print(count,' buoys moved')

#convert to radians
plot_bear *= pi/180.

#RG:  check conventions on direction
u = plot_dist * np.cos(plot_bear)
v = plot_dist * np.sin(plot_bear)

import matplotlib
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

matplotlib.use('Agg')

# NWP and then some
proj = ccrs.LambertConformal(central_longitude = -105,
                             central_latitude = 75., cutoff = 45.)
ax  = plt.axes(projection = proj)
fig = plt.figure(figsize=(11,8.5))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.set_extent((-170, -75, 60, 90), crs=ccrs.PlateCarree())

#ax.gridlines(crs=ccrs.PlateCarree(),
#      xlocs=[-225, -210, -195, -180, -165, -150, -135., -120, -105, -90,
#              -75, -60, -45, -30, -15, 0, 15],
#      ylocs=[60, 65, 70, 75, 80, 85, 90] )
ax.gridlines(crs=ccrs.PlateCarree() )
ax.add_feature(cfeature.GSHHSFeature(levels = [1,2,3,4], scale = "low") )

ax.quiver(init_lon[:count], init_lat[:count], u[:count], v[:count], transform = ccrs.PlateCarree())
fig.savefig("nwp.png")
plt.close()

# Hudson Bay and then some
proj = ccrs.LambertConformal(central_longitude = -85,
                             central_latitude = 55., cutoff = 40.)
ax  = plt.axes(projection = proj)
fig = plt.figure(figsize=(11,8.5))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.set_extent((-95, -65, 50, 70), crs=ccrs.PlateCarree())

ax.gridlines(crs=ccrs.PlateCarree() )
ax.add_feature(cfeature.GSHHSFeature(levels = [1,2,3,4], scale = "low") )

ax.quiver(init_lon[:count], init_lat[:count], u[:count], v[:count], transform = ccrs.PlateCarree())
fig.savefig("hudson.png")
plt.close()

# Seam around 255 longitude, -105
proj = ccrs.LambertConformal(central_longitude = -105,
                             central_latitude = 75., cutoff = 45.)
ax  = plt.axes(projection = proj)
fig = plt.figure(figsize=(11,8.5))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.set_extent((-145, -65, 70, 90), crs=ccrs.PlateCarree())

#ax.gridlines(crs=ccrs.PlateCarree(),
#      xlocs=[-225, -210, -195, -180, -165, -150, -135., -120, -105, -90,
#              -75, -60, -45, -30, -15, 0, 15],
#      ylocs=[60, 65, 70, 75, 80, 85, 90] )
ax.gridlines(crs=ccrs.PlateCarree() )
ax.add_feature(cfeature.GSHHSFeature(levels = [1,2,3,4], scale = "low") )

ax.quiver(init_lon[:count], init_lat[:count], u[:count], v[:count], transform = ccrs.PlateCarree())
fig.savefig("seam1.png")
plt.close()


# Seam around 75 longitude
proj = ccrs.LambertConformal(central_longitude = 75., 
                             central_latitude = 75., cutoff = 45.)
ax  = plt.axes(projection = proj)
fig = plt.figure(figsize=(11,8.5))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.set_extent((40, 110, 70, 90), crs=ccrs.PlateCarree())

ax.gridlines(crs=ccrs.PlateCarree() )
ax.add_feature(cfeature.GSHHSFeature(levels = [1,2,3,4], scale = "low") )

ax.quiver(init_lon[:count], init_lat[:count], u[:count], v[:count], transform = ccrs.PlateCarree())
fig.savefig("seam2.png")
plt.close()


# Antarctica -- Weddell
#proj = ccrs.LambertConformal(central_longitude = -30.,
#                             central_latitude = -65., cutoff = -45.)
proj = ccrs.PlateCarree()
ax  = plt.axes(projection = proj)
fig = plt.figure(figsize=(11,8.5))
ax  = fig.add_subplot(1, 1, 1, projection = proj)
ax.set_extent((-65, -10, -55, -80), crs=ccrs.PlateCarree())

ax.gridlines(crs=ccrs.PlateCarree() )
ax.add_feature(cfeature.GSHHSFeature(levels = [1,2,3,4], scale = "low") )

ax.quiver(init_lon[:count], init_lat[:count], u[:count], v[:count], transform = ccrs.PlateCarree())
fig.savefig("weddell.png")
plt.close()
