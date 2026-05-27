# newdrift
Replacement for Grumbine, 1998 drift model -- https://github.com/NOAA-EMC/drift_grumbine

# To be runnable as
* Standalone program ('offline')
* Subroutine to CICE ('inline')

# Inputs:
* .nc file of points to track
** reformat_points.py in ush takes text file of locations and create .nc file
** skiles.nc has just the skiles points
** merged.nc has skiles + ~25 km grid based on RTOFS points
* u, v of ice 
* metric (needs input lat,lon of points, computes dlatdi, etc.)

# Outputs:
* .nc file of buoy starting point, drift distance and direction for each time lead
* kml output option -- not yet in place

# Notes:
* Parallelized
     -- handle drifter passing out of domain of a processor (inline)
     -- -- > reference processor number that handles re-partitioning, knows tiling

#Auxiliary:
* Blender of ensemble outputs
* RG: ?how to represent uncertainty

