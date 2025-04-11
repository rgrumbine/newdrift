# newdrift
Replacement for Grumbine, 1998 drift model -- https://github.com/NOAA-EMC/drift_grumbine

# To be runnable as
* Standalone program ('offline')
* Subroutine to CICE ('inline')

# Inputs:
* .nc file of points to track
** reformat_points.py in ush to take text file of locations and create .nc file
** skiles.nc has just the skiles points
** merged.nc has skiles + ~25 km grid based on RTOFS points
* u, v of ice -- RG: rotation? grid?
* metric (needs input lat,lon of points, computes dlatdi, etc.; add rotation info here?)

# Outputs:
* .nc file of buoy starting point, drift distance and direction w.r.t. time lead
* kml output option

# Notes:
* Parallelized
     -- handle drifter passing out of domain of a processor (inline)
     -- -- > reference processor number that handles re-partitioning, knows tiling


#Auxiliary:
* Blender of ensemble outputs
* RG: ?how to represent uncertainty
* Sidfex-compatible output?



