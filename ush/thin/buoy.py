''' simple buoy class for buoys that only have location 
Robert Grumbine 8 September 2025
'''
from math import sin, cos, pi, atan2, sqrt

def arcdis(lat1, lon1, lat2, lon2):
  '''
    Return haversine distance between lat1,lon1, lat2, lon2 -- in meters
  '''
  # Radius of the Earth in meters
  earth_radius = 6371.e3

  # Convert latitude and longitude from degrees to radians
  tlat1 = pi/180.*(lat1)
  tlat2 = pi/180.*(lat2)
  #tlon1 = pi/180.*(lon1) #tlon1, tlon2 only used in finding dlon
  #tlon2 = pi/180.*(lon2)

  # Haversine formula
  dlon = pi/180.*(lon2 - lon1)
  dlat = tlat2 - tlat1
  a = sin(dlat / 2) ** 2 + cos(tlat1) * cos(tlat2) * sin(dlon / 2) ** 2
  c = 2 * atan2(sqrt(a), sqrt(1 - a))

  # Calculate the distance
  distance = earth_radius * c

  return distance

#-------------------------------------------------------------------------------------
class buoy:
  ''' simple buoy class for buoys that only have location '''

  def __init__(self, flat = 0., flon = 0.):
    ''' initialize buoy location to 0,0 '''
    self.latitude = flat
    self.longitude = flon

  def mindist(self, buoylist, toler):
    '''
      Return True if no buoy in list is within toler of the trial buoy (self) 
    '''
    #RG: can do a thinning test in latitude -- if greater than (toler/111.2e3) not near
    #RG: only need compute distances up to point of finding one closer than toler

    nb = len(buoylist)
    delta_lat = toler/111.2e3

    for fi in range(-1,-nb,-1):
      if (abs(buoylist[fi].latitude - self.latitude) < delta_lat):
        d = arcdis(self.latitude, self.longitude, buoylist[fi].latitude, buoylist[fi].longitude)
        if (d < toler):
          return False
    return True

  def may_be_ice(self, posteriori):
    ''' may_be_ice(posteriori) -- posteriori is the 1/12th degree posteriori ice mask ''' 
    if (self.latitude > -40. and self.latitude < 30.):
        return False
    delta = 1./12.
    tlon = self.longitude
    fi = int( (tlon - delta/2.)/delta + 0.5)
    fj = int( (90. - delta/2. - self.latitude)/delta + 0.5)
    if (posteriori[fj,fi] == 165):
      return True
    return False
