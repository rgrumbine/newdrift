MODULE constants
  PUBLIC
  REAL pi, rpd
  PARAMETER(pi = 3.141592653589793 )
  PARAMETER(rpd = pi/180.)

  REAL nansen_ampl, nansen_rotation
  PARAMETER (nansen_ampl = 1.468e-2)
  PARAMETER (nansen_rotation = 28.0)

  REAL kmtonm
  PARAMETER (kmtonm = 1. /  1.852 )

END module constants

!haversine arcdis
!  http://www.movable-type.co.uk/scripts/gis-faq-5.1.html
!assumes lat lon in degrees, distance in km
REAL FUNCTION rearth(lat)
  USE constants
  IMPLICIT none
  REAL, intent(in) :: lat
  rearth = (6378.137 - 21.385*sin(lat*rpd) )
  RETURN
END function rearth

REAL FUNCTION harcdis(lat1, lon1, lat2, lon2)
  USE constants
  IMPLICIT none
  REAL lat1, lon1, lat2, lon2
  REAL dlat, dlon, mlat
  REAL a, c

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  mlat = (lat1 + lat2)/2.

  a = sin(dlat*rpd/2)**2 + cos(lat1*rpd)*cos(lat2*rpd)*sin(dlon*rpd/2)**2
  c = 2.*asin(min(1.,sqrt(a)))

! approximating ellipsoidal flattening RG WGS84
  harcdis = c * (6378.137 - 21.385*sin(mlat*rpd) )
  RETURN 
END function harcdis

!bearing
!  convert_bw (bearing to weather direction convention)
!  convert_wb (converse)

! Bearing/unbearing from movable type:
!  bearing    (convert lat-lon pair to distance/bearing)
!  unbearing  (convert distance,bearing and starting point to final point)
SUBROUTINE bearing(lat1, lon1, lat2, lon2, dist, dir)
  USE constants
  IMPLICIT none
  REAL, intent(in) :: lat1, lon1, lat2, lon2
  REAL, intent(out) :: dist, dir
  REAL harcdis

  dist = harcdis(lat1, lon1, lat2, lon2)
  dir  = atan2(sin((lon1-lon2)*rpd)*cos(lat2*rpd) , &
               cos(lat1*rpd)*sin(lat2*rpd) -        &
               sin(lat1*rpd)*cos(lat2*rpd)*cos((lon1-lon2)*rpd) )
  dir = dir / rpd
  IF (dir < 0.) THEN
    dir = dir + 360.
  ENDIF

RETURN
END subroutine bearing

!dist in km, rearth in km
SUBROUTINE unbearing(lat1, lon1, dist, dir, lat2, lon2)
  USE constants
  IMPLICIT none
  REAL, intent(in) :: lat1, lon1, dist, dir
  REAL, intent(out) :: lat2, lon2
  REAL theta, R, rearth

  theta = dir*rpd
  R = rearth(lat1)
  lat2 = asin( sin(lat1*rpd)*cos(dist/R) + &
               cos(lat1*rpd)*sin(dist/R)*cos(dir*rpd) ) / rpd
  lon2 = lon1 + atan2(sin(theta)*sin(dist/R)*cos(lat1*rpd) ,  &
                      cos(dist/R)-sin(lat1*rpd)*sin(lat2*rpd) ) / rpd
RETURN
END subroutine unbearing
!wdir


!grid rotation

!lat, lon in degrees
!dx, dy in meters
SUBROUTINE local_metric(ulat, ulon, dx, dy, rot, dlatdi, dlondi, dlatdj, dlondj, nx, ny)
  IMPLICIT none
  INTEGER, intent(in) :: nx, ny
  REAL, intent(in)   :: ulat(nx, ny), ulon(nx, ny)
  REAL, intent(out)  :: dx(nx, ny), dy(nx, ny), rot(nx, ny)
  REAL, intent(out)  :: dlatdi(nx, ny), dlatdj(nx, ny)
  REAL, intent(out)  :: dlondi(nx, ny), dlondj(nx, ny)
  INTEGER i, j

  rot = 0.
  CALL local_cartesian(ulat, dx, dy, nx, ny)

  DO j = 1, ny-1
  DO i = 1, nx-1
    dlatdi(i,j) = ulat(i+1,j) - ulat(i,j)
    dlondi(i,j) = ulon(i+1,j) - ulon(i,j)
    dlatdj(i,j) = ulat(i,j+1) - ulat(i,j)
    dlondj(i,j) = ulon(i,j+1) - ulon(i,j)
  ENDDO
  ENDDO
  !rot = atan2(?,?)

  RETURN
END subroutine local_metric


SUBROUTINE local_cartesian(ulat, dx, dy, nx, ny)
  USE constants
  IMPLICIT none
  INTEGER, intent(in) :: nx, ny
  REAL, intent(in)    :: ulat(nx, ny)
  REAL, intent(out)   :: dx(nx, ny), dy(nx, ny)

  INTEGER i, j

! From WGS84 via Amy Solomon, ESRL
! Meters per degree
  DO j = 1, ny
  DO i = 1, nx
    dy(i,j) = 111132.92 - 559.82*cos(2*ulat(i,j)*rpd) + &
                           1.175*cos(4*ulat(i,j)*rpd) - &
                          0.0023*cos(6*ulat(i,j)*rpd)
    dx(i,j) = 111412.84*cos(  ulat(i,j)*rpd) - &
                   93.5*cos(3*ulat(i,j)*rpd) + &
                  0.118*cos(5*ulat(i,j)*rpd)
  ENDDO
  ENDDO

  RETURN
END subroutine local_cartesian

! Convert from buoy's lat-lon location to its ij coordinate (x,y in buoy member)
SUBROUTINE ll_to_xy(lat, lon, ulat, ulon, x, y, nx, ny)

  IMPLICIT none
  INTEGER, intent(in) :: nx, ny
  REAL, intent(in)    :: lat, lon, ulat(nx, ny), ulon(nx, ny)
  REAL, intent(inout)   :: x, y

  INTEGER :: AR2(2), AR1(2)
  REAL dlon(nx, ny), dlat(nx, ny)

  dlon = ABS(ulon - lon)
  dlat = ABS(ulat - lat)
  AR1 = MINLOC(dlon)
  AR2 = MINLOC(dlat)
  PRINT *,'ar1 ar2 ',AR1, AR2
  x = FLOAT(AR1(1))
  y = FLOAT(AR1(2))
  RETURN
END SUBROUTINE ll_to_xy 


!Bilinear interpolation to buoy.x,y
