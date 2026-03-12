MODULE constants
  USE iso_fortran_env, only : real64, real32
  IMPLICIT none

  PUBLIC
  REAL(real64), PARAMETER :: pi  = 3.14159265358979323846_real64
  REAL(real64), PARAMETER :: rpd = pi / 180.0_real64

  REAL(kind = real64) :: nansen_ampl, nansen_rotation
  PARAMETER (nansen_ampl = 1.468e-2)
  PARAMETER (nansen_rotation = 28.0)

  REAL(kind = real64) :: kmtonm
  PARAMETER (kmtonm = 1. /  1.852 )

  REAL(real64), PARAMETER :: flag = 1.0e9_real64

END module constants

REAL(kind=real64) FUNCTION rearth(lat)
  USE constants
  IMPLICIT none
  REAL(kind=real64), intent(in) :: lat
  rearth = (6378.137_real64 - 21.385_real64*sin(lat*rpd) )
  RETURN
END function rearth

!haversine arcdis
!  http://www.movable-type.co.uk/scripts/gis-faq-5.1.html
!assumes lat lon in degrees, distance in km
REAL(kind=real64) FUNCTION harcdis(lat1, lon1, lat2, lon2)
  USE constants
  IMPLICIT none
  REAL(kind=real64), intent(in) :: lat1, lon1, lat2, lon2
  REAL(kind=real64) :: dlat, dlon, mlat
  REAL(kind=real64) :: a, c, rearth

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  mlat = (lat1 + lat2)/2.

  a = sin(dlat*rpd/2)**2 + cos(lat1*rpd)*cos(lat2*rpd)*sin(dlon*rpd/2)**2
  !equivalent c = 2.*asin(min(1._real64,sqrt(a)))
  c = 2.*atan2(sqrt(a), sqrt(1._real64 - a) )

! approximating ellipsoidal flattening RG WGS84
  harcdis = c * rearth(mlat)
  !harcdis = c * 6371.2_real64 !spherical earth
!  RETURN 
END function harcdis

!bearing
!  convert_bw (bearing to weather direction convention)
!  convert_wb (converse)

! Bearing/unbearing from movable type:
!  bearing    (convert lat-lon pair to distance/bearing)
!  unbearing  (convert distance,bearing and starting point to final point)
!dist in km, rearth in km
SUBROUTINE bearing(lat1, lon1, lat2, lon2, dist, dir)
  USE constants
  IMPLICIT none
  REAL(kind=real64), intent(in) :: lat1, lon1, lat2, lon2
  REAL(kind=real64), intent(out) :: dist, dir
  REAL(kind=real64) :: harcdis

  IF (lat1 >= flag .or. lon1 >= flag .or. lat2 >= flag .or. lon2 >= flag) THEN
    dist = flag
    dir  = flag
    RETURN
  ENDIF

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
  REAL(kind=real64), intent(in) :: lat1, lon1, dist, dir
  REAL(kind=real64), intent(out) :: lat2, lon2
  REAL(kind=real64) :: theta, R
  REAL(kind=real64) :: rearth

  theta = dir*rpd
  R = rearth(lat1)
  lat2 = asin( sin(lat1*rpd)*cos(dist/R) + &
               cos(lat1*rpd)*sin(dist/R)*cos(dir*rpd) ) / rpd
  lon2 = lon1 + atan2(sin(theta)*sin(dist/R)*cos(lat1*rpd) ,  &
                      cos(dist/R)-sin(lat1*rpd)*sin(lat2*rpd) ) / rpd
RETURN
END subroutine unbearing

!Wrap longitude to be in range [0,360]
REAL(kind=real64) FUNCTION wrap(y)
  USE constants
  IMPLICIT none
  REAL(kind=real64), intent(in) :: y
  REAL(kind=real64) :: x

  x = y
  IF (x > 360.) x = x - 360.
  IF (x > 360.) x = x - 360.
  IF (x > 360.) x = x - 360.
  IF (x < 0) x = x + 360.
  !debug: PRINT *,'wrap ',y,x
  wrap = x
  RETURN
END FUNCTION wrap


! local_metric now in class metric
! local_cartesian now in class metric

! Convert from buoy's lat-lon location to its ij coordinate (x,y in buoy member)
! --> in class metric

! ---------- To add ----------
!wdir
!grid rotation
!Bilinear interpolation of a field to buoy.x,y

