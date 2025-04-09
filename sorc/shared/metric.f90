MODULE metric_mod
    USE constants

    !module
    TYPE, public :: metric
      REAL, allocatable :: dlatdi(:,:), dlondi(:,:), dlatdj(:,:), dlondj(:,:)
      REAL, allocatable :: dx(:,:), dy(:,:), rot(:,:)
      REAL, allocatable :: ulat(:,:), ulon(:,:)
      INTEGER nx, ny
    CONTAINS
      PROCEDURE :: set, local_metric, local_cartesian, ll_to_xy, xy_to_ll
    END TYPE metric
CONTAINS

SUBROUTINE set(this, nx, ny)
  IMPLICIT none
  INTEGER nx, ny
  CLASS(metric) :: this

  this%nx = nx
  this%ny = ny

  allocate(this%dx(nx, ny))
  allocate(this%dy(nx, ny))
  allocate(this%rot(nx, ny))
  allocate(this%ulat(nx, ny))
  allocate(this%ulon(nx, ny))
  allocate(this%dlatdi(nx, ny))
  allocate(this%dlatdj(nx, ny))
  allocate(this%dlondi(nx, ny))
  allocate(this%dlondj(nx, ny))

  RETURN
END SUBROUTINE set

!lat, lon in degrees
!dx, dy in meters
SUBROUTINE local_metric(this)
  IMPLICIT none
  CLASS(metric) :: this
  INTEGER i, j

  this%rot = 0.
  CALL this%local_cartesian()

  DO j = 1, this%ny-1
  DO i = 1, this%nx-1
    this%dlatdi(i,j) = this%ulat(i+1,j) - this%ulat(i,j)
    this%dlondi(i,j) = this%ulon(i+1,j) - this%ulon(i,j)
    this%dlatdj(i,j) = this%ulat(i,j+1) - this%ulat(i,j)
    this%dlondj(i,j) = this%ulon(i,j+1) - this%ulon(i,j)
  ENDDO
  ENDDO
  !rot = atan2(?,?)

  RETURN
END subroutine local_metric


SUBROUTINE local_cartesian(this)
  USE constants
  IMPLICIT none
  CLASS(metric) :: this
  INTEGER k
  INTEGER i, j

! From WGS84 via Amy Solomon, ESRL
! Meters per degree
  DO j = 1, this%ny
  DO i = 1, this%nx
    this%dy(i,j) = 111132.92 - 559.82*cos(2*this%ulat(i,j)*rpd) + &
                           1.175*cos(4*this%ulat(i,j)*rpd) - &
                          0.0023*cos(6*this%ulat(i,j)*rpd)
    this%dx(i,j) = 111412.84*cos(  this%ulat(i,j)*rpd) - &
                   93.5*cos(3*this%ulat(i,j)*rpd) + &
                  0.118*cos(5*this%ulat(i,j)*rpd)
  ENDDO
  ENDDO

  RETURN
END subroutine local_cartesian

! Convert from buoy's lat-lon location to its ij coordinate (x,y in buoy member)
SUBROUTINE ll_to_xy(this, lat, lon, x, y)

  IMPLICIT none
  CLASS(metric) :: this
  REAL, intent(in)    :: lat, lon
  REAL, intent(inout)   :: x, y

  INTEGER :: AR2(2), AR1(2)
  REAL dlon(this%nx, this%ny), dlat(this%nx, this%ny)

  dlon = ABS(this%ulon - lon)
  dlat = ABS(this%ulat - lat)
  AR1 = MINLOC(dlon)
  AR2 = MINLOC(dlat)
  PRINT *,'ar1 ar2 ',AR1, AR2
  x = FLOAT(AR1(1))
  y = FLOAT(AR1(2))
  RETURN
END SUBROUTINE ll_to_xy 

SUBROUTINE xy_to_ll(this, lat, lon, x, y)

  IMPLICIT none
  CLASS(metric) :: this
  REAL, intent(in)    :: lat, lon
  REAL, intent(inout)   :: x, y

  INTEGER :: AR2(2), AR1(2)
  REAL dlon(this%nx, this%ny), dlat(this%nx, this%ny)

  dlon = ABS(this%ulon - lon)
  dlat = ABS(this%ulat - lat)
  AR1 = MINLOC(dlon)
  AR2 = MINLOC(dlat)
  PRINT *,'ar1 ar2 ',AR1, AR2
  x = FLOAT(AR1(1))
  y = FLOAT(AR1(2))
  RETURN
END SUBROUTINE xy_to_ll 

END MODULE metric_mod

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
SUBROUTINE elephant(lat, lon, x, y)
    IMPLICIT none
    REAL, intent(in) :: lat, lon
    REAL, intent(inout) :: x, y
    PRINT *,'elephant',lat,lon,x,y
END SUBROUTINE elephant
