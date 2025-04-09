! Module/class for drifting bodies
! expects a set of lat-lon locations and a velocity field with velocity 
!    grid type (B/C/...) and lat-longs of velocity points.
!    also time step to extrapolate over
MODULE drifter_mod
  USE metric_mod

  IMPLICIT none
  TYPE, public :: drifter
    REAL x, y        ! current i,j location
    REAL ilat, ilon  ! initial latitude-longitude
    REAL clat, clon  ! current latitude-longitude
  CONTAINS
    PROCEDURE, pass :: init, move, zero
  END TYPE drifter


CONTAINS
  SUBROUTINE init(buoy, tlon, tlat, xmetric)
    USE metric_mod
    IMPLICIT none
    REAL, intent(in) :: tlon, tlat
    TYPE(metric) :: xmetric
    REAL x, y
    CLASS(drifter), intent(inout) :: buoy

    buoy%ilat = tlat
    buoy%ilon = tlon
    buoy%clat = tlat
    buoy%clon = tlon

    CALL xmetric%ll_to_xy(tlat, tlon, x, y)
    buoy%x = x
    buoy%y = y

    RETURN
  END SUBROUTINE init

  SUBROUTINE zero(buoy)
    IMPLICIT none
    CLASS(drifter), intent(inout) :: buoy
    buoy%x = 0.
    buoy%y = 0.
    buoy%ilat = 0.
    buoy%ilon = 0.
    buoy%clat = 0.
    buoy%clon = 0.
    RETURN
  END SUBROUTINE zero 

  SUBROUTINE move(buoy, u, v, dx, dy, dt, nx, ny)
!note dx, dy are the mesh variables
    IMPLICIT none

    CLASS(drifter), intent(inout) :: buoy
    INTEGER, intent(in) :: nx, ny
    REAL, intent(in) ::  u(nx, ny), v(nx, ny)
    REAL, intent(in) :: dx(nx, ny), dy(nx, ny)
    REAL, intent(in) :: dt

    REAL deltax, deltay
    INTEGER ti, tj

    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    !RG:  These could be interpolated (bilinear, ...)
    deltax = u(ti, tj) * dt
    deltay = v(ti, tj) * dt

!where deltax < dx, deltax > -x-int(x):
    !IF (deltax < dx) THEN               !note, deltas could be negative
    buoy%x = buoy%x + deltax/dx(ti, tj)

!where deltay < dy:
    buoy%y = buoy%y + deltay/dy(ti, tj)

  RETURN
  END subroutine move

!----------------------------------------------------------------
SUBROUTINE run(buoys, nbuoy, u, v, xmetric, dt, dtout)
  USE metric_mod
  IMPLICIT none

  TYPE(metric), intent(in) :: xmetric
  INTEGER, intent(in) :: nbuoy
  REAL, intent(in) :: u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
  REAL, intent(in) :: dt, dtout

  TYPE(drifter), intent(inout) ::  buoys(nbuoy)

  INTEGER k

  DO k = 1, nbuoy
    !c-like (object-like) 
    CALL buoys(k)%move(u, v, xmetric%dx, xmetric%dy, dt, xmetric%nx, xmetric%ny)

  ENDDO

END SUBROUTINE run

END MODULE drifter_mod
