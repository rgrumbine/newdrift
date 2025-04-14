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
    TYPE(metric), intent(in) :: xmetric
    CLASS(drifter), intent(inout) :: buoy
    REAL x, y

    buoy%ilat = tlat
    buoy%ilon = tlon
    buoy%clat = tlat
    buoy%clon = tlon

    CALL xmetric%ll_to_xy(tlat, tlon, x, y)
    buoy%x = x
    buoy%y = y
    !debug: WRITE(*,9001) tlat, tlon, x, y
 9001 FORMAT(4F10.3)

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

  SUBROUTINE move(buoy, u, v, xmetric, dt)
!note dx, dy are the mesh variables
    USE metric_mod
    IMPLICIT none

    CLASS(drifter), intent(inout) :: buoy
    TYPE(metric), intent(in) :: xmetric
    !INTEGER, intent(in) :: nx, ny
    !REAL, intent(in) :: dx(nx, ny), dy(nx, ny)
    REAL, intent(in) ::  u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
    REAL, intent(in) :: dt

    REAL tu, tv, deltax, deltay
    INTEGER ti, tj

    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    tu = u(ti, tj)
    tv = v(ti, tj)
    !flag values
    if (tu > 1.e30 .or. tv > 1.e30) RETURN

    !RG:  These could be interpolated (bilinear, ...)
    deltax = tu * dt
    deltay = tv * dt
    IF ((abs(deltax) > 3.*3600) .or. (abs(deltay) > 3.*3600) ) THEN
      PRINT *,'fast ',deltax, deltay, deltax/xmetric%dx(ti, tj), deltay/xmetric%dy(ti, tj)
    ENDIF

    !RG: beware of seams
    !RG: beware of running outside (1,1),(nx,ny)
    buoy%x = buoy%x + deltax/xmetric%dx(ti, tj)
    buoy%y = buoy%y + deltay/xmetric%dy(ti, tj)
    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    !RG: refine by adding for dlat/di, dlon/dj for residual ti-x, tj-y
    buoy%clat = xmetric%ulat(ti, tj) + (tj-buoy%y)*xmetric%dlatdj(ti,tj)
    buoy%clon = xmetric%ulon(ti, tj) + (ti-buoy%x)*xmetric%dlondi(ti,tj) 
    !debug: PRINT *,'move ',buoy%x, buoy%y

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
    CALL buoys(k)%move(u, v, xmetric, dt)
  ENDDO

END SUBROUTINE run

END MODULE drifter_mod
