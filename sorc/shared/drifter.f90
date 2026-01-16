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
  SUBROUTINE init(buoy, tlon, tlat, clon, clat, xmetric)
    USE constants
    USE metric_mod
    IMPLICIT none
    REAL, intent(in) :: tlon, tlat, clon, clat
    TYPE(metric), intent(in) :: xmetric
    CLASS(drifter), intent(inout) :: buoy
    REAL x, y

    buoy%ilat = tlat
    buoy%ilon = tlon
    if (clat < flag .and. clon < flag) THEN
      buoy%clat = clat
      buoy%clon = clon
    else
      buoy%clat = tlat
      buoy%clon = tlon
    endif

    CALL xmetric%ll_to_xy(clat, clon, x, y)
    buoy%x = x
    buoy%y = y
    if (x >= flag .or. y >= flag) THEN
      !debug: PRINT *,'skipping buoy at ',clat, clon
      buoy%clat = flag
      buoy%clon = flag
    endif
    !RG: initialize to skip if dlon/d(ij) too large

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
    REAL, intent(in) ::  u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
    REAL, intent(in) :: dt

    REAL tu, tv, deltax, deltay, di, dj
    INTEGER ti, tj, nti, ntj

    if (buoy%x >= flag .or. buoy%y >= flag) RETURN
    if (buoy%clat >= flag .or. buoy%clon >= flag) RETURN
    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    tu = u(ti, tj)
    tv = v(ti, tj)
    !flag values
    if (tu >= flag .or. tv >= flag) RETURN

    !RG:  These could be interpolated (bilinear, ...)
    deltax = tu * dt  !deltax, deltay are meters
    deltay = tv * dt
    IF ((abs(deltax) > 3.*3600) .or. (abs(deltay) > 3.*3600) ) THEN
!debug: 
      PRINT *,'fast ',deltax, deltay, deltax/xmetric%dx(ti, tj), deltay/xmetric%dy(ti, tj)
    ENDIF

    !RG: beware of seams
    di = deltax/xmetric%dx(ti, tj) !di,dj are degrees
    dj = deltay/xmetric%dy(ti, tj)
    buoy%x = buoy%x + di
    buoy%y = buoy%y + dj

    nti = NINT(buoy%x)
    ntj = NINT(buoy%y)
    ! beware of running outside (1,1),(nx,ny)
    IF (nti < 1 .or. nti > xmetric%nx .or. ntj < 1 .or. ntj > xmetric%ny) THEN
      PRINT *,'buoy out of bounds ',ti,tj,nti, ntj
      buoy%clat = flag
      buoy%clon = flag
      buoy%x    = flag
      buoy%y    = flag
    ELSE IF (abs(xmetric%dlondj(nti,ntj)) > 80 .or. abs(xmetric%dlondi(nti,ntj)) > 80 ) THEN
      !debug: PRINT *,'near seam ',ti,tj,nti, ntj
      buoy%clat = flag
      buoy%clon = flag
      buoy%x    = flag
      buoy%y    = flag
    ELSE 
      buoy%clat = xmetric%ulat(nti, ntj) + (ntj-buoy%y)*xmetric%dlatdj(nti,ntj)
      buoy%clon = xmetric%ulon(nti, ntj) + (nti-buoy%x)*xmetric%dlondi(nti,ntj) 
    !debug: PRINT *,'move ',buoy%x, buoy%y, ti, tj, nti, ntj
    ENDIF

  RETURN
  END subroutine move

!----------------------------------------------------------------
SUBROUTINE run(buoys, nbuoy, u, v, xmetric, dt)
  USE metric_mod
  IMPLICIT none

  TYPE(metric), intent(in) :: xmetric
  INTEGER, intent(in) :: nbuoy
  REAL, intent(in) :: u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
  REAL, intent(in) :: dt

  TYPE(drifter), intent(inout) ::  buoys(nbuoy)

  INTEGER k, track

  track = INT(0.5+nbuoy*5./6.)

  DO k = 1, nbuoy
    !c-like (object-like) 
    IF (k .eq. track) THEN
      PRINT *,buoys(k)%clat, buoys(k)%clon, buoys(k)%x, buoys(k)%y
    ENDIF
    CALL buoys(k)%move(u, v, xmetric, dt)
    IF (k .eq. track) THEN
      PRINT *,'  ',buoys(k)%clat, buoys(k)%clon, buoys(k)%x, buoys(k)%y, dt
    ENDIF
  ENDDO

END SUBROUTINE run

END MODULE drifter_mod
