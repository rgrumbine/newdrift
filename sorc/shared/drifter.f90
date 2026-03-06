! Module/class for drifting bodies
! expects a set of lat-lon locations and a velocity field with velocity 
!    grid type (B/C/...) and lat-longs of velocity points.
!    also time step to extrapolate over
MODULE drifter_mod
  USE metric_mod
  USE constants

  IMPLICIT none
  TYPE, public :: drifter
    REAL(kind=real64) :: x, y        ! current i,j location
    REAL(kind=real64) :: ilat, ilon  ! initial latitude-longitude
    REAL(kind=real64) :: clat, clon  ! current latitude-longitude
  CONTAINS
    PROCEDURE, pass :: init, move, zero
  END TYPE drifter


CONTAINS
  SUBROUTINE init(buoy, tlon, tlat, clon, clat, xmetric)
    USE constants
    USE metric_mod
    IMPLICIT none
    REAL(kind=real64), intent(in) :: tlon, tlat, clon, clat
    TYPE(metric), intent(in) :: xmetric
    CLASS(drifter), intent(inout) :: buoy
    REAL(kind=real64) :: x, y

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
    USE metric_mod, only: xy_to_ll
    IMPLICIT none

    CLASS(drifter), intent(inout) :: buoy
    TYPE(metric), intent(in) :: xmetric
    REAL(kind=real64), intent(in) ::  u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
    REAL(kind=real64), intent(in) :: dt

    REAL(kind=real64) :: tu, tv, deltax, deltay, di, dj
    INTEGER ti, tj, nti, ntj
    REAL(kind=real64) :: a, b, toler, flat, flon

    if (buoy%x >= flag .or. buoy%y >= flag) RETURN
    if (buoy%clat >= flag .or. buoy%clon >= flag) RETURN
    if (buoy%ilat >= flag .or. buoy%ilon >= flag) RETURN
    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    tu = u(ti, tj)
    tv = v(ti, tj)
    !flag values
    if (tu >= flag .or. tv >= flag) RETURN

    !RG:  These could be interpolated (bilinear, ...)
    deltax = tu * dt  !deltax, deltay are meters
    deltay = tv * dt
    a = deltax / xmetric%dx(ti, tj)
    b = deltay / xmetric%dy(ti, tj)
    IF ((abs(deltax) > 3.*dt) .or. (abs(deltay) > 3.*dt) ) THEN
!debug: 
      PRINT *,'fast ',xmetric%ulat(ti, tj), xmetric%ulon(ti, tj), ti, tj, deltax, deltay
    ENDIF

    !RG: beware of seams
    IF (xmetric%dlondi(ti,tj) .ne. 0 .and. xmetric%dlatdj(ti,tj) .ne. 0 .and. xmetric%area(ti, tj) .ne. 0) THEN
      ! these two are true only if x,y are parallel to i,j
      !di = a/xmetric%dlondi(ti,tj) !di,dj are grid points
      !dj = b/xmetric%dlatdj(ti,tj)
      di = (a*xmetric%dlatdj(ti,tj) - b*xmetric%dlondj(ti,tj))/xmetric%area(ti,tj)
      dj = (b*xmetric%dlondi(ti,tj) - a*xmetric%dlatdi(ti,tj))/xmetric%area(ti,tj)
      buoy%x = buoy%x + di
      buoy%y = buoy%y + dj
    ELSE
      buoy%x = flag
      buoy%y = flag
    ENDIF

    nti = NINT(buoy%x)
    ntj = NINT(buoy%y)
    ! beware of running outside (1,1),(nx,ny)
    IF (nti < 1 .or. nti > xmetric%nx .or. ntj < 1 .or. ntj > xmetric%ny) THEN
      !realbug: 
      PRINT *,'buoy out of bounds ',ti,tj,di, dj, nti, ntj,xmetric%dlondi(ti,tj), xmetric%dlatdj(ti,tj) 
 9001 FORMAT(2I5,2F10.3,2I5,2E13.6)
      buoy%ilat = flag
      buoy%ilon = flag
      buoy%clat = flag
      buoy%clon = flag
      buoy%x    = flag
      buoy%y    = flag
    ELSE IF (abs(xmetric%dlondj(nti,ntj)) > 40 .or. abs(xmetric%dlondi(nti,ntj)) > 40 ) THEN
      !debug: 
      PRINT *,'near seam ',ti,tj,nti, ntj
      buoy%clat = flag
      buoy%clon = flag
      buoy%x    = flag
      buoy%y    = flag
    ELSE 
      !buoy%clat = xmetric%ulat(nti, ntj) + (ntj-buoy%y)*xmetric%dlatdj(nti,ntj)
      !buoy%clon = xmetric%ulon(nti, ntj) + (nti-buoy%x)*xmetric%dlondi(nti,ntj) 
      CALL xy_to_ll(xmetric, flat, flon, buoy%x, buoy%y)
      buoy%clat = flat
      buoy%clon = flon
    !debug2: PRINT *,'move ',buoy%x, buoy%y, ti, tj, nti, ntj
    ENDIF

  RETURN
  END subroutine move

!----------------------------------------------------------------
SUBROUTINE run(buoys, nbuoy, u, v, xmetric, dt)
  USE metric_mod
  IMPLICIT none

  TYPE(metric), intent(in) :: xmetric
  INTEGER, intent(in) :: nbuoy
  REAL(kind=real64), intent(in) :: u(xmetric%nx, xmetric%ny), v(xmetric%nx, xmetric%ny)
  REAL(kind=real64), intent(in) :: dt

  TYPE(drifter), intent(inout) ::  buoys(nbuoy)

  INTEGER k, track

  !debug: 
  track = INT(0.5+nbuoy*6./7.)
  !track = 100

  DO k = 1, nbuoy
    !c-like (object-like) 
    !debug: 
    IF (k .eq. track) THEN
      WRITE(*,9001) buoys(k)%clat, buoys(k)%clon, buoys(k)%x, buoys(k)%y
    ENDIF

    CALL buoys(k)%move(u, v, xmetric, dt)
    !debug: 
    IF (k .eq. track) THEN
      WRITE(*,9002) buoys(k)%clat, buoys(k)%clon, buoys(k)%x, buoys(k)%y, dt
    ENDIF
  ENDDO
 9001 FORMAT('track',4F10.4)
 9002 FORMAT('track',4F10.4, 1x,F7.1)

END SUBROUTINE run

END MODULE drifter_mod
