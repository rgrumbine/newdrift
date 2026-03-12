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
    PROCEDURE, pass :: init, move, zero, invalidate_buoy
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

    IF (any([buoy%x, buoy%y, buoy%clat, buoy%clon, buoy%ilat, buoy%ilon] >= flag)) THEN
      CALL invalidate_buoy(buoy)
      RETURN
    ENDIF
    ti = NINT(buoy%x)
    tj = NINT(buoy%y)
    tu = u(ti, tj)
    tv = v(ti, tj)
    !flag values
    if (tu >= flag .or. tv >= flag) RETURN

    !RG:  These could be interpolated (bilinear, ...)
    a = (tu * dt) / xmetric%dx(ti, tj)
    b = (tv * dt) / xmetric%dy(ti, tj)

    !RG: beware of seams
    IF (xmetric%dlondi(ti,tj) .ne. 0 .and. xmetric%dlatdj(ti,tj) .ne. 0 .and. xmetric%area(ti, tj) .ne. 0) THEN
      ! these two are true only if x,y are parallel to i,j
      !di = a/xmetric%dlondi(ti,tj) !di,dj are units of grid points
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
      WRITE(*,9001) ti,tj,di, dj, nti, ntj, xmetric%dlondi(ti,tj), xmetric%dlatdj(ti,tj) 
 9001 FORMAT('buoy out of bounds ',2I5,2F10.3,2I11,2E14.6)
      CALL invalidate_buoy(buoy)
    ELSE IF (abs(xmetric%dlondj(nti,ntj)) > 80 .or. abs(xmetric%dlondi(nti,ntj)) > 80 ) THEN
      !debug: 
      WRITE(*,9002) ti,tj,nti, ntj,xmetric%dlondj(nti,ntj), xmetric%dlondi(nti,ntj)
 9002 FORMAT('near seam ',4I5,2E14.6)
      CALL invalidate_buoy(buoy)
    ELSE 
      CALL xy_to_ll(xmetric, flat, flon, buoy%x, buoy%y)
      buoy%clat = flat
      buoy%clon = flon
    !debug2: PRINT *,'move ',buoy%x, buoy%y, ti, tj, nti, ntj
    ENDIF

  RETURN
  END subroutine move

  ! Helper to reset all buoy values to flag
  SUBROUTINE invalidate_buoy(b)
     CLASS(drifter), INTENT(inout) :: b
     b%x = flag; b%y = flag
     b%clat = flag; b%clon = flag
     b%ilat = flag; b%ilon = flag
  END SUBROUTINE invalidate_buoy

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

  track = 0
  !debug: 
  !track = INT(0.5+nbuoy*5./6.)
  track = 118703

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
