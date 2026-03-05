MODULE metric_mod
    USE constants

    !module
    TYPE, public :: metric
      REAL(kind=real64), allocatable :: dlatdi(:,:), dlondi(:,:), dlatdj(:,:), dlondj(:,:)
      REAL(kind=real64), allocatable :: dxdi(:,:), dydi(:,:), dxdj(:,:), dydj(:,:)
      REAL(kind=real64), allocatable :: dx(:,:), dy(:,:), rot(:,:), deg_area(:,:), area(:,:)
      REAL(kind=real64), allocatable :: ulat(:,:), ulon(:,:)
      INTEGER :: nx = 0
      INTEGER :: ny = 0
    CONTAINS
      PROCEDURE :: set, local_metric, local_cartesian, ll_to_xy, ll_to_xy_brute
      PROCEDURE :: xy_to_ll
    END TYPE metric
CONTAINS

SUBROUTINE set(this, nx, ny)
  IMPLICIT none
  INTEGER, intent(in) :: nx, ny
  CLASS(metric), intent(inout) :: this

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
  allocate(this%dxdi(nx, ny))
  allocate(this%dxdj(nx, ny))
  allocate(this%dydi(nx, ny))
  allocate(this%dydj(nx, ny))
  allocate(this%area(nx, ny))
  allocate(this%deg_area(nx, ny))

  RETURN
END SUBROUTINE set

!lat, lon in degrees
!dx, dy in meters
SUBROUTINE local_metric(this)
  IMPLICIT none
  CLASS(metric), intent(inout) :: this
  INTEGER i, j
  REAL(kind=real64) :: toler, harcdis

  this%rot = 0.
  CALL this%local_cartesian()

  !RG: Add loops for j = ny or i = nx, and then the corner nx, ny
  DO j = 1, this%ny-1
  DO i = 1, this%nx-1
    this%dlatdi(i,j) = this%ulat(i+1,j) - this%ulat(i,j)
    this%dlondi(i,j) = this%ulon(i+1,j) - this%ulon(i,j)
    this%dlatdj(i,j) = this%ulat(i,j+1) - this%ulat(i,j)
    this%dlondj(i,j) = this%ulon(i,j+1) - this%ulon(i,j)
  ENDDO
  ENDDO
  j = this%ny
  DO i = 2, this%nx
    this%dlatdi(i,j) = this%ulat(i,j) - this%ulat(i-1,j)
    this%dlondi(i,j) = this%ulon(i,j) - this%ulon(i-1,j)
    this%dlatdj(i,j) = this%ulat(i,j) - this%ulat(i,j-1)
    this%dlondj(i,j) = this%ulon(i,j) - this%ulon(i,j-1)

  ENDDO
  i = this%nx
  DO j = 2, this%ny
    this%dlatdi(i,j) = this%ulat(i,j) - this%ulat(i-1,j)
    this%dlondi(i,j) = this%ulon(i,j) - this%ulon(i-1,j)
    this%dlatdj(i,j) = this%ulat(i,j) - this%ulat(i,j-1)
    this%dlondj(i,j) = this%ulon(i,j) - this%ulon(i,j-1)
  ENDDO
  ! Corners
  this%dlatdi(this%nx, this%ny) = this%dlatdi(this%nx - 1, this%ny-1)
  this%dlatdj(this%nx, this%ny) = this%dlatdj(this%nx - 1, this%ny-1)
  this%dlondi(this%nx, this%ny) = this%dlondi(this%nx - 1, this%ny-1)
  this%dlondj(this%nx, this%ny) = this%dlondj(this%nx - 1, this%ny-1)
  this%dlatdi(this%nx,1) = this%dlatdi(this%nx - 1,2)
  this%dlatdj(this%nx,1) = this%dlatdj(this%nx - 1,2)
  this%dlondi(this%nx,1) = this%dlondi(this%nx - 1,2)
  this%dlondj(this%nx,1) = this%dlondj(this%nx - 1,2)
  this%dlatdi(1,this%ny) = this%dlatdi(2, this%ny)
  this%dlatdj(1,this%ny) = this%dlatdj(2, this%ny)
  this%dlondi(1,this%ny) = this%dlondi(2, this%ny)
  this%dlondj(1,this%ny) = this%dlondj(2, this%ny)

  this%area = this%dlondi*this%dlatdj - this%dlondj*this%dlatdi
  !debug: PRINT *,'area max min',MAXVAL(this%area), MINVAL(this%area)
  !debug: PRINT *,'dlatdj max min',MAXVAL(this%dlatdj), MINVAL(this%dlatdj)
  !debug: PRINT *,'dlondj max min',MAXVAL(this%dlondj), MINVAL(this%dlondj)
  !debug: PRINT *,'dlatdi max min',MAXVAL(this%dlatdi), MINVAL(this%dlatdi)
  !debug: PRINT *,'dlondi max min',MAXVAL(this%dlondi), MINVAL(this%dlondi)
  
  !rot = atan2(?,?)

  RETURN
END subroutine local_metric


SUBROUTINE local_cartesian(this)
  USE constants
  IMPLICIT none
  CLASS(metric), intent(inout) :: this
  INTEGER k
  INTEGER i, j

! From WGS84 via Amy Solomon, ESRL
! Meters per degree
  DO j = 1, this%ny
  DO i = 1, this%nx
    this%dy(i,j) = 111132.92 - 559.82*cos(2*this%ulat(i,j)*rpd) + &
                           1.175 *cos(4*this%ulat(i,j)*rpd) - &
                           0.0023*cos(6*this%ulat(i,j)*rpd)
    this%dx(i,j) = 111412.84 *cos(  this%ulat(i,j)*rpd) - &
                       93.5  *cos(3*this%ulat(i,j)*rpd) + &
                        0.118*cos(5*this%ulat(i,j)*rpd)
  ENDDO
  ENDDO

  RETURN
END subroutine local_cartesian

! Convert from buoy's lat-lon location to its ij coordinate (x,y in buoy member)
SUBROUTINE ll_to_xy(this, lat, lon, x, y)
  USE constants
  IMPLICIT none
  CLASS(metric), intent(in) :: this
  REAL(kind=real64), intent(in)    :: lat, lon
  REAL(kind=real64), intent(inout)   :: x, y

  INTEGER :: iter, itmax, ii, ij
  REAL(kind=real64) :: toler, tlat, tlon, delta, dlat, dlon, fi, fj, dfi, dfj
  REAL(kind=real64) :: ratio, wrap
  REAL(kind=real64) :: flat, flon

! if flag values (lat or lon >= flag) skip, assign xy to flag, flag
  IF (lat >= flag .or. lon >= flag) THEN
    x = flag
    y = flag
    RETURN
  ENDIF

! Use something like Newton method with starting point as if grid were linear
  itmax = 90
  iter  = 0
  toler = 1./1024. ! degrees
  ratio = 1.0

  tlon = lon
  tlat = lat
  
  fi = 1.0
  fj = 1.0
  !fi = this%nx / 2
  !fj = this%ny / 2
  !fj = (tlat+78.64)*this%ny/(90+78.64)

  if (fi <= 0.5) fi = 1

  if (lat == 90) tlat = lat - 0.05

  ii    = nint(fi)
  ij    = nint(fj)

  dlat = lat - this%ulat(ii,ij)
  tlon = lon
  if (tlon > 360. .or. tlon < 0) tlon = wrap(tlon)
  dlon = this%ulon(ii,ij) - wrap(lon)

!debug: WRITE (*,9002) fi, fj, dlat, dlon, tlat, tlon, this%ulat(ii,ij), tlon
 9002 FORMAT('init ',8F10.3)

  newton : do
    iter = iter + 1
    delta = this%area(ii,ij)
    if (delta == 0) THEN
      !debug: 
      PRINT *,'delta == 0 ',delta
      delta = 3.e-3
    endif

    !debug: PRINT *,'dlat dlon ',dlat, dlon
    dfi   = ratio * ((dlat*this%dlondj(ii,ij)) - (dlon*this%dlatdj(ii,ij)) ) / delta
    dfj   = ratio * ((this%dlatdi(ii,ij)*dlon) - (dlat*this%dlondi(ii,ij)) ) / delta
    fi    = fi - dfi
    fj    = fj - dfj

    ! keep fi, fj, inside (1,1),(nx,ny)
    if (fi > this%nx + 0.5) THEN
      !debug: 
      PRINT *,'fi > nx ',fi
      fi = mod(fi, REAL(this%nx,kind=real64) )
      if (fi .eq. 0) fi = this%nx
    endif
    if (fi < 0) THEN !Assuming that grid wraps around in i
      fi = fi + this%nx
      !debug: 
      PRINT *,'fi < 0',fi
    endif
    if (fi <= 0.5) THEN
      fi = 1
      !debug: 
      PRINT *,'fi < 0.5'
    endif

    if (fj > this%ny) THEN
      fj = 0.75*this%ny
      !debug: 
      PRINT *,'fj > ny'
    endif
    if (fj < 0) THEN
      fj = 0.25*this%ny
      !debug: 
      PRINT *,'fj < 0'
    endif
    if (fj <= 0.5) THEN
      fj = 1
      !debug: 
      PRINT *,'fj < 0.5'
    endif
    ii    = nint(fi)
    ij    = nint(fj)

    CALL xy_to_ll(this, flat, flon, fi, fj)
    IF (flat .ge. flag .or. flon .ge. flag) THEN
      iter = itmax
      dlat = flag
      dlon = flag
    ELSE
      dlat = tlat - flat
      dlon = tlon - flon
      !debug: PRINT *,'dlat dlon ',dlat, dlon
    ENDIF
    !dlat = tlat - this%ulat(ii,ij)
    !tlon = this%ulon(ii,ij)
    !if (tlon > 360. .or. tlon < 0) tlon = wrap(tlon)
    !dlon = lon - tlon

!debug2:     WRITE(*,9001) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, flat, flon
    IF (iter > nint(0.6*itmax) ) THEN
      ratio = 0.25
    ELSE IF (iter > nint(itmax/3._real64) ) THEN
      ratio = 0.5
    ENDIF
    IF (iter >= itmax .or. (abs(dlat) < toler .and. abs(dlon) < toler)) exit newton
  end do newton
  !debug: PRINT *,'iterations ',iter
 9001 FORMAT(I3,6F10.3,4F10.3)

  IF (iter .eq. itmax) THEN  ! need brute force or something to cross seam
    fi = flag 
    fj = flag 
  ENDIF
 9004 FORMAT('itmax ',I3,6F10.3,4F10.3)

    !debug: WRITE(*,9003) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, flat, flon
 9003 FORMAT('final ',I3,6F10.3,4F10.3)

  ! x,y = i,j (floating) of floe
  x = fi
  y = fj
  RETURN
END SUBROUTINE ll_to_xy 
SUBROUTINE ll_to_xy_brute(this, lat, lon, fi, fj)
  USE constants
  CLASS(metric), intent(in) :: this
  REAL(kind=real64), intent(in)    :: lat, lon
  REAL(kind=real64), intent(inout) :: fi, fj
  INTEGER i,j, jmin, bi, bj 
  REAL(kind=real64) :: dlat, dlon, tlon, dbest, wrap
!RG: This is very slow. The newton search fails (mostly) because of the tripolar seam.
!       should be able to take advantage of that. Lats > 45.
!    Sometimes also encounter difficulty along 0 E 
  !debug: PRINT *,'brute ',lat,lon
  !fi = flag
  !fj = flag
  !RETURN

  dbest = 999.
  bi = 1
  bj = 1
  IF (lat < 45.) THEN
    jmin = 1
  ELSE
    jmin = int(0.5 + 0.65*this%ny)
  ENDIF
  DO j = jmin, this%ny
  DO i = 1, this%nx
    dlat = lat - this%ulat(i,j)
    tlon = this%ulon(i,j)
    if (tlon > 360. .or. tlon < 0.) tlon = wrap(tlon)
    dlon = lon - tlon
    if ((dlat*dlat + dlon*dlon) < dbest) THEN
      bi = i
      bj = j
      dbest = (dlat*dlat + dlon*dlon)
    endif
  END DO 
  END DO 
  fi = bi
  fj = bj
  !debug: PRINT *,'dbest ',dbest

  RETURN
END SUBROUTINE ll_to_xy_brute


SUBROUTINE xy_to_ll(this, lat, lon, x, y)
  USE constants
  IMPLICIT none
  CLASS(metric), intent(in) :: this
  REAL(kind=real64), intent(out)    :: lat, lon
  REAL(kind=real64), intent(in)     :: x, y
  INTEGER ii, ij
  REAL(kind=real64) :: di, dj

  if (x < 0.5 .or. x > this%nx+0.5 .or. y < 0.5 .or. y > this%ny+0.5) then
    PRINT *,'off grid in xytoll', x, y
    lat = flag
    lon = flag
    RETURN
  endif
  ii = NINT(x)
  ij = NINT(y)
  di = x - ii
  dj = y - ij
  lat = this%ulat(ii,ij)
  lon = this%ulon(ii,ij)
  lat = lat + di*this%dlatdi(ii,ij) + dj*this%dlatdj(ii,ij)
  lon = lon + di*this%dlondi(ii,ij) + dj*this%dlondj(ii,ij)
  !debug: WRITE (*,9001) lat, lon, x, y
 9001 FORMAT('xytoll ',2F10.3,2F10.3)

  RETURN
END SUBROUTINE xy_to_ll 

END MODULE metric_mod
