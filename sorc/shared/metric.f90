MODULE metric_mod
    USE constants

    !module
    TYPE, public :: metric
      REAL, allocatable :: dlatdi(:,:), dlondi(:,:), dlatdj(:,:), dlondj(:,:)
      REAL, allocatable :: dx(:,:), dy(:,:), rot(:,:)
      REAL, allocatable :: ulat(:,:), ulon(:,:)
      INTEGER nx, ny
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

  RETURN
END SUBROUTINE set

!lat, lon in degrees
!dx, dy in meters
SUBROUTINE local_metric(this)
  IMPLICIT none
  CLASS(metric), intent(inout) :: this
  INTEGER i, j
  REAL toler

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
  ENDDO
  i = this%nx
  DO j = 2, this%ny
    this%dlatdj(i,j) = this%ulat(i,j) - this%ulat(i,j-1)
    this%dlondj(i,j) = this%ulon(i,j) - this%ulon(i,j-1)
  ENDDO
  this%dlatdi(this%nx, this%ny) = this%dlatdi(this%nx - 1, this%ny-1)
  this%dlatdj(this%nx, this%ny) = this%dlatdj(this%nx - 1, this%ny-1)
  this%dlondi(this%nx, this%ny) = this%dlondi(this%nx - 1, this%ny-1)
  this%dlondj(this%nx, this%ny) = this%dlondj(this%nx - 1, this%ny-1)

  !rot = atan2(?,?)

  !debug: PRINT *,'lat metrics',MAXVAL(this%dlatdi), MINVAL(this%dlatdi), MAXVAL(this%dlatdj), MINVAL(this%dlatdj)
  !PRINT *,'lon metrics',MAXVAL(this%dlondi), MINVAL(this%dlondi), MAXVAL(this%dlondj), MINVAL(this%dlondj)
!  toler = 1.0
!  DO j = 1, this%ny
!  DO i = 1, this%nx
!    IF (abs(this%dlondj(i,j) ) > toler ) THEN
!      PRINT *,'dj ',i,j,this%dlondj(i,j), this%ulon(i,j)
!    ENDIF
!    IF (abs(this%dlondi(i,j) ) > toler ) THEN
!      PRINT *,'di ',i,j,this%dlondi(i,j), this%ulon(i,j)
!    ENDIF
!  ENDDO
!  ENDDO


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
  USE constants
  IMPLICIT none
  CLASS(metric), intent(in) :: this
  REAL, intent(in)    :: lat, lon
  REAL, intent(inout)   :: x, y

  INTEGER :: iter, itmax, ii, ij
  REAL toler, tlat, tlon, delta, dlat, dlon, fi, fj, dfi, dfj
  REAL ratio
  REAL wrap

! if flag values (lat or lon >= 1.e30) skip, assign xy to 1,1
  IF (lat >= flag .or. lon >= flag) THEN
    x = 1.
    y = 1.
    RETURN
  ENDIF

! Use something like Newton method with starting point as if grid were linear
  itmax = 60
  iter  = 0
  toler = 0.05 ! degrees
  ratio = 1.

  tlon = lon
  !RG: i = 1 at ~74 longitude
  fi = (tlon/360)*this%nx
  if (fi <= 0.5) fi = 1
  ii    = int(fi+0.5)

  tlat = lat
  if (lat == 90) tlat = lat - 0.05
  fj = (tlat+78.64)*this%ny/(90+78.64)
  if (fj <= 0) THEN
    fj = 1.
  ENDIF
  IF (fj >= this%ny) THEN
    PRINT *,'fj overrunning grid',fj, lat, tlat, this%ny
    STOP
  ENDIF
  ij    = int(fj+0.5)

  dlat = tlat - this%ulat(ii,ij)
  tlon = this%ulon(ii,ij)
  if (tlon > 360. .or. tlon < 0) tlon = wrap(tlon)
  dlon = lon - tlon


!debug: WRITE (*,9002) fi, fj, dlat, dlon, tlat, lon, this%ulat(ii,ij), tlon
  9002 FORMAT('init ',8F10.3)

  newton : do
    iter = iter + 1
    delta = (this%dlatdi(ii,ij)*this%dlondj(ii,ij) - this%dlatdj(ii,ij)*this%dlondi(ii,ij))
    if (delta == 0) THEN
      !debug: PRINT *,'delta == 0 ',delta
      delta = 3.e-3
    endif

    dfi   = ratio * ((dlat*this%dlondj(ii,ij)) - (dlon*this%dlatdj(ii,ij)) ) / delta
    dfj   = ratio * ((this%dlatdi(ii,ij)*dlon) - (dlat*this%dlondi(ii,ij)) ) / delta
    fi    = fi + dfi
    fj    = fj + dfj

    ! keep fi, fj, inside (1,1),(nx,ny)
    if (fi > this%nx) THEN
      !debug: PRINT *,'fi > nx ',fi
      fi = mod(fi, float(this%nx) )
    endif
    if (fi < 0) THEN !Assuming that grid wraps around in i
      fi = fi + this%nx
      !debug: PRINT *,'fi < 0',fi
    endif
    if (fi <= 0.5) THEN
      fi = 1
      !debug: PRINT *,'fi < 0.5'
    endif

    if (fj > this%ny) THEN
      fj = 0.75*this%ny
      !debug: PRINT *,'fj > ny'
    endif
    if (fj < 0) THEN
      fj = 0.25*this%ny
      !debug: PRINT *,'fj < 0'
    endif
    if (fj <= 0.5) THEN
      fj = 1
      !debug: PRINT *,'fj < 0.5'
    endif
    ii    = int(fi+0.5)
    ij    = int(fj+0.5)

    dlat = tlat - this%ulat(ii,ij)
    tlon = this%ulon(ii,ij)
    if (tlon > 360. .or. tlon < 0) tlon = wrap(tlon)
    dlon = lon - tlon

!debug:     WRITE(*,9001) iter, dfi, dfj, fi, fj, dlat, dlon, tlat, lon, this%ulat(ii,ij), tlon
    IF (iter > 35 ) THEN
      ratio = 0.25
    ELSE IF (iter > 20) THEN
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

    !debug: WRITE(*,9003) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, this%ulat(ii,ij), this%ulon(ii,ij)
 9003 FORMAT('final ',I3,6F10.3,4F10.3)

  ! x,y = i,j (floating) of floe
  x = fi
  y = fj
  RETURN
END SUBROUTINE ll_to_xy 
SUBROUTINE ll_to_xy_brute(this, lat, lon, fi, fj)
  USE constants
  CLASS(metric), intent(in) :: this
  REAL, intent(in)    :: lat, lon
  REAL, intent(inout) :: fi, fj
  INTEGER i,j, jmin, bi, bj 
  REAL dlat, dlon, tlon, dbest, wrap
!RG: This is very slow. The newton search fails (mostly) because of the tripolar seam.
!       should be able to take advantage of that. Lats > 45.
!    Sometimes also encounter difficulty along 0 E 
  !debug: PRINT *,'brute ',lat,lon
  fi = flag
  fj = flag
  RETURN

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

  IMPLICIT none
  CLASS(metric), intent(in) :: this
  REAL, intent(in)    :: lat, lon
  REAL, intent(inout)   :: x, y

  RETURN
END SUBROUTINE xy_to_ll 

END MODULE metric_mod
