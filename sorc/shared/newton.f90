! Convert from buoy's lat-lon location to its ij coordinate (x,y in buoy member)
! given a first guess x,y
SUBROUTINE newton(this, lat, lon, x, y)
  USE constants
  IMPLICIT none
  CLASS(metric), intent(in)        :: this
  REAL(kind=real64), intent(in)    :: lat, lon
  REAL(kind=real64), intent(inout) :: x, y

  INTEGER :: iter, itmax, ii, ij
  REAL(kind=real64) :: toler, tlat, tlon, dlat, dlon, fi, fj, dfi, dfj
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
  toler = 1./1024./1024. ! degrees
  ratio = 1.0

  tlon = lon
  tlat = lat
  
  ! First guess for location
  fi = x
  fj = y

  if (fi <= 0.5) fi = 1

  ii    = nint(fi)
  ij    = nint(fj)

  dlat = lat - this%ulat(ii,ij)
  if (tlon > 360. .or. tlon < 0) tlon = wrap(tlon)
  dlon = this%ulon(ii,ij) - wrap(lon)

  newton : do
    iter = iter + 1

    dfi   = ratio * ((dlat*this%dlondj(ii,ij)) - (dlon*this%dlatdj(ii,ij)) ) / this%area(ii,ij)
    dfj   = ratio * ((this%dlatdi(ii,ij)*dlon) - (dlat*this%dlondi(ii,ij)) ) / this%area(ii,ij)
    fi    = fi - dfi
    fj    = fj - dfj

    ! keep fi, fj, inside (1,1),(nx,ny)
    if (fi > this%nx + 0.5) THEN
      !debug: PRINT *,'fi > nx ',fi
      fi = mod(fi, REAL(this%nx,kind=real64) ) !RG: Assumes grid wraps in i
      if (fi .eq. 0) fi = this%nx
    endif
    if (fi < 0) THEN !Assuming that grid wraps around in i
      fi = fi + this%nx
      !debug: PRINT *,'fi < 0',fi
    endif
    if (fi <= 0.5) THEN
      fi = 1
      !debug: PRINT *,'fi < 0.5'
    endif

    if (fj > this%ny+0.5) THEN
      !debug: PRINT *,'fj > ny',fj
      fj = 0.75*this%ny
    endif
    if (fj < 0) THEN
      !debug: PRINT *,'fj < 0',fj
      fj = 0.25*this%ny
    endif
    if (fj <= 0.5) THEN
      fj = 1
      !debug: PRINT *,'fj < 0.5'
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

!debug2:     WRITE(*,9001) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, flat, flon
 9001 FORMAT(I3,6F10.3,4F10.3)

    IF (iter > nint(0.6*itmax) ) THEN
      ratio = 0.25
    ELSE IF (iter > nint(itmax/3._real64) ) THEN
      ratio = 0.5
    ENDIF
    IF (iter >= itmax .or. (abs(dlat) < toler .and. abs(dlon) < toler)) exit newton
  end do newton
  o
  !debug: PRINT *,'iterations ',iter

  IF (iter .eq. itmax) THEN  ! need brute force or something to cross seam
    fi = flag 
    fj = flag 
    !debug: WRITE(*,9004) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, flat, flon
  ENDIF
 9004 FORMAT('itmax ',I3,6(F10.3,1x),4F10.3)

    !debug: WRITE(*,9003) iter, dfi, dfj, fi, fj, dlat, dlon, lat, lon, flat, flon
 9003 FORMAT('final ',I3,6F10.3,4F10.3)

  ! x,y = i,j (floating) of floe
  x = fi
  y = fj
  RETURN
END SUBROUTINE newton 
