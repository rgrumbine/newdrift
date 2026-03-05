MODULE io
  USE netcdf
  USE drifter_mod
  USE metric_mod
  PUBLIC

CONTAINS

SUBROUTINE initialize_in(nvar, fname, ncid, varid, varnames, xname, yname, nx, ny)
!  USE metric_mod
  IMPLICIT none

  INTEGER, intent(in) :: nvar
  CHARACTER(*), intent(in) :: fname
  INTEGER, intent(out) :: nx, ny
  INTEGER, intent(out) :: ncid
  INTEGER, intent(out) :: varid(nvar)
! Names from rtofs output
  CHARACTER(len=40), intent(in) :: varnames(nvar)
  CHARACTER(len=50), intent(in) :: xname, yname
  CHARACTER(len=50) :: nxname, nyname
  
! For Netcdf processing
  INTEGER i
  INTEGER retcode

  !debug: PRINT *,'entered initialize_in'
!! This is the rtofs cice_inst variable set -- much more extensive
!  varnames(1) = "TLON"
!  varnames(2) = "TLAT"
!  varnames(3) = "ULON"
!  varnames(4) = "ULAT"
!  varnames(5) = "hi"
!  varnames(6) = "hs"
!  varnames(7) = "Tsfc"
!  varnames(8) = "aice"
!  varnames(9) = "uvel"
!  varnames(10) = "vvel"
!  varnames(11) = "fswdn"
!  varnames(12) = "flwdn"
!  varnames(13) = "snow"
!  varnames(14) = "snow_ai"
!  varnames(15) = "rain_ai"
!  varnames(16) = "sst"
!  varnames(17) = "uocn"
!  varnames(18) = "vocn"
!  varnames(19) = "meltt"
!  varnames(20) = "meltb"
!  varnames(21) = "meltl"
!  varnames(22) = "strength"
!  varnames(23) = "divu"
!  varnames(24) = "Tair"
!RG: make file type an input parameter 2ds_ice:
  !RTOFS
  !varnames(1) = "Longitude"
  !varnames(2) = "Latitude"
  !varnames(3) = "ice_coverage"
  !varnames(4) = "ice_temperature"
  !varnames(5) = "ice_thickness"
  !varnames(6) = "ice_uvelocity"
  !varnames(7) = "ice_vvelocity"
  !UFS
  !varnames(1) = "TLON"
  !varnames(2) = "TLAT"
  !varnames(3) = "aice_h"
  !varnames(4) = "Tsfc_h"
  !varnames(5) = "hi_h"
  !varnames(6) = "uvel_h"
  !varnames(7) = "vvel_h"
  
  !debug: PRINT *,'trying to open ',fname, len(fname)
  retcode = nf90_open(fname, NF90_NOWRITE, ncid)
  CALL check(retcode)

  DO i = 1, nvar
    !debug: PRINT *,i,' reading varname ',varnames(i)
    retcode = nf90_inq_varid(ncid, varnames(i), varid(i))
    CALL check(retcode)
  ENDDO
  !debug: PRINT *,'done reading varnames'

  nxname = xname
  nyname = yname
  !debug: PRINT *,'xname, yname',xname, yname
  !RG: Read nx, ny in from netcdf file -- 
  ! rtofs
  retcode = nf90_inquire_dimension(ncid, 3, nxname, nx)
  CALL check(retcode)
  retcode = nf90_inquire_dimension(ncid, 2, nyname, ny)
  CALL check(retcode)

  !!ufs
  !retcode = nf90_inquire_dimension(ncid, 2, nxname, nx)
  !CALL check(retcode)
  !retcode = nf90_inquire_dimension(ncid, 3, nyname, ny)
  !CALL check(retcode)

  !debug: PRINT *,'leaving initialize_in', nx, ny

RETURN
END subroutine initialize_in

!----------------------------------------------------------------
SUBROUTINE initial_read(fname, nvar, ncid, varid, &
                        allvars, xmetric)
  USE metric_mod
  IMPLICIT none

  CHARACTER(90), intent(in) :: fname
  INTEGER, intent(in)  :: nvar, ncid, varid(nvar)

  TYPE(metric),intent(inout) :: xmetric
  REAL(kind=real64), intent(inout) :: allvars(xmetric%nx, xmetric%ny, nvar)
  REAL start_time, end_time

  REAL(kind=real64) :: dlat, dlon
  INTEGER i, j

! Forcing / velocities
  !Get first set of data and construct the local metric for drifting
  !timing CALL cpu_time(start_time)
  CALL readin(xmetric%nx, xmetric%ny, nvar, ncid, varid, allvars)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'sub initial_read time for readin',end_time - start_time

!2ds_ice:
  xmetric%ulat = allvars(:,:,2)
  xmetric%ulon = allvars(:,:,1)
  !debug: PRINT *,'ulon nc readin ',MAXVAL(xmetric%ulon), MINVAL(xmetric%ulon)
  WHERE(xmetric%ulon > 720) xmetric%ulon = xmetric%ulon - 720
  WHERE(xmetric%ulon > 360) xmetric%ulon = xmetric%ulon - 360
  !debug: PRINT *,'ulat nc readin ',MAXVAL(xmetric%ulat), MINVAL(xmetric%ulat)
  !debug: PRINT *,'ulon nc readin after cleanup',MAXVAL(xmetric%ulon), MINVAL(xmetric%ulon)
  !debug -- regular latlon grid:
  !dlat = 180./xmetric%ny
  !dlon = 360./xmetric%nx
  !DO j = 1, xmetric%ny
  !  xmetric%ulat(:,j) = j*dlat - 90.0 - dlat/2.
  !ENDDO
  !DO i = 1, xmetric%nx
  !  xmetric%ulon(i,:) = i*dlon - dlon/2.
  !ENDDO
  !end debug

  !timing CALL cpu_time(start_time)
  CALL xmetric%local_metric()
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'sub initial_read time for local_metric',end_time - start_time

END SUBROUTINE initial_read
!----------------------------------------------------------------
SUBROUTINE initialize_drifters(nvar_drift, drift_name, ncid_drift, &
                              varid_drift, nbuoy, restart)
  IMPLICIT none
  INTEGER nvar_drift, ncid_drift, varid_drift(nvar_drift)
  INTEGER nbuoy
  CHARACTER(90), intent(in) :: drift_name
  LOGICAL, intent(in) :: restart

  INTEGER i, retcode
  CHARACTER(50) varnames(nvar_drift), dimname

  !debug: PRINT *,'entered drifter initialize'
  dimname = 'nbuoy'
  retcode = nf90_open(drift_name, NF90_NOWRITE, ncid_drift)
  CALL check(retcode)

  IF (restart) THEN
    !debug: PRINT *,'running from warm start'
    varnames(1) = 'Initial_Latitude'
    varnames(2) = 'Initial_Longitude'
    varnames(3) = 'Final_Latitude'
    varnames(4) = 'Final_Longitude'
  ELSE
    !debug: PRINT *,'running from cold start'
    varnames(1) = 'Initial_Latitude'
    varnames(2) = 'Initial_Longitude'
  ENDIF

  DO i = 1, nvar_drift
    retcode = nf90_inq_varid(ncid_drift, varnames(i), varid_drift(i))
    CALL check(retcode)
    !debug: PRINT *,'initialize_drifters',retcode, i, varnames(i), varid_drift(i)
  ENDDO
    
  retcode = nf90_inquire_dimension(ncid_drift, 1, dimname, nbuoy)
  CALL check(retcode)
  !debug: PRINT *,'initialize -- nbuoy ', nbuoy

  RETURN
END SUBROUTINE initialize_drifters

SUBROUTINE readin_drifters(nbuoy, nvar_drift, ncid_drift, varid_drift, buoylist, xmetric, restart)
  USE constants
  USE metric_mod
  USE, intrinsic :: ieee_arithmetic
  IMPLICIT none

  INTEGER nvar_drift, ncid_drift, varid_drift(nvar_drift)
  INTEGER nbuoy
  CLASS(drifter) :: buoylist(nbuoy)
  TYPE(metric)  :: xmetric
  LOGICAL, intent(in) :: restart

  INTEGER i, retcode, bad_count, very_bad
  REAL(kind=real64) :: tlon(nbuoy), tlat(nbuoy)
  REAL(kind=real64) :: clon(nbuoy), clat(nbuoy)
  REAL(kind=real64), allocatable :: bad_lat(:), bad_lon(:)
  REAL(kind=real64), allocatable :: bad_fi(:), bad_fj(:)
  INTEGER, allocatable :: bad_index(:)

  REAL start_time, end_time

  !debug: PRINT *,' entered drifter read in'
  !debug: PRINT *,ncid_drift, varid_drift
  retcode = nf90_get_var(ncid_drift, varid_drift(1), tlat)
  CALL check(retcode)
  !debug: PRINT *,'lat ',MAXVAL(tlat), MINVAL(tlat)

  retcode = nf90_get_var(ncid_drift, varid_drift(2), tlon)
  CALL check(retcode)
  !debug: PRINT *,'lon ',MAXVAL(tlon), MINVAL(tlon)

  IF (.not. restart) THEN
    clat = tlat
    clon = tlon
  ELSE
    retcode = nf90_get_var(ncid_drift, varid_drift(3), clat)
    CALL check(retcode)
    !debug: PRINT *,'clat ',MAXVAL(clat), MINVAL(clat)
    retcode = nf90_get_var(ncid_drift, varid_drift(4), clon)
    CALL check(retcode)
    !debug: PRINT *,'clon ',MAXVAL(clon), MINVAL(clon)
  ENDIF

  !debug: PRINT *,' about to create buoys '
  bad_count = 0
  !timing CALL cpu_time(start_time)
  DO i = 1, nbuoy
    !debug2: PRINT *,'initializing buoy ',i,tlat(i), tlon(i)
    CALL buoylist(i)%init(tlon(i), tlat(i), clon(i), clat(i), xmetric)
    if (clat(i) >= flag .or. clon(i) >= flag .or. &
        buoylist(i)%x >= flag .or. buoylist(i)%y >= flag ) THEN
      bad_count = bad_count + 1
    endif
    !debug2: WRITE(*,9001) i, tlon(i), tlat(i), clon(i), clat(i)
  ENDDO
 9001 FORMAT(I6,4F10.3)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'sub readin_drifters buoy list timing ',end_time - start_time

! Processing en masse all buoys which have bad locations 
! RG: make this its own routine for general use
  !debug:
  PRINT *,'count of bad locations: ',bad_count
  ALLOCATE(bad_index(bad_count), bad_lat(bad_count), bad_lon(bad_count))
  ALLOCATE(bad_fi(bad_count), bad_fj(bad_count))
  bad_count = 1
  DO i = 1, nbuoy
    if (clat(i) >= flag .or. clon(i) >= flag .or. &
        buoylist(i)%x >= flag .or. buoylist(i)%y >= flag ) THEN
      bad_index(bad_count) = i
      IF (tlat(i) >= flag .or. tlon(i) >= flag) THEN
        tlat(i) = 0.
        tlon(i) = 0.
      ENDIF
      bad_lat(bad_count)   = tlat(i)
      bad_lon(bad_count)   = tlon(i)
      bad_count = bad_count + 1
    endif
  ENDDO
  bad_count = bad_count - 1 
  !RG: now call mass searcher irreg_
  !timing CALL cpu_time(start_time)
!  CALL irreg_ll2ij_cice(xmetric%nx, xmetric%ny, xmetric%ulat, xmetric%ulon, &
!          bad_count, bad_lat, bad_lon, bad_fi, bad_fj)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'sub readin_drifters irreg timing ',end_time - start_time

  very_bad = 0
  DO i = 1, bad_count
    !debug: PRINT *,'retry ',i,bad_fi(i), bad_fj(i), bad_lat(i), bad_lon(i)
    IF (bad_fi(i) < 1 .or. bad_fi(i) >= flag .or. ieee_is_nan(bad_fi(i)) .or. &
        bad_fj(i) < 1 .or. bad_fj(i) >= flag .or. ieee_is_nan(bad_fj(i)) ) THEN
      PRINT *,'could not place ',buoylist(bad_index(i))%ilat, buoylist(bad_index(i))%ilon
      buoylist(bad_index(i))%x = flag
      buoylist(bad_index(i))%y = flag
      buoylist(bad_index(i))%clat = flag
      buoylist(bad_index(i))%clon = flag
      very_bad = very_bad + 1
    ELSE
      buoylist(bad_index(i))%x = bad_fi(i)
      buoylist(bad_index(i))%y = bad_fj(i)
      buoylist(bad_index(i))%clat = bad_lat(i)
      buoylist(bad_index(i))%clon = bad_lon(i)
    ENDIF
  ENDDO
  !debug: PRINT *,'could not place points: ',very_bad
  
  !debug: PRINT *,' leaving drifter read in'
  RETURN
END SUBROUTINE readin_drifters
!----------------------------------------------------------------
SUBROUTINE readin(nx, ny, nvars, ncid, varid, allvars)
  IMPLICIT none
  INTEGER, intent(in)    :: nvars
  INTEGER, intent(in)    :: nx, ny, ncid, varid(nvars)
  REAL(kind=real64), intent(inout)    :: allvars(nx, ny, nvars)
  
  INTEGER i, retcode
  
  !got nx, ny from the .nc file, in initialize_in
  DO i = 1, nvars
    retcode = nf90_get_var(ncid, varid(i), allvars(:,:,i) )
    CALL check(retcode)
    !debug: PRINT *,'readin',retcode, i, MAXVAL(allvars(:,:,i)), MINVAL(allvars(:,:,i))
    if (retcode < 0) then
      PRINT *,'error on ',ncid, varid(i)
      STOP
    endif
  ENDDO

  RETURN
END

!----------------------------------------------------------------
!utility for netcdf
SUBROUTINE check(status)
  IMPLICIT none
  INTEGER, intent(in) :: status
  IF (status /= nf90_noerr) THEN
    PRINT *,nf90_strerror(status)
    !debug: STOP "erredout"
    PRINT *, "erredout"
  ENDIF
  RETURN
END subroutine check

!----------------------------------------------------------------
SUBROUTINE initialize_out(fname, ncid, varid, nvar, nbuoy, dimids)
  IMPLICIT none
  CHARACTER(*), intent(in) :: fname
  INTEGER, intent(in) :: nvar
  INTEGER, intent(out) :: varid(nvar)
  INTEGER ncid, nbuoy

  CHARACTER(20) varnames(nvar)
  INTEGER, intent(inout) :: dimids(1)

  INTEGER i, retcode
  INTEGER x_dimid
!names: 
  varnames(1) = "Initial_Latitude"
  varnames(2) = "Initial_Longitude"
  varnames(3) = "Final_Latitude"
  varnames(4) = "Final_Longitude"
  varnames(5) = "Drift_Distance"
  varnames(6) = "Drift_Bearing"
  
! open
  retcode = nf90_create(fname, NF90_NOCLOBBER, ncid)
  CALL check(retcode)
  
! dimensionalize
  retcode = nf90_def_dim(ncid, "nbuoy", nbuoy, x_dimid)
  CALL check(retcode)
  dimids = (/ x_dimid /)

! assign varid to varnames
  DO i = 1, nvar
    retcode = nf90_def_var(ncid, trim(varnames(i)), NF90_DOUBLE, dimids, varid(i))
    CALL check(retcode)
  ENDDO

  retcode = nf90_enddef(ncid)
  CALL check(retcode)

  RETURN
END subroutine initialize_out

!----------------------------------------------------------------
SUBROUTINE outvars(ncid, varid, nvar, buoys, nbuoy)
  IMPLICIT none

  INTEGER, intent(in) :: ncid, nvar
  INTEGER, intent(in) :: varid(nvar)
  INTEGER, intent(in) :: nbuoy
  TYPE(drifter), intent(in) :: buoys(nbuoy)

  INTEGER retcode
  REAL(kind=real64), allocatable :: var(:,:)
  INTEGER i, k
  REAL(kind=real64) :: wrap, distance, bear

!Note that netcdf dimensions are in C order, not fortran
  !debug: PRINT *,'entered outvars'
  IF (ALLOCATED(var)) THEN
    PRINT *,'outvars deallocating var'
    DEALLOCATE(var)
  ENDIF
  ALLOCATE(var(nbuoy, nvar))
  distance = 0.
  bear = 0.

  DO k = 1, nbuoy
    var(k,1) = buoys(k)%ilat
    var(k,2) = buoys(k)%ilon
    var(k,3) = buoys(k)%clat
    if (buoys(k)%clon >= 360. .or. buoys(k)%clon < 0) THEN
      var(k,4) = wrap(buoys(k)%clon)
    else
      var(k,4) = buoys(k)%clon
    endif
    CALL bearing(var(k,1), var(k,2), var(k,3), var(k,4), distance, bear)
    var(k,5) = distance
    var(k,6) = bear
  ENDDO
  !debug: PRINT *,'ilat ',MAXVAL(var(:,1))
  !debug: PRINT *,'ilon ',MAXVAL(var(:,2))
  !debug: PRINT *,'clat ',MAXVAL(var(:,3))
  !debug: PRINT *,'clon ',MAXVAL(var(:,4))
  !debug: PRINT *,'dist ',MAXVAL(var(:,5))
  !debug: PRINT *,'bear ',MAXVAL(var(:,6))

  !RG: separate initial write -- just ilat, ilon, from later writes, clat, clon, distance, bear
  DO i = 1, nvar
    retcode = nf90_put_var(ncid, varid(i), var(:,i) )
    CALL check(retcode)
  ENDDO

  DEALLOCATE(var)
  !debug: PRINT *,'leaving outvars'
  RETURN
END subroutine outvars

!----------------------------------------------------------------
SUBROUTINE close_out(ncid)
  IMPLICIT none
  INTEGER ncid
  INTEGER retcode

  retcode = nf90_close(ncid)
  CALL check(retcode)

  RETURN
END subroutine close_out

!----------------------------------------------------------------
! Write out results -- drift distance and direction
SUBROUTINE writeout(ncid_out, varid_out, nvar_out, buoys, nbuoy, closeout)
  USE drifter_mod

  IMPLICIT none

  INTEGER, intent(in) ::  ncid_out, nvar_out, nbuoy
  INTEGER, intent(in) ::  varid_out(nvar_out)
  LOGICAL, intent(in) :: closeout
  TYPE(drifter), intent(in) :: buoys(nbuoy)

  CALL outvars(ncid_out, varid_out, nvar_out, buoys, nbuoy )
  IF (closeout) THEN
    !debug: PRINT *,'calling close_out'
    CALL close_out(ncid_out)
  ENDIF

END SUBROUTINE writeout

END module io
