MODULE io
  USE netcdf
  USE drifter_mod
  PUBLIC

CONTAINS

SUBROUTINE initialize_in(nvar, fname, ncid, varid, nx, ny)
  IMPLICIT none

  INTEGER ncid
  INTEGER, intent(in) :: nvar
  CHARACTER(*), intent(in) :: fname
  INTEGER, intent(out) :: nx, ny

! Names from rtofs output
  CHARACTER(len=40) :: varnames(nvar)
  INTEGER varid(nvar)
  
! For Netcdf processing
  INTEGER i
  INTEGER retcode
  CHARACTER(len=50) :: xname, yname

!debug:
  PRINT *,'entered initialize_in'
!! This is the cice_inst variable set -- much more extensive
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
  varnames(1) = "Longitude"
  varnames(2) = "Latitude"
  varnames(3) = "ice_coverage"
  varnames(4) = "ice_temperature"
  varnames(5) = "ice_thickness"
  varnames(6) = "ice_uvelocity"
  varnames(7) = "ice_vvelocity"
  
  PRINT *,'trying to open ',fname, len(fname)
  retcode = nf90_open(fname, NF90_NOWRITE, ncid)
  CALL check(retcode)

  DO i = 1, nvar
    PRINT *,i,' reading varname ',varnames(i)
    retcode = nf90_inq_varid(ncid, varnames(i), varid(i))
    CALL check(retcode)
  ENDDO
  PRINT *,'done reading varnames'

  !RG: Read this in from netcdf file
  !debug: nx = 4500
  !debug: ny = 3298
  xname = "X"
  yname = "Y"
  retcode = nf90_inquire_dimension(ncid, 3, xname, nx)
  CALL check(retcode)
  retcode = nf90_inquire_dimension(ncid, 2, yname, ny)
  CALL check(retcode)

  PRINT *,'leaving initialize_in', nx, ny

RETURN
END subroutine initialize_in

!----------------------------------------------------------------
SUBROUTINE initial_read(fname, outname, nx, ny, nvar, ncid, varid, &
                        allvars, ulon, ulat, dx, dy, rot, &
                        dlatdi, dlatdj, dlondi, dlondj,   &
                        outdimids, ncid_out, varid_out, nvar_out, nbuoy)

  USE drifter_mod
  IMPLICIT none

  INTEGER, intent(in) :: nx, ny
  CHARACTER(90), intent(in) :: fname, outname
  INTEGER, intent(in)  :: nvar, nvar_out, ncid, varid(nvar)
  INTEGER, intent(out) :: outdimids(1)
  INTEGER, intent(out) :: ncid_out, varid_out(nvar_out)

  INTEGER ncid_drift, nbuoy
  INTEGER varid_driftic(nvar_out)
  
  REAL, intent(inout) :: allvars(nx, ny, nvar)
  REAL, intent(inout) :: ulat(nx, ny), ulon(nx, ny)
  REAL, intent(out)   :: dx(nx, ny), dy(nx, ny), rot(nx, ny)
  REAL, intent(out)   :: dlatdi(nx, ny), dlatdj(nx, ny), dlondi(nx, ny), dlondj(nx, ny)

! Locals:
!  INTEGER i, j, k

! Allocate space for variables and initialize the netcdf reading

! Forcing / velocities
  !Get first set of data and construct the local metric for drifting
  CALL readin(nx, ny, nvar, ncid, varid, allvars)

!cice_inst:
!  ulat = allvars(:,:,4)
!  ulon = allvars(:,:,3)
!2ds_ice:
  ulat = allvars(:,:,2)
  ulon = allvars(:,:,1)
  !old CALL local_metric(ulat, ulon, dx, dy, rot, nx, ny)
  CALL local_metric(ulat, ulon, dx, dy, rot, dlatdi, dlondi, dlatdj, dlondj, nx, ny)


!  !----------- Initialize buoys, this should be a read in -- separate function

! Initialize Output -- need definite sizes
  CALL initialize_out(outname, ncid_out, varid_out, nvar_out, nbuoy, outdimids)

END SUBROUTINE initial_read
!----------------------------------------------------------------
SUBROUTINE initialize_drifters(nvar_drift, drift_name, ncid_drift, varid_drift, nbuoy)
  IMPLICIT none
  INTEGER nvar_drift, ncid_drift, varid_drift(nvar_drift)
  INTEGER nbuoy
  CHARACTER(90) drift_name

  INTEGER i, retcode
  CHARACTER(50) varnames(nvar_drift), dimname

  PRINT *,'entered drifter initialize'
  varnames(1) = 'Initial_Latitude'
  varnames(2) = 'Initial_Longitude'
  dimname = 'nbuoy'

  retcode = nf90_open(drift_name, NF90_NOWRITE, ncid_drift)
  CALL check(retcode)

  DO i = 1, nvar_drift
    retcode = nf90_inq_varid(ncid_drift, varnames(i), varid_drift(i))
    CALL check(retcode)
  ENDDO

  retcode = nf90_inquire_dimension(ncid_drift, 1, dimname, nbuoy)
  CALL check(retcode)
  PRINT *,'initialize -- nbuoy ', nbuoy

  RETURN
END SUBROUTINE initialize_drifters

SUBROUTINE readin_drifters(nbuoy, nvar_drift, ncid_drift, varid_drift, buoylist)
  IMPLICIT none
  INTEGER nvar_drift, ncid_drift, varid_drift(nvar_drift)
  INTEGER nbuoy
  CLASS(drifter) :: buoylist(nbuoy)
  INTEGER i, retcode
  REAL tlon(nbuoy), tlat(nbuoy)

  PRINT *,' entered drifter read in'
  PRINT *,ncid_drift, varid_drift
  retcode = nf90_get_var(ncid_drift, varid_drift(1), tlat)
  CALL check(retcode)
  PRINT *,'lat ',MAXVAL(tlat), MINVAL(tlat)

  retcode = nf90_get_var(ncid_drift, varid_drift(2), tlon)
  CALL check(retcode)
  PRINT *,'lon ',MAXVAL(tlon), MINVAL(tlon)

  PRINT *,' about to create buoys '
  DO i = 1, nbuoy
    CALL buoylist(i)%init(tlon(i), tlat(i))
    PRINT *,i, buoylist(i)%ilat, buoylist(i)%ilon
  ENDDO
  
  PRINT *,' leaving drifter read in'
  RETURN
END SUBROUTINE readin_drifters
!----------------------------------------------------------------
SUBROUTINE readin(nx, ny, nvars, ncid, varid, allvars)
  IMPLICIT none
  INTEGER, intent(in)    :: nvars
  INTEGER, intent(in)    :: nx, ny, ncid, varid(nvars)
  REAL, intent(inout)    :: allvars(nx, ny, nvars)
  
  INTEGER i, retcode
  
  !got nx, ny from the .nc file, in initialize_in

  DO i = 1, nvars
    retcode = nf90_get_var(ncid, varid(i), allvars(:,:,i) )
    CALL check(retcode)
    PRINT *,i, MAXVAL(allvars(:,:,i)), MINVAL(allvars(:,:,i))
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
    STOP "erredout"
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
    retcode = nf90_def_var(ncid, trim(varnames(i)), NF90_REAL, dimids, varid(i))
    CALL check(retcode)
  ENDDO

  retcode = nf90_enddef(ncid)
  CALL check(retcode)


  RETURN
END subroutine initialize_out

!----------------------------------------------------------------
SUBROUTINE outvars(ncid, varid, nvar, buoys, nbuoy)
  IMPLICIT none

  INTEGER ncid, nvar
  INTEGER varid(nvar)
  INTEGER nbuoy
  TYPE(drifter) :: buoys(nbuoy)

  INTEGER retcode
  REAL, allocatable :: var(:,:)
  INTEGER i, k
  REAL distance, bear

!Note that netcdf dimensions are in C order, not fortran
  ALLOCATE(var(nbuoy, nvar))
  distance = 0.
  bear = 0.

  DO k = 1, nbuoy
    var(k,1) = buoys(k)%ilat
    var(k,2) = buoys(k)%ilon
    var(k,3) = buoys(k)%clat
    var(k,4) = buoys(k)%clon
    CALL bearing(var(k,1), var(k,2), var(k,3), var(k,4), distance, bear)
    var(k,5) = distance
    var(k,6) = bear
  ENDDO

  !RG: separate initial write -- just ilat, ilon, from later writes, clat, clon, distance, bear
  DO i = 1, nvar
    retcode = nf90_put_var(ncid, varid(i), var(:,i) )
    CALL check(retcode)
  ENDDO

  DEALLOCATE(var)

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

  TYPE(drifter) :: buoys(nbuoy)
  LOGICAL, intent(in) :: closeout

  IF (closeout) THEN
    CALL close_out(ncid_out)
  ELSE
    CALL outvars(ncid_out, varid_out, nvar_out, buoys, nbuoy )
  ENDIF

END SUBROUTINE writeout

END module io
