PROGRAM newdrift

  USE drifter_mod
  USE io
  USE metric_mod

  IMPLICIT none

! Netcdf-related
  INTEGER nvar, nvar_out, nvar_drift
  !cice_inst: PARAMETER (nvar = 24)
  PARAMETER (nvar = 7)
  PARAMETER (nvar_out = 6)
  PARAMETER (nvar_drift = 2)
  INTEGER ncid, varid(nvar)
  INTEGER ncid_out, ncid_drift

  INTEGER varid_out(nvar_out), varid_drift(nvar_drift)
  INTEGER dimids(1)
  
! Read from input (or argument to main)
  REAL dt, dtout
  INTEGER outfreq
  LOGICAL restart

  TYPE(metric) :: xmetric

  REAL, allocatable  :: allvars(:,:,:)
  REAL, allocatable  :: u(:,:), v(:,:)
  REAL, allocatable  :: aice(:,:)

! Utilities for main
  INTEGER i, j
  INTEGER n, nstep
!  REAL x, y

!For drifter 
  CLASS(drifter), allocatable :: buoys(:)
  INTEGER nbuoys
  CHARACTER(300) drift_name

! Read from .nc file (or argument to main)
  INTEGER nx, ny
! Names
  CHARACTER(300) fname, outname, tmp, tmp2
  LOGICAL closeout


! -- Begin main for offline ----
  READ (*,*) tmp
  fname = trim(tmp)
  !debug: PRINT *,"name of file to get inputs from",fname
!debug:  STOP

  OPEN(10, FILE=fname, FORM='FORMATTED', STATUS='OLD')

  READ (10,*) tmp2
  fname = trim(tmp2)
  !debug: PRINT *,fname,len(fname), len(trim(tmp2))
  READ (10,*) tmp2
  drift_name = trim(tmp2)
  !debug: PRINT *,drift_name,len(drift_name)
  READ (10,*) tmp
  outname = trim(tmp)
  !debug: PRINT *,outname
  !debug: PRINT *,fname, drift_name, outname

! Read in run parameters:
  READ (10,*) dt
  READ (10,*) nstep
  READ (10,*) outfreq
  READ (10,*) restart
  PRINT *,'dt, nstep, outfreq, restart = ',dt, nstep, outfreq, restart

! Initialize input Forcing / velocities
  CALL initialize_in(nvar, trim(fname), ncid, varid, nx, ny, xmetric)
!debug: PRINT *,"main ",nvar, trim(fname), ncid, varid, nx, ny

! Initialize Buoy points
  !debug: PRINT *,'calling initialize_drifters'
  CALL initialize_drifters(nvar_drift, drift_name, ncid_drift, varid_drift, nbuoys, restart)
  ALLOCATE(buoys(nbuoys))
  !debug: PRINT *,'back from initialize_drifters, nbuoys = ',nbuoys

! Initialize Output -- need definite sizes
  !debug: PRINT *,'allocating input variables'
  ALLOCATE(allvars(nx, ny, nvar))
  ALLOCATE(aice(nx, ny))

  !RG: really initialize_io
  !Get first set of data and construct the local metric for drifting
  !debug: PRINT *,'calling initial read '
  CALL initial_read(trim(fname), outname, nx, ny, nvar, ncid, varid, &
                    allvars, xmetric, &
                    dimids, ncid_out, varid_out, nvar_out, nbuoys)
  PRINT *,'initial read results'
!debug: PRINT *,trim(fname), drift_name, outname, nx, ny, nvar, ncid, varid
  PRINT *,MAXVAL(xmetric%ulon), MAXVAL(xmetric%ulat), MAXVAL(xmetric%dx), MAXVAL(xmetric%dy)
  PRINT *,dimids, ncid_out, varid_out, nvar_out
  PRINT *,'nbuoys = ', nbuoys

!cice_inst:
!  aice = allvars(:,:,8)
!  u    = allvars(:,:,9)
!  v    = allvars(:,:,10)

!2ds_ice (not _prog or _diag)
  aice = allvars(:,:,3)
  u    = allvars(:,:,6)
  v    = allvars(:,:,7)

  !debug: PRINT *,'in main, parcelling out aice, u, v'
  !debug: PRINT *,"aice",MAXVAL(aice), MINVAL(aice)
  !debug: PRINT *,"u",MAXVAL(u), MINVAL(u)
  !debug: PRINT *,"v",MAXVAL(v), MINVAL(v)

  CALL readin_drifters(nbuoys, nvar_drift, ncid_drift, varid_drift, buoys, xmetric, restart)
  !debug: PRINT *,'nbuoys after readin_drifters',nbuoys

  !debug STOP
!---------------------------------------------------------
! RUN

! First/only time step (u,v, etc. in hand):
  CALL run(buoys, nbuoys, u, v, xmetric, dt, dtout)
  CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)
  closeout = .TRUE.
  CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)

! Iterate as needed:
!  !RG: need to set up multi-level time outputs in io first
!  DO n = 2, nstep
!    CALL readin(nx, ny, nvar, ncid, varid, allvars)
!    u = allvars(:,:,6)
!    v = allvars(:,:,7)
!    CALL run(buoys, nbuoys, u, v, xmetric, dt, dtout)
!    CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)
!  ENDDO
!  closeout = .TRUE.
!  CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)

!----------------------------------------------------------------
! WRITE Write out results -- drift distance and direction

END program newdrift
