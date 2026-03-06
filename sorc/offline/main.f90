PROGRAM newdrift

  USE constants
  USE drifter_mod
  USE io
  USE metric_mod

  IMPLICIT none

! Netcdf-related
  INTEGER nvar, nvar_out, nvar_drift
  !cice_inst: PARAMETER (nvar = 24)
  PARAMETER (nvar = 7)
  PARAMETER (nvar_out = 6)
  INTEGER ncid, varid(nvar)
  INTEGER ncid_out, ncid_drift

  INTEGER varid_out(nvar_out)
  
  INTEGER, allocatable :: varid_drift(:)
  INTEGER dimids(1)

  CHARACTER(len=40) :: varnames(nvar)
  CHARACTER(len=50) :: xname, yname
  
! Read from input (or argument to main)
  REAL(kind=real64) dt, dtout
  INTEGER outfreq
  LOGICAL restart

  TYPE(metric) :: xmetric

  REAL(kind=real64), allocatable  :: allvars(:,:,:)
  REAL(kind=real64), allocatable  :: u(:,:), v(:,:), aice(:,:) ! u,v,aice are extracted from allvars

! Utilities for main
  INTEGER i, j
  INTEGER n, nstep
  REAL start_time, end_time

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
  !timing CALL cpu_time(start_time)
  READ (*,*) tmp
  fname = trim(tmp)
  OPEN(10, FILE=fname, FORM='FORMATTED', STATUS='OLD')

  READ (10,*) tmp2
  fname = trim(tmp2)
  READ (10,*) tmp2
  drift_name = trim(tmp2)
  READ (10,*) tmp
  outname = trim(tmp)

! Read in run parameters:
  READ (10,*) dt
  READ (10,*) nstep
  READ (10,*) outfreq
  READ (10,*) restart
  !debug: PRINT *,'dt, nstep, outfreq, restart = ',dt, nstep, outfreq, restart

! RG: Read in .nc variable names
  DO i = 1, nvar
    READ (10,*) varnames(i)
  ENDDO
  READ (10,*) xname  ! x,y dimensions
  READ (10,*) yname
  !debug: PRINT *,'xname, yname',xname, yname


! RTOFS et al. files, not inlineable --------------------------------
! Initialize input Forcing / velocities
  !debug: PRINT *,'calling initialize_in'
! RG: add varnames to arg list
  CALL initialize_in(nvar, trim(fname), ncid, varid, varnames, xname, yname, nx, ny)
!RG: really initialize_io
  !Get first set of data and construct the local metric for drifting
  ! also constructs xmetric
  !debug: PRINT *,'allocating input variables'
  ALLOCATE(allvars(nx, ny, nvar))
  ALLOCATE(aice(nx, ny), u(nx, ny), v(nx, ny))
  CALL xmetric%set(nx, ny)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'time through xmetric set ',end_time - start_time
  !timing start_time = end_time

  !debug: PRINT *,'calling initial read '
  CALL initial_read(trim(fname), nvar, ncid, varid, &
                    allvars, xmetric )
!2ds_ice (not _prog or _diag)
  aice = allvars(:,:,3)
  u    = allvars(:,:,6)
  v    = allvars(:,:,7)
  !debug: PRINT *,'returned from initial read '
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'time for initial read ',end_time - start_time
  !timing start_time = end_time

  !debug: STOP
!-------------------------- Buoys, inlineable ---------------------
! Initialize Buoy -- points, input file, output file
! Allocate varid_drift, initialize nvar_drift based on whether this is restart
  IF (restart) THEN
    nvar_drift = 4
    ELSE
    nvar_drift = 2
  ENDIF

  ALLOCATE(varid_drift(nvar_drift))
  !debug: PRINT *,'calling initialize_drifters, nbuoys = ',nbuoys
  CALL initialize_drifters(nvar_drift, drift_name, ncid_drift, varid_drift, nbuoys, restart)
  ALLOCATE(buoys(nbuoys))
  !debug: PRINT *,'back from initialize_drifters, nbuoys = ',nbuoys
  !STOP
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'time for initialize_drifters ',end_time - start_time
  !timing start_time = end_time

  ! For buoy output -- inlineable
  CALL initialize_out(outname, ncid_out, varid_out, nvar_out, nbuoys, dimids)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'time for initialize_out ',end_time - start_time
  !timing start_time = end_time
  !debug: PRINT *,'back from initialize out '
  !debug: PRINT *,'nbuoys = ', nbuoys

  CALL readin_drifters(nbuoys, nvar_drift, ncid_drift, varid_drift, buoys, xmetric, restart)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'readin_drifters time ',end_time - start_time
  !debug: PRINT *,'back from readin_drifters'
  !debug: STOP
!---------------------------------------------------------
! RUN

! First/only time step (u,v, etc. in hand):
  !timing CALL cpu_time(start_time)
  !debug 0.2778 ~= 1 km/hr: 
  !u = 0.2778
  !v = 0.2778 
  !debug: PRINT *,'calling run'
  CALL run(buoys, nbuoys, u, v, xmetric, dt)
  !debug: PRINT *,'back from run'
  closeout = .TRUE.
  CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)
  !timing CALL cpu_time(end_time)
  !timing PRINT *,'run, write time ',end_time - start_time

!----------------------------------------------------------------
! WRITE Write out results -- drift distance and direction

END program newdrift
