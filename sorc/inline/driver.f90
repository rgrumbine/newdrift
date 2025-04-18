PROGRAM driver

  USE drifter_mod
  USE io
  USE metric_mod
  IMPLICIT none

  INTEGER nvar
  PARAMETER (nvar = 7)

  INTEGER ncid, varid(nvar), nx, ny
  REAL, allocatable  :: allvars(:,:,:)
  TYPE(metric) :: xmetric

  CHARACTER(300) fname, tmp, tmp2
  REAL dt
  INTEGER nstep, outfreq
  LOGICAL restart

  REAL dtout
  REAL, allocatable :: u(:,:), v(:,:), aice(:,:)
  INTEGER phase
  LOGICAL closeout
  CHARACTER(300) drift_in, drift_out

! RTOFS et al. files, not inlineable --------------------------------
  READ (*,*) tmp
  fname = trim(tmp)
  OPEN(10, FILE=fname, FORM='FORMATTED', STATUS='OLD')
  READ (10,*) tmp2
  fname = trim(tmp2)
  READ (10,*) tmp2
  drift_in = trim(tmp2)
  READ (10,*) tmp2
  drift_out = trim(tmp2)
! Read in run parameters:
  READ (10,*) dt
  READ (10,*) nstep
  READ (10,*) outfreq
  READ (10,*) restart

! Initialize input Forcing / velocities
  CALL initialize_in(nvar, trim(fname), ncid, varid, nx, ny)
  PRINT *,'initial nx ny ',nx, ny

  ALLOCATE(allvars(nx, ny, nvar))
  ALLOCATE(u(nx, ny), v(nx, ny), aice(nx, ny) )
  CALL initial_read(trim(fname), nx, ny, nvar, ncid, varid, &
                    allvars, xmetric )
  u = allvars(:,:,6)
  v = allvars(:,:,7)
  aice = allvars(:,:,3)
  PRINT *,'xmetric nx ny ', xmetric%nx, xmetric%ny
  nx = xmetric%nx
  ny = xmetric%ny

! now execute calls for inline usage -------------------------------------
  closeout = .FALSE.
  phase = 1
  CALL driftmain(dt, dtout, outfreq, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, nx, ny,  &
        phase )
  PRINT *,'nx ny after drift init ',nx, ny
  
  phase = 2
  CALL driftmain(dt, dtout, outfreq, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, nx, ny,  &
        phase )
  PRINT *,'nx ny after drift run ',nx, ny

  closeout = .TRUE.
  phase = 3
  CALL driftmain(dt, dtout, outfreq, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, nx, ny,  &
        phase )


END PROGRAM driver
