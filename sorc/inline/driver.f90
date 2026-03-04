PROGRAM driver

  USE drifter_mod
  USE io
  USE metric_mod
  IMPLICIT none

! Not needed by drifter inline
  INTEGER nvar
  PARAMETER (nvar = 7)
  INTEGER ncid, varid(nvar), nx, ny
  CHARACTER(len=40) :: varnames(nvar)
  CHARACTER(len=50) :: xname, yname
  REAL(kind=real64), allocatable  :: allvars(:,:,:)
  CHARACTER(300) fname, tmp, tmp2
  INTEGER i

! Arguments to inline drifters
  TYPE(metric) :: xmetric
  CHARACTER(300) drift_in, drift_out
  REAL(kind=real64) :: dt
  LOGICAL restart, closeout
  REAL(kind=real64), allocatable :: u(:,:), v(:,:), aice(:,:)
  INTEGER phase

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
  READ (10,*) restart

! Read in .nc variable names
  DO i = 1, nvar
    READ (10,*) varnames(i)
    ENDDO
  READ (10,*) xname  ! x,y dimensions
  READ (10,*) yname

! Initialize input Forcing / velocities
  CALL initialize_in(nvar, trim(fname), ncid, varid, varnames, xname, yname, nx, ny)

  ALLOCATE(allvars(nx, ny, nvar))
  ALLOCATE(u(nx, ny), v(nx, ny), aice(nx, ny) )
  CALL xmetric%set(nx, ny)
  CALL initial_read(trim(fname), nvar, ncid, varid, &
                    allvars, xmetric)
  u = allvars(:,:,6)
  v = allvars(:,:,7)
  aice = allvars(:,:,3)

! now execute calls for inline usage -------------------------------------
  closeout = .FALSE.
  phase = 1
  CALL driftmain(dt, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, &
        phase )
  
  phase = 2
  CALL driftmain(dt, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, &
        phase )

  closeout = .TRUE.
  phase = 3
  CALL driftmain(dt, restart, closeout, &
        trim(drift_in), trim(drift_out), &
        xmetric, u, v, aice, &
        phase )


END PROGRAM driver
