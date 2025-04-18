SUBROUTINE driftmain(dt, dtout, outfreq, restart, closeout, &
        drift_in, drift_out, &
        u, v, aice, nx, ny,  &
        phase )
  USE drifter_mod
  USE io
  USE metric_mod

  IMPLICIT none
!arguments:
  REAL, intent(in)    :: dt, dtout
  INTEGER, intent(in) :: outfreq
  LOGICAL, intent(in) :: restart
  LOGICAL, intent(in) :: closeout
  INTEGER, intent(in) :: nx, ny
  REAL, intent(in)    :: u(nx, ny), v(nx, ny)
  REAL, intent(in)    :: aice(nx, ny)
  INTEGER,intent(in)  :: phase
  CHARACTER(*), intent(in) :: drift_in, drift_out

! Netcdf-related
  INTEGER nvar_out, nvar_drift
  PARAMETER (nvar_out = 6)
  INTEGER ncid_out, ncid_drift
  INTEGER varid_out(nvar_out)
  INTEGER, allocatable :: varid_drift(:)
  INTEGER dimids(1)
  
! Read from input (or argument to main)
  TYPE(metric), save :: xmetric

!For drifter 
  CLASS(drifter), allocatable, save :: buoys(:)
  INTEGER, save        :: nbuoys
  CHARACTER(300), save :: drift_name, outname

! Utilities for main
  INTEGER i, j
  INTEGER n

! -- Begin main for inline usage ----

!-------------------------------------------------------
select case(phase)
  case(1)
!    INIT
  drift_name = trim(drift_in)
  outname = trim(drift_out)

! Run parameters by argument

! Initialize input Forcing / velocities
  CALL initialize_in(nvar, varid, nx, ny, xmetric)

! Initialize Buoy points
! Allocate varid_drift, initialize nvar_drift based on whether this is restart
  IF (restart) THEN
    nvar_drift = 4
    ELSE
    nvar_drift = 2
  ENDIF
  ALLOCATE(varid_drift(nvar_drift))

  CALL initialize_drifters(nvar_drift, drift_name, ncid_drift, varid_drift, nbuoys, restart)
  ALLOCATE(buoys(nbuoys))

! Initialize Output -- need definite sizes
  !RG: really initialize_io
  !Get first set of data and construct the local metric for drifting
  CALL initial_read(outname, nx, ny, nvar, varid, &
                    xmetric, &
                    dimids, ncid_out, varid_out, nvar_out, nbuoys)
  PRINT *,'nbuoys = ', nbuoys

  CALL readin_drifters(nbuoys, nvar_drift, ncid_drift, varid_drift, buoys, xmetric, restart)
!---------------------------------------------------------
! RUN
  case(2)
! First/only time step (u,v, etc. in hand):
  CALL run(buoys, nbuoys, u, v, xmetric, dt, dtout)

!---------------------------------------------------------
! Write
  case(3)
  CALL writeout(ncid_out, varid_out, nvar_out, buoys, nbuoys, closeout)

  end select
END subroutine driftmain
