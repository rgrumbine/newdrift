PROGRAM driver


  SUBROUTINE readin(nx, ny, nvars, ncid, varid, allvars)
  IMPLICIT none
  INTEGER, intent(in)    :: nvars
  INTEGER, intent(in)    :: nx, ny, ncid, varid(nvars)
  REAL, intent(inout)    :: allvars(nx, ny, nvars)
  REAL ulon(nx, ny), ulat(nx, ny)


  TYPE(metric) :: xmetric

  CALL xmetric%set(nx, ny)
  xmetric%ulat = allvars(:,:,2)
  xmetric%ulon = allvars(:,:,1)
  CALL xmetric%local_metric()

