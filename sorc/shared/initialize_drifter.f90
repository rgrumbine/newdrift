! Use the drifter class
!SUBROUTINE zero_buoys(buoys, nbuoy)
!  USE drifter_mod
!
!  IMPLICIT none
!
!  INTEGER, intent(in) :: nbuoy
!  TYPE(drifter), intent(inout) :: buoys(nbuoy)
!
!  INTEGER k
!
!  !Zero the buoys:
!  DO k = 1, nbuoy
!    CALL buoys(k)%zero(k)
!  ENDDO
!  RETURN
!END SUBROUTINE zero_buoys

SUBROUTINE dummy_buoys(aice, dx, dy, u, v, ulat, ulon, nx, ny, ratio, buoys, nbuoys, ngood)
  USE drifter_mod
  IMPLICIT none
  INTEGER, intent(in) :: nx, ny, ratio
  REAL, intent(in) :: aice(nx, ny), dx(nx, ny), dy(nx, ny)
  REAL, intent(in) :: u(nx, ny), v(nx, ny)
  REAL, intent(in) :: ulat(nx, ny), ulon(nx, ny)
  INTEGER, intent(inout) :: nbuoys, ngood
  TYPE(drifter), intent(inout) :: buoys(nbuoys)

  INTEGER i, j, k, imax, jmax

  PRINT *,'entered dummy_buoys'
 
! Dummy for testing
  imax = INT(nx/ratio)
  jmax = INT(ny/ratio)

  IF (nbuoys .NE. imax*jmax) THEN
    PRINT *,'mismatch in sizes ',nbuoys,' vs ',imax*jmax
    STOP
  ENDIF
 
  k = 0
  DO j = 1, jmax
  DO i = 1, imax
    !IF (aice(i*ratio, j*ratio) > 0. .AND. aice(i*ratio,j*ratio) <= 1.0  .AND. &
    !    ABS(u(i*ratio, j*ratio)) < 100. .AND. ABS(v(i*ratio, j*ratio)) < 100. ) THEN
    IF ( &
        dx(i*ratio, j*ratio) .NE. 0. .AND. dy(i*ratio, j*ratio) .NE. 0. &
       ) THEN 

      k = k + 1
      buoys(k)%ilat = ulat(i*ratio, j*ratio)
      buoys(k)%ilon = ulon(i*ratio, j*ratio)
      buoys(k)%clat = ulat(i*ratio, j*ratio)
      buoys(k)%clon = ulon(i*ratio, j*ratio)

      !RG: Must compute physical value here, location vs. grid mesh function
      !CALL ll_to_xy(buoys(k)%ilat, buoys(k)%ilon, ulat, ulon, x, y, nx, ny)
      buoys(k)%x = i*ratio
      buoys(k)%y = j*ratio

    ENDIF
  ENDDO
  ENDDO

  ngood = k
  RETURN
END SUBROUTINE dummy_buoys
