program test_newton_inversion
    use metric_mod
    implicit none
    
    interface
        subroutine newton(this, lat, lon, x, y)
            use metric_mod
            use constants
            class(metric), intent(in)        :: this
            real(kind=real64), intent(in)    :: lat, lon
            real(kind=real64), intent(inout) :: x, y
        end subroutine newton
    end interface

    type(metric) :: spatial_mesh
    real(kind=real64) :: target_lat, target_lon, test_x, test_y
    integer :: i, j

    print *, "=== RUNNING REGRESSION TEST: newton.f90 ==="

    ! 1. Instantiate testing landscape grid framework
    call spatial_mesh%set(10, 10)
    do j = 1, 10
        do i = 1, 10
            spatial_mesh%ulat(i, j) = 10.0_real64 + real(j - 1, 8) * 1.0_real64
            spatial_mesh%ulon(i, j) = -80.0_real64 + real(i - 1, 8) * 1.0_real64
        end do
    end do
    call spatial_mesh%local_metric()

    ! 2. Target coordinates known to match exact location tracking node index i=4, j=4
    target_lat = 13.0_real64
    target_lon = -77.0_real64
    
    ! Provide initial rough approximate guess
    test_x = 3.2_real64
    test_y = 3.1_real64

    ! 3. Call execution
    call newton(spatial_mesh, target_lat, target_lon, test_x, test_y)

    ! 4. Evaluate convergence matching conditions
    if (abs(test_x - 4.0_real64) > 1.0e-3 .or. abs(test_y - 4.0_real64) > 1.0e-3) then
        print *, "[FAIL] Newton optimization failed convergence targets. Obtained: ", test_x, test_y
        stop 1
    else
        print *, "[PASS] Newton solver accurately converged to true coordinates."
    end if

end program test_newton_inversion