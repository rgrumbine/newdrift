program test_metric_matrices
    use metric_mod
    implicit none

    type(metric) :: test_mesh
    integer, parameter :: mx = 5, my = 5
    integer :: i, j

    print *, "=== RUNNING REGRESSION TEST: metric.f90 ==="

    ! 1. Memory instantiation and mapping
    call test_mesh%set(mx, my)
    
    ! Populate linear spacing coordinates
    do j = 1, my
        do i = 1, mx
            test_mesh%ulat(i, j) = real(j - 1, 8) * 2.0
            test_mesh%ulon(i, j) = real(i - 1, 8) * 2.0
        end do
    end do

    ! 2. Execute local transformation derivation pipelines
    call test_mesh%local_metric()

    ! 3. Assert interior matrices maintain accurate spatial delta parameters
    if (test_mesh%dlatdi(2, 2) /= 0.0_real64) then
        print *, "[FAIL] Latitude shouldn't shift across index modifications in i-direction."
        stop 1
    end if
    if (test_mesh%dlatdj(2, 2) <= 0.0_real64 .or. test_mesh%area(2, 2) <= 0.0_real64) then
        print *, "[FAIL] Matrix area or horizontal transformation steps returned zero values."
        stop 1
    end if

    print *, "[PASS] Spatial metric calculations match structure invariants cleanly."
end program test_metric_matrices