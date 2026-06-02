program test_drifter_kinematics
    use drifter_mod
    use metric_mod
    implicit none

    type(metric) :: dynamic_mesh
    type(drifter) :: test_buoy
    real(kind=real64), allocatable :: u(:,:), v(:,:)
    real(kind=real64) :: simulation_dt
    integer :: i, j

    print *, "=== RUNNING REGRESSION TEST: drifter.f90 ==="

    ! 1. Construct standard metrics tracking matrix
    call dynamic_mesh%set(10, 10)
    do j = 1, 10
        do i = 1, 10
            dynamic_mesh%ulat(i, j) = real(j - 1, 8) * 0.1_real64
            dynamic_mesh%ulon(i, j) = real(i - 1, 8) * 0.1_real64
        end do
    end do
    call dynamic_mesh%local_metric()
    !debug: PRINT *,'dx = ',dynamic_mesh%dx

    
    allocate(u(10, 10), v(10, 10))
    u = 10.0_real64                  ! Constant eastward drift force (10 m/s)
    v = 0.0_real64

    ! 2. Initialize the buoy object at specific grid indexing starting points
    test_buoy%x = 5.0_real64
    test_buoy%y = 5.0_real64
    test_buoy%clat = 0.4_real64
    test_buoy%clon = 0.4_real64
    test_buoy%ilat = 0.4_real64
    test_buoy%ilon = 0.4_real64

    ! 3. Propagate buoy object via kinematics step function
    simulation_dt = 100.0_real64 ! 100 seconds tracking window
    call test_buoy%move(u, v, dynamic_mesh, simulation_dt)

    ! Expected position displacement calculation: 
    ! di delta = (u * dt) / dx = (10 * 100) / 1000 = 1.0 complete grid node step forward
    if (abs((test_buoy%clon - test_buoy%ilon)*dynamic_mesh%dx(5,5) - 1000.0_real64) > 1.0e-2) then
        print *, "[FAIL] Kinetic movement updates missed tracking displacement targets: ", test_buoy%x
        PRINT *,test_buoy%x, test_buoy%y, test_buoy%clat, test_buoy%clon
        PRINT *,test_buoy%clon, test_buoy%ilon, dynamic_mesh%dx(5,5)
        PRINT *,(test_buoy%clon - test_buoy%ilon)*dynamic_mesh%dx(5,5)
        stop 1
    else
        print *, "[PASS] Buoy object updated coordinate displacements cleanly."
    end if

    deallocate(u, v)
end program test_drifter_kinematics
