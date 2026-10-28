program test_kd_mod
    use kd_tree_mod
    implicit none

    integer, parameter :: N = 5
    type(GlobalGridNode) :: nodes(N)
    integer :: indices(N)
    integer :: root_idx, best_idx, i
    real :: target_xyz(3)

    print *, "=== RUNNING REGRESSION TEST: kd.f90 ==="

    ! 1. Initialize manual nodes across different coordinate quadrants
    nodes(1)%lat_lon = [0.0, 0.0]     ! Equator, Prime Meridian
    nodes(2)%lat_lon = [45.0, 45.0]   ! Northern / Eastern Hemisphere
    nodes(3)%lat_lon = [-45.0, -45.0] ! Southern / Western Hemisphere
    nodes(4)%lat_lon = [89.0, 0.0]    ! Near North Pole
    nodes(5)%lat_lon = [0.0, 180.0]   ! Equator, Date Line

    do i = 1, N
        call lat_lon_to_3d(nodes(i)%lat_lon(1), nodes(i)%lat_lon(2), nodes(i)%coord)
        nodes(i)%payload%grid_id = i
        nodes(i)%payload%value = real(i)
        indices(i) = i
    end do

    ! 2. Verify structural tree composition
    root_idx = build_global_kd_tree(nodes, indices, 1, N, 0)
    if (root_idx <= 0 .or. root_idx > N) then
        print *, "[FAIL] Invalid tree root index returned: ", root_idx
        stop 1
    end if
    print *, "[PASS] Tree structural build completed. Root node index: ", root_idx

    ! 3. Verify exact Match Nearest Neighbor
    call lat_lon_to_3d(44.0, 46.0, target_xyz) ! Close to nodes(2)
    best_idx = 0
    call find_nearest_global(nodes, root_idx, target_xyz, 0, best_idx)
    
    if (best_idx == 2) then
        print *, "[PASS] Nearest neighbor test match successful. Found index: ", best_idx
    else
        print *, "[FAIL] Mismatched neighbor lookup. Expected index 2, got: ", best_idx
        stop 1
    end if

end program test_kd_mod