program test_kdfind_interface
    use kd_tree_mod
    implicit none

    integer, parameter :: nx = 10, ny = 10
    real :: lats(nx, ny), lons(nx, ny)
    integer, parameter :: bad_count = 2
    real :: bad_lat(bad_count), bad_lon(bad_count)
    real :: bad_fi(bad_count), bad_fj(bad_count)
    integer :: i, j

    print *, "=== RUNNING REGRESSION TEST: kdfind.f90 ==="

    ! 1. Synthesize a mock linear geographic field grid
    do j = 1, ny
        do i = 1, nx
            lats(i, j) = real(j - 1) * 5.0  ! Latitude maps across columns
            lons(i, j) = real(i - 1) * 5.0  ! Longitude maps across rows
        end do
    end do

    ! 2. Set explicit queries directly matching known index combinations
    ! Target 1: Exact center matching grid element i=5, j=5 (lat=20.0, lon=20.0)
    bad_lat(1) = 20.0
    bad_lon(1) = 20.0
    ! Target 2: Near edge elements matching index i=2, j=3 (lat=10.0, lon=5.0)
    bad_lat(2) = 10.1
    bad_lon(2) = 4.9

    call kdfind(nx, ny, lats, lons, bad_count, bad_lat, bad_lon, bad_fi, bad_fj)

    ! 3. Enforce positional matching assertions
    if (nint(bad_fi(1)) == 5 .and. nint(bad_fj(1)) == 5) then
        print *, "[PASS] Target 1 unwrapped spatial grid coordinates match: (5.0, 5.0)"
    else
        print *, "[FAIL] Target 1 grid mapping incorrect: ", bad_fi(1), bad_fj(1)
        stop 1
    end if

    if (nint(bad_fi(2)) == 2 .and. nint(bad_fj(2)) == 3) then
        print *, "[PASS] Target 2 proximity spatial grid coordinates match: (2.0, 3.0)"
    else
        print *, "[FAIL] Target 2 grid mapping incorrect: ", bad_fi(2), bad_fj(2)
        stop 1
    end if

end program test_kdfind_interface