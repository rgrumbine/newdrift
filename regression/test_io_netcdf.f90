program test_io_netcdf_layer
    use io
    use metric_mod
    implicit none

    character(len=*), parameter :: test_filename = "nonexistent_verification_file.nc"
    character(len=40) :: mock_vars(1)
    integer :: test_ncid, out_nx, out_ny
    integer :: mock_varids(1)
    
    print *, "=== RUNNING REGRESSION TEST: io_nc.f90 ==="
    print *, "Testing runtime exception handling of missing tracking assets..."

    mock_vars(1) = "ssh"
    
    ! This execution path should gracefully throw error signals or safely abort execution,
    ! proving the NetCDF validation hooks inside the subroutine check functions work.
    call initialize_in(1, test_filename, test_ncid, mock_varids, mock_vars, "lon", "lat", out_nx, out_ny, 1)

    print *, "[FAIL] Code execution failed to identify missing physical data dependencies."
    stop 1

end program test_io_netcdf_layer