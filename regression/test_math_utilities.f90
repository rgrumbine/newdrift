program test_math_utilities
    use constants
    implicit none

    ! External function signatures inside math.f90
    real(kind=real64) :: wrap, harcdis
    external :: wrap, harcdis, bearing, unbearing
    
    real(kind=real64) :: t_lat1, t_lon1, t_lat2, t_lon2
    real(kind=real64) :: distance, heading, out_lat, out_lon

    print *, "=== RUNNING REGRESSION TEST: math.f90 ==="

    ! 1. Assertive validation of Longitude Boundary Wrapping
    if (abs(wrap(370.0_real64) - 10.0_real64) > 1.0e-5) then
        print *, "[FAIL] Positive out-of-bounds wrap failure: ", wrap(370.0_real64)
        stop 1
    end if
    if (abs(wrap(-10.0_real64) - 350.0_real64) > 1.0e-5) then
        print *, "[FAIL] Negative out-of-bounds wrap failure: ", wrap(-10.0_real64)
        stop 1
    end if
    print *, "[PASS] Longitude wrap handling boundary rules correct."

    ! 2. Verification of Distance & Coordinate Projection Loops
    t_lat1 = 0.0_real64;   t_lon1 = 0.0_real64
    t_lat2 = 0.0_real64;   t_lon2 = 1.0_real64 ! Moving 1 degree east along equator

    distance = harcdis(t_lat1, t_lon1, t_lat2, t_lon2)
    ! 1 degree along the equator should roughly equal ~111.3 km
    if (distance < 110.0 .or. distance > 112.0) then
        print *, "[FAIL] Haversine distance out of realistic range: ", distance
        stop 1
    end if
    print *, "[PASS] Distance validation check output: ", distance, " km"

    ! 3. Assert forward and reverse projection parity matches closure
    call bearing(t_lat1, t_lon1, t_lat2, t_lon2, distance, heading)
    call unbearing(t_lat1, t_lon1, distance, heading, out_lat, out_lon)

    if (abs(out_lat - t_lat2) > 1.0e-4 .or. abs(out_lon - t_lon2) > 1.0e-4) then
        print *, "[FAIL] Closure roundtrip projection parity error. Expected: ", t_lat2, t_lon2, " Got: ", out_lat, out_lon
        stop 1
    else
        print *, "[PASS] Bearing/Unbearing closure checks run cleanly."
    end if

end program test_math_utilities