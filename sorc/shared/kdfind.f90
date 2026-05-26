
! Code base from Gemini implementing a k-d tree and searching it in fortran
! 21 May 2026

SUBROUTINE kdfind( nx, ny, lats, lons, bad_count, bad_lat, bad_lon, bad_fi, bad_fj)
    use kd_tree_mod
    implicit none
    INTEGER, intent(in) :: nx, ny, bad_count
    REAL, intent(in) :: lons(nx, ny), lats(nx, ny)
    REAL, intent(inout) :: bad_lat(bad_count), bad_lon(bad_count)
    REAL, intent(inout) :: bad_fi(bad_count), bad_fj(bad_count)

    ! Note: building and throwing away the kdtree
    !Gemini
    integer :: N 
    type(GlobalGridNode), allocatable :: tree_nodes(:)
    integer, allocatable :: indices(:)
    integer :: root_node, best_node, i, j, k
    real :: query_lat, query_lon, query_xyz(3)
    real :: temp_xyz(3)

    N = nx*ny
    allocate(tree_nodes(N))
    allocate(indices(N))

    ! RG
    k = 0
    DO j = 1, ny
    DO i = 1, nx
      k = k + 1
      tree_nodes(k)%lat_lon(1) = lats(i,j)
      tree_nodes(k)%lat_lon(2) = lons(i,j)
      indices(k) = k

      ! Pre-calculate 3D spatial points
      call lat_lon_to_3d(tree_nodes(k)%lat_lon(1), tree_nodes(k)%lat_lon(2), tree_nodes(k)%coord)
     
      tree_nodes(k)%payload%grid_id = k
      tree_nodes(k)%payload%value =  k ! Sample payload data
    ENDDO
    ENDDO
    !END RG


    ! 2. Build the Tree
    !debug: print *, "Building K-D Tree "
    root_node = build_global_kd_tree(tree_nodes, indices, 1, N, 0)
    print *, "Tree structural build finished! Root node is index: ", root_node

    ! 3. Search Example
    do i = 1, bad_count
      query_lat = bad_lat(i)
      query_lon = bad_lon(i)
      call lat_lon_to_3d(query_lat, query_lon, query_xyz)
  
      best_node = 0
      !debug2: print *, "Querying closest grid point to Lat, lon",query_lat, query_lon
      call find_nearest_global(tree_nodes, root_node, query_xyz, 0, best_node)
  
      !print *, "--- NEAREST NEIGHBOR FOUND ---"
      !print *, "Node Array Index: ", best_node
      !print *, "Grid ID: ", tree_nodes(best_node)%payload%grid_id
      !print *, "Actual Lat/Lon: ", tree_nodes(best_node)%lat_lon
      !print *, "Data Value: ", tree_nodes(best_node)%payload%value
      bad_fi(i) =  1. + MOD(best_node-1, nx)
      bad_fj(i) =  1. + INT((best_node-1) / nx)
      !debug WRITE(*,9001) i, query_lat, query_lon, bad_fi(i), bad_fj(i), &
      !debug         tree_nodes(best_node)%lat_lon(1), tree_nodes(best_node)%lat_lon(2), &
      !debug         best_node
 9001 FORMAT(I5,6F9.3,I8)

    enddo !i

    deallocate(tree_nodes)
    deallocate(indices)
end subroutine kdfind
