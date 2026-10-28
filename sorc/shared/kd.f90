! Code from Gemini to implement a k-d tree on the surface of a sphere
! 21 May 2026

module kd_tree_mod
    implicit none
    private
    public :: GlobalGridNode, build_global_kd_tree, find_nearest_global, lat_lon_to_3d

    ! A structure to hold the payload data for your nonuniform grid
    type :: GridData
        integer :: grid_id
        real :: value        ! e.g., temperature, elevation, pressure
    end type GridData

    type :: GlobalGridNode
        real :: coord(3)     ! [X, Y, Z] on unit sphere
        real :: lat_lon(2)   ! [Lat, Lon] in degrees
        type(GridData) :: payload
        integer :: left = 0  ! Index of left child (0 means null)
        integer :: right = 0 ! Index of right child (0 means null)
    end type GlobalGridNode

contains

    subroutine lat_lon_to_3d(lat, lon, xyz)
        real, intent(in) :: lat, lon
        real, intent(out) :: xyz(3)
        real :: lat_rad, lon_rad
        real, parameter :: PI = 3.141592653589793

        lat_rad = lat * PI / 180.0
        lon_rad = lon * PI / 180.0

        xyz(1) = cos(lat_rad) * cos(lon_rad)
        xyz(2) = cos(lat_rad) * sin(lon_rad)
        xyz(3) = sin(lat_rad)
    end subroutine lat_lon_to_3d

    real function distance_3d(p1, p2)
        real, intent(in) :: p1(3), p2(3)
        distance_3d = sqrt((p1(1)-p2(1))**2 + (p1(2)-p2(2))**2 + (p1(3)-p2(3))**2)
    end function distance_3d


recursive function build_global_kd_tree(nodes, indices, start_idx, end_idx, depth) result(root_idx)
        type(GlobalGridNode), intent(inout) :: nodes(:)
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: start_idx, end_idx, depth
        integer :: root_idx
        integer :: axis, median

        if (start_idx > end_idx) then
            root_idx = 0
            return
        end if

        ! Fortran arrays are 1-indexed. Cycle axis between 1, 2, and 3
        axis = modulo(depth, 3) + 1

        ! Sort the indices array slice based on the current axis
        call sort_indices_by_axis(nodes, indices, start_idx, end_idx, axis)
        
        median = start_idx + (end_idx - start_idx) / 2
        root_idx = indices(median)

        ! Recursively build left and right subtrees
        nodes(root_idx)%left = build_global_kd_tree(nodes, indices, start_idx, median - 1, depth + 1)
        nodes(root_idx)%right = build_global_kd_tree(nodes, indices, median + 1, end_idx, depth + 1)
    end function build_global_kd_tree

    ! Simple QuickSort helper for sorting the index pointer array
    recursive subroutine sort_indices_by_axis(nodes, indices, left, right, axis)
        type(GlobalGridNode), intent(in) :: nodes(:)
        integer, intent(inout) :: indices(:)
        integer, intent(in) :: left, right, axis
        integer :: i, j, pivot_idx, temp
        real :: pivot_val

        if (left >= right) return
        pivot_idx = indices(left + (right - left) / 2)
        pivot_val = nodes(pivot_idx)%coord(axis)
        i = left
        j = right

        do while (i <= j)
            do while (nodes(indices(i))%coord(axis) < pivot_val)
                i = i + 1
            end do
            do while (nodes(indices(j))%coord(axis) > pivot_val)
                j = j - 1
            end do
            if (i <= j) then
                temp = indices(i)
                indices(i) = indices(j)
                indices(j) = temp
                i = i + 1
                j = j - 1
            end if
        end do
        if (left < j) call sort_indices_by_axis(nodes, indices, left, j, axis)
        if (i < right) call sort_indices_by_axis(nodes, indices, i, right, axis)
    end subroutine sort_indices_by_axis


recursive subroutine find_nearest_global(nodes, curr_idx, target_xyz, depth, best_idx)
        type(GlobalGridNode), intent(in) :: nodes(:)
        integer, intent(in) :: curr_idx
        real, intent(in) :: target_xyz(3)
        integer, intent(in) :: depth
        integer, intent(inout) :: best_idx
        
        integer :: axis, next_branch, other_branch
        real :: current_dist, best_dist

        if (curr_idx == 0) return

        ! Calculate current minimum distance
        if (best_idx == 0) then
            best_idx = curr_idx
        else
            if (distance_3d(target_xyz, nodes(curr_idx)%coord) < distance_3d(target_xyz, nodes(best_idx)%coord)) then
                best_idx = curr_idx
            end if
        end if

        axis = modulo(depth, 3) + 1

        ! Determine branch ordering
        if (target_xyz(axis) < nodes(curr_idx)%coord(axis)) then
            next_branch = nodes(curr_idx)%left
            other_branch = nodes(curr_idx)%right
        else
            next_branch = nodes(curr_idx)%right
            other_branch = nodes(curr_idx)%left
        end if

        ! Search down the promising side
        call find_nearest_global(nodes, next_branch, target_xyz, depth + 1, best_idx)

        ! Check if we must cross the boundary plane
        best_dist = distance_3d(target_xyz, nodes(best_idx)%coord)
        if (abs(target_xyz(axis) - nodes(curr_idx)%coord(axis)) < best_dist) then
            call find_nearest_global(nodes, other_branch, target_xyz, depth + 1, best_idx)
        end if
    end subroutine find_nearest_global

end module kd_tree_mod
