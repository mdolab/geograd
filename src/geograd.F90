module geograd
    implicit none

    contains
    subroutine compare_and_swap_minimum(current_min, value)
        implicit none
        real(kind=8), intent(inout) :: current_min
        real(kind=8), intent(in) :: value
        if (value < current_min) then
            current_min = value
        end if
    end subroutine compare_and_swap_minimum

    subroutine minval_and_loc(A, length_A, min_value, min_index)
        implicit none
        real(kind=8), dimension(length_A), intent(in) :: A
        integer, intent(in) :: length_A
        real(kind=8), intent(out) :: min_value
        integer, intent(out) :: min_index
        integer :: count
        min_value = A(1)
        min_index = 1
        do count = 2, length_A
            if (A(count) < min_value) then
                min_value = A(count)
                min_index = count
            end if
        end do
        RETURN
    end subroutine minval_and_loc

    subroutine compute(KS, intersect_length, mindist, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho)
        use triangles
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho ! known global minimum distance (used for second pass, computing KS)
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        integer :: tri_ind_1, tri_ind_2, minloc_index ! loop indices
        real(kind=8) :: d, cur_min_dist, intersect_temp, ks_accumulator
        real(kind=8), dimension(15) :: distance_vec
        real(kind=8), PARAMETER :: second_pass_flag = -1.0

        intersect_length = 0.0_8
        cur_min_dist = 9.9e10
        ks_accumulator = 0.0_8

        do tri_ind_1 = 1, n1
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests
                call line_line(A1(:,tri_ind_1), B1(:,tri_ind_1), A2(:,tri_ind_2), B2(:,tri_ind_2), distance_vec(1))
                call line_line(A1(:,tri_ind_1), B1(:,tri_ind_1), B2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(2))
                call line_line(A1(:,tri_ind_1), B1(:,tri_ind_1), A2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(3))
                call line_line(B1(:,tri_ind_1), C1(:,tri_ind_1), A2(:,tri_ind_2), B2(:,tri_ind_2), distance_vec(4))
                call line_line(B1(:,tri_ind_1), C1(:,tri_ind_1), B2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(5))
                call line_line(B1(:,tri_ind_1), C1(:,tri_ind_1), A2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(6))
                call line_line(A1(:,tri_ind_1), C1(:,tri_ind_1), A2(:,tri_ind_2), B2(:,tri_ind_2), distance_vec(7))
                call line_line(A1(:,tri_ind_1), C1(:,tri_ind_1), B2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(8))
                call line_line(A1(:,tri_ind_1), C1(:,tri_ind_1), A2(:,tri_ind_2), C2(:,tri_ind_2), distance_vec(9))

                ! do 6 point-triangle comparison tests
                call point_tri(A1(:,tri_ind_1), B1(:,tri_ind_1), C1(:,tri_ind_1), A2(:,tri_ind_2), distance_vec(10))
                call point_tri(A1(:,tri_ind_1), B1(:,tri_ind_1), C1(:,tri_ind_1), B2(:,tri_ind_2), distance_vec(11))
                call point_tri(A1(:,tri_ind_1), B1(:,tri_ind_1), C1(:,tri_ind_1), C2(:,tri_ind_2), distance_vec(12))
                call point_tri(A2(:,tri_ind_2), B2(:,tri_ind_2), C2(:,tri_ind_2), A1(:,tri_ind_1), distance_vec(13))
                call point_tri(A2(:,tri_ind_2), B2(:,tri_ind_2), C2(:,tri_ind_2), B1(:,tri_ind_1), distance_vec(14))
                call point_tri(A2(:,tri_ind_2), B2(:,tri_ind_2), C2(:,tri_ind_2), C1(:,tri_ind_1), distance_vec(15))

                call minval_and_loc(distance_vec, 15, d, minloc_index)
                call compare_and_swap_minimum(cur_min_dist, d)
                
                if (mindist_in /= second_pass_flag) then
                ! compute KS function on second pass only
                    ks_accumulator = ks_accumulator + sum(exp((mindist_in - distance_vec)*rho))
                ! compute the intersection on second pass only
                    call intersect(A1(:,tri_ind_1), B1(:,tri_ind_1), C1(:,tri_ind_1), &
                    A2(:,tri_ind_2), B2(:,tri_ind_2), C2(:,tri_ind_2), intersect_temp)
                    intersect_length = intersect_temp + intersect_length
                end if
            end do
        end do
        if (mindist_in /= second_pass_flag) then
            KS = (1 / rho) * log(ks_accumulator) - mindist_in
        else
            KS = zero
        end if
        mindist = cur_min_dist
    end subroutine compute

end module geograd
