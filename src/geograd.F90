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

    subroutine compute(KS, intersect_length, mindist, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in)
        use triangles
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in ! known global minimum distance (used for second pass, computing KS)
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        integer :: tri_ind_1, tri_ind_2 ! loop indices
        real(kind=8) :: dsquared, cur_min_dist, intersect_temp
        real(kind=8), dimension(15) :: distance_vec
        real(kind=8), PARAMETER :: second_pass_flag = -1.0

        intersect_length = 0.0_8
        cur_min_dist = 9.9e10
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

                dsquared = minval(distance_vec)
                call compare_and_swap_minimum(cur_min_dist, dsquared)
                
                if (mindist_in /= second_pass_flag) then
                ! compute the intersection on second pass only
                    call intersect(A1(:,tri_ind_1), B1(:,tri_ind_1), C1(:,tri_ind_1), &
                    A2(:,tri_ind_2), B2(:,tri_ind_2), C2(:,tri_ind_2), intersect_temp)
                    intersect_length = intersect_temp + intersect_length
                end if
            end do
        end do
        KS = cur_min_dist
        mindist = cur_min_dist
    end subroutine compute

end module geograd
