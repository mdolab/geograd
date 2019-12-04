module geograd
    implicit none

    contains
    subroutine compare_and_swap_minimum(current_min, value, flag)
        implicit none
        real(kind=8), intent(inout) :: current_min
        real(kind=8), intent(in) :: value
        integer, intent(out), optional :: flag
        if(present(flag)) then
            flag = -1
        end if
        if (value < current_min) then
            current_min = value
            if(present(flag)) then
                flag = 1
            end if
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

#ifndef USE_COMPLEX
    subroutine compute_derivs(KS, intersect_length, mindist, dKSdA1, dKSdB1, dKSdC1, dKSdA2, dKSdB2, dKSdC2, &
         A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho)
        use triangles_db
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho ! known global minimum distance (used for second pass, computing KS
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        real(kind=8), intent(out), dimension(3,n1) :: dKSdA1, dKSdB1, dKSdC1
        real(kind=8), intent(out), dimension(3,n2) :: dKSdA2, dKSdB2, dKSdC2
        integer :: tri_ind_1, tri_ind_2, minloc_index, count, is_less_flag ! loop indices
        real(kind=8) :: d, cur_min_dist, garbage, base_exp_accumulator, deriv_exp_accumulator
        real(kind=8), dimension(3) :: sumdA1, sumdB1, sumdC1, sumdA2, sumdB2, sumdC2
        real(kind=8), dimension(15) :: distance_vec ! holds the pairwise distances from the batch
        real(kind=8), dimension(3,15) :: dA1, dB1, dC1, dA2, dB2, dC2, deriv_exp_temp3 ! holds the reverse mode derivatives for the batch
        real(kind=8), dimension(15) :: dist_subtract_temp, base_exp_temp, deriv_exp_temp ! these hold some intermediate quantities that need to get accumulated
        real(kind=8) :: rev_seed
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch
        real(kind=8), dimension(3) :: ddminA1, ddminB1, ddminC1, ddminA2, ddminB2, ddminC2
        integer :: dmin_index_1, dmin_index_2
        dKSdA1 = 0.0
        dKSdB1 = 0.0
        dKSdC1 = 0.0
        dKSdA2 = 0.0
        dKSdB2 = 0.0
        dKSdC2 = 0.0
        base_exp_accumulator = 0.0
        deriv_exp_accumulator = 0.0
        cur_min_dist = 9.9e10


        do tri_ind_1 = 1, n1
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests and derivatives
                A1batch = A1(:,tri_ind_1)
                B1batch = B1(:,tri_ind_1)
                C1batch = C1(:,tri_ind_1)
                A2batch = A2(:,tri_ind_2)
                B2batch = B2(:,tri_ind_2)
                C2batch = C2(:,tri_ind_2) 
                dA1 = 0.0
                dB1 = 0.0
                dC1 = 0.0
                dA2 = 0.0
                dB2 = 0.0
                dC2 = 0.0
                rev_seed = 1.0
                call line_line(A1batch, B1batch, A2batch, B2batch, distance_vec(1))
                call line_line_b(A1batch, dA1(:,1), B1batch, dB1(:,1), &
                                 A2batch, dA2(:,1), B2batch, dB2(:,1), garbage, rev_seed)
                rev_seed = 1.0
                call line_line(A1batch, B1batch, B2batch, C2batch, distance_vec(2))
                call line_line_b(A1batch, dA1(:,2), B1batch, dB1(:,2), &
                                 B2batch, dB2(:,2), C2batch, dC2(:,2), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(A1batch, B1batch, A2batch, C2batch, distance_vec(3))
                call line_line_b(A1batch, dA1(:,3), B1batch, dB1(:,3), &
                                 A2batch, dA2(:,3), C2batch, dC2(:,3), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(B1batch, C1batch, A2batch, B2batch, distance_vec(4))
                call line_line_b(B1batch, dB1(:,4), C1batch, dC1(:,4), &
                                 A2batch, dA2(:,4), B2batch, dB2(:,4), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(B1batch, C1batch, B2batch, C2batch, distance_vec(5))
                call line_line_b(B1batch, dB1(:,5), C1batch, dC1(:,5), &
                                 B2batch, dB2(:,5), C2batch, dC2(:,5), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(B1batch, C1batch, A2batch, C2batch, distance_vec(6))
                call line_line_b(B1batch, dB1(:,6), C1batch, dC1(:,6), &
                                 A2batch, dA2(:,6), C2batch, dC2(:,6), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(A1batch, C1batch, A2batch, B2batch, distance_vec(7))
                call line_line_b(A1batch, dA1(:,7), C1batch, dC1(:,7), &
                                 A2batch, dA2(:,7), B2batch, dB2(:,7), garbage, rev_seed)


                rev_seed = 1.0
                call line_line(A1batch, C1batch, B2batch, C2batch, distance_vec(8))
                call line_line_b(A1batch, dA1(:,8), C1batch, dC1(:,8), &
                                 B2batch, dB2(:,8), C2batch, dC2(:,8), garbage, rev_seed)

                rev_seed = 1.0
                call line_line(A1batch, C1batch, A2batch, C2batch, distance_vec(9))
                call line_line_b(A1batch, dA1(:,9), C1batch, dC1(:,9), &
                                 A2batch, dA2(:,9), C2batch, dC2(:,9), garbage, rev_seed)

                ! do 6 point-triangle comparison tests
                rev_seed = 1.0
                call point_tri(A1batch, B1batch, C1batch, A2batch, distance_vec(10))
                call point_tri_b(A1batch, dA1(:,10), B1batch, dB1(:,10), C1batch, dC1(:,10), &
                                 A2batch, dA2(:,10), garbage, rev_seed)    

                rev_seed = 1.0
                call point_tri(A1batch, B1batch, C1batch, B2batch, distance_vec(11))
                call point_tri_b(A1batch, dA1(:,11), B1batch, dB1(:,11), C1batch, dC1(:,11), &
                                 B2batch, dB2(:,11), garbage, rev_seed)

                rev_seed = 1.0
                call point_tri(A1batch, B1batch, C1batch, C2batch, distance_vec(12))
                call point_tri_b(A1batch, dA1(:,12), B1batch, dB1(:,12), C1batch, dC1(:,12), &
                                 C2batch, dC2(:,12), garbage, rev_seed)

                rev_seed = 1.0
                call point_tri(A2batch, B2batch, C2batch, A1batch, distance_vec(13))
                call point_tri_b(A2batch, dA2(:,13), B2batch, dB2(:,13), C2batch, dC2(:,13), &
                                 A1batch, dA1(:,13), garbage, rev_seed)

                rev_seed = 1.0
                call point_tri(A2batch, B2batch, C2batch, B1batch, distance_vec(14))
                call point_tri_b(A2batch, dA2(:,14), B2batch, dB2(:,14), C2batch, dC2(:,14), &
                                 B1batch, dB1(:,14), garbage, rev_seed)

                rev_seed = 1.0
                call point_tri(A2batch, B2batch, C2batch, C1batch, distance_vec(15))
                call point_tri_b(A2batch, dA2(:,15), B2batch, dB2(:,15), C2batch, dC2(:,15), &
                                 C1batch, dC1(:,15), garbage, rev_seed)

                call minval_and_loc(distance_vec, 15, d, minloc_index)
                call compare_and_swap_minimum(cur_min_dist, d, is_less_flag)
                if (is_less_flag == 1) then
                    ! If the new minimum, set the global dmin derivs and both indices
                    dmin_index_1 = tri_ind_1
                    dmin_index_2 = tri_ind_2
                    ddminA1 = dA1(:,minloc_index)
                    ddminB1 = dB1(:,minloc_index)
                    ddminC1 = dC1(:,minloc_index)
                    ddminA2 = dA2(:,minloc_index)
                    ddminB2 = dB2(:,minloc_index)
                    ddminC2 = dC2(:,minloc_index)
                end if

                dist_subtract_temp = mindist_in - distance_vec
                base_exp_temp = exp(rho*dist_subtract_temp)
                !deriv_exp_temp = dist_subtract_temp * base_exp_temp
                ! accumulate exponentials
                base_exp_accumulator = base_exp_accumulator + sum(base_exp_temp)
                !deriv_exp_accumulator = deriv_exp_accumulator + sum(deriv_exp_temp)

                ! deriv_exp_temp3 = spread(deriv_exp_temp, 1, 3)
                
                sumdA1 = 0.0
                sumdB1 = 0.0
                sumdC1 = 0.0
                sumdA2 = 0.0
                sumdB2 = 0.0
                sumdC2 = 0.0
                do count = 1,15
                    sumdA1 = sumdA1 - dA1(:,count)*base_exp_temp(count)
                    sumdB1 = sumdB1 - dB1(:,count)*base_exp_temp(count)
                    sumdC1 = sumdC1 - dC1(:,count)*base_exp_temp(count)
                    sumdA2 = sumdA2 - dA2(:,count)*base_exp_temp(count)
                    sumdB2 = sumdB2 - dB2(:,count)*base_exp_temp(count)
                    sumdC2 = sumdC2 - dC2(:,count)*base_exp_temp(count)
                    
                end do
                ! dA1 = dA1 * deriv_exp_temp3
                ! dB1 = dB1 * deriv_exp_temp3
                ! dC1 = dC1 * deriv_exp_temp3
                ! dA2 = dA2 * deriv_exp_temp3
                ! dB2 = dB2 * deriv_exp_temp3
                ! dC2 = dC2 * deriv_exp_temp3
                dKSdA1(:, tri_ind_1) = dKSdA1(:, tri_ind_1) + sumdA1
                dKSdB1(:, tri_ind_1) = dKSdB1(:, tri_ind_1) + sumdB1
                dKSdC1(:, tri_ind_1) = dKSdC1(:, tri_ind_1) + sumdC1
                dKSdA2(:, tri_ind_2) = dKSdA2(:, tri_ind_2) + sumdA2
                dKSdB2(:, tri_ind_2) = dKSdB2(:, tri_ind_2) + sumdB2
                dKSdC2(:, tri_ind_2) = dKSdC2(:, tri_ind_2) + sumdC2
            end do
        end do
        ! include the dmin contribution to the accumulated derivatives
        ! dKSdA1(:, dmin_index_1) = dKSdA1(:, dmin_index_1) + ddminA1 * base_exp_accumulator
        ! dKSdB1(:, dmin_index_1) = dKSdB1(:, dmin_index_1) + ddminB1 * base_exp_accumulator
        ! dKSdC1(:, dmin_index_1) = dKSdC1(:, dmin_index_1) + ddminC1 * base_exp_accumulator
        ! dKSdA2(:, dmin_index_2) = dKSdA2(:, dmin_index_2) + ddminA2 * base_exp_accumulator
        ! dKSdB2(:, dmin_index_2) = dKSdB2(:, dmin_index_2) + ddminB2 * base_exp_accumulator
        ! dKSdC2(:, dmin_index_2) = dKSdC2(:, dmin_index_2) + ddminC2 * base_exp_accumulator
        dKSdA1 = (1 / base_exp_accumulator) * dKSdA1
        dKSdB1 = (1 / base_exp_accumulator) * dKSdB1
        dKSdC1 = (1 / base_exp_accumulator) * dKSdC1
        dKSdA2 = (1 / base_exp_accumulator) * dKSdA2
        dKSdB2 = (1 / base_exp_accumulator) * dKSdB2
        dKSdC2 = (1 / base_exp_accumulator) * dKSdC2
        ! dKSdA1(:, dmin_index_1) = dKSdA1(:, dmin_index_1) - ddminA1
        ! dKSdB1(:, dmin_index_1) = dKSdB1(:, dmin_index_1) - ddminB1
        ! dKSdC1(:, dmin_index_1) = dKSdC1(:, dmin_index_1) - ddminC1
        ! dKSdA2(:, dmin_index_2) = dKSdA2(:, dmin_index_2) - ddminA2
        ! dKSdB2(:, dmin_index_2) = dKSdB2(:, dmin_index_2) - ddminB2
        ! dKSdC2(:, dmin_index_2) = dKSdC2(:, dmin_index_2) - ddminC2

        ! print *, dmin_index_1, dmin_index_2

        KS = (1/rho) * log(base_exp_accumulator) - mindist_in
        mindist = cur_min_dist
        if (cur_min_dist /= mindist) then
            print *, 'These should match'
        end if
        intersect_length = 0
        return

    end subroutine compute_derivs
#endif

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
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch

        intersect_length = 0.0_8
        cur_min_dist = 9.9e10
        ks_accumulator = 0.0_8

        do tri_ind_1 = 1, n1
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests
                A1batch = A1(:,tri_ind_1)
                B1batch = B1(:,tri_ind_1)
                C1batch = C1(:,tri_ind_1)
                A2batch = A2(:,tri_ind_2)
                B2batch = B2(:,tri_ind_2)
                C2batch = C2(:,tri_ind_2)             
                call line_line(A1batch, B1batch, A2batch, B2batch, distance_vec(1))
                call line_line(A1batch, B1batch, B2batch, C2batch, distance_vec(2))
                call line_line(A1batch, B1batch, A2batch, C2batch, distance_vec(3))
                call line_line(B1batch, C1batch, A2batch, B2batch, distance_vec(4))
                call line_line(B1batch, C1batch, B2batch, C2batch, distance_vec(5))
                call line_line(B1batch, C1batch, A2batch, C2batch, distance_vec(6))
                call line_line(A1batch, C1batch, A2batch, B2batch, distance_vec(7))
                call line_line(A1batch, C1batch, B2batch, C2batch, distance_vec(8))
                call line_line(A1batch, C1batch, A2batch, C2batch, distance_vec(9))

                ! do 6 point-triangle comparison tests
                call point_tri(A1batch, B1batch, C1batch, A2batch, distance_vec(10))
                call point_tri(A1batch, B1batch, C1batch, B2batch, distance_vec(11))
                call point_tri(A1batch, B1batch, C1batch, C2batch, distance_vec(12))
                call point_tri(A2batch, B2batch, C2batch, A1batch, distance_vec(13))
                call point_tri(A2batch, B2batch, C2batch, B1batch, distance_vec(14))
                call point_tri(A2batch, B2batch, C2batch, C1batch, distance_vec(15))

                call minval_and_loc(distance_vec, 15, d, minloc_index)
                call compare_and_swap_minimum(cur_min_dist, d)
                
                if (mindist_in /= second_pass_flag) then
                ! compute KS function on second pass only
                    ks_accumulator = ks_accumulator + sum(exp((mindist_in - distance_vec)*rho))
                ! compute the intersection on second pass only
                    call intersect(A1batch, B1batch, C1batch, &
                    A2batch, B2batch, C2batch, intersect_temp)
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
