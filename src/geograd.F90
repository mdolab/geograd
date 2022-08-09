module geograd
    implicit none

contains
    subroutine compare_and_swap_minimum(current_min, value, flag)
        implicit none
        real(kind=8), intent(inout) :: current_min
        real(kind=8), intent(in) :: value
        integer, intent(out), optional :: flag
        if (present(flag)) then
            flag = -1
        end if
        if (value < current_min) then
            current_min = value
            if (present(flag)) then
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
                              dPdA1, dPdB1, dPdC1, dPdA2, dPdB2, dPdC2, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho)
        use triangles_db
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3, n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3, n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho ! known global minimum distance (used for second pass, computing KS
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        real(kind=8), intent(out), dimension(3, n1) :: dKSdA1, dKSdB1, dKSdC1
        real(kind=8), intent(out), dimension(3, n2) :: dKSdA2, dKSdB2, dKSdC2
        real(kind=8), intent(out), dimension(3, n1) :: dPdA1, dPdB1, dPdC1
        real(kind=8), intent(out), dimension(3, n2) :: dPdA2, dPdB2, dPdC2
        integer :: tri_ind_1, tri_ind_2, minloc_index, count, is_less_flag ! loop indices
        real(kind=8) :: d, cur_min_dist, garbage, base_exp_accumulator, intersect_accumulator, intersect_batch
        real(kind=8), dimension(3) :: sumdA1, sumdB1, sumdC1, sumdA2, sumdB2, sumdC2
        real(kind=8), dimension(15) :: distance_vec ! holds the pairwise distances from the batch
        real(kind=8), dimension(3) :: dKSdA1batch, dKSdB1batch, dKSdC1batch, dKSdA2batch, dKSdB2batch, dKSdC2batch ! holds the reverse mode derivatives for the batch
        real(kind=8), dimension(3) :: dPdA1batch, dPdB1batch, dPdC1batch, dPdA2batch, dPdB2batch, dPdC2batch ! holds the reverse mode derivatives for the batch
        real(kind=8) :: base_exp_temp ! these hold some intermediate quantities that need to get accumulated
        real(kind=8) :: rev_seed
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch
        dKSdA1 = 0.0
        dKSdB1 = 0.0
        dKSdC1 = 0.0
        dKSdA2 = 0.0
        dKSdB2 = 0.0
        dKSdC2 = 0.0
        dPdA1 = 0.0
        dPdB1 = 0.0
        dPdC1 = 0.0
        dPdA2 = 0.0
        dPdB2 = 0.0
        dPdC2 = 0.0
        base_exp_accumulator = 0.0
        intersect_accumulator = 0.0
        cur_min_dist = 9.9e10

        do tri_ind_1 = 1, n1
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests and derivatives
                A1batch = A1(:, tri_ind_1)
                B1batch = B1(:, tri_ind_1)
                C1batch = C1(:, tri_ind_1)
                A2batch = A2(:, tri_ind_2)
                B2batch = B2(:, tri_ind_2)
                C2batch = C2(:, tri_ind_2)
                intersect_batch = 0.0
                dKSdA1batch = 0.0
                dKSdB1batch = 0.0
                dKSdC1batch = 0.0
                dKSdA2batch = 0.0
                dKSdB2batch = 0.0
                dKSdC2batch = 0.0
                dPdA1batch = 0.0
                dPdB1batch = 0.0
                dPdC1batch = 0.0
                dPdA2batch = 0.0
                dPdB2batch = 0.0
                dPdC2batch = 0.0
                ! find the minimum distance test of the 15
                call line_line(A1batch, B1batch, A2batch, B2batch, distance_vec(1))
                call line_line(A1batch, B1batch, B2batch, C2batch, distance_vec(2))
                call line_line(A1batch, B1batch, A2batch, C2batch, distance_vec(3))
                call line_line(B1batch, C1batch, A2batch, B2batch, distance_vec(4))
                call line_line(B1batch, C1batch, B2batch, C2batch, distance_vec(5))
                call line_line(B1batch, C1batch, A2batch, C2batch, distance_vec(6))
                call line_line(A1batch, C1batch, A2batch, B2batch, distance_vec(7))
                call line_line(A1batch, C1batch, B2batch, C2batch, distance_vec(8))
                call line_line(A1batch, C1batch, A2batch, C2batch, distance_vec(9))
                call point_tri(A1batch, B1batch, C1batch, A2batch, distance_vec(10))
                call point_tri(A1batch, B1batch, C1batch, B2batch, distance_vec(11))
                call point_tri(A1batch, B1batch, C1batch, C2batch, distance_vec(12))
                call point_tri(A2batch, B2batch, C2batch, A1batch, distance_vec(13))
                call point_tri(A2batch, B2batch, C2batch, B1batch, distance_vec(14))
                call point_tri(A2batch, B2batch, C2batch, C1batch, distance_vec(15))

                call minval_and_loc(distance_vec, 15, d, minloc_index)
                call compare_and_swap_minimum(cur_min_dist, d, is_less_flag)

                rev_seed = 1.0
                select case (minloc_index)
                case (1)
                    call line_line_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, &
                                     A2batch, dKSdA2batch, B2batch, dKSdB2batch, garbage, rev_seed)
                case (2)
                    call line_line_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, &
                                     B2batch, dKSdB2batch, C2batch, dKSdC2batch, garbage, rev_seed)
                case (3)
                    call line_line_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, &
                                     A2batch, dKSdA2batch, C2batch, dKSdC2batch, garbage, rev_seed)
                case (4)
                    call line_line_b(B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     A2batch, dKSdA2batch, B2batch, dKSdB2batch, garbage, rev_seed)
                case (5)
                    call line_line_b(B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     B2batch, dKSdB2batch, C2batch, dKSdC2batch, garbage, rev_seed)
                case (6)
                    call line_line_b(B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     A2batch, dKSdA2batch, C2batch, dKSdC2batch, garbage, rev_seed)
                case (7)
                    call line_line_b(A1batch, dKSdA1batch, C1batch, dKSdC1batch, &
                                     A2batch, dKSdA2batch, B2batch, dKSdB2batch, garbage, rev_seed)
                case (8)
                    call line_line_b(A1batch, dKSdA1batch, C1batch, dKSdC1batch, &
                                     B2batch, dKSdB2batch, C2batch, dKSdC2batch, garbage, rev_seed)
                case (9)
                    call line_line_b(A1batch, dKSdA1batch, C1batch, dKSdC1batch, &
                                     A2batch, dKSdA2batch, C2batch, dKSdC2batch, garbage, rev_seed)

                    ! do 6 point-triangle comparison tests
                case (10)
                    call point_tri_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     A2batch, dKSdA2batch, garbage, rev_seed)
                case (11)
                    call point_tri_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     B2batch, dKSdB2batch, garbage, rev_seed)
                case (12)
                    call point_tri_b(A1batch, dKSdA1batch, B1batch, dKSdB1batch, C1batch, dKSdC1batch, &
                                     C2batch, dKSdC2batch, garbage, rev_seed)
                case (13)
                    call point_tri_b(A2batch, dKSdA2batch, B2batch, dKSdB2batch, C2batch, dKSdC2batch, &
                                     A1batch, dKSdA1batch, garbage, rev_seed)
                case (14)
                    call point_tri_b(A2batch, dKSdA2batch, B2batch, dKSdB2batch, C2batch, dKSdC2batch, &
                                     B1batch, dKSdB1batch, garbage, rev_seed)
                case (15)
                    call point_tri_b(A2batch, dKSdA2batch, B2batch, dKSdB2batch, C2batch, dKSdC2batch, &
                                     C1batch, dKSdC1batch, garbage, rev_seed)
                END select
                rev_seed = 1.0
                call intersect(A1batch, B1batch, C1batch, A2batch, B2batch, C2batch, intersect_batch)
                if (intersect_batch > 0.0) then
                    call intersect_b(A1batch, dPdA1batch, B1batch, dPdB1batch, C1batch, dPdC1batch, &
                                     A2batch, dPdA2batch, B2batch, dPdB2batch, C2batch, dPdC2batch, &
                                     garbage, rev_seed)
                end if
                base_exp_temp = exp((mindist_in - d) * rho)
                ! accumulate exponentials
                base_exp_accumulator = base_exp_accumulator + base_exp_temp
                intersect_accumulator = intersect_accumulator + intersect_batch
                sumdA1 = -dKSdA1batch * base_exp_temp
                sumdB1 = -dKSdB1batch * base_exp_temp
                sumdC1 = -dKSdC1batch * base_exp_temp
                sumdA2 = -dKSdA2batch * base_exp_temp
                sumdB2 = -dKSdB2batch * base_exp_temp
                sumdC2 = -dKSdC2batch * base_exp_temp

                dKSdA1(:, tri_ind_1) = dKSdA1(:, tri_ind_1) + sumdA1
                dKSdB1(:, tri_ind_1) = dKSdB1(:, tri_ind_1) + sumdB1
                dKSdC1(:, tri_ind_1) = dKSdC1(:, tri_ind_1) + sumdC1
                dKSdA2(:, tri_ind_2) = dKSdA2(:, tri_ind_2) + sumdA2
                dKSdB2(:, tri_ind_2) = dKSdB2(:, tri_ind_2) + sumdB2
                dKSdC2(:, tri_ind_2) = dKSdC2(:, tri_ind_2) + sumdC2
                dPdA1(:, tri_ind_1) = dPdA1(:, tri_ind_1) + dPdA1batch
                dPdB1(:, tri_ind_1) = dPdB1(:, tri_ind_1) + dPdB1batch
                dPdC1(:, tri_ind_1) = dPdC1(:, tri_ind_1) + dPdC1batch
                dPdA2(:, tri_ind_2) = dPdA2(:, tri_ind_2) + dPdA2batch
                dPdB2(:, tri_ind_2) = dPdB2(:, tri_ind_2) + dPdB2batch
                dPdC2(:, tri_ind_2) = dPdC2(:, tri_ind_2) + dPdC2batch
            end do
        end do
        ! the derivative with respect to minimum distance drops out
        ! it's just a conditioning step - it doesn't affect the result
        ! except for finite precision arithmetic
        dKSdA1 = (1 / base_exp_accumulator) * dKSdA1
        dKSdB1 = (1 / base_exp_accumulator) * dKSdB1
        dKSdC1 = (1 / base_exp_accumulator) * dKSdC1
        dKSdA2 = (1 / base_exp_accumulator) * dKSdA2
        dKSdB2 = (1 / base_exp_accumulator) * dKSdB2
        dKSdC2 = (1 / base_exp_accumulator) * dKSdC2

        KS = (1 / rho) * log(base_exp_accumulator) - mindist_in
        intersect_length = intersect_accumulator
        return

    end subroutine compute_derivs
#endif

    subroutine compute(KS, intersect_length, mindist, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho)
        use triangles
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3, n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3, n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
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
                A1batch = A1(:, tri_ind_1)
                B1batch = B1(:, tri_ind_1)
                C1batch = C1(:, tri_ind_1)
                A2batch = A2(:, tri_ind_2)
                B2batch = B2(:, tri_ind_2)
                C2batch = C2(:, tri_ind_2)
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
                    ks_accumulator = ks_accumulator + exp((mindist_in - d) * rho)
                    ! compute the intersection on second pass only
                    call intersect(A1batch, B1batch, C1batch, &
                                   A2batch, B2batch, C2batch, intersect_temp)
                    if (intersect_temp /= intersect_temp) then
                        print *, 'NaN detected intersection routine'
                        print *, 'Triangle pairs (1-index): ', tri_ind_1, tri_ind_2
                    else
                        intersect_length = intersect_temp + intersect_length
                    end if
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
