module geograd_parallel
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
        use mpi
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho ! known global minimum distance (used for second pass, computing KS
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        real(kind=8), intent(out), dimension(3,n1) :: dKSdA1, dKSdB1, dKSdC1
        real(kind=8), intent(out), dimension(3,n2) :: dKSdA2, dKSdB2, dKSdC2
        integer :: tri_ind_1_local, tri_ind_2, minloc_index, count, is_less_flag ! loop indices
        real(kind=8) :: d, cur_min_dist, garbage, base_exp_accumulator, base_exp_accumulator_local
        real(kind=8), dimension(3) :: sumdA1, sumdB1, sumdC1, sumdA2, sumdB2, sumdC2
        real(kind=8), dimension(15) :: distance_vec ! holds the pairwise distances from the batch
        real(kind=8), dimension(3,15) :: dA1, dB1, dC1, dA2, dB2, dC2, deriv_exp_temp3 ! holds the reverse mode derivatives for the batch
        real(kind=8), dimension(15) :: dist_subtract_temp, base_exp_temp, deriv_exp_temp ! these hold some intermediate quantities that need to get accumulated
        real(kind=8) :: rev_seed
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch
        real(kind=8), allocatable :: dKSdA1_local(:,:), dKSdB1_local(:,:), dKSdC1_local(:,:)
        real(kind=8), dimension(3,n2) :: dKSdA2_local, dKSdB2_local, dKSdC2_local
        
        integer :: error, id, n_procs
        integer, allocatable :: proc_split(:), proc_disp(:)
        integer :: n_tris_per_proc_base, n_tris_remaining, proc_idx, displ, slice_start, slice_end
        integer, parameter :: n_dim = 3
 
        real(kind=8), allocatable :: A1_local(:,:), B1_local(:,:), C1_local(:,:)
 
        call MPI_Comm_size ( MPI_COMM_WORLD, n_procs, error )
        call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
 
        allocate(proc_split(0:(n_procs-1)), proc_disp(0:(n_procs-1)))
        ! compute the even processor split
        n_tris_per_proc_base = n1 / n_procs
        n_tris_remaining = mod(n1, n_procs)
        displ = 0
        do proc_idx = 0, n_procs-1
            if (proc_idx < n_tris_remaining) then
                proc_split(proc_idx) = (n_tris_per_proc_base + 1)
            else
                proc_split(proc_idx) = (n_tris_per_proc_base)
            end if
            proc_disp(proc_idx) = displ
            displ = displ + proc_split(proc_idx)
        end do
        ! get array slices
        allocate(A1_local(3, proc_split(id)), B1_local(3, proc_split(id)), C1_local(3, proc_split(id)))
        slice_start=proc_disp(id) + 1
        slice_end=proc_split(id) + slice_start - 1

        A1_local = A1(:,slice_start:slice_end)
        B1_local = B1(:,slice_start:slice_end)
        C1_local = C1(:,slice_start:slice_end)

        allocate(dKSdA1_local(3, proc_split(id)), &
                 dKSdB1_local(3, proc_split(id)), &
                 dKSdC1_local(3, proc_split(id)))

        dKSdA1_local = 0.0
        dKSdB1_local = 0.0
        dKSdC1_local = 0.0
        dKSdA2_local = 0.0
        dKSdB2_local = 0.0
        dKSdC2_local = 0.0
        base_exp_accumulator_local = 0.0
        cur_min_dist = 9.9e10


        do tri_ind_1_local = 1, proc_split(id)
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests and derivatives
                A1batch = A1_local(:,tri_ind_1_local)
                B1batch = B1_local(:,tri_ind_1_local)
                C1batch = C1_local(:,tri_ind_1_local)
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

                dist_subtract_temp = mindist_in - distance_vec
                base_exp_temp = exp(rho*dist_subtract_temp)
                base_exp_accumulator_local = base_exp_accumulator_local + sum(base_exp_temp)
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
                dKSdA1_local(:, tri_ind_1_local) = dKSdA1_local(:, tri_ind_1_local) + sumdA1
                dKSdB1_local(:, tri_ind_1_local) = dKSdB1_local(:, tri_ind_1_local) + sumdB1
                dKSdC1_local(:, tri_ind_1_local) = dKSdC1_local(:, tri_ind_1_local) + sumdC1
                dKSdA2_local(:, tri_ind_2) = dKSdA2_local(:, tri_ind_2) + sumdA2
                dKSdB2_local(:, tri_ind_2) = dKSdB2_local(:, tri_ind_2) + sumdB2
                dKSdC2_local(:, tri_ind_2) = dKSdC2_local(:, tri_ind_2) + sumdC2
            end do
        end do
        ! allreduce the base exp
        call MPI_Allreduce(base_exp_accumulator_local, base_exp_accumulator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)

        dKSdA1_local = (1 / base_exp_accumulator) * dKSdA1_local
        dKSdB1_local = (1 / base_exp_accumulator) * dKSdB1_local
        dKSdC1_local = (1 / base_exp_accumulator) * dKSdC1_local
        dKSdA2_local = (1 / base_exp_accumulator) * dKSdA2_local
        dKSdB2_local = (1 / base_exp_accumulator) * dKSdB2_local
        dKSdC2_local = (1 / base_exp_accumulator) * dKSdC2_local

        ! allgatherv the ABC1 derivatives
        call MPI_Allgatherv(dKSdA1_local, proc_split(id)*n_dim, MPI_DOUBLE, dKSdA1, proc_split*n_dim, proc_disp*n_dim, &
        MPI_DOUBLE, MPI_COMM_WORLD, error)

        call MPI_Allgatherv(dKSdB1_local, proc_split(id)*n_dim, MPI_DOUBLE, dKSdB1, proc_split*n_dim, proc_disp*n_dim, &
        MPI_DOUBLE, MPI_COMM_WORLD, error)        

        call MPI_Allgatherv(dKSdC1_local, proc_split(id)*n_dim, MPI_DOUBLE, dKSdC1, proc_split*n_dim, proc_disp*n_dim, &
        MPI_DOUBLE, MPI_COMM_WORLD, error)

        ! allreduce the ABC2 derivatives

        call MPI_Allreduce(dKSdA2_local, dKSdA2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)
        call MPI_Allreduce(dKSdB2_local, dKSdB2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)
        call MPI_Allreduce(dKSdC2_local, dKSdC2, 3*n2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)

        KS = (1/rho) * log(base_exp_accumulator) - mindist_in
        ! allreduce the mindist
        call MPI_Allreduce(cur_min_dist, mindist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, error)
        if (mindist_in /= mindist) then
            print *, 'These should match'
        end if
        intersect_length = 0
        ! TODO put intersection back in here
        return

    end subroutine compute_derivs
#endif

    subroutine compute(KS, intersect_length, mindist, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho)
        use triangles
        use mpi
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho ! known global minimum distance (used for second pass, computing KS)
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        integer :: tri_ind_1, tri_ind_1_local, tri_ind_2, minloc_index ! loop indices
        real(kind=8) :: d, cur_min_dist, intersect_temp, intersect_length_local, ks_accumulator_local, ks_accumulator
        real(kind=8), dimension(15) :: distance_vec
        real(kind=8), PARAMETER :: second_pass_flag = -1.0
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch

        integer :: error, id, n_procs
        integer, allocatable :: proc_split(:), proc_disp(:)
        integer :: n_tris_per_proc_base, n_tris_remaining, proc_idx, displ, slice_start, slice_end
        integer, parameter :: n_dim = 3
 
        real(kind=8), allocatable :: A1_local(:,:), B1_local(:,:), C1_local(:,:)
 
        call MPI_Comm_size ( MPI_COMM_WORLD, n_procs, error )
        call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
 
        allocate(proc_split(0:(n_procs-1)), proc_disp(0:(n_procs-1)))
        ! compute the even processor split
        n_tris_per_proc_base = n1 / n_procs
        n_tris_remaining = mod(n1, n_procs)
        displ = 0
        do proc_idx = 0, n_procs-1
            if (proc_idx < n_tris_remaining) then
                proc_split(proc_idx) = (n_tris_per_proc_base + 1)
            else
                proc_split(proc_idx) = (n_tris_per_proc_base)
            end if
            proc_disp(proc_idx) = displ
            displ = displ + proc_split(proc_idx)
        end do
        ! get array slices
        allocate(A1_local(3, proc_split(id)), B1_local(3, proc_split(id)), C1_local(3, proc_split(id)))
        slice_start=proc_disp(id) + 1
        slice_end=proc_split(id) + slice_start - 1

        A1_local = A1(:,slice_start:slice_end)
        B1_local = B1(:,slice_start:slice_end)
        C1_local = C1(:,slice_start:slice_end)

        intersect_length_local = 0.0_8
        cur_min_dist = 9.9e10
        ks_accumulator_local = 0.0_8

        do tri_ind_1_local = 1, proc_split(id)
            do tri_ind_2 = 1, n2
                ! do 9 line-line comparison tests
                A1batch = A1_local(:,tri_ind_1_local)
                B1batch = B1_local(:,tri_ind_1_local)
                C1batch = C1_local(:,tri_ind_1_local)
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
                    ks_accumulator_local = ks_accumulator_local + sum(exp((mindist_in - distance_vec)*rho))
                ! compute the intersection on second pass only
                    call intersect(A1batch, B1batch, C1batch, &
                    A2batch, B2batch, C2batch, intersect_temp)
                    intersect_length_local = intersect_temp + intersect_length_local
                end if
            end do
        end do
        ! do the reductions
        call MPI_Allreduce(ks_accumulator_local, ks_accumulator, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)
        call MPI_Allreduce(intersect_length_local, intersect_length, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)
        call MPI_Allreduce(cur_min_dist, mindist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, error)
        if (mindist_in /= second_pass_flag) then
            KS = (1 / rho) * log(ks_accumulator) - mindist_in
        else
            KS = zero
        end if
    end subroutine compute

end module geograd_parallel
