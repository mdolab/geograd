module geograd_parallel
    implicit none

    contains
    subroutine even_proc_split(proc_split, proc_disp, n_tris, n_procs)
        implicit none
        integer, intent(in) :: n_tris, n_procs
        integer, intent(out), dimension(0:(n_procs-1)) :: proc_split, proc_disp
        integer :: n_tris_per_proc_base, n_tris_remaining, displ, proc_idx
        n_tris_per_proc_base = n_tris / n_procs
        n_tris_remaining = mod(n_tris, n_procs)
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
    end subroutine even_proc_split

    subroutine load_balance_split(proc_split, proc_disp, bb_test_result, n_tris, n_procs)
        implicit none
        integer, intent(in) :: n_tris, n_procs
        integer, intent(in), dimension(n_tris) :: bb_test_result
        integer, intent(out), dimension(0:(n_procs-1)) :: proc_split, proc_disp
        integer :: proc_idx, tri_idx, cumsum, bb_active_count, start_val_this_proc
        integer, dimension(0:(n_procs-1)) :: active_split, active_disp

        bb_active_count = sum(bb_test_result) ! total number of triangles that need expensive computation
        if (bb_active_count < n_procs) then
            ! degenerate case where there is less than one triangle to compute per processor
            return
        end if

        call even_proc_split(active_split, active_disp, bb_active_count, n_procs) ! split the active triangles evenly across procs

        tri_idx = 1
        cumsum = bb_test_result(tri_idx) !keep a running total
        proc_disp(0) = 0 ! first processor always starts at first entry in the array

        do proc_idx = 1, n_procs-1
            start_val_this_proc = active_disp(proc_idx) + 1
            do while (cumsum /= start_val_this_proc)
                
                tri_idx = tri_idx + 1
                cumsum = cumsum + bb_test_result(tri_idx)
            end do
            proc_disp(proc_idx) = tri_idx - 1
            proc_split(proc_idx - 1) = proc_disp(proc_idx) - proc_disp(proc_idx - 1)
        end do
        proc_split(n_procs-1) = n_tris - proc_disp(n_procs-1)
    end subroutine load_balance_split

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

    subroutine check_bb_tol(obj_mins, obj_maxs, A2, B2, C2, n2, obj_tol, id)
        real(kind=8), dimension(3,n2), intent(in) :: A2, B2, C2
        real(kind=8), dimension(3), intent(out) :: obj_mins, obj_maxs
        real(kind=8), intent(in) :: obj_tol
        integer, intent(in) :: id, n2
        real(kind=8), dimension(3) :: obj_mins_A, obj_mins_B, obj_mins_C, &
                                      obj_maxs_A, obj_maxs_B, obj_maxs_C        
        real(kind=8) :: obj_dx, obj_dy, obj_dz, obj_max_d
        
        ! compute the maximum and minimum extent of the second mesh in the x direction
        integer :: tri_idx
        obj_mins = 9e10
        obj_maxs = -9e10
        do tri_idx = 1,n2
            obj_maxs(1) = max(obj_maxs(1), A2(1,tri_idx), B2(1, tri_idx), C2(1, tri_idx))
            obj_maxs(2) = max(obj_maxs(2), A2(2,tri_idx), B2(2, tri_idx), C2(2, tri_idx))
            obj_maxs(3) = max(obj_maxs(3), A2(3,tri_idx), B2(3, tri_idx), C2(3, tri_idx))
            obj_mins(1) = min(obj_mins(1), A2(1,tri_idx), B2(1, tri_idx), C2(1, tri_idx))
            obj_mins(2) = min(obj_mins(2), A2(2,tri_idx), B2(2, tri_idx), C2(2, tri_idx))
            obj_mins(3) = min(obj_mins(3), A2(3,tri_idx), B2(3, tri_idx), C2(3, tri_idx))
        end do
        obj_dx = obj_maxs(1) - obj_mins(1)
        obj_dy = obj_maxs(2) - obj_mins(2)
        obj_dz = obj_maxs(3) - obj_mins(3)
        obj_max_d = max(obj_dx, obj_dy, obj_dz)
        if ((obj_max_d > obj_tol) .AND. (id == 0)) then
            print *,'Object bounding box tol ',obj_tol, 'is smaller than object max dimension ',obj_max_d
        end if
    end subroutine check_bb_tol

#ifndef USE_COMPLEX
subroutine compute_derivs(KS, intersect_length, mindist, dKSdA1, dKSdB1, dKSdC1, dKSdA2, dKSdB2, dKSdC2, &
    A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho, obj_tol_in)
   use triangles_db
   use mpi
   implicit none
   integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
   real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
   real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
   real(kind=8), INTENT(in) :: mindist_in, rho, obj_tol_in ! known global minimum distance (used for second pass, computing KS)
   real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
   real(kind=8), intent(out), dimension(3,n1) :: dKSdA1, dKSdB1, dKSdC1
   real(kind=8), intent(out), dimension(3,n2) :: dKSdA2, dKSdB2, dKSdC2
   integer :: tri_ind_1_local, tri_ind_2, minloc_index, inner_count ! loop indices
   real(kind=8) :: d, cur_min_dist, intersect_temp, intersect_length_local, &
                   garbage, base_exp_accumulator, base_exp_accumulator_local, obj_tol
   real(kind=8), dimension(3) :: sumdA1, sumdB1, sumdC1, sumdA2, sumdB2, sumdC2
   real(kind=8), dimension(15) :: distance_vec ! holds the pairwise distances from the batch
   real(kind=8), dimension(3,15) :: dA1, dB1, dC1, dA2, dB2, dC2 ! holds the reverse mode derivatives for the batch
   real(kind=8), dimension(15) :: dist_subtract_temp, base_exp_temp ! these hold some intermediate quantities that need to get accumulated
   real(kind=8) :: rev_seed
   real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch
   real(kind=8), allocatable :: dKSdA1_local(:,:), dKSdB1_local(:,:), dKSdC1_local(:,:)
   real(kind=8), dimension(3,n2) :: dKSdA2_local, dKSdB2_local, dKSdC2_local
   
   integer :: error, id, n_procs, outside_bb_flag
   integer, allocatable :: proc_split(:), proc_disp(:), bb_flag_vec_local(:), bb_flag_vec(:)
   integer, parameter :: n_dim = 3

   real(kind=8), dimension(3) :: obj_mins, obj_maxs, this_tri_mins, this_tri_maxs
   real(kind=8) :: obj_bb_xmin, obj_bb_xmax, obj_bb_ymin, obj_bb_ymax, &
                   obj_bb_zmin, obj_bb_zmax

   real(kind=8) :: start_time, end_time, elapsed_time, max_time, min_time
   call MPI_Comm_size ( MPI_COMM_WORLD, n_procs, error )
   call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
   allocate(proc_split(0:(n_procs-1)), proc_disp(0:(n_procs-1)))

   ! get array slices
   mindist = 9.9e10
   obj_tol = obj_tol_in

   call even_proc_split(proc_split, proc_disp, n1, n_procs)
   allocate(bb_flag_vec_local(proc_split(id)), bb_flag_vec(n1))

   do while (real(mindist) > 9e10)
       ! if the bounding box margin is too small, the minimum distance test can fail to find a point.
       ! if it fails, expand the bounding box margin

       ! split all the triangles evenly across the processors to compute the bounding box tests
       call even_proc_split(proc_split, proc_disp, n1, n_procs)


       cur_min_dist = 9.9e10
       intersect_length_local = 0.0_8
       base_exp_accumulator_local = 0.0
       call check_bb_tol(obj_mins, obj_maxs, A2, B2, C2, n2, obj_tol, id)
       ! compute the bounds of the box (add the width in both directions)
       obj_bb_xmin = obj_mins(1) - obj_tol
       obj_bb_xmax = obj_maxs(1) + obj_tol 
       obj_bb_ymin = obj_mins(2) - obj_tol
       obj_bb_ymax = obj_maxs(2) + obj_tol   
       obj_bb_zmin = obj_mins(3) - obj_tol
       obj_bb_zmax = obj_maxs(3) + obj_tol

       bb_flag_vec_local = 1 ! initialize to 1
       do tri_ind_1_local = 1, proc_split(id)
           ! compute the bounding box test here for the load balance

           A1batch = A1(:,tri_ind_1_local+proc_disp(id))
           B1batch = B1(:,tri_ind_1_local+proc_disp(id))
           C1batch = C1(:,tri_ind_1_local+proc_disp(id))
           ! do a cheap bounding box check and potentially skip the n2 loop
           ! find the extents of A1
           this_tri_mins = min(A1batch, B1batch, C1batch)
           this_tri_maxs = max(A1batch, B1batch, C1batch)
           
           ! check spanwise direction first
           outside_bb_flag = -1
           if (this_tri_mins(3) > obj_bb_zmax) then
               bb_flag_vec_local(tri_ind_1_local) = 0
           else if (this_tri_maxs(3) < obj_bb_zmin) then
               bb_flag_vec_local(tri_ind_1_local) = 0
           else if (this_tri_mins(1) > obj_bb_xmax) then
               bb_flag_vec_local(tri_ind_1_local) = 0
           else if (this_tri_maxs(1) < obj_bb_xmin) then
               bb_flag_vec_local(tri_ind_1_local) = 0
           else if (this_tri_mins(2) > obj_bb_ymax) then 
               bb_flag_vec_local(tri_ind_1_local) = 0
           else if (this_tri_maxs(2) < obj_bb_ymin) then
               bb_flag_vec_local(tri_ind_1_local) = 0
           end if
       end do

       call MPI_Allgatherv(bb_flag_vec_local, proc_split(id), MPI_INTEGER, bb_flag_vec, proc_split, proc_disp, &
                           MPI_INTEGER, MPI_COMM_WORLD, error)
       
       call load_balance_split(proc_split, proc_disp, bb_flag_vec, n1, n_procs)
       ! CALL CPU_TIME(real(start_time))
       allocate(dKSdA1_local(3, proc_split(id)), &
               dKSdB1_local(3, proc_split(id)), &
               dKSdC1_local(3, proc_split(id)))

       dKSdA1_local = 0.0
       dKSdB1_local = 0.0
       dKSdC1_local = 0.0
       dKSdA2_local = 0.0
       dKSdB2_local = 0.0
       dKSdC2_local = 0.0
       do tri_ind_1_local = 1, proc_split(id)
           ! check this triangle is getting computed
           ! if no, compute a KS contribution equal to the width of the component times the number of facets
           ! zero out the gradient entries
           ! if yes, do the loop below
           A1batch = A1(:,tri_ind_1_local+proc_disp(id))
           B1batch = B1(:,tri_ind_1_local+proc_disp(id))
           C1batch = C1(:,tri_ind_1_local+proc_disp(id))
          
           if (bb_flag_vec(tri_ind_1_local+proc_disp(id)) == 0) then
               ! do not bother computing the actual pairwise tests. Add a conservative estimate
               base_exp_accumulator_local = base_exp_accumulator_local + exp(-obj_tol*rho)*15
           else                
               do tri_ind_2 = 1, n2
                   ! do 9 line-line comparison tests and derivatives
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
                   call compare_and_swap_minimum(cur_min_dist, d)

                   dist_subtract_temp = mindist_in - distance_vec
                   base_exp_temp = exp(rho*dist_subtract_temp)
                   base_exp_accumulator_local = base_exp_accumulator_local + sum(base_exp_temp)
                   sumdA1 = 0.0
                   sumdB1 = 0.0
                   sumdC1 = 0.0
                   sumdA2 = 0.0
                   sumdB2 = 0.0
                   sumdC2 = 0.0
                   do inner_count = 1,15
                       sumdA1 = sumdA1 - dA1(:,inner_count)*base_exp_temp(inner_count)
                       sumdB1 = sumdB1 - dB1(:,inner_count)*base_exp_temp(inner_count)
                       sumdC1 = sumdC1 - dC1(:,inner_count)*base_exp_temp(inner_count)
                       sumdA2 = sumdA2 - dA2(:,inner_count)*base_exp_temp(inner_count)
                       sumdB2 = sumdB2 - dB2(:,inner_count)*base_exp_temp(inner_count)
                       sumdC2 = sumdC2 - dC2(:,inner_count)*base_exp_temp(inner_count)
                       
                   end do
                   dKSdA1_local(:, tri_ind_1_local) = dKSdA1_local(:, tri_ind_1_local) + sumdA1
                   dKSdB1_local(:, tri_ind_1_local) = dKSdB1_local(:, tri_ind_1_local) + sumdB1
                   dKSdC1_local(:, tri_ind_1_local) = dKSdC1_local(:, tri_ind_1_local) + sumdC1
                   dKSdA2_local(:, tri_ind_2) = dKSdA2_local(:, tri_ind_2) + sumdA2
                   dKSdB2_local(:, tri_ind_2) = dKSdB2_local(:, tri_ind_2) + sumdB2
                   dKSdC2_local(:, tri_ind_2) = dKSdC2_local(:, tri_ind_2) + sumdC2
               end do
           end if
           ! CALL CPU_TIME(real(end_time))
           !elapsed_time = end_time - start_time
       end do
       ! allreduce the base exp
       call MPI_Allreduce(base_exp_accumulator_local, base_exp_accumulator, 1, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)

       dKSdA1_local = (1 / base_exp_accumulator) * dKSdA1_local
       dKSdB1_local = (1 / base_exp_accumulator) * dKSdB1_local
       dKSdC1_local = (1 / base_exp_accumulator) * dKSdC1_local
       dKSdA2_local = (1 / base_exp_accumulator) * dKSdA2_local
       dKSdB2_local = (1 / base_exp_accumulator) * dKSdB2_local
       dKSdC2_local = (1 / base_exp_accumulator) * dKSdC2_local

       ! allgatherv the ABC1 derivatives
       call MPI_Allgatherv(dKSdA1_local, proc_split(id)*n_dim, MPI_DOUBLE_PRECISION, dKSdA1, &
       proc_split*n_dim, proc_disp*n_dim, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, error)

       call MPI_Allgatherv(dKSdB1_local, proc_split(id)*n_dim, MPI_DOUBLE_PRECISION, dKSdB1, &
       proc_split*n_dim, proc_disp*n_dim, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, error)        

       call MPI_Allgatherv(dKSdC1_local, proc_split(id)*n_dim, MPI_DOUBLE_PRECISION, dKSdC1, &
       proc_split*n_dim, proc_disp*n_dim, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, error)

       ! allreduce the ABC2 derivatives

       call MPI_Allreduce(dKSdA2_local, dKSdA2, 3*n2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
       call MPI_Allreduce(dKSdB2_local, dKSdB2, 3*n2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
       call MPI_Allreduce(dKSdC2_local, dKSdC2, 3*n2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)

       KS = (1/rho) * log(base_exp_accumulator) - mindist_in
       ! allreduce the mindist
       call MPI_Allreduce(cur_min_dist, mindist, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, error)
       ! compute the imbalance
       ! call MPI_Allreduce(elapsed_time, max_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, error)
       ! call MPI_Allreduce(elapsed_time, min_time, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, error)
    !    if (id == 0) then
    !        print*,'Imbalance: ',(max_time - min_time) * 100 / max_time
    !    end if
      

       intersect_length = 0
       ! TODO put intersection back in here
       obj_tol = obj_tol * 2.0       
       deallocate(dKSdA1_local,dKSdB1_local,dKSdC1_local)
    end do
    if (mindist_in /= mindist) then
        print *, 'These should match: in: ',mindist_in, 'calc: ', mindist
    end if
end subroutine compute_derivs
#endif

    subroutine compute(KS, intersect_length, mindist, A1, B1, C1, A2, B2, C2, n1, n2, mindist_in, rho, obj_tol_in)
        use triangles
        use mpi
        implicit none
        integer, intent(in) :: n1, n2 ! array dimensions of first and second triangulated surfaces
        real(kind=8), dimension(3,n1), INTENT(in) :: A1, B1, C1 ! first triangulated surface vertices
        real(kind=8), dimension(3,n2), INTENT(in) :: A2, B2, C2 ! second triangulated surface vertices
        real(kind=8), INTENT(in) :: mindist_in, rho, obj_tol_in ! known global minimum distance (used for second pass, computing KS)
        real(kind=8), intent(out) :: KS, intersect_length, mindist ! results
        integer :: tri_ind_1_local, tri_ind_2, minloc_index ! loop indices
        real(kind=8) :: d, cur_min_dist, intersect_temp, intersect_length_local, ks_accumulator_local, ks_accumulator, obj_tol
        real(kind=8), dimension(15) :: distance_vec
        real(kind=8), PARAMETER :: second_pass_flag = -1.0
        real(kind=8), dimension(3) :: A1batch, B1batch, C1batch, A2batch, B2batch, C2batch

        integer :: error, id, n_procs, outside_bb_flag
        integer, allocatable :: proc_split(:), proc_disp(:), bb_flag_vec_local(:), bb_flag_vec(:)
        integer, parameter :: n_dim = 3
 
        real(kind=8), dimension(3) :: obj_mins, obj_maxs, this_tri_mins, this_tri_maxs
        real(kind=8) :: obj_bb_xmin, obj_bb_xmax, obj_bb_ymin, obj_bb_ymax, &
                        obj_bb_zmin, obj_bb_zmax

        real(kind=8) :: start_time, end_time, elapsed_time, max_time, min_time
        call MPI_Comm_size ( MPI_COMM_WORLD, n_procs, error )
        call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
        allocate(proc_split(0:(n_procs-1)), proc_disp(0:(n_procs-1)))
 


        ! idea? fixed tol in all dirs vs box extents

        ! get array slices
        mindist = 9.9e10
        obj_tol = obj_tol_in

        call even_proc_split(proc_split, proc_disp, n1, n_procs)
        allocate(bb_flag_vec_local(proc_split(id)), bb_flag_vec(n1))

        do while (real(mindist) > 9e10)
            ! if the bounding box margin is too small, the minimum distance test can fail to find a point.
            ! if it fails, expand the bounding box margin
            ! print *, 'Beginning the do-while obj tol: ',obj_tol
            ! split all the triangles evenly across the processors to compute the bounding box tests
            call even_proc_split(proc_split, proc_disp, n1, n_procs)

            cur_min_dist = 9.9e10
            intersect_length_local = 0.0_8
            ks_accumulator_local = 0.0_8
            call check_bb_tol(obj_mins, obj_maxs, A2, B2, C2, n2, obj_tol, id)
            ! compute the bounds of the box (add the width in both directions)
            obj_bb_xmin = obj_mins(1) - obj_tol
            obj_bb_xmax = obj_maxs(1) + obj_tol 
            obj_bb_ymin = obj_mins(2) - obj_tol
            obj_bb_ymax = obj_maxs(2) + obj_tol   
            obj_bb_zmin = obj_mins(3) - obj_tol
            obj_bb_zmax = obj_maxs(3) + obj_tol

            bb_flag_vec_local = 1 ! initialize to 1
            do tri_ind_1_local = 1, proc_split(id)
                ! compute the bounding box test here for the load balance

                A1batch = A1(:,tri_ind_1_local+proc_disp(id))
                B1batch = B1(:,tri_ind_1_local+proc_disp(id))
                C1batch = C1(:,tri_ind_1_local+proc_disp(id))
                ! do a cheap bounding box check and potentially skip the n2 loop
                ! find the extents of A1
                this_tri_mins = min(A1batch, B1batch, C1batch)
                this_tri_maxs = max(A1batch, B1batch, C1batch)
                
                ! check spanwise direction first
                outside_bb_flag = -1
                if (this_tri_mins(3) > obj_bb_zmax) then
                    bb_flag_vec_local(tri_ind_1_local) = 0
                else if (this_tri_maxs(3) < obj_bb_zmin) then
                    bb_flag_vec_local(tri_ind_1_local) = 0
                else if (this_tri_mins(1) > obj_bb_xmax) then
                    bb_flag_vec_local(tri_ind_1_local) = 0
                else if (this_tri_maxs(1) < obj_bb_xmin) then
                    bb_flag_vec_local(tri_ind_1_local) = 0
                else if (this_tri_mins(2) > obj_bb_ymax) then 
                    bb_flag_vec_local(tri_ind_1_local) = 0
                else if (this_tri_maxs(2) < obj_bb_ymin) then
                    bb_flag_vec_local(tri_ind_1_local) = 0
                end if
            end do
            
            call MPI_Allgatherv(bb_flag_vec_local, proc_split(id), MPI_INTEGER, bb_flag_vec, proc_split, proc_disp, &
                                MPI_INTEGER, MPI_COMM_WORLD, error)
            
            call load_balance_split(proc_split, proc_disp, bb_flag_vec, n1, n_procs)
            ! CALL CPU_TIME(real(start_time))
            do tri_ind_1_local = 1, proc_split(id)
                ! TODO implement this
                ! check this triangle is getting computed
                ! if no, compute a KS contribution equal to the width of the component times the number of facets
                ! zero out the gradient entries
                ! if yes, do the loop below
                A1batch = A1(:,tri_ind_1_local+proc_disp(id))
                B1batch = B1(:,tri_ind_1_local+proc_disp(id))
                C1batch = C1(:,tri_ind_1_local+proc_disp(id))
               
                if (bb_flag_vec(tri_ind_1_local+proc_disp(id)) == 0) then
                    ! do not bother computing the actual pairwise tests. Add a conservative estimate
                    ks_accumulator_local = ks_accumulator_local + exp(-obj_tol*rho)*15
                else
                    do tri_ind_2 = 1, n2
                        ! do 9 line-line comparison tests

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
                end if
                ! CALL CPU_TIME(real(end_time))
                ! elapsed_time = end_time - start_time
            end do
            ! do the reductions
            
            ! compute the imbalance
            ! call MPI_Allreduce(elapsed_time, max_time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, error)
            ! call MPI_Allreduce(elapsed_time, min_time, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, error)
            ! if (id == 0) then
            !     print*,'Imbalance: ',(max_time - min_time) * 100 / max_time
            ! end if
            call MPI_Allreduce(ks_accumulator_local, ks_accumulator, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
            call MPI_Allreduce(intersect_length_local, intersect_length, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
#ifdef USE_COMPLEX
            call MPI_Allreduce(real(cur_min_dist), mindist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, error)
#else
            call MPI_Allreduce(cur_min_dist, mindist, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, error)
#endif
            if (mindist_in /= second_pass_flag) then
                KS = (1 / rho) * log(ks_accumulator) - mindist_in
            else
                KS = zero
            end if
            obj_tol = obj_tol * 2.0
        end do
    end subroutine compute

end module geograd_parallel
