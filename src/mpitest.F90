module mpitest
    implicit none

    contains
    subroutine scat_gath_test(A, n, sum_out, A_out)
       use mpi
       implicit none
       integer, intent(in) :: n ! array dimension
       real(kind=8), intent(in), dimension(3,n) :: A
       real(kind=8), intent(out) :: sum_out
       real(kind=8), intent(out), dimension(3,n) :: A_out
       real(kind=8), dimension(3,n) :: A_out_temp

       integer :: error, id, n_procs
       integer, allocatable :: proc_split(:), proc_disp(:)
       integer :: n_tris_per_proc_base, n_tris_remaining, proc_idx, displ
       integer, parameter :: n_dim = 3

       real(kind=8), allocatable :: A_local(:,:)
       real(kind=8) :: sum_local

       call MPI_Comm_size ( MPI_COMM_WORLD, n_procs, error )
       call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )

       allocate(proc_split(0:(n_procs-1)), proc_disp(0:(n_procs-1)))
       ! compute the even processor split
        n_tris_per_proc_base = n / n_procs
        n_tris_remaining = mod(n, n_procs)
        displ = 0
        do proc_idx = 0, n_procs-1
            if (proc_idx < n_tris_remaining) then
                proc_split(proc_idx) = (n_tris_per_proc_base + 1) * n_dim
            else
                proc_split(proc_idx) = (n_tris_per_proc_base) * n_dim
            end if
            proc_disp(proc_idx) = displ
            displ = displ + proc_split(proc_idx)
        end do

        ! now we have a proc_split array
        ! allocate a local array to work from and recieve the scatter
        allocate(A_local(3, proc_split(id)/3))
        call MPI_Scatterv(A, proc_split, proc_disp, &
                 MPI_DOUBLE, A_local, proc_split(id), &
                 MPI_DOUBLE, &
                 0, MPI_COMM_WORLD, error)
        !print *,'Rank: ',id,' ',A_local(:,1)
        sum_local = sum(A_local)
        
       ! print *,'Rank: ',id,' Sum: ', A_local

        ! Sum many scalars into one scalar and distribute back to all processes
        call MPI_Allreduce(sum_local, sum_out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)
        A_local = A_local * 2
        ! Concatenate an array with variable stride and distribute back to all procs
        call MPI_Allgatherv(A_local, proc_split(id), MPI_DOUBLE, A_out_temp, proc_split, proc_disp, &
        MPI_DOUBLE, MPI_COMM_WORLD, error)
        A_out_temp = A_out_temp * (id + 1)
        ! Element-wise array sum and redistribute
        call MPI_Allreduce(A_out_temp, A_out, 3*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, error)


    end subroutine scat_gath_test

end module mpitest