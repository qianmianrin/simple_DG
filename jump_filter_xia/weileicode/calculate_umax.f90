
    subroutine calculate_umax

    use com

    umax = 0
    umax1 = 0

    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                if (abs(uh(i,j,k,1,1)) > umax) then
                    umax = abs(uh(i,j,k,1,1))
                end if
            end do
        end do
    end do

    do the_id = 2,N_process

    if (myid1 == the_id) then
        call MPI_SEND(umax,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == 1) then
        call MPI_RECV(umax1,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
        if (umax1 > umax) then
            umax = umax1
        end if
    end if

    end do

    end subroutine calculate_umax
