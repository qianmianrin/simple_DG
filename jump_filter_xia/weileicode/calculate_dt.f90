
    subroutine calculate_dt

    use com

    alphax = 0
    alphay = 0
    alphax0 = 0
    alphay0 = 0

    do i = 0,Nx1
        do j = 0,Ny1
            do k = 0,Nphi
                call eigenvalueMm(alpha1,alpha2,uh(i,j,k,1,1),uh(i,j,k,1,2),uh(i,j,k,1,3),uh(i,j,k,1,4),1,0)
                if (abs(alpha1) > alphax .or. abs(alpha2) > alphax) then
                    alphax = max(abs(alpha1),abs(alpha2))
                end if
                call eigenvalueMm(alpha1,alpha2,uh(i,j,k,1,1),uh(i,j,k,1,2),uh(i,j,k,1,3),uh(i,j,k,1,4),0,1)
                if (abs(alpha1) > alphay .or. abs(alpha2) > alphay) then
                    alphay = max(abs(alpha1),abs(alpha2))
                end if
            end do
        end do
    end do

    do the_id = 2,N_process

    if (myid1 == the_id) then
        call MPI_SEND(alphax,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == 1) then
        call MPI_RECV(alphax0,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
        if (alphax0 > alphax) then
            alphax = alphax0
        end if
    end if

    if (myid1 == the_id) then
        call MPI_SEND(alphay,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == 1) then
        call MPI_RECV(alphay0,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
        if (alphay0 > alphay) then
            alphay = alphay0
        end if
    end if

    end do

    do the_id = 2,N_process

    if (myid1 == 1) then
        call MPI_SEND(alphax,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == the_id) then
        call MPI_RECV(alphax,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
    end if

    if (myid1 == 1) then
        call MPI_SEND(alphay,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == the_id) then
        call MPI_RECV(alphay,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
    end if

    end do

    dt = CFL/(alphax/hx + alphay/hy)

    end subroutine calculate_dt