
    subroutine set_bc

    use com

    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny_process
                    do i = 1,Nx_process

                    the_id = i + Nx_process*(j - 1)

                    
                    if (i == Nx_process) then

                    if (bcR == 1) then  
                        the_idx = 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                        end if
                    else if (bcR == 2) then  
                        if (myid1 == the_id) then
                            uh(Nx1,0:Ny1,k,d,n) = uh(Nx,0:Ny1,k,d,n)
                        end if
                    else if (bcR == 5) then 
                        if (myid1 == the_id) then
                            
                        end if   
                    else if (bcR == 4) then  
                        the_idx = 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (n == 6 .or. n == 7) then
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1,0:Ny1,k,d,n)*xb/xa,Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                            end if
                        else
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if

                    end if

                    else

                    the_idx = i + 1
                    the_idy = j
                    the_id2 = the_idx + Nx_process*(the_idy - 1)

                    if (myid1 == the_id2) then
                        call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                    end if
                    if (myid1 == the_id) then
                        call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                    end if
                    end if

                    if (i == 1) then

                    if (bcL == 1) then
                        the_idx = Nx_process
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                        end if
                    else if (bcL == 2) then
                        if (myid1 == the_id) then
                            uh(0,0:Ny1,k,d,n) = uh(1,0:Ny1,k,d,n)
                        end if
                    else if (bcL == 5) then 
                        if (myid1 == the_id) then
                            
                        end if
                    else if (bcL == 3) then
                        if (myid1 == the_id) then
                            
                        end if
                    else if (bcL == 4) then
                        the_idx = Nx_process
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (n == 6 .or. n == 7) then
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(Nx,0:Ny1,k,d,n)*xa/xb,Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                            end if
                        else
                            if (myid1 == the_id2) then
                                call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                            end if
                            if (myid1 == the_id) then
                                call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                            end if
                        end if
                    end if
                    else
                        the_idx = i - 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)

                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(Nx,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,2,MPI_COMM_WORLD,status,ierr)
                        end if
                    end if

                    if (j == Ny_process) then

                    if (bcU == 1) then
                        the_idx = i
                        the_idy = 1
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,Ny1,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                        end if
                    else if (bcU == 2) then !  
                        if (myid1 == the_id) then
                            uh(0:Nx1,Ny1,k,d,n) = uh(0:Nx1,Ny,k,d,n)
                        end if
                    else if (bcU == 5) then 
                        if (myid1 == the_id) then
                             
                        end if  
                    else if (bcU == 3) then
                        if (myid1 == the_id) then
                            uh(:,Ny1,:,:,:) = 0
                            do ii = 1,Nx
                               
                            end do
                        end if
                    end if

                    else
                        the_idx = i
                        the_idy = j + 1
                        the_id2 = the_idx + Nx_process*(the_idy - 1)

                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(0:Nx1,1,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,3,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,Ny1,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,3,MPI_COMM_WORLD,status,ierr)
                        end if

                    end if

 
                    if (j == 1) then

                    if (bcD == 1) then
                        the_idx = i
                        the_idy = Ny_process
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,0,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                        end if
                    else if (bcD == 2) then
                        if (myid1 == the_id) then
                            uh(0:Nx1,0,k,d,n) = uh(0:Nx1,1,k,d,n)
                        end if
                    else if (bcD == 3) then 
                        if (myid1 == the_id) then
                             
                        end if
                    else if (bcD == 4) then 
                        if (myid1 == the_id) then
                           
                        end if
                    end if
                    else
                        the_idx = i
                        the_idy = j - 1
                        the_id2 = the_idx + Nx_process*(the_idy - 1)

                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(0:Nx1,Ny,k,d,n),Nx + 2,MPI_REAL8,the_id - 1,4,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(0:Nx1,0,k,d,n),Nx + 2,MPI_REAL8,the_id2 - 1,4,MPI_COMM_WORLD,status,ierr)
                        end if

                    end if

                    end do
                end do
            end do
        end do 
    end do

    end subroutine set_bc