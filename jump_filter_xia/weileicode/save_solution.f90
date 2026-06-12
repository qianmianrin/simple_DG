
    subroutine save_solution

    use com

    real uhsave(NumEq)
    real(kind=8) :: jumpa,omega11
    integer the_idx1,the_idy1
 
    if (myid1 == 1) then

    open(unit = 1,file = 'Q1.txt')
    open(unit = 2,file = 'Q2.txt')
    open(unit = 3,file = 'Q3.txt')
    open(unit = 4,file = 'Q4.txt')
    open(unit = 9,file = 'Xc.txt')
    open(unit = 10,file = 'Yc.txt')

    do i = 1,Nx0
        write(9,*) Xc0(i)
    end do

    do j = 1,Ny0
        write(10,*) Yc0(j)
    end do

    close(9)
    close(10)

    end if

    
    do j = 1,Ny0
        do i = 1,Nx0

        the_idx1 = mod(i,Nx)
        if (the_idx1 == 0) then
            the_idx1 = Nx
        end if
        the_idx = (i - the_idx1)/Nx + 1

        the_idy1 = mod(j,Ny)
        if (the_idy1 == 0) then
            the_idy1 = Ny
        end if
        the_idy = (j - the_idy1)/Ny + 1

        the_id = the_idx + Nx_process*(the_idy - 1)

        if (the_id /= 1) then
            if (myid1 == the_id) then
                call MPI_SEND(uh(the_idx1,the_idy1,0,1,:),NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
            end if
            if (myid1 == 1) then
                call MPI_RECV(uhsave,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
            end if
        else if (the_id == 1) then
            if (myid1 == 1) then
                uhsave = uh(the_idx1,the_idy1,0,1,:)
            end if
        end if

        if (myid1 == 1) then
            write(1,*) uhsave(1)
            write(2,*) uhsave(2)
            write(3,*) uhsave(3)
            write(4,*) uhsave(4)
        end if

        end do
    end do
   
    if (myid1 == 1) then
        open(unit = 1112,file = 'jumpa.txt')
    end if 

    do j = 1,Ny0
        do i = 1,Nx0

        the_idx1 = mod(i,Nx)
        if (the_idx1 == 0) then
            the_idx1 = Nx
        end if
        the_idx = (i - the_idx1)/Nx + 1

        the_idy1 = mod(j,Ny)
        if (the_idy1 == 0) then
            the_idy1 = Ny
        end if
        the_idy = (j - the_idy1)/Ny + 1

        the_id = the_idx + Nx_process*(the_idy - 1)

        if (the_id /= 1) then
            if (myid1 == the_id) then
                call MPI_SEND(jump(the_idx1,the_idy1),NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
            end if
            if (myid1 == 1) then
                call MPI_RECV(jumpa,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
            end if
        else if (the_id == 1) then
            if (myid1 == 1) then
                jumpa = jump(the_idx1,the_idy1)
            end if
        end if

        if (myid1 == 1) then
            write(1112,*) jumpa
        end if

        end do
    end do
    
    if (myid1 == 1) then
        open(unit = 11112,file = 'omega1.txt')
    end if 

    do j = 1,Ny0
        do i = 1,Nx0

        the_idx1 = mod(i,Nx)
        if (the_idx1 == 0) then
            the_idx1 = Nx
        end if
        the_idx = (i - the_idx1)/Nx + 1

        the_idy1 = mod(j,Ny)
        if (the_idy1 == 0) then
            the_idy1 = Ny
        end if
        the_idy = (j - the_idy1)/Ny + 1

        the_id = the_idx + Nx_process*(the_idy - 1)

        if (the_id /= 1) then
            if (myid1 == the_id) then
                call MPI_SEND(omega_i(the_idx1,the_idy1),NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
            end if
            if (myid1 == 1) then
                call MPI_RECV(omega11,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
            end if
        else if (the_id == 1) then
            if (myid1 == 1) then
                omega11 = omega_i(the_idx1,the_idy1)
            end if
        end if

        if (myid1 == 1) then
            write(11112,*) omega11 
        end if

        end do
    end do
    
    
    if (myid1 == 1) then
        close(1)
        close(2)
        close(3)
        close(4)
        close(5)
        close(6)
        close(7)
        close(8)
    end if

    end subroutine save_solution
