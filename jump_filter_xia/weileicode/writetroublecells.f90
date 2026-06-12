
     subroutine writetroubledcells
   
    use com

    integer :: i,j,d,ss
    real(8) Troubledcellsall 


    open(unit = 100,file = 'troubledcells.txt')

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
                call MPI_SEND(change_all(the_idx1,the_idy1),1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
            end if
            if (myid1 == 1) then
                call MPI_RECV(Troubledcellsall,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr) 
            end if
        else if (the_id == 1) then
            if (myid1 == 1) then
                Troubledcellsall = change_all(the_idx1,the_idy1)
            end if
        end if

        if (myid1 == 1) then
            write(100,*)  Troubledcellsall           
        end if

        end do
    end do

    end  subroutine writetroubledcells