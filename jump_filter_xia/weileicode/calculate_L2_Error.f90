
    subroutine calculate_L2_Error

    use com

    use init1

    real(8) U1
    U1(x,y,z) = rho(x,y,z)
    real(8) U2
    U2(x,y,z) = rho(x,y,z)*v1(x,y,z)
    real(8) U3
    U3(x,y,z) = rho(x,y,z)*v2(x,y,z)
    real(8) U4
    U4(x,y,z) = p(x,y,z)/gamma1 + 0.5d0*rho(x,y,z)*(v1(x,y,z)**2 + v2(x,y,z)**2)

    if (tend > 3) then
        tend = 0
    end if

    L2 = 0d0
 Linfty = 0d0
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi

            uGint = 0
            do d = 1,dimPk
                do n = 1,NumEq
                    uGint(:,:,n) = uGint(:,:,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                end do
            end do

            do i1 = 1,NumGLP
                do j1 = 1,NumGLP
                    L2(1) = L2(1) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,1) - U1(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                    L2(2) = L2(2) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,2) - U2(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                    L2(3) = L2(3) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,3) - U3(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                    L2(4) = L2(4) + 0.25*weight(i1)*weight(j1)*(uGint(i1,j1,4) - U4(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k)))**2
                end do
            end do

            end do
        end do
    end do

    do the_id = 2,N_process

    if (myid1 == the_id) then
        call MPI_SEND(L2,NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == 1) then
        call MPI_RECV(L2pre,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
        L2 = L2 + L2pre
    end if

    end do

    if (myid1 == 1) then
        L2 = (L2/(Nx0*Ny0*Nphi1))**0.5d0
    end if
     
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                uGint = 0
                do d = 1,dimPk
                    do n = 1,NumEq
                        uGint(:,:,n) = uGint(:,:,n) + uh(i,j,k,d,n)*phiG(:,:,d) 
                    end do
                end do

                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        Linfty(1) =  max(Linfty(1) ,abs( uGint(i1,j1,1) - U1(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k))) ) 
                        Linfty(2) =  max(Linfty(2) ,abs( uGint(i1,j1,2) - U2(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k))) ) 
                        Linfty(3) =  max(Linfty(3) ,abs( uGint(i1,j1,3) - U3(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k))) ) 
                        Linfty(4) =  max(Linfty(4) ,abs( uGint(i1,j1,4) - U4(Xc(i) + hx1*lambda(i1) - tend,Yc(j) + hy1*lambda(j1) - tend,Phi(k))) ) 
                     
                    end do
                end do
            end do  
        end do
    end do

    do the_id = 2,N_process

    if (myid1 == the_id) then
        call MPI_SEND(Linfty,NumEq,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
    end if

    if (myid1 == 1) then
        call MPI_RECV(Linftypre,NumEq,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
        do i=1,NumEq
        if (Linftypre(i)>Linfty(i)) then
        Linfty(i) = Linftypre(i)
        end if 
        end do
    end if

    end do
 
    
 
    if (myid1 == 1) then
    open(unit = 111,file = 'L2error.txt')
    open(unit = 112,file = 'Linftyerror.txt')
    write(111,*) L2
    write(112,*) Linfty
    end if 


    end subroutine calculate_L2_Error