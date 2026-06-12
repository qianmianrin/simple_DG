
    subroutine init_data

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

    call mesh

    hx = (xb - xa)/Nx0
    hy = (yb - ya)/Ny0
    hx1 = 0.5d0*hx
    hy1 = 0.5d0*hy
    hphi = 2d0*pi/Nphi1

    do i = 1,Nx0
        Xc0(i) = xa + (i - 0.5)*hx
    end do
    Xc = Xc0((myidx - 1)*Nx + 1:myidx*Nx)

    do j = 1,Ny0
        Yc0(j) = ya + (j - 0.5)*hy
    end do
    Yc = Yc0((myidy - 1)*Ny + 1:myidy*Ny)

    do k = 0,Nphi
        Phi(k) = k*hphi
    end do

    call get_basis

    uh = 0

    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                do d = 1,dimPk1
                    do i1 = 1,NumGLP
                        do j1 = 1,NumGLP
                            uh(i,j,k,d,1) = uh(i,j,k,d,1) + 0.25*weight(i1)*weight(j1)*U1(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,2) = uh(i,j,k,d,2) + 0.25*weight(i1)*weight(j1)*U2(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,3) = uh(i,j,k,d,3) + 0.25*weight(i1)*weight(j1)*U3(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                            uh(i,j,k,d,4) = uh(i,j,k,d,4) + 0.25*weight(i1)*weight(j1)*U4(Xc(i) + hx1*lambda(i1),Yc(j) + hy1*lambda(j1),Phi(k))*phiG(i1,j1,d)
                        end do
                    end do
                end do
            end do
        end do
    end do

    do d = 1,dimPk1
        uh(:,:,:,d,:) = uh(:,:,:,d,:)/mm(d)
    end do

    end subroutine init_data