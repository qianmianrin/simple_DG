
    subroutine Lh

    use com

    real(8) Fx(NumGLP,NumGLP,0:Nphi,NumEq), Fy(NumGLP,NumGLP,0:Nphi,NumEq), Fz(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) rhoij,uij,vij,wij,Eij,B1ij,B2ij,B3ij,pij,Sij,Tij,Kij,rB1ij,rB2ij,rB3ij

    !real(8),allocatable :: UR(:,:,:,:,:),UL(:,:,:,:,:),UU(:,:,:,:,:),UD(:,:,:,:,:)
    real(8),allocatable :: FR(:,:,:,:,:),FL(:,:,:,:,:),FU(:,:,:,:,:),FD(:,:,:,:,:)
    real(8),allocatable :: Fxhat(:,:,:,:,:), Fyhat(:,:,:,:,:)

    !allocate(UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    !allocate(UL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    !allocate(UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    !allocate(UD(Nx,Ny1,0:Nphi,NumGLP,NumEq))

    allocate(FR(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FL(Nx1,Ny,0:Nphi,NumGLP,NumEq))
    allocate(FU(Nx,0:Ny,0:Nphi,NumGLP,NumEq))
    allocate(FD(Nx,Ny1,0:Nphi,NumGLP,NumEq))

    !allocate(URU(0:Nx1,0:Ny1, 0:Nphi, NumEq))
    !allocate(ULU(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(URD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(ULD(0:Nx1,0:Ny1,0:Nphi,NumEq))
    !allocate(EzVertex(0:Nx,0:Ny,0:Nphi))

    allocate(Fxhat(0:Nx,Ny,0:Nphi,NumGLP,NumEq))
    allocate(Fyhat(Nx,0:Ny,0:Nphi,NumGLP,NumEq))

    call set_bc

    du = 0

    do j = 1,Ny
        do i = 1,Nx

        uGint3D = 0
        do n = 1,NumEq
            do d = 1,dimPk
                do k = 0,Nphi
                    uGint3D(:,:,k,n) = uGint3D(:,:,k,n) + uh(i,j,k,d,n)*phiG(:,:,d)
                end do
            end do
        end do

        do k = 0,Nphi
            do j1 = 1,NumGLP
                do i1 = 1,NumGLP
                    rhoij = uGint3D(i1,j1,k,1)
                    uij = uGint3D(i1,j1,k,2)/rhoij
                    vij = uGint3D(i1,j1,k,3)/rhoij
                    Eij = uGint3D(i1,j1,k,4)

                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2))

                    Fx(i1,j1,k,1) = rhoij*uij
                    Fx(i1,j1,k,2) = rhoij*uij**2 + pij
                    Fx(i1,j1,k,3) = rhoij*uij*vij
                    Fx(i1,j1,k,4) = uij*(Eij + pij)

                    Fy(i1,j1,k,1) = rhoij*vij
                    Fy(i1,j1,k,2) = rhoij*uij*vij
                    Fy(i1,j1,k,3) = rhoij*vij**2 + pij
                    Fy(i1,j1,k,4) = vij*(Eij + pij)
                end do
            end do
        end do

        do n = 1,NumEq
            do d = 1,dimPk1
                do k = 0,Nphi
                    do j1 = 1,NumGLP
                        do i1 = 1,NumGLP
                            if (d > 1) then
                                du(i,j,k,d,n) = du(i,j,k,d,n) + 0.25d0*weight(i1)*weight(j1)*(Fx(i1,j1,k,n)*phixG(i1,j1,d) + Fy(i1,j1,k,n)*phiyG(i1,j1,d))
                            end if
                        end do
                    end do
                end do
            end do
        end do

        end do
    end do

    UR = 0
    UL = 0

    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR(i,j,k,:,n) = UR(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR(:,d)
                        UL(i + 1,j,k,:,n) = UL(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL(:,d)
                    end do
                end do
            end do
        end do
    end do

    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    rhoij = uR(i,j,k,j1,1)
                    uij = uR(i,j,k,j1,2)/rhoij
                    vij = uR(i,j,k,j1,3)/rhoij
                    Eij = uR(i,j,k,j1,4)

                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2))

                    FR(i,j,k,j1,1) = rhoij*uij
                    FR(i,j,k,j1,2) = rhoij*uij**2 + pij
                    FR(i,j,k,j1,3) = rhoij*uij*vij
                    FR(i,j,k,j1,4) = uij*(Eij + pij)
                end do
            end do
        end do
    end do

    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 1,Nx1
                    rhoij = uL(i,j,k,j1,1)
                    uij = uL(i,j,k,j1,2)/rhoij
                    vij = uL(i,j,k,j1,3)/rhoij
                    Eij = uL(i,j,k,j1,4)

                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2))

                    FL(i,j,k,j1,1) = rhoij*uij
                    FL(i,j,k,j1,2) = rhoij*uij**2 + pij
                    FL(i,j,k,j1,3) = rhoij*uij*vij
                    FL(i,j,k,j1,4) = uij*(Eij + pij)
                end do
            end do
        end do
    end do

    UU = 0
    UD = 0

    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 0,Ny
                    do i = 1,Nx
                        UU(i,j,k,:,n) = UU(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU(:,d)
                        UD(i,j + 1,k,:,n) = UD(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD(:,d)
                    end do
                end do
            end do
        end do
    end do

    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    rhoij = UU(i,j,k,i1,1)
                    uij = UU(i,j,k,i1,2)/rhoij
                    vij = UU(i,j,k,i1,3)/rhoij
                    Eij = UU(i,j,k,i1,4)

                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2))

                    FU(i,j,k,i1,1) = rhoij*vij
                    FU(i,j,k,i1,2) = rhoij*uij*vij
                    FU(i,j,k,i1,3) = rhoij*vij**2 + pij
                    FU(i,j,k,i1,4) = vij*(Eij + pij)
                end do
            end do
        end do
    end do

    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny1
                do i = 1,Nx
                    rhoij = UD(i,j,k,i1,1)
                    uij = UD(i,j,k,i1,2)/rhoij
                    vij = UD(i,j,k,i1,3)/rhoij
                    Eij = UD(i,j,k,i1,4)

                    pij = gamma1*(Eij - 0.5d0*rhoij*(uij**2 + vij**2))

                    FD(i,j,k,i1,1) = rhoij*vij
                    FD(i,j,k,i1,2) = rhoij*uij*vij
                    FD(i,j,k,i1,3) = rhoij*vij**2 + pij
                    FD(i,j,k,i1,4) = vij*(Eij + pij)
                end do
            end do
        end do
    end do

    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    call eigenvalueMm(SRmax,SRmin,UR(i,j,k,j1,1),UR(i,j,k,j1,2),UR(i,j,k,j1,3),UR(i,j,k,j1,4),1,0)
                    call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,k,j1,1),UL(i + 1,j,k,j1,2),UL(i + 1,j,k,j1,3),UL(i + 1,j,k,j1,4),1,0)
                    
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FL(i + 1,j,k,j1,:)
                    FL1 = FR(i,j,k,j1,:)
                    UR1 = UL(i + 1,j,k,j1,:)
                    UL1 = UR(i,j,k,j1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    end if
                    Fxhat(i,j,k,j1,:) = Fhat1
                end do
            end do
        end do
    end do

   
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    call eigenvalueMm(SRmax,SRmin,UU(i,j,k,i1,1),UU(i,j,k,i1,2),UU(i,j,k,i1,3),UU(i,j,k,i1,4),0,1)
                    call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,k,i1,1),UD(i,j + 1,k,i1,2),UD(i,j + 1,k,i1,3),UD(i,j + 1,k,i1,4),0,1)
                    
                    SR = max(SRmax,SLmax)
                    SL = min(SRmin,SLmin)
                    FR1 = FD(i,j + 1,k,i1,:)
                    FL1 = FU(i,j,k,i1,:)
                    UR1 = UD(i,j + 1,k,i1,:)
                    UL1 = UU(i,j,k,i1,:)
                    if (flux_type == 1) then
                        call LF_Flux
                    else if (flux_type == 2) then
                        call HLL_Flux
                    end if
                    Fyhat(i,j,k,i1,:) = Fhat1
                end do
            end do
        end do
    end do

   
    do n = 1,NumEq
        do d = 1,dimPk1
            do j1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hx)*weight(j1)*(Fxhat(i,j,k,j1,n)*phiGR(j1,d) - Fxhat(i - 1,j,k,j1,n)*phiGL(j1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do

    do n = 1,NumEq
        do d = 1,dimPk1
            do i1 = 1,NumGLP
                do k = 0,Nphi
                    do j = 1,Ny
                        do i = 1,Nx
                            du(i,j,k,d,n) = du(i,j,k,d,n) - (0.5d0/hy)*weight(i1)*(Fyhat(i,j,k,i1,n)*phiGU(i1,d) - Fyhat(i,j - 1,k,i1,n)*phiGD(i1,d))
                        end do
                    end do
                end do
            end do
        end do
    end do

    do d = 1,dimPk1
        du(:,:,:,d,:) = du(:,:,:,d,:)/mm(d)
    end do

    end subroutine Lh