    module com

    include 'mpif.h'

    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP
    parameter(N_process = 16)
    parameter(Nx0 = 80, Ny0 = 80, kk = 3, NumGLP = 5, flux_type = 1)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)

    real(8) pi,gamma,gamma1
    parameter(Lphi = 0)
    parameter(dimPk = (kk + 1)*(kk + 2)/2)
    parameter(dimPk1 = (kk + 1)*(kk + 2)/2)
    parameter(Nphi = max(2*Lphi - 1,0))
    parameter(Nphi1 = Nphi + 1)
    parameter(gamma = 1.4d0) ! other
    !parameter(gamma = 5d0/3d0) ! jet
    parameter(gamma1 = gamma - 1)
    parameter(pi = 4*atan(1d0))
    parameter(NumEq=4)
    parameter(RKorder=4)
    parameter(limitertypeall = 5)
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,tRK,t1,t2,alphax,alphay,totaldiv,rij
    real(8) uh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),du(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uI(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),uII(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uh00(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) hx,hy,Xc(Nx),Yc(Ny),Xc0(Nx0),Yc0(Ny0),Phi(0:Nphi),hx1,hy1,hphi
    real(8) Bx(0:Nx,0:Ny1,0:Nphi,kk + 1),By(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) dBx(0:Nx,0:Ny1,0:Nphi,kk + 1),dBy(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxI(0:Nx,0:Ny1,0:Nphi,kk + 1),ByI(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxII(0:Nx,0:Ny1,0:Nphi,kk + 1),ByII(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) lambda(NumGLP),weight(NumGLP),sink(0:Nphi,Lphi),cosk(0:Nphi,Lphi)
    real(8) phiG(NumGLP,NumGLP,dimPk),phixG(NumGLP,NumGLP,dimPk),phiyG(NumGLP,NumGLP,dimPk),mm(dimPk)
    real(8) phiGLL(NumGLP,NumGLP,dimPk,2),lambdaL(NumGLP)
    real(8) phiGR(NumGLP,dimPk), phiGL(NumGLP,dimPk), phiGU(NumGLP,dimPk), phiGD(NumGLP,dimPk)
    real(8) phiRU(dimPk), phiLU(dimPk), phiRD(dimPk), phiLD(dimPk)
    real(8) EzG(NumGLP,kk + 1),EzxG(NumGLP,kk + 1),EzyG(NumGLP,kk + 1),mmE(kk + 1)
    real(8) EzR(kk + 1),EzL(kk + 1),EzU(kk + 1),EzD(kk + 1),omega1(Nx,Ny,0:Nphi)
    real(8) uGint3D(NumGLP,NumGLP,0:Nphi,NumEq),uGint(NumGLP,NumGLP,NumEq)
    real(8) RHSC(NumGLP,NumGLP,0:Nphi,NumEq),RHSCopen,RG(NumGLP,NumGLP,0:Nphi,NumEq)
    real(8) RHS(NumGLP,NumGLP,0:Nphi,NumEq),Fzsin(Lphi),Fzcos(Lphi),Fzzsin(Lphi),Fzzcos(Lphi)
    real(8) FR1(NumEq),FL1(NumEq),UR1(NumEq),UL1(NumEq),Fhat1(NumEq),SR,SL
    real(8) URstar(NumEq),ULstar(NumEq),Ustar(NumEq),UUstar(NumEq),UDstar(NumEq)
    real(8) URU1(NumEq),ULU1(NumEq),URD1(NumEq),ULD1(NumEq)
    real(8) URstarstar(NumEq),ULstarstar(NumEq),Ezhat
    real(8) EzVertex(0:Nx,0:Ny,0:Nphi)
    real(8) URU(0:Nx1,0:Ny1,0:Nphi,NumEq),ULU(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) URD(0:Nx1,0:Ny1,0:Nphi,NumEq),ULD(0:Nx1,0:Ny1,0:Nphi,NumEq)
    real(8) L2(NumEq),L2pre(NumEq),Linfty(NumEq),Linftypre(NumEq)
    real(8) UR(0:Nx,Ny,0:Nphi,NumGLP,NumEq),UL(Nx1,Ny,0:Nphi,NumGLP,NumEq),UU(Nx,0:Ny,0:Nphi,NumGLP,NumEq),UD(Nx,Ny1,0:Nphi,NumGLP,NumEq) 
    real(8) M,beta
    real(8) DeltaUR1(NumEq,1),DeltaUL1(NumEq,1),DeltaUU1(NumEq,1),DeltaUD1(NumEq,1),DeltaU1(NumEq,1),DeltaUmod1(NumEq,1)
    real(8) DeltaUR1mod(NumEq,1),DeltaUL1mod(NumEq,1),DeltaUU1mod(NumEq,1),DeltaUD1mod(NumEq,1)
    real(8) R(NumEq,NumEq),L(NumEq,NumEq)
    real(8) DeltaUR(NumEq,1),DeltaUL(NumEq,1),DeltaU(NumEq,1),DeltaUmod(NumEq,1)
    real(8) Is_trouble_cell(Nx,Ny,0:Nphi)
    integer change_all(Nx,Ny)
    real(8) densityave,momentum1ave,momentum2ave,Enerave
    real(8) phiGR_ver(2,10),phiGL_ver(2,10), phiGU_ver(2,10), phiGD_ver(2,10)
    real(8) phiGR_ver_derx(2,10),phiGL_ver_derx(2,10), phiGU_ver_derx(2,10), phiGD_ver_derx(2,10)
    real(8) phiGR_ver_dery(2,10),phiGL_ver_dery(2,10), phiGU_ver_dery(2,10), phiGD_ver_dery(2,10)
    real(8) phiGR_ver_derxx(2,10),phiGL_ver_derxx(2,10), phiGU_ver_derxx(2,10), phiGD_ver_derxx(2,10)
    real(8) phiGR_ver_derxy(2,10),phiGL_ver_derxy(2,10), phiGU_ver_derxy(2,10), phiGD_ver_derxy(2,10)
    real(8) phiGR_ver_deryy(2,10),phiGL_ver_deryy(2,10), phiGU_ver_deryy(2,10), phiGD_ver_deryy(2,10)
    real(8) phiGR_ver_derxxx(2,10),phiGL_ver_derxxx(2,10), phiGU_ver_derxxx(2,10), phiGD_ver_derxxx(2,10)
    real(8) phiGR_ver_derxxy(2,10),phiGL_ver_derxxy(2,10), phiGU_ver_derxxy(2,10), phiGD_ver_derxxy(2,10)
    real(8) phiGR_ver_derxyy(2,10),phiGL_ver_derxyy(2,10), phiGU_ver_derxyy(2,10), phiGD_ver_derxyy(2,10)
    real(8) phiGR_ver_deryyy(2,10),phiGL_ver_deryyy(2,10), phiGU_ver_deryyy(2,10), phiGD_ver_deryyy(2,10)
    real(kind=8) UR_ver(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver(Nx1,Ny,0:Nphi,2,NumEq),UU_ver(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver(Nx,Ny1,0:Nphi,2,NumEq) 
    real(kind=8) UR_ver_derx(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derx(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_derx(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_derx(Nx,Ny1,0:Nphi,2,NumEq)
    real(kind=8) UR_ver_dery(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_dery(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_dery(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_dery(Nx,Ny1,0:Nphi,2,NumEq) 
    real(kind=8) UR_ver_derxx(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derxx(Nx1,Ny,0:Nphi,2,NumEq),UR_ver_derxy(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derxy(Nx1,Ny,0:Nphi,2,NumEq), UR_ver_deryy(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_deryy(Nx1,Ny,0:Nphi,2,NumEq)
    real(kind=8) UU_ver_derxx(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_derxx(Nx,Ny1,0:Nphi,2,NumEq) ,UU_ver_derxy(Nx,0:Ny,0:Nphi,2,NumEq) ,UD_ver_derxy(Nx,Ny1,0:Nphi,2,NumEq) ,UU_ver_deryy(Nx,0:Ny,0:Nphi,2,NumEq) ,UD_ver_deryy(Nx,Ny1,0:Nphi,2,NumEq) 
    real(kind=8) UR_ver_derxxx(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derxxx(Nx1,Ny,0:Nphi,2,NumEq),UR_ver_derxxy(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derxxy(Nx1,Ny,0:Nphi,2,NumEq), UR_ver_derxyy(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derxyy(Nx1,Ny,0:Nphi,2,NumEq), UR_ver_deryyy(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_deryyy(Nx1,Ny,0:Nphi,2,NumEq)
    real(kind=8) UU_ver_derxxx(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_derxxx(Nx,Ny1,0:Nphi,2,NumEq) ,UU_ver_derxxy(Nx,0:Ny,0:Nphi,2,NumEq) ,UD_ver_derxxy(Nx,Ny1,0:Nphi,2,NumEq) ,UU_ver_derxyy(Nx,0:Ny,0:Nphi,2,NumEq) ,UD_ver_derxyy(Nx,Ny1,0:Nphi,2,NumEq) ,UU_ver_deryyy(Nx,0:Ny,0:Nphi,2,NumEq) ,UD_ver_deryyy(Nx,Ny1,0:Nphi,2,NumEq)

    real(kind=8) jump(Nx,Ny)
    real(kind=8) omega_i(Nx,Ny)
    
    integer bcR,bcL,bcU,bcD,direction
    integer myid,myid1,the_id,the_id2
    integer myidx,myidy,the_idx,the_idy
    integer numprocs, namelen, rc,ierr,status(MPI_STATUS_SIZE),myid0
    character * (MPI_MAX_PROCESSOR_NAME) processor_name

    end module com

    module init1

    use com

    contains

    function rr(x,y)
    real(8) x,y,rr
    rr = (x - 5)**2 + y**2
    end function rr

    function rho(x,y,z)
    real(8) x,y,rho,z
    beta = 5
    rho = (1 - gamma1/(16*gamma*pi**2)*beta**2*exp(2*(1 - rr(x,y))))**(1d0/gamma1)
    end function rho

    function p(x,y,z)
    real(8) x,y,p,z
    p = rho(x,y,z)**gamma
    end function p

    function v1(x,y,z)
    real(8) v1,x,y,z
    beta = 5
    v1 = 1 - beta*exp(1 - rr(x,y))*y/(2*pi)
    end function v1

    function v2(x,y,z)
    real(8) v2,x,y,z
    beta = 5
    v2 = beta*exp(1 - rr(x,y))*(x - 5)/(2*pi)
    end function v2

    subroutine mesh

    use com

    xa = 0
    xb = 10
    ya = -5
    yb = 5

    bcR = 1
    bcL = 1
    bcU = 1
    bcD = 1

    tend = 10

    M = 100000000
    beta = 1

    end subroutine mesh

    end module init1

     
    !*****************************************************************************************************

    program main

    use com

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
    call CPU_TIME(t1)

    myid1 = myid + 1

    myidx = mod(myid1,Nx_process)
    if (myidx == 0) then
        myidx = Nx_process
    end if

    myidy = (myid1 - myidx)/Nx_process + 1

    print *,"process",myid1,"is alive,the index is",myidx,myidy

    call init_data

    if (RKorder == 1) then
    else if (RKorder == 3) then
        call RK3
    else if (RKorder == 4) then
        call RK4
    end if

    call set_bc

    call save_solution
    
    call  writetroubledcells
    
    call calculate_L2_Error

    call CPU_TIME(t2)

    if (myid1 == 1) then
        open(unit = 1,file = 'time.txt')
        write(1,*) t2 - t1
        print *,"Run time is",t2 - t1,"second"
        close(1)
    end if

    call MPI_FINALIZE(rc)

    end program main

    !*****************************************************************************************************

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

    !*****************************************************************************************************

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

    !*****************************************************************************************************

    subroutine get_basis

    use com

    ! get Gauss points
    if (NumGLP == 2) then
        lambda(1) = -0.5773502691896257645091488
        lambda(2) = 0.5773502691896257645091488

        weight(1) = 1
        weight(2) = 1
    else if (NumGLP == 3) then
        lambda(1) = -0.7745966692414833770358531d0
        lambda(2) = 0
        lambda(3) = 0.7745966692414833770358531d0

        weight(1) = 0.5555555555555555555555556d0
        weight(2) = 0.8888888888888888888888889d0
        weight(3) = 0.5555555555555555555555556d0
    else if (NumGLP == 4) then
        lambda(1) = -0.8611363115940525752239465d0
        lambda(2) = -0.3399810435848562648026658d0
        lambda(3) = 0.3399810435848562648026658d0
        lambda(4) = 0.8611363115940525752239465d0

        weight(1) = 0.3478548451374538573730639d0
        weight(2) = 0.6521451548625461426269361d0
        weight(3) = 0.6521451548625461426269361d0   
        weight(4) = 0.3478548451374538573730639d0
    else if (NumGLP == 5) then
        lambda(1) = -0.9061798459386639927976269d0     
        lambda(2) = -0.5384693101056830910363144d0     
        lambda(3) = 0d0                                 
        lambda(4) = 0.5384693101056830910363144d0     
        lambda(5) = 0.9061798459386639927976269d0     

        weight(1) = 0.2369268850561890875142640d0
        weight(2) = 0.4786286704993664680412915d0
        weight(3) = 0.5688888888888888888888889d0
        weight(4) = 0.4786286704993664680412915d0
        weight(5) = 0.2369268850561890875142640d0

        lambdaL(1) = -1
        lambdaL(2) = -0.6546536707079771437983
        lambdaL(3) = 0
        lambdaL(4) = 0.654653670707977143798
        lambdaL(5) = 1
    else if (NumGLP == 6) then
        lambda(1) = -0.9324695142031520278123016d0     
        lambda(2) = -0.6612093864662645136613996d0    
        lambda(3) = -0.2386191860831969086305017d0     
        !lambda(4) = 0.2386191860831969086305017d0     
        !lambda(5) = 0.6612093864662645136613996d0     
        !lambda(6) = 0.9324695142031520278123016d0     

        weight(1) = 0.1713244923791703450402961d0
        weight(2) = 0.3607615730481386075698335d0
        weight(3) = 0.4679139345726910473898703d0
        !weight(4) = 0.4679139345726910473898703d0
        !weight(5) = 0.3607615730481386075698335d0
        !weight(6) = 0.1713244923791703450402961d0
    end if

    do i = 1,NumGLP
        do j = 1,NumGLP
            phiG(i,j,1) = 1
            phiGLL(i,j,1,1) = 1
            phiGLL(i,j,1,2) = 1
            phiGR(j,1) = 1
            phiGL(j,1) = 1
            phiGU(i,1) = 1
            phiGD(i,1) = 1
            phixG(i,j,1) = 0
            phiyG(i,j,1) = 0
            mm(1) = 1

            phiG(i,j,2) = lambda(i)
            phiGLL(i,j,2,1) = lambdaL(i)
            phiGLL(i,j,2,2) = lambda(i)
            phiGR(j,2) = 1
            phiGL(j,2) = -1
            phiGU(i,2) = lambda(i)
            phiGD(i,2) = lambda(i)
            phixG(i,j,2) = 1d0/hx1
            phiyG(i,j,2) = 0
            mm(2) = 1d0/3d0

            phiG(i,j,3) = lambda(j)
            phiGLL(i,j,3,1) = lambda(j)
            phiGLL(i,j,3,2) = lambdaL(j)
            phiGR(j,3) = lambda(j)
            phiGL(j,3) = lambda(j)
            phiGU(i,3) = 1
            phiGD(i,3) = -1
            phixG(i,j,3) = 0
            phiyG(i,j,3) = 1d0/hy1
            mm(3) = 1d0/3d0

            phiG(i,j,4) = lambda(i)**2 - 1d0/3d0
            phiGLL(i,j,4,1) = lambdaL(i)**2 - 1d0/3d0
            phiGLL(i,j,4,2) = lambda(i)**2 - 1d0/3d0
            phiGR(j,4) = 2d0/3d0
            phiGL(j,4) = 2d0/3d0
            phiGU(i,4) = lambda(i)**2 - 1d0/3d0
            phiGD(i,4) = lambda(i)**2 - 1d0/3d0
            phixG(i,j,4) = 2d0*lambda(i)/hx1
            phiyG(i,j,4) = 0
            mm(4) = 4d0/45d0

            phiG(i,j,5) = lambda(i)*lambda(j)
            phiGLL(i,j,5,1) = lambdaL(i)*lambda(j)
            phiGLL(i,j,5,2) = lambda(i)*lambdaL(j)
            phiGR(j,5) = lambda(j)
            phiGL(j,5) = -lambda(j)
            phiGU(i,5) = lambda(i)
            phiGD(i,5) = -lambda(i)
            phixG(i,j,5) = lambda(j)/hx1
            phiyG(i,j,5) = lambda(i)/hy1
            mm(5) = 1d0/9d0

            phiG(i,j,6) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,1) = lambda(j)**2 - 1d0/3d0
            phiGLL(i,j,6,2) = lambdaL(j)**2 - 1d0/3d0
            phiGR(j,6) = lambda(j)**2 - 1d0/3d0
            phiGL(j,6) = lambda(j)**2 - 1d0/3d0
            phiGU(i,6) = 2d0/3d0
            phiGD(i,6) = 2d0/3d0
            phixG(i,j,6) = 0
            phiyG(i,j,6) = 2d0*lambda(j)/hy1
            mm(6) = 4d0/45d0

            phiG(i,j,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGLL(i,j,7,1) = lambdaL(i)**3 - 3d0*lambdaL(i)/5d0
            phiGLL(i,j,7,2) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGR(j,7) = 2d0/5d0
            phiGL(j,7) = -2d0/5d0
            phiGU(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phiGD(i,7) = lambda(i)**3 - 3d0*lambda(i)/5d0
            phixG(i,j,7) = (3*lambda(i)**2 - 3d0/5d0)/hx1
            phiyG(i,j,7) = 0
            mm(7) = 4d0/175d0

            phiG(i,j,8) = (lambda(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,1) = (lambdaL(i)**2 - 1d0/3d0)*(lambda(j))
            phiGLL(i,j,8,2) = (lambda(i)**2 - 1d0/3d0)*(lambdaL(j))
            phiGR(j,8) = (2d0/3d0)*(lambda(j))
            phiGL(j,8) = (2d0/3d0)*(lambda(j))
            phiGU(i,8) = (lambda(i)**2 - 1d0/3d0)
            phiGD(i,8) = -(lambda(i)**2 - 1d0/3d0)
            phixG(i,j,8) = 2d0*lambda(i)*lambda(j)/hx1
            phiyG(i,j,8) = (lambda(i)**2 - 1d0/3d0)/hy1
            mm(8) = 4d0/135d0

            phiG(i,j,9) = (lambda(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,1) = (lambdaL(i))*(lambda(j)**2 - 1d0/3d0)
            phiGLL(i,j,9,2) = (lambda(i))*(lambdaL(j)**2 - 1d0/3d0)
            phiGR(j,9) = (lambda(j)**2 - 1d0/3d0)
            phiGL(j,9) = -(lambda(j)**2 - 1d0/3d0)
            phiGU(i,9) = lambda(i)*(2d0/3d0)
            phiGD(i,9) = lambda(i)*(2d0/3d0)
            phixG(i,j,9) = (lambda(j)**2 - 1d0/3d0)/hx1
            phiyG(i,j,9) = 2d0*lambda(i)*lambda(j)/hy1
            mm(9) = 4d0/135d0

            phiG(i,j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,1) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGLL(i,j,10,2) = lambdaL(j)**3 - 3d0*lambdaL(j)/5d0
            phiGR(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGL(j,10) = lambda(j)**3 - 3d0*lambda(j)/5d0
            phiGU(i,10) = 2d0/5d0
            phiGD(i,10) = -2d0/5d0
            phixG(i,j,10) = 0
            phiyG(i,j,10) = (3*lambda(j)**2 - 3d0/5d0)/hy1
            mm(10) = 4d0/175d0
        end do

    end do
 
    phiGR_ver(1,1) = 1d0
    phiGL_ver(1,1) = 1d0
    phiGU_ver(1,1) = 1d0
    phiGD_ver(1,1) = 1d0

    phiGR_ver(2,1) = 1d0
    phiGL_ver(2,1) = 1d0
    phiGU_ver(2,1) = 1d0
    phiGD_ver(2,1) = 1d0

    phiGR_ver_derx(1,1) = 0d0
    phiGL_ver_derx(1,1) = 0d0
    phiGU_ver_derx(1,1) = 0d0
    phiGD_ver_derx(1,1) = 0d0

    phiGR_ver_derx(2,1) = 0d0
    phiGL_ver_derx(2,1) = 0d0
    phiGU_ver_derx(2,1) = 0d0
    phiGD_ver_derx(2,1) = 0d0

    phiGR_ver_dery(1,1) = 0d0
    phiGL_ver_dery(1,1) = 0d0
    phiGU_ver_dery(1,1) = 0d0
    phiGD_ver_dery(1,1) = 0d0

    phiGR_ver_dery(2,1) = 0d0
    phiGL_ver_dery(2,1) = 0d0
    phiGU_ver_dery(2,1) = 0d0
    phiGD_ver_dery(2,1) = 0d0

    phiGR_ver_derxx(1,1) = 0d0 
    phiGL_ver_derxx(1,1) = 0d0
    phiGU_ver_derxx(1,1) = 0d0
    phiGD_ver_derxx(1,1) = 0d0

    phiGR_ver_derxx(2,1) = 0d0 
    phiGL_ver_derxx(2,1) = 0d0 
    phiGU_ver_derxx(2,1) = 0d0
    phiGD_ver_derxx(2,1) = 0d0

    phiGR_ver_derxy(1,1) = 0d0 
    phiGL_ver_derxy(1,1) = 0d0 
    phiGU_ver_derxy(1,1) = 0d0
    phiGD_ver_derxy(1,1) = 0d0

    phiGR_ver_derxy(2,1) = 0d0 
    phiGL_ver_derxy(2,1) = 0d0 
    phiGU_ver_derxy(2,1) = 0d0
    phiGD_ver_derxy(2,1) = 0d0

    phiGR_ver_deryy(1,1) = 0d0 
    phiGL_ver_deryy(1,1) = 0d0 
    phiGU_ver_deryy(1,1) = 0d0
    phiGD_ver_deryy(1,1) = 0d0

    phiGR_ver_deryy(2,1) = 0d0 
    phiGL_ver_deryy(2,1) = 0d0 
    phiGU_ver_deryy(2,1) = 0d0
    phiGD_ver_deryy(2,1) = 0d0

     
    phiGR_ver_derxxx = 0d0 
    phiGR_ver_derxxy = 0d0
    phiGR_ver_derxyy = 0d0
    phiGR_ver_deryyy = 0d0 
    phiGL_ver_derxxx = 0d0 
    phiGL_ver_derxxy = 0d0
    phiGL_ver_derxyy = 0d0
    phiGL_ver_deryyy = 0d0 
    phiGU_ver_derxxx = 0d0 
    phiGU_ver_derxxy = 0d0
    phiGU_ver_derxyy = 0d0
    phiGU_ver_deryyy = 0d0 
    phiGD_ver_derxxx = 0d0 
    phiGD_ver_derxxy = 0d0
    phiGD_ver_derxyy = 0d0
    phiGD_ver_deryyy = 0d0 
   
    phiGR_ver(1,2) = 1d0
    phiGL_ver(1,2) = -1d0
    phiGU_ver(1,2) = -1d0
    phiGD_ver(1,2) = -1d0

    phiGR_ver(2,2) = 1d0
    phiGL_ver(2,2) = -1d0
    phiGU_ver(2,2) = 1d0
    phiGD_ver(2,2) = 1d0

    phiGR_ver_derx(1,2) = 1d0
    phiGL_ver_derx(1,2) = 1d0
    phiGU_ver_derx(1,2) = 1d0
    phiGD_ver_derx(1,2) = 1d0

    phiGR_ver_derx(2,2) = 1d0
    phiGL_ver_derx(2,2) = 1d0
    phiGU_ver_derx(2,2) = 1d0
    phiGD_ver_derx(2,2) = 1d0

    phiGR_ver_dery(1,2) = 0d0
    phiGL_ver_dery(1,2) = 0d0
    phiGU_ver_dery(1,2) = 0d0
    phiGD_ver_dery(1,2) = 0d0

    phiGR_ver_dery(2,2) = 0d0
    phiGL_ver_dery(2,2) =  0d0
    phiGU_ver_dery(2,2) = 0d0
    phiGD_ver_dery(2,2) = 0d0

    phiGR_ver_derxx(1,2) = 0d0 
    phiGL_ver_derxx(1,2) = 0d0 
    phiGU_ver_derxx(1,2) = 0d0
    phiGD_ver_derxx(1,2) = 0d0

    phiGR_ver_derxx(2,2) = 0d0 
    phiGL_ver_derxx(2,2) = 0d0 
    phiGU_ver_derxx(2,2) = 0d0
    phiGD_ver_derxx(2,2) = 0d0

    phiGR_ver_derxy(1,2) = 0d0 
    phiGL_ver_derxy(1,2) = 0d0 
    phiGU_ver_derxy(1,2) = 0d0
    phiGD_ver_derxy(1,2) = 0d0

    phiGR_ver_derxy(2,2) = 0d0 
    phiGL_ver_derxy(2,2) = 0d0 
    phiGU_ver_derxy(2,2) = 0d0
    phiGD_ver_derxy(2,2) = 0d0

    phiGR_ver_deryy(1,2) = 0d0 
    phiGL_ver_deryy(1,2) = 0d0 
    phiGU_ver_deryy(1,2) = 0d0
    phiGD_ver_deryy(1,2) = 0d0

    phiGR_ver_deryy(2,2) = 0d0 
    phiGL_ver_deryy(2,2) = 0d0 
    phiGU_ver_deryy(2,2) = 0d0
    phiGD_ver_deryy(2,2) = 0d0
 
    phiGR_ver(1,3) = 1d0
    phiGL_ver(1,3) = 1d0
    phiGU_ver(1,3) = 1d0
    phiGD_ver(1,3) = -1d0

    phiGR_ver(2,3) = -1d0
    phiGL_ver(2,3) = -1d0
    phiGU_ver(2,3) = 1d0
    phiGD_ver(2,3) = -1d0

    phiGR_ver_derx(1,3) = 0d0
    phiGL_ver_derx(1,3) = 0d0
    phiGU_ver_derx(1,3) = 0d0
    phiGD_ver_derx(1,3) = 0d0

    phiGR_ver_derx(2,3) = 0d0
    phiGL_ver_derx(2,3) = 0d0
    phiGU_ver_derx(2,3) = 0d0
    phiGD_ver_derx(2,3) = 0d0

    phiGR_ver_dery(1,3) = 1d0
    phiGL_ver_dery(1,3) = 1d0
    phiGU_ver_dery(1,3) = 1d0
    phiGD_ver_dery(1,3) = 1d0

    phiGR_ver_dery(2,3) = 1d0
    phiGL_ver_dery(2,3) = 1d0
    phiGU_ver_dery(2,3) = 1d0
    phiGD_ver_dery(2,3) = 1d0

    phiGR_ver_derxx(1,3) = 0d0 
    phiGL_ver_derxx(1,3) = 0d0 
    phiGU_ver_derxx(1,3) = 0d0
    phiGD_ver_derxx(1,3) = 0d0

    phiGR_ver_derxx(2,3) = 0d0 
    phiGL_ver_derxx(2,3) = 0d0 
    phiGU_ver_derxx(2,3) = 0d0
    phiGD_ver_derxx(2,3) = 0d0

    phiGR_ver_derxy(1,3) = 0d0 
    phiGL_ver_derxy(1,3) = 0d0
    phiGU_ver_derxy(1,3) = 0d0
    phiGD_ver_derxy(1,3) = 0d0

    phiGR_ver_derxy(2,3) = 0d0 
    phiGL_ver_derxy(2,3) = 0d0 
    phiGU_ver_derxy(2,3) = 0d0
    phiGD_ver_derxy(2,3) = 0d0

    phiGR_ver_deryy(1,3) = 0d0 
    phiGL_ver_deryy(1,3) = 0d0 
    phiGU_ver_deryy(1,3) = 0d0
    phiGD_ver_deryy(1,3) = 0d0

    phiGR_ver_deryy(2,3) = 0d0 
    phiGL_ver_deryy(2,3) = 0d0 
    phiGU_ver_deryy(2,3) = 0d0
    phiGD_ver_deryy(2,3) = 0d0
 
    phiGR_ver(1,4) =  2d0/3d0
    phiGL_ver(1,4) =  2d0/3d0
    phiGU_ver(1,4) =  2d0/3d0
    phiGD_ver(1,4) =  2d0/3d0

    phiGR_ver(2,4) =  2d0/3d0
    phiGL_ver(2,4) =  2d0/3d0
    phiGU_ver(2,4) =  2d0/3d0
    phiGD_ver(2,4) =  2d0/3d0

    phiGR_ver_derx(1,4) = 2d0
    phiGL_ver_derx(1,4) = -2d0
    phiGU_ver_derx(1,4) = -2d0
    phiGD_ver_derx(1,4) = -2d0

    phiGR_ver_derx(2,4) = 2d0
    phiGL_ver_derx(2,4) = -2d0
    phiGU_ver_derx(2,4) = 2d0
    phiGD_ver_derx(2,4) = 2d0

    phiGR_ver_dery(1,4) = 0d0
    phiGL_ver_dery(1,4) = 0d0
    phiGU_ver_dery(1,4) = 0d0
    phiGD_ver_dery(1,4) = 0d0

    phiGR_ver_dery(2,4) = 0d0
    phiGL_ver_dery(2,4) = 0d0
    phiGU_ver_dery(2,4) = 0d0
    phiGD_ver_dery(2,4) = 0d0

    phiGR_ver_derxx(1,4) = 2d0
    phiGL_ver_derxx(1,4) = 2d0
    phiGU_ver_derxx(1,4) = 2d0
    phiGD_ver_derxx(1,4) = 2d0

    phiGR_ver_derxx(2,4) = 2d0 
    phiGL_ver_derxx(2,4) = 2d0 
    phiGU_ver_derxx(2,4) = 2d0
    phiGD_ver_derxx(2,4) = 2d0

    phiGR_ver_derxy(1,4) = 0d0 
    phiGL_ver_derxy(1,4) = 0d0
    phiGU_ver_derxy(1,4) = 0d0
    phiGD_ver_derxy(1,4) = 0d0

    phiGR_ver_derxy(2,4) = 0d0 
    phiGL_ver_derxy(2,4) = 0d0 
    phiGU_ver_derxy(2,4) = 0d0
    phiGD_ver_derxy(2,4) = 0d0

    phiGR_ver_deryy(1,4) = 0d0 
    phiGL_ver_deryy(1,4) = 0d0
    phiGU_ver_deryy(1,4) = 0d0
    phiGD_ver_deryy(1,4) = 0d0

    phiGR_ver_deryy(2,4) = 0d0 
    phiGL_ver_deryy(2,4) = 0d0 
    phiGU_ver_deryy(2,4) = 0d0
    phiGD_ver_deryy(2,4) = 0d0

    phiGR_ver(1,5) =  1d0
    phiGL_ver(1,5) =  -1d0
    phiGU_ver(1,5) =  -1d0
    phiGD_ver(1,5) =  1d0

    phiGR_ver(2,5) =  -1d0
    phiGL_ver(2,5) =  1d0
    phiGU_ver(2,5) =  1d0
    phiGD_ver(2,5) =  -1d0

    phiGR_ver_derx(1,5) = 1d0
    phiGL_ver_derx(1,5) = 1d0
    phiGU_ver_derx(1,5) = 1d0
    phiGD_ver_derx(1,5) = -1d0

    phiGR_ver_derx(2,5) = -1d0
    phiGL_ver_derx(2,5) = -1d0
    phiGU_ver_derx(2,5) = 1d0
    phiGD_ver_derx(2,5) = -1d0

    phiGR_ver_dery(1,5) = 1d0
    phiGL_ver_dery(1,5) = -1d0
    phiGU_ver_dery(1,5) = -1d0
    phiGD_ver_dery(1,5) = -1d0

    phiGR_ver_dery(2,5) = 1d0
    phiGL_ver_dery(2,5) = -1d0
    phiGU_ver_dery(2,5) = 1d0
    phiGD_ver_dery(2,5) = 1d0

    phiGR_ver_derxx(1,5) = 0d0
    phiGL_ver_derxx(1,5) = 0d0
    phiGU_ver_derxx(1,5) = 0d0
    phiGD_ver_derxx(1,5) = 0d0

    phiGR_ver_derxx(2,5) = 0d0 
    phiGL_ver_derxx(2,5) = 0d0 
    phiGU_ver_derxx(2,5) =0d0
    phiGD_ver_derxx(2,5) = 0d0


    phiGR_ver_derxy(1,5) = 1d0  
    phiGL_ver_derxy(1,5) = 1d0
    phiGU_ver_derxy(1,5) = 1d0
    phiGD_ver_derxy(1,5) = 1d0

    phiGR_ver_derxy(2,5) = 1d0 
    phiGL_ver_derxy(2,5) = 1d0 
    phiGU_ver_derxy(2,5) = 1d0
    phiGD_ver_derxy(2,5) = 1d0

    phiGR_ver_deryy(1,5) = 0d0
    phiGL_ver_deryy(1,5) = 0d0
    phiGU_ver_deryy(1,5) = 0d0
    phiGD_ver_deryy(1,5) = 0d0

    phiGR_ver_deryy(2,5) = 0d0 
    phiGL_ver_deryy(2,5) = 0d0
    phiGU_ver_deryy(2,5) =0d0
    phiGD_ver_deryy(2,5) = 0d0

    
    phiGR_ver(1,6) =   2d0/3d0
    phiGL_ver(1,6) =   2d0/3d0
    phiGU_ver(1,6) =   2d0/3d0
    phiGD_ver(1,6) =   2d0/3d0

    phiGR_ver(2,6) =   2d0/3d0
    phiGL_ver(2,6) =   2d0/3d0
    phiGU_ver(2,6) =   2d0/3d0
    phiGD_ver(2,6) =   2d0/3d0

    phiGR_ver_derx(1,6) = 0d0
    phiGL_ver_derx(1,6) = 0d0
    phiGU_ver_derx(1,6) = 0d0
    phiGD_ver_derx(1,6) = 0d0

    phiGR_ver_derx(2,6) = 0d0
    phiGL_ver_derx(2,6) = 0d0
    phiGU_ver_derx(2,6) = 0d0
    phiGD_ver_derx(2,6) = 0d0

    phiGR_ver_dery(1,6) = 2d0
    phiGL_ver_dery(1,6) = 2d0
    phiGU_ver_dery(1,6) = 2d0
    phiGD_ver_dery(1,6) = -2d0

    phiGR_ver_dery(2,6) = -2d0
    phiGL_ver_dery(2,6) = -2d0
    phiGU_ver_dery(2,6) = 2d0
    phiGD_ver_dery(2,6) = -2d0

    phiGR_ver_derxx(1,6) = 0d0
    phiGL_ver_derxx(1,6) = 0d0
    phiGU_ver_derxx(1,6) = 0d0
    phiGD_ver_derxx(1,6) = 0d0

    phiGR_ver_derxx(2,6) = 0d0 
    phiGL_ver_derxx(2,6) = 0d0 
    phiGU_ver_derxx(2,6) = 0d0
    phiGD_ver_derxx(2,6) = 0d0

    phiGR_ver_derxy(1,6) = 0d0 
    phiGL_ver_derxy(1,6) = 0d0
    phiGU_ver_derxy(1,6) = 0d0
    phiGD_ver_derxy(1,6) = 0d0

    phiGR_ver_derxy(2,6) = 0d0 
    phiGL_ver_derxy(2,6) = 0d0 
    phiGU_ver_derxy(2,6) = 0d0
    phiGD_ver_derxy(2,6) = 0d0

    phiGR_ver_deryy(1,6) = 2d0
    phiGL_ver_deryy(1,6) = 2d0
    phiGU_ver_deryy(1,6) = 2d0
    phiGD_ver_deryy(1,6) = 2d0

    phiGR_ver_deryy(2,6) = 2d0 
    phiGL_ver_deryy(2,6) = 2d0
    phiGU_ver_deryy(2,6) = 2d0
    phiGD_ver_deryy(2,6) = 2d0
 
    phiGR_ver(1,7) =   0.4d0 
    phiGL_ver(1,7) =   -0.4d0
    phiGU_ver(1,7) =   -0.4d0
    phiGD_ver(1,7) =   -0.4d0

    phiGR_ver(2,7) =   0.4d0
    phiGL_ver(2,7) =   -0.4d0
    phiGU_ver(2,7) =   0.4d0
    phiGD_ver(2,7) =   0.4d0

    phiGR_ver_derx(1,7) = 2.4d0
    phiGL_ver_derx(1,7) = 2.4d0
    phiGU_ver_derx(1,7) = 2.4d0
    phiGD_ver_derx(1,7) = 2.4d0

    phiGR_ver_derx(2,7) = 2.4d0
    phiGL_ver_derx(2,7) = 2.4d0
    phiGU_ver_derx(2,7) = 2.4d0
    phiGD_ver_derx(2,7) = 2.4d0

    phiGR_ver_dery(1,7) = 0d0
    phiGL_ver_dery(1,7) = 0d0
    phiGU_ver_dery(1,7) = 0d0
    phiGD_ver_dery(1,7) = 0d0

    phiGR_ver_dery(2,7) = 0d0
    phiGL_ver_dery(2,7) = 0d0
    phiGU_ver_dery(2,7) = 0d0
    phiGD_ver_dery(2,7) = 0d0

    phiGR_ver_derxx(1,7) = 6d0
    phiGL_ver_derxx(1,7) = -6d0
    phiGU_ver_derxx(1,7) = -6d0
    phiGD_ver_derxx(1,7) = -6d0

    phiGR_ver_derxx(2,7) = 6d0 
    phiGL_ver_derxx(2,7) = -6d0 
    phiGU_ver_derxx(2,7) = 6d0
    phiGD_ver_derxx(2,7) = 6d0

    phiGR_ver_derxy(1,7) = 0d0 
    phiGL_ver_derxy(1,7) = 0d0
    phiGU_ver_derxy(1,7) = 0d0
    phiGD_ver_derxy(1,7) = 0d0

    phiGR_ver_derxy(2,7) = 0d0 
    phiGL_ver_derxy(2,7) = 0d0 
    phiGU_ver_derxy(2,7) = 0d0
    phiGD_ver_derxy(2,7) = 0d0

    phiGR_ver_deryy(1,7) = 0d0
    phiGL_ver_deryy(1,7) = 0d0
    phiGU_ver_deryy(1,7) = 0d0
    phiGD_ver_deryy(1,7) = 0d0

    phiGR_ver_deryy(2,7) = 0d0 
    phiGL_ver_deryy(2,7) = 0d0
    phiGU_ver_deryy(2,7) = 0d0
    phiGD_ver_deryy(2,7) = 0d0

    phiGR_ver_derxxx(1,7) = 6d0  
    phiGL_ver_derxxx(1,7) = 6d0 
    phiGU_ver_derxxx(1,7) = 6d0
    phiGD_ver_derxxx(1,7) = 6d0

    phiGR_ver_derxxx(2,7) = 6d0  
    phiGL_ver_derxxx(2,7) = 6d0 
    phiGU_ver_derxxx(2,7) = 6d0
    phiGD_ver_derxxx(2,7) = 6d0

    phiGR_ver(1,8) =   2d0/3d0 
    phiGL_ver(1,8) =   2d0/3d0
    phiGU_ver(1,8) =   2d0/3d0
    phiGD_ver(1,8) =   -2d0/3d0

    phiGR_ver(2,8) =   -2d0/3d0
    phiGL_ver(2,8) =   -2d0/3d0
    phiGU_ver(2,8) =   2d0/3d0
    phiGD_ver(2,8) =   -2d0/3d0

    phiGR_ver_derx(1,8) = 2d0
    phiGL_ver_derx(1,8) = -2d0
    phiGU_ver_derx(1,8) = -2d0
    phiGD_ver_derx(1,8) = 2d0

    phiGR_ver_derx(2,8) = -2d0
    phiGL_ver_derx(2,8) = 2d0
    phiGU_ver_derx(2,8) = 2d0
    phiGD_ver_derx(2,8) = -2d0

    phiGR_ver_dery(1,8) = 2d0/3d0
    phiGL_ver_dery(1,8) = 2d0/3d0
    phiGU_ver_dery(1,8) = 2d0/3d0
    phiGD_ver_dery(1,8) = 2d0/3d0

    phiGR_ver_dery(2,8) = 2d0/3d0
    phiGL_ver_dery(2,8) = 2d0/3d0
    phiGU_ver_dery(2,8) = 2d0/3d0
    phiGD_ver_dery(2,8) = 2d0/3d0

    phiGR_ver_derxx(1,8) = 2d0  
    phiGL_ver_derxx(1,8) = 2d0
    phiGU_ver_derxx(1,8) = 2d0
    phiGD_ver_derxx(1,8) = -2d0

    phiGR_ver_derxx(2,8) = -2d0 
    phiGL_ver_derxx(2,8) = -2d0 
    phiGU_ver_derxx(2,8) = 2d0
    phiGD_ver_derxx(2,8) = -2d0

    phiGR_ver_derxy(1,8) = 2d0  
    phiGL_ver_derxy(1,8) = -2d0
    phiGU_ver_derxy(1,8) = -2d0
    phiGD_ver_derxy(1,8) = -2d0

    phiGR_ver_derxy(2,8) = 2d0 
    phiGL_ver_derxy(2,8) = -2d0 
    phiGU_ver_derxy(2,8) = 2d0
    phiGD_ver_derxy(2,8) = 2d0

    phiGR_ver_deryy(1,8) = 0d0!=0
    phiGL_ver_deryy(1,8) = 0d0
    phiGU_ver_deryy(1,8) = 0d0
    phiGD_ver_deryy(1,8) = 0d0

    phiGR_ver_deryy(2,8) = 0d0 
    phiGL_ver_deryy(2,8) = 0d0
    phiGU_ver_deryy(2,8) = 0d0
    phiGD_ver_deryy(2,8) = 0d0

    phiGR_ver_derxxy(1,8) = 2d0 
    phiGL_ver_derxxy(1,8) = 2d0 
    phiGU_ver_derxxy(1,8) = 2d0
    phiGD_ver_derxxy(1,8) = 2d0

    phiGR_ver_derxxy(2,8) = 2d0  
    phiGL_ver_derxxy(2,8) = 2d0 
    phiGU_ver_derxxy(2,8) = 2d0
    phiGD_ver_derxxy(2,8) = 2d0
    phiGR_ver(1,9) =   2d0/3d0 
    phiGL_ver(1,9) =   -2d0/3d0
    phiGU_ver(1,9) =   -2d0/3d0
    phiGD_ver(1,9) =   -2d0/3d0

    phiGR_ver(2,9) =   2d0/3d0
    phiGL_ver(2,9) =   -2d0/3d0
    phiGU_ver(2,9) =   2d0/3d0
    phiGD_ver(2,9) =   2d0/3d0

    phiGR_ver_derx(1,9) =  2d0/3d0
    phiGL_ver_derx(1,9) =  2d0/3d0
    phiGU_ver_derx(1,9) =  2d0/3d0
    phiGD_ver_derx(1,9) =  2d0/3d0

    phiGR_ver_derx(2,9) =  2d0/3d0
    phiGL_ver_derx(2,9) =  2d0/3d0
    phiGU_ver_derx(2,9) =  2d0/3d0
    phiGD_ver_derx(2,9) =  2d0/3d0

    phiGR_ver_dery(1,9) = 2d0 
    phiGL_ver_dery(1,9) = -2d0
    phiGU_ver_dery(1,9) = -2d0
    phiGD_ver_dery(1,9) = 2d0

    phiGR_ver_dery(2,9) = -2d0
    phiGL_ver_dery(2,9) = 2d0
    phiGU_ver_dery(2,9) = 2d0
    phiGD_ver_dery(2,9) = -2d0

    phiGR_ver_derxx(1,9) = 0d0  
    phiGL_ver_derxx(1,9) = 0d0
    phiGU_ver_derxx(1,9) = 0d0
    phiGD_ver_derxx(1,9) = 0d0

    phiGR_ver_derxx(2,9) = 0d0 
    phiGL_ver_derxx(2,9) = 0d0 
    phiGU_ver_derxx(2,9) = 0d0
    phiGD_ver_derxx(2,9) = 0d0

    phiGR_ver_derxy(1,9) = 2d0  
    phiGL_ver_derxy(1,9) = 2d0
    phiGU_ver_derxy(1,9) = 2d0
    phiGD_ver_derxy(1,9) = -2d0

    phiGR_ver_derxy(2,9) = -2d0 
    phiGL_ver_derxy(2,9) = -2d0 
    phiGU_ver_derxy(2,9) = 2d0
    phiGD_ver_derxy(2,9) = -2d0

    phiGR_ver_deryy(1,9) = 2d0 
    phiGL_ver_deryy(1,9) = -2d0
    phiGU_ver_deryy(1,9) = -2d0
    phiGD_ver_deryy(1,9) = -2d0

    phiGR_ver_deryy(2,9) = 2d0 
    phiGL_ver_deryy(2,9) = -2d0
    phiGU_ver_deryy(2,9) = 2d0
    phiGD_ver_deryy(2,9) = 2d0

    phiGR_ver_derxyy(1,9) = 2d0  
    phiGL_ver_derxyy(1,9) = 2d0 
    phiGU_ver_derxyy(1,9) = 2d0
    phiGD_ver_derxyy(1,9) = 2d0

    phiGR_ver_derxyy(2,9) = 2d0  
    phiGL_ver_derxyy(2,9) = 2d0 
    phiGU_ver_derxyy(2,9) = 2d0
    phiGD_ver_derxyy(2,9) = 2d0

    
    phiGR_ver(1,10) =   0.4d0 
    phiGL_ver(1,10) =   0.4d0
    phiGU_ver(1,10) =   0.4d0
    phiGD_ver(1,10) =   -0.4d0

    phiGR_ver(2,10) =   -0.4d0
    phiGL_ver(2,10) =   -0.4d0
    phiGU_ver(2,10) =   0.4d0
    phiGD_ver(2,10) =   -0.4d0

    phiGR_ver_derx(1,10) = 0d0
    phiGL_ver_derx(1,10) = 0d0
    phiGU_ver_derx(1,10) = 0d0
    phiGD_ver_derx(1,10) = 0d0

    phiGR_ver_derx(2,10) = 0d0
    phiGL_ver_derx(2,10) = 0d0
    phiGU_ver_derx(2,10) = 0d0
    phiGD_ver_derx(2,10) = 0d0

    phiGR_ver_dery(1,10) = 2.4d0
    phiGL_ver_dery(1,10) = 2.4d0
    phiGU_ver_dery(1,10) = 2.4d0
    phiGD_ver_dery(1,10) = 2.4d0

    phiGR_ver_dery(2,10) = 2.4d0
    phiGL_ver_dery(2,10) = 2.4d0
    phiGU_ver_dery(2,10) = 2.4d0
    phiGD_ver_dery(2,10) = 2.4d0

    phiGR_ver_derxx(1,10) = 0d0
    phiGL_ver_derxx(1,10) = 0d0
    phiGU_ver_derxx(1,10) = 0d0
    phiGD_ver_derxx(1,10) = 0d0

    phiGR_ver_derxx(2,10) = 0d0 
    phiGL_ver_derxx(2,10) = 0d0 
    phiGU_ver_derxx(2,10) = 0d0
    phiGD_ver_derxx(2,10) = 0d0

    phiGR_ver_derxy(1,10) = 0d0 
    phiGL_ver_derxy(1,10) = 0d0
    phiGU_ver_derxy(1,10) = 0d0
    phiGD_ver_derxy(1,10) = 0d0

    phiGR_ver_derxy(2,10) = 0d0 
    phiGL_ver_derxy(2,10) = 0d0 
    phiGU_ver_derxy(2,10) = 0d0
    phiGD_ver_derxy(2,10) = 0d0

    phiGR_ver_deryy(1,10) = 6d0
    phiGL_ver_deryy(1,10) = 6d0
    phiGU_ver_deryy(1,10) = 6d0
    phiGD_ver_deryy(1,10) = -6d0

    phiGR_ver_deryy(2,10) = -6d0 
    phiGL_ver_deryy(2,10) = -6d0
    phiGU_ver_deryy(2,10) = 6d0
    phiGD_ver_deryy(2,10) = -6d0

    phiGR_ver_deryyy(1,10) = 6d0  
    phiGL_ver_deryyy(1,10) = 6d0 
    phiGU_ver_deryyy(1,10) = 6d0
    phiGD_ver_deryyy(1,10) = 6d0

    phiGR_ver_deryyy(2,10) = 6d0  
    phiGL_ver_deryyy(2,10) = 6d0 
    phiGU_ver_deryyy(2,10) = 6d0
    phiGD_ver_deryyy(2,10) = 6d0




    end subroutine get_basis

    !*****************************************************************************************************

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

     
    subroutine RK3

    use com

    CFL = 0.10
    t = 0
    
   
    call calculate_umax

    if (myid1 == 1) then
        open(unit = 12,file = 'Latest_result.txt')
        print *,t,umax
        write(12,*) t,umax
    end if

     
    do while (t < tend)

    call calculate_dt

    tRK = t

    if (t + dt > tend) then
        dt = tend - t
        t = tend
    else
        t = t + dt
    end if

    ! Stage I
    call Lh

    uh00 = uh

    uI = uh + dt*du

    uh = uI
    call jumpfilter
    ! Stage II
    tRK = tRK + dt
    call Lh

    uII = (3d0/4d0)*uh00 + (1d0/4d0)*uh + (1d0/4d0)*dt*du

    uh = uII

    call jumpfilter
         
    ! Stage III
    tRK = tRK - 0.5*dt
    call Lh

    uh = (1d0/3d0)*uh00 + (2d0/3d0)*uh + (2d0/3d0)*dt*du
 
    call jumpfilter
          
    call calculate_umax

    if (myid1 == 1) then
!        print *,t,umax
    write(12,*) t,umax

    end if

    end do

    end subroutine RK3


    !*****************************************************************************************************

    subroutine RK4

    use com

    CFL = 0.75
    t = 0
    
    
    call calculate_umax

    if (myid1 == 1) then
        open(unit = 12,file = 'Latest_result.txt')
        print *,t,umax
        write(12,*) t,umax
    end if

    do while (t < tend)

    call calculate_dt
    tRK = t

    if (t + dt > tend) then
        dt = tend - t
        t = tend
    else
        t = t + dt
    end if
    uI = uh
    uII = uh
    do i = 1,5
    call Lh
    uI = uh + (dt/6d0)*du
    tRK = tRK + (dt/6d0)
    uh = uI
    call jumpfilter
    end do
    uII = 0.04d0*uII + 0.36d0*uI
    uI = 15*uII - 5*uI
    uh = uI
    tRK = tRK - 0.5*dt

    do i = 6,9
    call Lh
    uI = uh + (dt/6d0)*du
    tRK = tRK + dt/6d0
    uh = uI
    call jumpfilter
    end do
    call Lh
    uh = uII + 0.6d0*uI + (dt/10d0)*du

    call jumpfilter
    call calculate_umax

    if (myid1 == 1) then
        !print *,t,umax,sum(Is_trouble_cell)
       ! print *,t,umax,dt
        write(12,*) t,umax
    end if

    end do

    end subroutine RK4

    !*****************************************************************************************************

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

    !*****************************************************************************************************

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

    !*****************************************************************************************************

    subroutine eigenvalueMm(Amax,Amin,rho,rhou,rhov,E,n1,n2)

    use com

    real(8) u,v,w,p,c,BP,Bn,un,Amax,Amin

    u = rhou/rho
    v = rhov/rho

    un = u*n1 + v*n2

    p = gamma1*(E - 0.5d0*rho*(u**2 + v**2))

    c = sqrt(abs(gamma*p/rho))

    Amax = un + c
    Amin = un - c

    end subroutine eigenvalueMm

    !*****************************************************************************************************

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

    !*****************************************************************************************************

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

    !*****************************************************************************************************

    subroutine HLL_Flux

    use com

    if (SR < 0) then
        Fhat1 = FR1
    else if (SL > 0) then
        Fhat1 = FL1
    else
        Fhat1 = ( SR*FL1 - SL*FR1 + SL*SR*(UR1 - UL1) )/(SR - SL)
    end if

    end subroutine HLL_Flux

    !*****************************************************************************************************

    subroutine LF_Flux

    use com

    Fhat1 = 0.5d0*(FR1 + FL1 - max(abs(SR),abs(SL))*(UR1 - UL1))

    end subroutine LF_Flux

     
    !*****************************************************************************************************

    subroutine compute_Rinv(Rmat,Rinv,rho,u,v,E,n1,n2)

    real nf(2),n1,n2
    real Rmat(4,4)
    real Rinv(4,4)
    real c,rho,rhou,rhov,u,v,pr,eH,ek,unf,magnorm

    gam = 1.4

    nf(1) = n1
    nf(2) = n2

    magnorm = dsqrt(nf(1)**2+nf(2)**2)
    nf(1) = nf(1)/magnorm
    nf(2) = nf(2)/magnorm

    pr  = (E-0.5d0*rho*(u**2+v**2))*(gam-1)
    c   = dsqrt(gam*pr/rho) ! speed of sound
    eH  = (E+pr)/rho ! specific enthalpy
    unf = u*nf(1)+v*nf(2)
    ek = 0.5d0*(u**2+v**2)

    Rmat(1,1) = 1.d0
    Rmat(2,1) = u-c*nf(1)
    Rmat(3,1) = v-c*nf(2)
    Rmat(4,1) = eH-c*unf

    Rmat(1,2) = 1.d0
    Rmat(2,2) = u
    Rmat(3,2) = v
    Rmat(4,2) = ek

    Rmat(1,3) = 1.d0
    Rmat(2,3) = u+c*nf(1)
    Rmat(3,3) = v+c*nf(2)
    Rmat(4,3) = eH+c*unf

    Rmat(1,4) = 0.d0
    Rmat(2,4) = nf(2)
    Rmat(3,4) = -nf(1)
    Rmat(4,4) = u*nf(2)-v*nf(1)

    ! ===========================
    Rinv(1,1) = ((gam-1)*ek+c*unf)*0.5d0/(c**2.d0)
    Rinv(2,1) = (c**2-(gam-1)*ek)/(c**2.d0)
    Rinv(3,1) = ((gam-1)*ek-c*unf)*0.5d0/(c**2.d0)
    Rinv(4,1) = v*nf(1)-u*nf(2)

    Rinv(1,2) = ((1-gam)*u-c*nf(1))*0.5d0/(c**2.d0)
    Rinv(2,2) = (gam-1)*u/(c**2.d0)
    Rinv(3,2) = ((1-gam)*u+c*nf(1))*0.5d0/(c**2.d0)
    Rinv(4,2) = nf(2)

    Rinv(1,3) = ((1-gam)*v-c*nf(2))*0.5d0/(c**2.d0)
    Rinv(2,3) = (gam-1)*v/(c**2.d0)
    Rinv(3,3) = ((1-gam)*v+c*nf(2))*0.5d0/(c**2.d0)
    Rinv(4,3) = -nf(1)

    Rinv(1,4) = (gam-1)*0.5d0/(c**2.d0)
    Rinv(2,4) = (1-gam)/(c**2.d0)
    Rinv(3,4) = (gam-1)*0.5d0/(c**2.d0)
    Rinv(4,4) = 0.d0
 

    end subroutine compute_Rinv    

    !*****************************************************************************************************

    subroutine minmod

    use com

    if (direction == 1) then
        hd = hx
    else if (direction == 2) then
        hd = hy
    end if

    do i = 1,NumEq
        if (abs(DeltaU(i,1)) <= M*hd**2) then
            DeltaUmod(i,1) = DeltaU(i,1)
        else
            a = sign(1d0,DeltaU(i,1))
            b = sign(1d0,DeltaUR(i,1))
            c = sign(1d0,DeltaUL(i,1))
            s = (a + b + c)/3d0
            if (abs(s) == 1) then
                DeltaUmod(i,1) = s*min(abs(DeltaU(i,1)),beta*abs(DeltaUR(i,1)),beta*abs(DeltaUL(i,1)))
            else
                DeltaUmod(i,1) = 0
            end if
        end if

    end do

    end subroutine minmod

  
    subroutine norm(x,d)

    real x(4),d

    d = 0

    do i = 1,4
        d = d + x(i)**2
    end do

    d = d**0.5

    end subroutine norm


    subroutine calculate_tq(ubar,uq,tq,gamma)

    real ubar(4),uq(4),tq,ta,tb,ut(4),gamma
    integer count

    ta = 0
    tb = 1
    count = 0

    do while (tb - ta > 1e-14)
        tq = 0.5*(ta + tb)
        ut = tq*uq + (1 - tq)*ubar
        if (pressure(ut(1),ut(2),ut(3),ut(4),gamma) < 1e-13) then
            !ta = ta
            tb = tq
        else
            ta = tq
            !tb = tb
        end if
    end do

    tq = ta
    ut = tq*uq + (1 - tq)*ubar
    !print *,pressure(ut(1),ut(2),ut(3),ut(4),ut(5),ut(6),ut(7),ut(8),gamma),ta,tb

    end subroutine calculate_tq

    !*****************************************************************************************************

    function pressure(rho,rhou,rhov,E,gamma)

    real(8) rho,rhou,rhov,rhow,E,B1,B2,B3,gamma
    real(8) pressure

    pressure = (gamma - 1)*(E - 0.5*(rhou**2 + rhov**2)/rho)

    end

    !*****************************************************************************************************

    subroutine evenex_y(a,b)

    real a(10),b(10)

    a(1) = b(1)

    !a(2) = b(2)
    !a(3) = -b(3)
    !
    !a(4) = b(4)
    !a(5) = -b(5)
    !a(6) = b(6)
    !
    !a(7) = b(7)
    !a(8) = -b(8)
    !a(9) = b(9)
    !a(10) = -b(10)
    a(2) = 0d0
    a(3) = 0d0
    
    a(4) = 0d0
    a(5) = 0d0
    a(6) = 0d0
    
    a(7) = 0d0
    a(8) = 0d0
    a(9) = 0d0
    a(10) = 0d0
    end subroutine evenex_y

    !*****************************************************************************************************

    subroutine oddex_y(a,b)

    real a(10),b(10)

    a(1) = -b(1)

    !a(2) = -b(2)
    !a(3) = b(3)
    !
    !a(4) = -b(4)
    !a(5) = b(5)
    !a(6) = -b(6)
    !
    !a(7) = -b(7)
    !a(8) = b(8)
    !a(9) = -b(9)
    !a(10) = b(10)
    a(2) = 0d0
    a(3) = 0d0
    
    a(4) = 0d0
    a(5) = 0d0
    a(6) = 0d0
    
    a(7) = 0d0
    a(8) = 0d0
    a(9) = 0d0
    a(10) = 0d0
    end subroutine oddex_y
 
     

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
    
    subroutine  jumpfilter

    use com

    real(kind=8) rhoLinfty,rhouLinfty,rhovLinfty,EnerLinfty,rhoij,label_density(Nx,Ny),rhouij,label_rhou(Nx,Ny),rhovij,label_rhov(Nx,Ny),Enerij,label_Ener(Nx,Ny)
    real(kind=8) rhoLinfty1,rhouLinfty1,rhovLinfty1,EnerLinfty1
    real(kind=8)  delta0max , delta1max , delta2max , delta3max ,damping
    real(kind=8)  damping1,damping2,damping3 
    real(kind=8)  deltadensity1,deltarhou1,deltarhov1,deltaEner1
    real(kind=8)  deltadensity2,deltarhou2,deltarhov2,deltaEner2 
    real(kind=8)  deltadensity3,deltarhou3,deltarhov3,deltaEner3

    integer ploydeg ,ss,m_oe 
    real(kind=8) scal,u1ij,u2ij,Eij,pressureij,cij,betai,betaj,entropyij,enthayij  
    real(kind=8) left_edge_jump_density ,right_edge_jump_density,bottom_edge_jump_density,top_edge_jump_density, deltadensity0 
    real(kind=8) left_edge_jump_rhou,right_edge_jump_rhou,bottom_edge_jump_rhou,top_edge_jump_rhou,deltarhou0
    real(kind=8) left_edge_jump_rhov,right_edge_jump_rhov,bottom_edge_jump_rhov,top_edge_jump_rhov,deltarhov0            
    real(kind=8) left_edge_jump_Ener,right_edge_jump_Ener,bottom_edge_jump_Ener,top_edge_jump_Ener,deltaEner0  
 
    real(kind=8)  left_edge_jumpderx_density, right_edge_jumpderx_density, bottom_edge_jumpderx_density, top_edge_jumpderx_density 
    real(kind=8)  left_edge_jumpdery_density, right_edge_jumpdery_density, bottom_edge_jumpdery_density, top_edge_jumpdery_density 
    real(kind=8) left_edge_jumpderx_rhou,right_edge_jumpderx_rhou,bottom_edge_jumpderx_rhou, top_edge_jumpderx_rhou 
    real(kind=8) left_edge_jumpdery_rhou,right_edge_jumpdery_rhou,bottom_edge_jumpdery_rhou,top_edge_jumpdery_rhou 
    real(kind=8) left_edge_jumpderx_rhov,right_edge_jumpderx_rhov,bottom_edge_jumpderx_rhov, top_edge_jumpderx_rhov 
    real(kind=8) left_edge_jumpdery_rhov,right_edge_jumpdery_rhov,bottom_edge_jumpdery_rhov,top_edge_jumpdery_rhov       
    real(kind=8) left_edge_jumpderx_Ener,right_edge_jumpderx_Ener,bottom_edge_jumpderx_Ener, top_edge_jumpderx_Ener 
    real(kind=8) left_edge_jumpdery_Ener,right_edge_jumpdery_Ener,bottom_edge_jumpdery_Ener,top_edge_jumpdery_Ener     
    real(kind=8) left_edge_jumpderxx_density,right_edge_jumpderxx_density,bottom_edge_jumpderxx_density,top_edge_jumpderxx_density 
    real(kind=8) left_edge_jumpderxy_density,right_edge_jumpderxy_density,bottom_edge_jumpderxy_density,top_edge_jumpderxy_density 
    real(kind=8) left_edge_jumpderyy_density,right_edge_jumpderyy_density,bottom_edge_jumpderyy_density,top_edge_jumpderyy_density 
    real(kind=8) left_edge_jumpderxx_rhou,right_edge_jumpderxx_rhou,bottom_edge_jumpderxx_rhou,top_edge_jumpderxx_rhou 
    real(kind=8) left_edge_jumpderxy_rhou,right_edge_jumpderxy_rhou,bottom_edge_jumpderxy_rhou,top_edge_jumpderxy_rhou
    real(kind=8) left_edge_jumpderyy_rhou,right_edge_jumpderyy_rhou,bottom_edge_jumpderyy_rhou,top_edge_jumpderyy_rhou
    real(kind=8) left_edge_jumpderxx_rhov,right_edge_jumpderxx_rhov,bottom_edge_jumpderxx_rhov,top_edge_jumpderxx_rhov
    real(kind=8) left_edge_jumpderxy_rhov,right_edge_jumpderxy_rhov,bottom_edge_jumpderxy_rhov,top_edge_jumpderxy_rhov 
    real(kind=8) left_edge_jumpderyy_rhov,right_edge_jumpderyy_rhov,bottom_edge_jumpderyy_rhov,top_edge_jumpderyy_rhov
    real(kind=8) left_edge_jumpderxx_Ener,right_edge_jumpderxx_Ener,bottom_edge_jumpderxx_Ener,top_edge_jumpderxx_Ener
    real(kind=8) left_edge_jumpderxy_Ener,right_edge_jumpderxy_Ener,bottom_edge_jumpderxy_Ener,top_edge_jumpderxy_Ener
    real(kind=8) left_edge_jumpderyy_Ener,right_edge_jumpderyy_Ener,bottom_edge_jumpderyy_Ener,top_edge_jumpderyy_Ener
    real(kind=8) left_edge_jumpderxxx_density,right_edge_jumpderxxx_density,bottom_edge_jumpderxxx_density,top_edge_jumpderxxx_density 
    real(kind=8) left_edge_jumpderxxy_density,right_edge_jumpderxxy_density,bottom_edge_jumpderxxy_density,top_edge_jumpderxxy_density 
    real(kind=8) left_edge_jumpderxyy_density,right_edge_jumpderxyy_density,bottom_edge_jumpderxyy_density,top_edge_jumpderxyy_density 
    real(kind=8) left_edge_jumpderyyy_density,right_edge_jumpderyyy_density,bottom_edge_jumpderyyy_density,top_edge_jumpderyyy_density 
    real(kind=8) left_edge_jumpderxxx_rhou,right_edge_jumpderxxx_rhou,bottom_edge_jumpderxxx_rhou,top_edge_jumpderxxx_rhou 
    real(kind=8) left_edge_jumpderxxy_rhou,right_edge_jumpderxxy_rhou,bottom_edge_jumpderxxy_rhou,top_edge_jumpderxxy_rhou
    real(kind=8) left_edge_jumpderxyy_rhou,right_edge_jumpderxyy_rhou,bottom_edge_jumpderxyy_rhou,top_edge_jumpderxyy_rhou
    real(kind=8) left_edge_jumpderyyy_rhou,right_edge_jumpderyyy_rhou,bottom_edge_jumpderyyy_rhou,top_edge_jumpderyyy_rhou
    real(kind=8) left_edge_jumpderxxx_rhov,right_edge_jumpderxxx_rhov,bottom_edge_jumpderxxx_rhov,top_edge_jumpderxxx_rhov
    real(kind=8) left_edge_jumpderxxy_rhov,right_edge_jumpderxxy_rhov,bottom_edge_jumpderxxy_rhov,top_edge_jumpderxxy_rhov 
    real(kind=8) left_edge_jumpderxyy_rhov,right_edge_jumpderxyy_rhov,bottom_edge_jumpderxyy_rhov,top_edge_jumpderxyy_rhov
    real(kind=8) left_edge_jumpderyyy_rhov,right_edge_jumpderyyy_rhov,bottom_edge_jumpderyyy_rhov,top_edge_jumpderyyy_rhov
    real(kind=8) left_edge_jumpderxxx_Ener,right_edge_jumpderxxx_Ener,bottom_edge_jumpderxxx_Ener,top_edge_jumpderxxx_Ener
    real(kind=8) left_edge_jumpderxxy_Ener,right_edge_jumpderxxy_Ener,bottom_edge_jumpderxxy_Ener,top_edge_jumpderxxy_Ener
    real(kind=8) left_edge_jumpderxyy_Ener,right_edge_jumpderxyy_Ener,bottom_edge_jumpderxyy_Ener,top_edge_jumpderxyy_Ener
    real(kind=8) left_edge_jumpderyyy_Ener,right_edge_jumpderyyy_Ener,bottom_edge_jumpderyyy_Ener,top_edge_jumpderyyy_Ener
    real(kind=8) rhofor_x,rhoback_x,rhofor_y,rhoback_y,Linftyrhox,Linftyrhoy
    real(kind=8) rhoufor_x,rhouback_x,rhoufor_y,rhouback_y,Linftyrhoux,Linftyrhouy
    real(kind=8) rhovfor_x,rhovback_x,rhovfor_y,rhovback_y,Linftyrhovx,Linftyrhovy
    real(kind=8) Enerfor_x,Enerback_x,Enerfor_y,Enerback_y,LinftyEnerx,LinftyEnery
    real(kind=8) rhomax,rhomin,rhoumax,rhoumin,rhovmax,rhovmin,Enermax,Enermin
    real(kind=8) Machij
    real(kind=8) uhmod(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(kind=8) dampingfprint
    dampingfprint = 10d0
    call set_bc
    uhmod = uh
    ! The x-direction
    UR_ver = 0d0
    UL_ver = 0d0
    UR_ver_derx = 0d0
    UL_ver_derx = 0d0
    UR_ver_dery = 0d0
    UL_ver_dery = 0d0

    UR_ver_derxx = 0d0
    UL_ver_derxx = 0d0
    UR_ver_derxy = 0d0
    UL_ver_derxy = 0d0
    UR_ver_deryy = 0d0
    UL_ver_deryy = 0d0

    
    UR_ver_derxxx = 0d0
    UL_ver_derxxx = 0d0
    UR_ver_derxxy = 0d0
    UL_ver_derxxy = 0d0
    UR_ver_derxyy = 0d0
    UL_ver_derxyy = 0d0
    UR_ver_deryyy = 0d0
    UL_ver_deryyy = 0d0
    do n = 1,NumEq  
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR_ver(i,j,k,:,n) = UR_ver(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver(:,d)
                        UL_ver(i + 1,j,k,:,n) = UL_ver(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver(:,d)
                       
                        UR_ver_derx(i,j,k,:,n) = UR_ver_derx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derx(:,d)
                        UL_ver_derx(i + 1,j,k,:,n) = UL_ver_derx(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derx(:,d)

                        UR_ver_dery(i,j,k,:,n) = UR_ver_dery(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_dery(:,d)
                        UL_ver_dery(i + 1,j,k,:,n) = UL_ver_dery(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_dery(:,d)

                        UR_ver_derxx(i,j,k,:,n) = UR_ver_derxx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derxx(:,d)
                        UL_ver_derxx(i + 1,j,k,:,n) = UL_ver_derxx(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derxx(:,d)

                        UR_ver_derxy(i,j,k,:,n) = UR_ver_derxy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derxy(:,d)
                        UL_ver_derxy(i + 1,j,k,:,n) = UL_ver_derxy(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derxy(:,d)

                        UR_ver_deryy(i,j,k,:,n) = UR_ver_deryy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_deryy(:,d)
                        UL_ver_deryy(i + 1,j,k,:,n) = UL_ver_deryy(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_deryy(:,d)
                        
                        UR_ver_derxxx(i,j,k,:,n) = UR_ver_derxxx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derxxx(:,d)
                        UL_ver_derxxx(i + 1,j,k,:,n) = UL_ver_derxxx(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derxxx(:,d)

                        UR_ver_derxxy(i,j,k,:,n) = UR_ver_derxxy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derxxy(:,d)
                        UL_ver_derxxy(i + 1,j,k,:,n) = UL_ver_derxxy(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derxxy(:,d)

                        UR_ver_derxyy(i,j,k,:,n) = UR_ver_derxyy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derxyy(:,d)
                        UL_ver_derxyy(i + 1,j,k,:,n) = UL_ver_derxyy(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derxyy(:,d)

                        UR_ver_deryyy(i,j,k,:,n) = UR_ver_deryyy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_deryyy(:,d)
                        UL_ver_deryyy(i + 1,j,k,:,n) = UL_ver_deryyy(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_deryyy(:,d)
                    end do
                end do
            end do
        end do
    end do

    UR_ver_derx =  UR_ver_derx*2/hx
    UL_ver_derx =  UL_ver_derx*2/hx
    UR_ver_dery =  UR_ver_dery*2/hy
    UL_ver_dery =  UL_ver_dery*2/hy

    UR_ver_derxx =  UR_ver_derxx*2/hx*2/hx
    UL_ver_derxx =  UL_ver_derxx*2/hx*2/hx
    UR_ver_derxy =  UR_ver_derxy*2/hx*2/hy
    UL_ver_derxy =  UL_ver_derxy*2/hx*2/hy
    UR_ver_deryy =  UR_ver_deryy*2/hy*2/hy
    UL_ver_deryy =  UL_ver_deryy*2/hy*2/hy

    UR_ver_derxxx =  UR_ver_derxxx*2/hx*2/hx*2/hx
    UL_ver_derxxx =  UL_ver_derxxx*2/hx*2/hx*2/hx
    UR_ver_derxxy =  UR_ver_derxxy*2/hx*2/hy*2/hx
    UL_ver_derxxy =  UL_ver_derxxy*2/hx*2/hy*2/hx
    UR_ver_derxyy =  UR_ver_derxyy*2/hx*2/hy*2/hy
    UL_ver_derxyy =  UL_ver_derxyy*2/hx*2/hy*2/hy
    UR_ver_deryyy =  UR_ver_deryyy*2/hy*2/hy*2/hy
    UL_ver_deryyy =  UL_ver_deryyy*2/hy*2/hy*2/hy

    ! The y-direction
    UU_ver = 0d0
    UD_ver = 0d0
    UU_ver_derx = 0d0
    UD_ver_derx = 0d0
    UU_ver_dery = 0d0
    UD_ver_dery = 0d0

    UU_ver_derxx = 0d0
    UD_ver_derxx = 0d0
    UU_ver_derxy = 0d0
    UD_ver_derxy = 0d0
    UU_ver_deryy = 0d0
    UD_ver_deryy = 0d0

    UU_ver_derxxx = 0d0
    UD_ver_derxxx = 0d0
    UU_ver_derxxy = 0d0
    UD_ver_derxxy = 0d0
    UU_ver_derxyy = 0d0
    UD_ver_derxyy = 0d0
    UU_ver_deryyy = 0d0
    UD_ver_deryyy = 0d0
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 0,Ny
                    do i = 1,Nx
                        UU_ver(i,j,k,:,n) = UU_ver(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver(:,d)
                        UD_ver(i,j + 1,k,:,n) = UD_ver(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver(:,d)

                        UU_ver_derx(i,j,k,:,n) = UU_ver_derx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derx(:,d)
                        UD_ver_derx(i,j + 1,k,:,n) = UD_ver_derx(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derx(:,d)

                        UU_ver_dery(i,j,k,:,n) = UU_ver_dery(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_dery(:,d)
                        UD_ver_dery(i,j + 1,k,:,n) = UD_ver_dery(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_dery(:,d)

                        UU_ver_derxx(i,j,k,:,n) = UU_ver_derxx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derxx(:,d)
                        UD_ver_derxx(i,j + 1,k,:,n) = UD_ver_derxx(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derxx(:,d)

                        UU_ver_derxy(i,j,k,:,n) = UU_ver_derxy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derxy(:,d)
                        UD_ver_derxy(i,j + 1,k,:,n) = UD_ver_derxy(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derxy(:,d)

                        UU_ver_deryy(i,j,k,:,n) = UU_ver_deryy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_deryy(:,d)
                        UD_ver_deryy(i,j + 1,k,:,n) = UD_ver_deryy(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_deryy(:,d)
                         
                        UU_ver_derxxx(i,j,k,:,n) = UU_ver_derxxx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derxxx(:,d)
                        UD_ver_derxxx(i,j + 1,k,:,n) = UD_ver_derxxx(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derxxx(:,d)

                        UU_ver_derxxy(i,j,k,:,n) = UU_ver_derxxy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derxxy(:,d)
                        UD_ver_derxxy(i,j + 1,k,:,n) = UD_ver_derxxy(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derxxy(:,d)

                        UU_ver_derxyy(i,j,k,:,n) = UU_ver_derxyy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_derxyy(:,d)
                        UD_ver_derxyy(i,j + 1,k,:,n) = UD_ver_derxyy(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_derxyy(:,d)

                        UU_ver_deryyy(i,j,k,:,n) = UU_ver_deryyy(i,j,k,:,n) + uh(i,j,k,d,n)*phiGU_ver_deryyy(:,d)
                        UD_ver_deryyy(i,j + 1,k,:,n) = UD_ver_deryyy(i,j + 1,k,:,n) + uh(i,j + 1,k,d,n)*phiGD_ver_deryyy(:,d)
                    end do
                end do
            end do
        end do
    end do
    UU_ver_derx =  UU_ver_derx*2/hx
    UD_ver_derx =  UD_ver_derx*2/hx
    UU_ver_dery =  UU_ver_dery*2/hy
    UD_ver_dery =  UD_ver_dery*2/hy

    UU_ver_derxx =  UU_ver_derxx*2/hx*2/hx
    UD_ver_derxx =  UD_ver_derxx*2/hx*2/hx
    UU_ver_derxy =  UU_ver_derxy*2/hx*2/hy
    UD_ver_derxy =  UD_ver_derxy*2/hx*2/hy
    UU_ver_deryy =  UU_ver_deryy*2/hy*2/hy
    UD_ver_deryy =  UD_ver_deryy*2/hy*2/hy

    UU_ver_derxxx =  UU_ver_derxxx*2/hx*2/hx*2/hx
    UD_ver_derxxx =  UD_ver_derxxx*2/hx*2/hx*2/hx
    UU_ver_derxxy =  UU_ver_derxxy*2/hx*2/hy*2/hx
    UD_ver_derxxy =  UD_ver_derxxy*2/hx*2/hy*2/hx
    UU_ver_derxyy =  UU_ver_derxyy*2/hy*2/hy*2/hy
    UD_ver_derxyy =  UD_ver_derxyy*2/hy*2/hy*2/hy
    UU_ver_deryyy =  UU_ver_deryyy*2/hy*2/hy*2/hy
    UD_ver_deryyy =  UD_ver_deryyy*2/hy*2/hy*2/hy
   
    ploydeg = 3 
    outerj: do i = 1,Nx
        outeri : do j = 1,Ny
            outerk :  do k = 0,Nphi
                scal = 0.5d0
                damping = 0.0d0
                rhoij = uh(i,j,k,1,1)
                u1ij = uh(i,j,k,1,2)/rhoij
                u2ij = uh(i,j,k,1,3)/rhoij
                Eij = uh(i,j,k,1,4)
                pressureij = gamma1*(Eij - 0.5d0*rhoij*(u1ij**2+u2ij**2))
                enthayij  =   ( Eij + pressureij ) /rhoij
                cij = sqrt(abs(gamma*pressureij/rhoij))
                betai = abs(u1ij) + cij 
                betaj = abs(u2ij) + cij
                Machij  = sqrt(u1ij**2 + u2ij**2)/cij
                scal = scal*1d0/enthayij
                left_edge_jump_density=0
                right_edge_jump_density =0
                bottom_edge_jump_density=0
                top_edge_jump_density=0
                do ss = 1,2
                    left_edge_jump_density = left_edge_jump_density  + abs(UL_ver(i,j,k,ss,1)-UR_ver(i-1,j,k,ss,1)) 
                    right_edge_jump_density = right_edge_jump_density  + abs(UL_ver(i+1,j,k,ss,1)-UR_ver(i,j,k,ss,1)) 
                    bottom_edge_jump_density = bottom_edge_jump_density + abs(UU_ver(i,j-1,k,ss,1) - UD_ver(i,j,k,ss,1))
                    top_edge_jump_density = top_edge_jump_density + abs(UU_ver(i,j,k,ss,1) - UD_ver(i,j+1,k,ss,1))
                end do
                deltadensity0 =  betai*(left_edge_jump_density+right_edge_jump_density)  + &
                betaj*(bottom_edge_jump_density+top_edge_jump_density)   
                left_edge_jump_rhou=0
                right_edge_jump_rhou =0
                bottom_edge_jump_rhou=0
                top_edge_jump_rhou=0
                do ss = 1,2
                    left_edge_jump_rhou = left_edge_jump_rhou  + abs(UL_ver(i,j,k,ss,2)-UR_ver(i-1,j,k,ss,2)) 
                    right_edge_jump_rhou = right_edge_jump_rhou  + abs(UL_ver(i+1,j,k,ss,2)-UR_ver(i,j,k,ss,2)) 
                    bottom_edge_jump_rhou = bottom_edge_jump_rhou + abs(UU_ver(i,j-1,k,ss,2) - UD_ver(i,j,k,ss,2))
                    top_edge_jump_rhou = top_edge_jump_rhou+ abs(UU_ver(i,j,k,ss,2) - UD_ver(i,j+1,k,ss,2))
                end do
                deltarhou0 =  betai*(left_edge_jump_rhou+right_edge_jump_rhou)  + &
                betaj*(bottom_edge_jump_rhou+top_edge_jump_rhou) 
                left_edge_jump_rhov=0
                right_edge_jump_rhov =0
                bottom_edge_jump_rhov=0
                top_edge_jump_rhov=0
                do ss = 1,2
                    left_edge_jump_rhov = left_edge_jump_rhov  + abs(UL_ver(i,j,k,ss,3)- UR_ver(i-1,j,k,ss,3)) 
                    right_edge_jump_rhov = right_edge_jump_rhov  + abs(UL_ver(i+1,j,k,ss,3)-UR_ver(i,j,k,ss,3)) 
                    bottom_edge_jump_rhov = bottom_edge_jump_rhov + abs(UU_ver(i,j-1,k,ss,3) - UD_ver(i,j,k,ss,3))
                    top_edge_jump_rhov = top_edge_jump_rhov + abs(UU_ver(i,j,k,ss,3) - UD_ver(i,j+1,k,ss,3))
                end do
                deltarhov0 =  betai*(left_edge_jump_rhov+right_edge_jump_rhov) + &
                betaj*(bottom_edge_jump_rhov+top_edge_jump_rhov) 
                left_edge_jump_Ener=0
                right_edge_jump_Ener =0
                bottom_edge_jump_Ener=0
                top_edge_jump_Ener=0
                do ss = 1,2
                    left_edge_jump_Ener= left_edge_jump_Ener  + abs(UL_ver(i,j,k,ss,4)-UR_ver(i-1,j,k,ss,4)) 
                    right_edge_jump_Ener= right_edge_jump_Ener  + abs(UL_ver(i+1,j,k,ss,4)-UR_ver(i,j,k,ss,4)) 
                    bottom_edge_jump_Ener= bottom_edge_jump_Ener + abs(UU_ver(i,j-1,k,ss,4) - UD_ver(i,j,k,ss,4))
                    top_edge_jump_Ener = top_edge_jump_Ener+ abs(UU_ver(i,j,k,ss,4) - UD_ver(i,j+1,k,ss,4))
                end do
                deltaEner0 =  betai*(left_edge_jump_Ener+right_edge_jump_Ener)  + &
                betaj*(bottom_edge_jump_Ener+top_edge_jump_Ener) 
                delta0max = max( deltadensity0 ,  deltarhou0 ,  deltarhov0 ,  deltaEner0)
                left_edge_jumpderx_density=0
                right_edge_jumpderx_density =0
                bottom_edge_jumpderx_density=0
                top_edge_jumpderx_density=0
                left_edge_jumpdery_density=0
                right_edge_jumpdery_density =0
                bottom_edge_jumpdery_density=0
                top_edge_jumpdery_density=0
                do ss = 1,2
                    left_edge_jumpderx_density = left_edge_jumpderx_density  + abs(UL_ver_derx(i,j,k,ss,1)-UR_ver_derx(i-1,j,k,ss,1)) 
                    right_edge_jumpderx_density = right_edge_jumpderx_density  + abs(UL_ver_derx(i+1,j,k,ss,1)-UR_ver_derx(i,j,k,ss,1)) 
                    bottom_edge_jumpderx_density = bottom_edge_jumpderx_density + abs(UU_ver_derx(i,j-1,k,ss,1) - UD_ver_derx(i,j,k,ss,1))
                    top_edge_jumpderx_density = top_edge_jumpderx_density + abs(UU_ver_derx(i,j,k,ss,1) - UD_ver_derx(i,j+1,k,ss,1))
                    left_edge_jumpdery_density = left_edge_jumpdery_density  + abs(UL_ver_dery(i,j,k,ss,1)-UR_ver_dery(i-1,j,k,ss,1)) 
                    right_edge_jumpdery_density = right_edge_jumpdery_density  + abs(UL_ver_dery(i+1,j,k,ss,1)-UR_ver_dery(i,j,k,ss,1)) 
                    bottom_edge_jumpdery_density = bottom_edge_jumpdery_density + abs(UU_ver_dery(i,j-1,k,ss,1) - UD_ver_dery(i,j,k,ss,1))
                    top_edge_jumpdery_density = top_edge_jumpdery_density + abs(UU_ver_dery(i,j,k,ss,1) - UD_ver_dery(i,j+1,k,ss,1))
                end do
                m_oe = 1
                deltadensity1 =  betai*(left_edge_jumpderx_density+right_edge_jumpderx_density + left_edge_jumpdery_density + right_edge_jumpdery_density)*2*hx + &
                betaj*(bottom_edge_jumpderx_density + top_edge_jumpdery_density + bottom_edge_jumpderx_density +   top_edge_jumpdery_density )*2*hy
                left_edge_jumpderx_rhou=0
                right_edge_jumpderx_rhou =0
                bottom_edge_jumpderx_rhou=0
                top_edge_jumpderx_rhou=0

                left_edge_jumpdery_rhou=0
                right_edge_jumpdery_rhou =0
                bottom_edge_jumpdery_rhou=0
                top_edge_jumpdery_rhou=0
                do ss = 1,2
                    left_edge_jumpderx_rhou = left_edge_jumpderx_rhou  + abs(UL_ver_derx(i,j,k,ss,2)-UR_ver_derx(i-1,j,k,ss,2)) 
                    right_edge_jumpderx_rhou = right_edge_jumpderx_rhou  + abs(UL_ver_derx(i+1,j,k,ss,2)-UR_ver_derx(i,j,k,ss,2)) 
                    bottom_edge_jumpderx_rhou = bottom_edge_jumpderx_rhou + abs(UU_ver_derx(i,j-1,k,ss,2) - UD_ver_derx(i,j,k,ss,2))
                    top_edge_jumpderx_rhou = top_edge_jumpderx_rhou + abs(UU_ver_derx(i,j,k,ss,2) - UD_ver_derx(i,j+1,k,ss,2))
                    left_edge_jumpdery_rhou = left_edge_jumpdery_rhou  + abs(UL_ver_dery(i,j,k,ss,2)-UR_ver_dery(i-1,j,k,ss,2)) 
                    right_edge_jumpdery_rhou = right_edge_jumpdery_rhou  + abs(UL_ver_dery(i+1,j,k,ss,2)-UR_ver_dery(i,j,k,ss,2)) 
                    bottom_edge_jumpdery_rhou = bottom_edge_jumpdery_rhou + abs(UU_ver_dery(i,j-1,k,ss,2) - UD_ver_dery(i,j,k,ss,2))
                    top_edge_jumpdery_rhou = top_edge_jumpdery_rhou + abs(UU_ver_dery(i,j,k,ss,2) - UD_ver_dery(i,j+1,k,ss,2))
                end do
                m_oe = 1
                deltarhou1 =  betai*(left_edge_jumpderx_rhou +right_edge_jumpderx_rhou + left_edge_jumpdery_rhou+right_edge_jumpdery_rhou)*2*hx  + &
                betaj*(bottom_edge_jumpderx_rhou + top_edge_jumpderx_rhou + bottom_edge_jumpdery_rhou + top_edge_jumpdery_rhou )*2*hy 
                
                left_edge_jumpderx_rhov=0
                right_edge_jumpderx_rhov =0
                bottom_edge_jumpderx_rhov=0
                top_edge_jumpderx_rhov=0
                left_edge_jumpdery_rhov=0
                right_edge_jumpdery_rhov =0
                bottom_edge_jumpdery_rhov=0
                top_edge_jumpdery_rhov=0
                do ss = 1,2
                    left_edge_jumpderx_rhov = left_edge_jumpderx_rhov  + abs(UL_ver_derx(i,j,k,ss,3)-UR_ver_derx(i-1,j,k,ss,3)) 
                    right_edge_jumpderx_rhov = right_edge_jumpderx_rhov  + abs(UL_ver_derx(i+1,j,k,ss,3)-UR_ver_derx(i,j,k,ss,3)) 
                    bottom_edge_jumpderx_rhov = bottom_edge_jumpderx_rhov + abs(UU_ver_derx(i,j-1,k,ss,3) - UD_ver_derx(i,j,k,ss,3))
                    top_edge_jumpderx_rhov = top_edge_jumpderx_rhov + abs(UU_ver_derx(i,j,k,ss,3) - UD_ver_derx(i,j+1,k,ss,3))
                    left_edge_jumpdery_rhov = left_edge_jumpdery_rhov+ abs(UL_ver_dery(i,j,k,ss,3)-UR_ver_dery(i-1,j,k,ss,3)) 
                    right_edge_jumpdery_rhov = right_edge_jumpdery_rhov + abs(UL_ver_dery(i+1,j,k,ss,3)-UR_ver_dery(i,j,k,ss,3)) 
                    bottom_edge_jumpdery_rhov= bottom_edge_jumpdery_rhov + abs(UU_ver_dery(i,j-1,k,ss,3) - UD_ver_dery(i,j,k,ss,3))
                    top_edge_jumpdery_rhov = top_edge_jumpdery_rhov + abs(UU_ver_dery(i,j,k,ss,3) - UD_ver_dery(i,j+1,k,ss,3))
                end do
                deltarhov1 =  betai*(left_edge_jumpderx_rhov +right_edge_jumpderx_rhov + left_edge_jumpdery_rhov+right_edge_jumpdery_rhov)*2*hx  + &
                betaj*(bottom_edge_jumpderx_rhov + top_edge_jumpderx_rhov + bottom_edge_jumpdery_rhov + top_edge_jumpdery_rhov )*2*hy  
                left_edge_jumpderx_Ener=0
                right_edge_jumpderx_Ener =0
                bottom_edge_jumpderx_Ener =0
                top_edge_jumpderx_Ener =0
                left_edge_jumpdery_Ener =0
                right_edge_jumpdery_Ener =0
                bottom_edge_jumpdery_Ener=0
                top_edge_jumpdery_Ener=0
                do ss = 1,2 
                    left_edge_jumpderx_Ener = left_edge_jumpderx_Ener  + abs(UL_ver_derx(i,j,k,ss,4)-UR_ver_derx(i-1,j,k,ss,4)) 
                    right_edge_jumpderx_Ener = right_edge_jumpderx_Ener  + abs(UL_ver_derx(i+1,j,k,ss,4)-UR_ver_derx(i,j,k,ss,4)) 
                    bottom_edge_jumpderx_Ener = bottom_edge_jumpderx_Ener + abs(UU_ver_derx(i,j-1,k,ss,4) - UD_ver_derx(i,j,k,ss,4))
                    top_edge_jumpderx_Ener = top_edge_jumpderx_Ener + abs(UU_ver_derx(i,j,k,ss,4) - UD_ver_derx(i,j+1,k,ss,4))
                    left_edge_jumpdery_Ener = left_edge_jumpdery_Ener + abs(UL_ver_dery(i,j,k,ss,4)-UR_ver_dery(i-1,j,k,ss,4)) 
                    right_edge_jumpdery_Ener = right_edge_jumpdery_Ener + abs(UL_ver_dery(i+1,j,k,ss,4)-UR_ver_dery(i,j,k,ss,4)) 
                    bottom_edge_jumpdery_Ener = bottom_edge_jumpdery_Ener + abs(UU_ver_dery(i,j-1,k,ss,4) - UD_ver_dery(i,j,k,ss,4))
                    top_edge_jumpdery_Ener = top_edge_jumpdery_Ener + abs(UU_ver_dery(i,j,k,ss,4) - UD_ver_dery(i,j+1,k,ss,4))
                end do
                deltaEner1 =  betai*(left_edge_jumpderx_Ener +right_edge_jumpderx_Ener + left_edge_jumpdery_Ener+right_edge_jumpdery_Ener)*2*hx  + &
                betaj*(bottom_edge_jumpderx_Ener + top_edge_jumpderx_Ener + bottom_edge_jumpdery_Ener + top_edge_jumpdery_Ener)*2*hy  
                delta1max = max( deltadensity1 ,  deltarhou1 ,  deltarhov1 ,  deltaEner1)
                damping1 =   delta0max + delta1max  
                damping =   scal*hx*damping1/hx
                uhmod(i,j,k,2:3,1:4) = exp(-dt*damping)*uh(i,j,k,2:3,1:4) 
                left_edge_jumpderxx_density=0
                right_edge_jumpderxx_density =0
                bottom_edge_jumpderxx_density=0
                top_edge_jumpderxx_density=0
                left_edge_jumpderxy_density=0
                right_edge_jumpderxy_density =0
                bottom_edge_jumpderxy_density=0
                top_edge_jumpderxy_density=0
                left_edge_jumpderyy_density=0
                right_edge_jumpderyy_density =0
                bottom_edge_jumpderyy_density=0
                top_edge_jumpderyy_density=0
                do ss = 1,2
                    left_edge_jumpderxx_density = left_edge_jumpderxx_density  + abs(UL_ver_derxx(i,j,k,ss,1)-UR_ver_derxx(i-1,j,k,ss,1)) 
                    right_edge_jumpderxx_density = right_edge_jumpderxx_density  + abs(UL_ver_derxx(i+1,j,k,ss,1)-UR_ver_derxx(i,j,k,ss,1)) 
                    bottom_edge_jumpderxx_density = bottom_edge_jumpderxx_density + abs(UU_ver_derxx(i,j-1,k,ss,1) - UD_ver_derxx(i,j,k,ss,1))
                    top_edge_jumpderxx_density = top_edge_jumpderxx_density + abs(UU_ver_derxx(i,j,k,ss,1) - UD_ver_derxx(i,j+1,k,ss,1))
                    left_edge_jumpderxy_density = left_edge_jumpderxy_density  + abs(UL_ver_derxy(i,j,k,ss,1)-UR_ver_derxy(i-1,j,k,ss,1)) 
                    right_edge_jumpderxy_density = right_edge_jumpderxy_density  + abs(UL_ver_derxy(i+1,j,k,ss,1)-UR_ver_derxy(i,j,k,ss,1)) 
                    bottom_edge_jumpderxy_density = bottom_edge_jumpderxy_density + abs(UU_ver_derxy(i,j-1,k,ss,1) - UD_ver_derxy(i,j,k,ss,1))
                    top_edge_jumpderxy_density = top_edge_jumpderxy_density + abs(UU_ver_derxy(i,j,k,ss,1) - UD_ver_derxy(i,j+1,k,ss,1))
                    left_edge_jumpderyy_density = left_edge_jumpderyy_density  + abs(UL_ver_deryy(i,j,k,ss,1)-UR_ver_deryy(i-1,j,k,ss,1)) 
                    right_edge_jumpderyy_density = right_edge_jumpderyy_density  + abs(UL_ver_deryy(i+1,j,k,ss,1)-UR_ver_deryy(i,j,k,ss,1)) 
                    bottom_edge_jumpderyy_density = bottom_edge_jumpderyy_density + abs(UU_ver_deryy(i,j-1,k,ss,1) - UD_ver_deryy(i,j,k,ss,1))
                    top_edge_jumpderyy_density = top_edge_jumpderyy_density + abs(UU_ver_deryy(i,j,k,ss,1) - UD_ver_deryy(i,j+1,k,ss,1))
                end do
                m_oe = 2
                deltadensity2 =  betai*(left_edge_jumpderxx_density+right_edge_jumpderxx_density + left_edge_jumpderxy_density+right_edge_jumpderxy_density + left_edge_jumpderyy_density+right_edge_jumpderyy_density )*2*3*hx**2  + &
                betaj*(bottom_edge_jumpderxx_density + top_edge_jumpderxx_density + bottom_edge_jumpderxy_density + top_edge_jumpderxy_density + bottom_edge_jumpderyy_density + top_edge_jumpderyy_density)*2*3*hy**2  
                left_edge_jumpderxx_rhou=0
                right_edge_jumpderxx_rhou =0
                bottom_edge_jumpderxx_rhou=0
                top_edge_jumpderxx_rhou=0
                left_edge_jumpderxy_rhou=0
                right_edge_jumpderxy_rhou =0
                bottom_edge_jumpderxy_rhou=0
                top_edge_jumpderxy_rhou=0
                left_edge_jumpderyy_rhou=0
                right_edge_jumpderyy_rhou =0
                bottom_edge_jumpderyy_rhou=0
                top_edge_jumpderyy_rhou=0
                do ss = 1,2
                    left_edge_jumpderxx_rhou = left_edge_jumpderxx_rhou  + abs(UL_ver_derxx(i,j,k,ss,2)-UR_ver_derxx(i-1,j,k,ss,2)) 
                    right_edge_jumpderxx_rhou = right_edge_jumpderxx_rhou  + abs(UL_ver_derxx(i+1,j,k,ss,2)-UR_ver_derxx(i,j,k,ss,2)) 
                    bottom_edge_jumpderxx_rhou = bottom_edge_jumpderxx_rhou + abs(UU_ver_derxx(i,j-1,k,ss,2) - UD_ver_derxx(i,j,k,ss,2))
                    top_edge_jumpderxx_rhou = top_edge_jumpderxx_rhou + abs(UU_ver_derxx(i,j,k,ss,2) - UD_ver_derxx(i,j+1,k,ss,2))
                    left_edge_jumpderxy_rhou = left_edge_jumpderxy_rhou  + abs(UL_ver_derxy(i,j,k,ss,2)-UR_ver_derxy(i-1,j,k,ss,2)) 
                    right_edge_jumpderxy_rhou = right_edge_jumpderxy_rhou  + abs(UL_ver_derxy(i+1,j,k,ss,2)-UR_ver_derxy(i,j,k,ss,2)) 
                    bottom_edge_jumpderxy_rhou = bottom_edge_jumpderxy_rhou + abs(UU_ver_derxy(i,j-1,k,ss,2) - UD_ver_derxy(i,j,k,ss,2))
                    top_edge_jumpderxy_rhou = top_edge_jumpderxy_rhou + abs(UU_ver_derxy(i,j,k,ss,2) - UD_ver_derxy(i,j+1,k,ss,2))
                    left_edge_jumpderyy_rhou = left_edge_jumpderyy_rhou  + abs(UL_ver_deryy(i,j,k,ss,2)-UR_ver_deryy(i-1,j,k,ss,2)) 
                    right_edge_jumpderyy_rhou = right_edge_jumpderyy_rhou  + abs(UL_ver_deryy(i+1,j,k,ss,2)-UR_ver_deryy(i,j,k,ss,2)) 
                    bottom_edge_jumpderyy_rhou = bottom_edge_jumpderyy_rhou + abs(UU_ver_deryy(i,j-1,k,ss,2) - UD_ver_deryy(i,j,k,ss,2))
                    top_edge_jumpderyy_rhou = top_edge_jumpderyy_rhou + abs(UU_ver_deryy(i,j,k,ss,2) - UD_ver_deryy(i,j+1,k,ss,2))
                end do
                m_oe = 2
                deltarhou2 =  betai*(left_edge_jumpderxx_rhou +right_edge_jumpderxx_rhou  + left_edge_jumpderxy_rhou +right_edge_jumpderxy_rhou  + left_edge_jumpderyy_rhou +right_edge_jumpderyy_rhou )*2*3*hx**2   + &
                betaj*(bottom_edge_jumpderxx_rhou  + top_edge_jumpderxx_rhou  + bottom_edge_jumpderxy_rhou  + top_edge_jumpderxy_rhou  + bottom_edge_jumpderyy_rhou  + top_edge_jumpderyy_rhou )*2*3*hy**2 
                left_edge_jumpderxx_rhov=0
                right_edge_jumpderxx_rhov =0
                bottom_edge_jumpderxx_rhov=0
                top_edge_jumpderxx_rhov=0
                left_edge_jumpderxy_rhov=0
                right_edge_jumpderxy_rhov =0
                bottom_edge_jumpderxy_rhov=0
                top_edge_jumpderxy_rhov=0
                left_edge_jumpderyy_rhov=0
                right_edge_jumpderyy_rhov =0
                bottom_edge_jumpderyy_rhov=0
                top_edge_jumpderyy_rhov=0
                do ss = 1,2
                    left_edge_jumpderxx_rhov = left_edge_jumpderxx_rhov  + abs(UL_ver_derxx(i,j,k,ss,3)-UR_ver_derxx(i-1,j,k,ss,3)) 
                    right_edge_jumpderxx_rhov = right_edge_jumpderxx_rhov  + abs(UL_ver_derxx(i+1,j,k,ss,3)-UR_ver_derxx(i,j,k,ss,3)) 
                    bottom_edge_jumpderxx_rhov = bottom_edge_jumpderxx_rhov + abs(UU_ver_derxx(i,j-1,k,ss,3) - UD_ver_derxx(i,j,k,ss,3))
                    top_edge_jumpderxx_rhov = top_edge_jumpderxx_rhov + abs(UU_ver_derxx(i,j,k,ss,3) - UD_ver_derxx(i,j+1,k,ss,3))
                    left_edge_jumpderxy_rhov = left_edge_jumpderxy_rhov  + abs(UL_ver_derxy(i,j,k,ss,3)-UR_ver_derxy(i-1,j,k,ss,3)) 
                    right_edge_jumpderxy_rhov = right_edge_jumpderxy_rhov + abs(UL_ver_derxy(i+1,j,k,ss,3)-UR_ver_derxy(i,j,k,ss,3)) 
                    bottom_edge_jumpderxy_rhov = bottom_edge_jumpderxy_rhov + abs(UU_ver_derxy(i,j-1,k,ss,3) - UD_ver_derxy(i,j,k,ss,3))
                    top_edge_jumpderxy_rhov = top_edge_jumpderxy_rhov + abs(UU_ver_derxy(i,j,k,ss,3) - UD_ver_derxy(i,j+1,k,ss,3))
                    left_edge_jumpderyy_rhov = left_edge_jumpderyy_rhov  + abs(UL_ver_deryy(i,j,k,ss,3)-UR_ver_deryy(i-1,j,k,ss,3)) 
                    right_edge_jumpderyy_rhov = right_edge_jumpderyy_rhov  + abs(UL_ver_deryy(i+1,j,k,ss,3)-UR_ver_deryy(i,j,k,ss,3)) 
                    bottom_edge_jumpderyy_rhov = bottom_edge_jumpderyy_rhov + abs(UU_ver_deryy(i,j-1,k,ss,3) - UD_ver_deryy(i,j,k,ss,3))
                    top_edge_jumpderyy_rhov = top_edge_jumpderyy_rhov + abs(UU_ver_deryy(i,j,k,ss,3) - UD_ver_deryy(i,j+1,k,ss,3))
                end do
                m_oe = 2
                deltarhov2 =  betai*(left_edge_jumpderxx_rhov +right_edge_jumpderxx_rhov  + left_edge_jumpderxy_rhov +right_edge_jumpderxy_rhov  + left_edge_jumpderyy_rhov +right_edge_jumpderyy_rhov  )*2*3*hx**2  + &
                betaj*(bottom_edge_jumpderxx_rhov  + top_edge_jumpderxx_rhov  + bottom_edge_jumpderxy_rhov + top_edge_jumpderxy_rhov  + bottom_edge_jumpderyy_rhov  + top_edge_jumpderyy_rhov)*2*3*hy**2 
                left_edge_jumpderxx_Ener=0
                right_edge_jumpderxx_Ener =0
                bottom_edge_jumpderxx_Ener=0
                top_edge_jumpderxx_Ener=0
                left_edge_jumpderxy_Ener=0
                right_edge_jumpderxy_Ener =0
                bottom_edge_jumpderxy_Ener=0
                top_edge_jumpderxy_Ener=0
                left_edge_jumpderyy_Ener=0
                right_edge_jumpderyy_Ener =0
                bottom_edge_jumpderyy_Ener=0
                top_edge_jumpderyy_Ener=0
                do ss = 1,2
                    left_edge_jumpderxx_Ener = left_edge_jumpderxx_Ener  + abs(UL_ver_derxx(i,j,k,ss,4)-UR_ver_derxx(i-1,j,k,ss,4)) 
                    right_edge_jumpderxx_Ener = right_edge_jumpderxx_Ener  + abs(UL_ver_derxx(i+1,j,k,ss,4)-UR_ver_derxx(i,j,k,ss,4)) 
                    bottom_edge_jumpderxx_Ener = bottom_edge_jumpderxx_Ener + abs(UU_ver_derxx(i,j-1,k,ss,4) - UD_ver_derxx(i,j,k,ss,4))
                    top_edge_jumpderxx_Ener = top_edge_jumpderxx_Ener + abs(UU_ver_derxx(i,j,k,ss,4) - UD_ver_derxx(i,j+1,k,ss,4))
                    left_edge_jumpderxy_Ener = left_edge_jumpderxy_Ener  + abs(UL_ver_derxy(i,j,k,ss,4)-UR_ver_derxy(i-1,j,k,ss,4)) 
                    right_edge_jumpderxy_Ener = right_edge_jumpderxy_Ener + abs(UL_ver_derxy(i+1,j,k,ss,4)-UR_ver_derxy(i,j,k,ss,4)) 
                    bottom_edge_jumpderxy_Ener = bottom_edge_jumpderxy_Ener + abs(UU_ver_derxy(i,j-1,k,ss,4) - UD_ver_derxy(i,j,k,ss,4))
                    top_edge_jumpderxy_Ener = top_edge_jumpderxy_Ener+ abs(UU_ver_derxy(i,j,k,ss,4) - UD_ver_derxy(i,j+1,k,ss,4))
                    left_edge_jumpderyy_Ener = left_edge_jumpderyy_Ener  + abs(UL_ver_deryy(i,j,k,ss,4)-UR_ver_deryy(i-1,j,k,ss,4)) 
                    right_edge_jumpderyy_Ener = right_edge_jumpderyy_Ener  + abs(UL_ver_deryy(i+1,j,k,ss,4)-UR_ver_deryy(i,j,k,ss,4)) 
                    bottom_edge_jumpderyy_Ener = bottom_edge_jumpderyy_Ener + abs(UU_ver_deryy(i,j-1,k,ss,4) - UD_ver_deryy(i,j,k,ss,4))
                    top_edge_jumpderyy_Ener = top_edge_jumpderyy_Ener + abs(UU_ver_deryy(i,j,k,ss,4) - UD_ver_deryy(i,j+1,k,ss,4))
                end do
                m_oe = 2
                deltaEner2 =  betai*(left_edge_jumpderxx_Ener +right_edge_jumpderxx_Ener  + left_edge_jumpderxy_Ener +right_edge_jumpderxy_Ener  + left_edge_jumpderyy_Ener +right_edge_jumpderyy_Ener  )*2*3*hx**2  + &
                betaj*(bottom_edge_jumpderxx_Ener  + top_edge_jumpderxx_Ener + bottom_edge_jumpderxy_Ener + top_edge_jumpderxy_Ener  + bottom_edge_jumpderyy_Ener  + top_edge_jumpderyy_Ener)*2*3*hy**2 
                delta2max = max( deltadensity2 ,  deltarhou2 ,  deltarhov2 ,  deltaEner2)
                damping2 =   damping1 + delta2max   
                damping =  scal*hx*damping2/hx
                uhmod(i,j,k,4:6,1:4) = exp(-dt*damping)*uh(i,j,k,4:6,1:4) 
                left_edge_jumpderxxx_density=0
                right_edge_jumpderxxx_density =0
                bottom_edge_jumpderxxx_density=0
                top_edge_jumpderxxx_density=0
                left_edge_jumpderxxy_density=0
                right_edge_jumpderxxy_density =0
                bottom_edge_jumpderxxy_density=0
                top_edge_jumpderxxy_density=0
                left_edge_jumpderxyy_density=0
                right_edge_jumpderxyy_density =0
                bottom_edge_jumpderxyy_density=0
                top_edge_jumpderxyy_density=0
                left_edge_jumpderyyy_density=0
                right_edge_jumpderyyy_density =0
                bottom_edge_jumpderyyy_density=0
                top_edge_jumpderyyy_density=0
                do ss = 1,2
                    left_edge_jumpderxxx_density = left_edge_jumpderxxx_density  + abs(UL_ver_derxxx(i,j,k,ss,1)-UR_ver_derxxx(i-1,j,k,ss,1)) 
                    right_edge_jumpderxxx_density = right_edge_jumpderxxx_density  + abs(UL_ver_derxxx(i+1,j,k,ss,1)-UR_ver_derxxx(i,j,k,ss,1)) 
                    bottom_edge_jumpderxxx_density = bottom_edge_jumpderxxx_density + abs(UU_ver_derxxx(i,j-1,k,ss,1) - UD_ver_derxxx(i,j,k,ss,1))
                    top_edge_jumpderxxx_density = top_edge_jumpderxxx_density + abs(UU_ver_derxxx(i,j,k,ss,1) - UD_ver_derxxx(i,j+1,k,ss,1))
                    left_edge_jumpderxxy_density = left_edge_jumpderxxy_density  + abs(UL_ver_derxxy(i,j,k,ss,1)-UR_ver_derxxy(i-1,j,k,ss,1)) 
                    right_edge_jumpderxxy_density = right_edge_jumpderxxy_density  + abs(UL_ver_derxxy(i+1,j,k,ss,1)-UR_ver_derxxy(i,j,k,ss,1)) 
                    bottom_edge_jumpderxxy_density = bottom_edge_jumpderxxy_density + abs(UU_ver_derxxy(i,j-1,k,ss,1) - UD_ver_derxxy(i,j,k,ss,1))
                    top_edge_jumpderxxy_density = top_edge_jumpderxxy_density + abs(UU_ver_derxxy(i,j,k,ss,1) - UD_ver_derxxy(i,j+1,k,ss,1))
                    left_edge_jumpderxyy_density = left_edge_jumpderxyy_density  + abs(UL_ver_derxyy(i,j,k,ss,1)-UR_ver_derxyy(i-1,j,k,ss,1)) 
                    right_edge_jumpderxyy_density = right_edge_jumpderxyy_density  + abs(UL_ver_derxyy(i+1,j,k,ss,1)-UR_ver_derxyy(i,j,k,ss,1)) 
                    bottom_edge_jumpderxyy_density = bottom_edge_jumpderxyy_density + abs(UU_ver_derxyy(i,j-1,k,ss,1) - UD_ver_derxyy(i,j,k,ss,1))
                    top_edge_jumpderxyy_density = top_edge_jumpderxyy_density + abs(UU_ver_derxyy(i,j,k,ss,1) - UD_ver_derxyy(i,j+1,k,ss,1))
                    left_edge_jumpderyyy_density = left_edge_jumpderyyy_density  + abs(UL_ver_deryyy(i,j,k,ss,1)-UR_ver_deryyy(i-1,j,k,ss,1)) 
                    right_edge_jumpderyyy_density = right_edge_jumpderyyy_density  + abs(UL_ver_deryyy(i+1,j,k,ss,1)-UR_ver_deryyy(i,j,k,ss,1)) 
                    bottom_edge_jumpderyyy_density = bottom_edge_jumpderyyy_density + abs(UU_ver_deryyy(i,j-1,k,ss,1) - UD_ver_deryyy(i,j,k,ss,1))
                    top_edge_jumpderyyy_density = top_edge_jumpderyyy_density + abs(UU_ver_deryyy(i,j,k,ss,1) - UD_ver_deryyy(i,j+1,k,ss,1))
                end do
                m_oe = 3
                deltadensity3 =  betai*(left_edge_jumpderxxx_density+right_edge_jumpderxxx_density + left_edge_jumpderxxy_density+right_edge_jumpderxxy_density + left_edge_jumpderxyy_density+right_edge_jumpderxyy_density + left_edge_jumpderyyy_density+right_edge_jumpderyyy_density )*3*4*hx**3  + &
                betaj*(bottom_edge_jumpderxxx_density + top_edge_jumpderxxx_density + bottom_edge_jumpderxxy_density + top_edge_jumpderxxy_density + bottom_edge_jumpderxyy_density + top_edge_jumpderxyy_density+ bottom_edge_jumpderyyy_density + top_edge_jumpderyyy_density)*3*4*hy**3  
 
                left_edge_jumpderxxx_rhou=0
                right_edge_jumpderxxx_rhou =0
                bottom_edge_jumpderxxx_rhou=0
                top_edge_jumpderxxx_rhou=0
                left_edge_jumpderxxy_rhou=0
                right_edge_jumpderxxy_rhou =0
                bottom_edge_jumpderxxy_rhou=0
                top_edge_jumpderxxy_rhou=0
                left_edge_jumpderxyy_rhou=0
                right_edge_jumpderxyy_rhou =0
                bottom_edge_jumpderxyy_rhou=0
                top_edge_jumpderxyy_rhou=0
                left_edge_jumpderyyy_rhou=0
                right_edge_jumpderyyy_rhou =0
                bottom_edge_jumpderyyy_rhou=0
                top_edge_jumpderyyy_rhou=0
                do ss = 1,2
                    left_edge_jumpderxxx_rhou = left_edge_jumpderxxx_rhou  + abs(UL_ver_derxxx(i,j,k,ss,2)-UR_ver_derxxx(i-1,j,k,ss,2)) 
                    right_edge_jumpderxxx_rhou = right_edge_jumpderxxx_rhou  + abs(UL_ver_derxxx(i+1,j,k,ss,2)-UR_ver_derxxx(i,j,k,ss,2)) 
                    bottom_edge_jumpderxxx_rhou = bottom_edge_jumpderxxx_rhou + abs(UU_ver_derxxx(i,j-1,k,ss,2) - UD_ver_derxxx(i,j,k,ss,2))
                    top_edge_jumpderxxx_rhou = top_edge_jumpderxxx_rhou + abs(UU_ver_derxxx(i,j,k,ss,2) - UD_ver_derxxx(i,j+1,k,ss,2)) 
                    left_edge_jumpderxxy_rhou = left_edge_jumpderxxy_rhou  + abs(UL_ver_derxxy(i,j,k,ss,2)-UR_ver_derxxy(i-1,j,k,ss,2)) 
                    right_edge_jumpderxxy_rhou = right_edge_jumpderxxy_rhou  + abs(UL_ver_derxxy(i+1,j,k,ss,2)-UR_ver_derxxy(i,j,k,ss,2)) 
                    bottom_edge_jumpderxxy_rhou = bottom_edge_jumpderxxy_rhou + abs(UU_ver_derxxy(i,j-1,k,ss,2) - UD_ver_derxxy(i,j,k,ss,2))
                    top_edge_jumpderxxy_rhou = top_edge_jumpderxxy_rhou + abs(UU_ver_derxxy(i,j,k,ss,2) - UD_ver_derxxy(i,j+1,k,ss,2))
                    left_edge_jumpderxyy_rhou = left_edge_jumpderxyy_rhou  + abs(UL_ver_derxyy(i,j,k,ss,2)-UR_ver_derxyy(i-1,j,k,ss,2)) 
                    right_edge_jumpderxyy_rhou = right_edge_jumpderxyy_rhou  + abs(UL_ver_derxyy(i+1,j,k,ss,2)-UR_ver_derxyy(i,j,k,ss,2)) 
                    bottom_edge_jumpderxyy_rhou = bottom_edge_jumpderxyy_rhou + abs(UU_ver_derxyy(i,j-1,k,ss,2) - UD_ver_derxyy(i,j,k,ss,2))
                    top_edge_jumpderxyy_rhou = top_edge_jumpderxyy_rhou + abs(UU_ver_derxyy(i,j,k,ss,2) - UD_ver_derxyy(i,j+1,k,ss,2))
                    left_edge_jumpderyyy_rhou = left_edge_jumpderyyy_rhou  + abs(UL_ver_deryyy(i,j,k,ss,2)-UR_ver_deryyy(i-1,j,k,ss,2)) 
                    right_edge_jumpderyyy_rhou = right_edge_jumpderyyy_rhou  + abs(UL_ver_deryyy(i+1,j,k,ss,2)-UR_ver_deryyy(i,j,k,ss,2)) 
                    bottom_edge_jumpderyyy_rhou = bottom_edge_jumpderyyy_rhou + abs(UU_ver_deryyy(i,j-1,k,ss,2) - UD_ver_deryyy(i,j,k,ss,2))
                    top_edge_jumpderyyy_rhou = top_edge_jumpderyyy_rhou + abs(UU_ver_deryyy(i,j,k,ss,2) - UD_ver_deryyy(i,j+1,k,ss,2))
                end do
                m_oe = 3
                deltarhou3 =  betai*(left_edge_jumpderxxx_rhou +right_edge_jumpderxxx_rhou  + left_edge_jumpderxxy_rhou +right_edge_jumpderxxy_rhou  + left_edge_jumpderxyy_rhou +right_edge_jumpderxyy_rhou + left_edge_jumpderyyy_rhou +right_edge_jumpderyyy_rhou)*3*4*hx**3   + &
                betaj*(bottom_edge_jumpderxxx_rhou  + top_edge_jumpderxxx_rhou  + bottom_edge_jumpderxxy_rhou  + top_edge_jumpderxxy_rhou  + bottom_edge_jumpderxyy_rhou  + top_edge_jumpderxyy_rhou + bottom_edge_jumpderyyy_rhou  + top_edge_jumpderyyy_rhou)*3*4*hy**3 
                ! momentumn2 --rhov
                left_edge_jumpderxxx_rhov=0
                right_edge_jumpderxxx_rhov =0
                bottom_edge_jumpderxxx_rhov=0
                top_edge_jumpderxxx_rhov=0
                left_edge_jumpderxxy_rhov=0
                right_edge_jumpderxxy_rhov =0
                bottom_edge_jumpderxxy_rhov=0
                top_edge_jumpderxxy_rhov=0
                left_edge_jumpderxyy_rhov=0
                right_edge_jumpderxyy_rhov =0
                bottom_edge_jumpderxyy_rhov=0
                top_edge_jumpderxyy_rhov=0
                left_edge_jumpderyyy_rhov=0
                right_edge_jumpderyyy_rhov =0
                bottom_edge_jumpderyyy_rhov=0
                top_edge_jumpderyyy_rhov=0
                do ss = 1,2
                    left_edge_jumpderxxx_rhov = left_edge_jumpderxxx_rhov  + abs(UL_ver_derxxx(i,j,k,ss,3)-UR_ver_derxxx(i-1,j,k,ss,3)) 
                    right_edge_jumpderxxx_rhov = right_edge_jumpderxxx_rhov  + abs(UL_ver_derxxx(i+1,j,k,ss,3)-UR_ver_derxxx(i,j,k,ss,3)) 
                    bottom_edge_jumpderxxx_rhov = bottom_edge_jumpderxxx_rhov + abs(UU_ver_derxxx(i,j-1,k,ss,3) - UD_ver_derxxx(i,j,k,ss,3))
                    top_edge_jumpderxxx_rhov = top_edge_jumpderxxx_rhov + abs(UU_ver_derxxx(i,j,k,ss,3) - UD_ver_derxxx(i,j+1,k,ss,3))
                    left_edge_jumpderxxy_rhov = left_edge_jumpderxxy_rhov  + abs(UL_ver_derxxy(i,j,k,ss,3)-UR_ver_derxxy(i-1,j,k,ss,3)) 
                    right_edge_jumpderxxy_rhov = right_edge_jumpderxxy_rhov + abs(UL_ver_derxxy(i+1,j,k,ss,3)-UR_ver_derxxy(i,j,k,ss,3)) 
                    bottom_edge_jumpderxxy_rhov = bottom_edge_jumpderxxy_rhov + abs(UU_ver_derxxy(i,j-1,k,ss,3) - UD_ver_derxxy(i,j,k,ss,3))
                    top_edge_jumpderxxy_rhov = top_edge_jumpderxxy_rhov + abs(UU_ver_derxxy(i,j,k,ss,3) - UD_ver_derxxy(i,j+1,k,ss,3))
                    left_edge_jumpderxyy_rhov = left_edge_jumpderxyy_rhov  + abs(UL_ver_derxyy(i,j,k,ss,3)-UR_ver_derxyy(i-1,j,k,ss,3)) 
                    right_edge_jumpderxyy_rhov = right_edge_jumpderxyy_rhov  + abs(UL_ver_derxyy(i+1,j,k,ss,3)-UR_ver_derxyy(i,j,k,ss,3)) 
                    bottom_edge_jumpderxyy_rhov = bottom_edge_jumpderxyy_rhov + abs(UU_ver_derxyy(i,j-1,k,ss,3) - UD_ver_derxyy(i,j,k,ss,3))
                    top_edge_jumpderxyy_rhov = top_edge_jumpderxyy_rhov + abs(UU_ver_derxyy(i,j,k,ss,3) - UD_ver_derxyy(i,j+1,k,ss,3))
                    left_edge_jumpderyyy_rhov = left_edge_jumpderyyy_rhov  + abs(UL_ver_deryyy(i,j,k,ss,3)-UR_ver_deryyy(i-1,j,k,ss,3)) 
                    right_edge_jumpderyyy_rhov = right_edge_jumpderyyy_rhov  + abs(UL_ver_deryyy(i+1,j,k,ss,3)-UR_ver_deryyy(i,j,k,ss,3)) 
                    bottom_edge_jumpderyyy_rhov = bottom_edge_jumpderyyy_rhov + abs(UU_ver_deryyy(i,j-1,k,ss,3) - UD_ver_deryyy(i,j,k,ss,3))
                    top_edge_jumpderyyy_rhov = top_edge_jumpderyyy_rhov + abs(UU_ver_deryyy(i,j,k,ss,3) - UD_ver_deryyy(i,j+1,k,ss,3))
                end do
                m_oe = 3

                deltarhov3 =  betai*(left_edge_jumpderxxx_rhov +right_edge_jumpderxxx_rhov  + left_edge_jumpderxxy_rhov +right_edge_jumpderxxy_rhov  + left_edge_jumpderxyy_rhov +right_edge_jumpderxyy_rhov  + left_edge_jumpderyyy_rhov +right_edge_jumpderyyy_rhov )*3*4*hx**3  + &
                betaj*(bottom_edge_jumpderxxx_rhov  + top_edge_jumpderxxx_rhov  + bottom_edge_jumpderxxy_rhov + top_edge_jumpderxxy_rhov  + bottom_edge_jumpderxyy_rhov  + top_edge_jumpderxyy_rhov+bottom_edge_jumpderyyy_rhov  + top_edge_jumpderyyy_rhov)*3*4*hy**3
                left_edge_jumpderxxx_Ener=0
                right_edge_jumpderxxx_Ener =0
                bottom_edge_jumpderxxx_Ener=0
                top_edge_jumpderxxx_Ener=0
                left_edge_jumpderxxy_Ener=0
                right_edge_jumpderxxy_Ener =0
                bottom_edge_jumpderxxy_Ener=0
                top_edge_jumpderxxy_Ener=0
                left_edge_jumpderxyy_Ener=0
                right_edge_jumpderxyy_Ener =0
                bottom_edge_jumpderxyy_Ener=0
                top_edge_jumpderxyy_Ener=0
                left_edge_jumpderyyy_Ener=0
                right_edge_jumpderyyy_Ener =0
                bottom_edge_jumpderyyy_Ener=0
                top_edge_jumpderyyy_Ener=0
                do ss = 1,2
                    left_edge_jumpderxxx_Ener = left_edge_jumpderxxx_Ener  + abs(UL_ver_derxxx(i,j,k,ss,4)-UR_ver_derxxx(i-1,j,k,ss,4)) 
                    right_edge_jumpderxxx_Ener = right_edge_jumpderxxx_Ener  + abs(UL_ver_derxxx(i+1,j,k,ss,4)-UR_ver_derxxx(i,j,k,ss,4)) 
                    bottom_edge_jumpderxxx_Ener = bottom_edge_jumpderxxx_Ener + abs(UU_ver_derxxx(i,j-1,k,ss,4) - UD_ver_derxxx(i,j,k,ss,4))
                    top_edge_jumpderxxx_Ener = top_edge_jumpderxxx_Ener + abs(UU_ver_derxxx(i,j,k,ss,4) - UD_ver_derxxx(i,j+1,k,ss,4))
                    left_edge_jumpderxxy_Ener = left_edge_jumpderxxy_Ener  + abs(UL_ver_derxxy(i,j,k,ss,4)-UR_ver_derxxy(i-1,j,k,ss,4)) 
                    right_edge_jumpderxxy_Ener = right_edge_jumpderxxy_Ener + abs(UL_ver_derxxy(i+1,j,k,ss,4)-UR_ver_derxxy(i,j,k,ss,4)) 
                    bottom_edge_jumpderxxy_Ener = bottom_edge_jumpderxxy_Ener + abs(UU_ver_derxxy(i,j-1,k,ss,4) - UD_ver_derxxy(i,j,k,ss,4))
                    top_edge_jumpderxxy_Ener = top_edge_jumpderxxy_Ener+ abs(UU_ver_derxxy(i,j,k,ss,4) - UD_ver_derxxy(i,j+1,k,ss,4))
                    left_edge_jumpderxyy_Ener = left_edge_jumpderxyy_Ener  + abs(UL_ver_derxyy(i,j,k,ss,4)-UR_ver_derxyy(i-1,j,k,ss,4)) 
                    right_edge_jumpderxyy_Ener = right_edge_jumpderxyy_Ener  + abs(UL_ver_derxyy(i+1,j,k,ss,4)-UR_ver_derxyy(i,j,k,ss,4)) 
                    bottom_edge_jumpderxyy_Ener = bottom_edge_jumpderxyy_Ener + abs(UU_ver_derxyy(i,j-1,k,ss,4) - UD_ver_derxyy(i,j,k,ss,4))
                    top_edge_jumpderxyy_Ener = top_edge_jumpderxyy_Ener + abs(UU_ver_derxyy(i,j,k,ss,4) - UD_ver_derxyy(i,j+1,k,ss,4))
                    left_edge_jumpderyyy_Ener = left_edge_jumpderyyy_Ener  + abs(UL_ver_deryyy(i,j,k,ss,4)-UR_ver_deryyy(i-1,j,k,ss,4)) 
                    right_edge_jumpderyyy_Ener = right_edge_jumpderyyy_Ener  + abs(UL_ver_deryyy(i+1,j,k,ss,4)-UR_ver_deryyy(i,j,k,ss,4)) 
                    bottom_edge_jumpderyyy_Ener = bottom_edge_jumpderyyy_Ener + abs(UU_ver_deryyy(i,j-1,k,ss,4) - UD_ver_deryyy(i,j,k,ss,4))
                    top_edge_jumpderyyy_Ener = top_edge_jumpderyyy_Ener + abs(UU_ver_deryyy(i,j,k,ss,4) - UD_ver_deryyy(i,j+1,k,ss,4))
                end do
                m_oe = 3
                deltaEner3 =  betai*(left_edge_jumpderxxx_Ener +right_edge_jumpderxxx_Ener  + left_edge_jumpderxxy_Ener +right_edge_jumpderxxy_Ener  + left_edge_jumpderxyy_Ener +right_edge_jumpderxyy_Ener + left_edge_jumpderyyy_Ener +right_edge_jumpderyyy_Ener )*3*4*hx**3  + &
                betaj*(bottom_edge_jumpderxxx_Ener  + top_edge_jumpderxxx_Ener + bottom_edge_jumpderxxy_Ener + top_edge_jumpderxxy_Ener  + bottom_edge_jumpderxyy_Ener  + top_edge_jumpderxyy_Ener+bottom_edge_jumpderyyy_Ener  + top_edge_jumpderyyy_Ener)*3*4*hy**3 
                delta3max = max( deltadensity3 ,  deltarhou3 ,  deltarhov3 ,  deltaEner3)
                damping3 =   damping2 + delta3max    
                damping =  scal*hx*damping3/hx 
                uhmod(i,j,k,7:10,1:4) = exp(-dt*damping)*uh(i,j,k,7:10,1:4) 
            end do outerk
        end do outeri 
    end do outerj 
    uh = uhmod
    end  subroutine  jumpfilter
    
    
    
