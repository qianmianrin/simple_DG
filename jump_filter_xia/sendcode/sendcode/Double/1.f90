     
    module com
    ! P1 Double Mach reflection case
    ! ????? 640 ?? 200??768 ?? 240??1408 ?? 440
    include 'mpif.h'

    integer Nx0,Ny0
    integer N_process,Nx_process,Ny_process
    integer Nx,Ny,kk,NumEq,NumGLP
    parameter(N_process = 16)
    parameter(Nx0 = 768, Ny0 = 240 , Lphi = 0, kk = 1, NumEq = 4, NumGLP = 5, RKorder = 4, flux_type = 1)
    parameter(Nx_process = sqrt(1.0*N_process), Ny_process = sqrt(1.0*N_process))
    parameter(Nx = Nx0/Nx_process, Ny = Ny0/Ny_process)
    parameter(Nx1 = Nx + 1,Ny1 = Ny + 1)

    real(8) pi,gamma,gamma1
    parameter(dimPk = (kk + 1)*(kk + 2)/2)
    parameter(dimPk1 = (kk + 1)*(kk + 2)/2)
    parameter(Nphi = max(2*Lphi - 1,0))
    parameter(Nphi1 = Nphi + 1)
    parameter(gamma = 1.4d0)
    parameter(gamma1 = gamma - 1)
    parameter(pi = 4*atan(1d0))

    ! Limiter type: 1-->jumpfilter; 2-->hybrid jump filter
    parameter(limitertypeall = 1)

    ! The numerical solution and mesh
    real(8) xa,xb,ya,yb,t,dt,tend,CFL,umax,umax1,tRK,t1,t2,alphax,alphay,totaldiv,rij
    real(8) uh(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),du(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uI(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq),uII(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) uh00(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    real(8) hx,hy,Xc(Nx),Yc(Ny),Xc0(Nx0),Yc0(Ny0),Phi(0:Nphi),hx1,hy1,hphi
    real(8) Bx(0:Nx,0:Ny1,0:Nphi,kk + 1),By(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) dBx(0:Nx,0:Ny1,0:Nphi,kk + 1),dBy(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxI(0:Nx,0:Ny1,0:Nphi,kk + 1),ByI(0:Nx1,0:Ny,0:Nphi,kk + 1)
    real(8) BxII(0:Nx,0:Ny1,0:Nphi,kk + 1),ByII(0:Nx1,0:Ny,0:Nphi,kk + 1)

    ! The basis
    real(8) lambda(NumGLP),weight(NumGLP),sink(0:Nphi,Lphi),cosk(0:Nphi,Lphi)
    real(8) phiG(NumGLP,NumGLP,dimPk),phixG(NumGLP,NumGLP,dimPk),phiyG(NumGLP,NumGLP,dimPk),mm(dimPk)
    real(8) phiGLL(NumGLP,NumGLP,dimPk,2),lambdaL(NumGLP)
    real(8) phiGR(NumGLP,dimPk), phiGL(NumGLP,dimPk), phiGU(NumGLP,dimPk), phiGD(NumGLP,dimPk)
    real(8) phiRU(dimPk), phiLU(dimPk), phiRD(dimPk), phiLD(dimPk)
    real(8) EzG(NumGLP,kk + 1),EzxG(NumGLP,kk + 1),EzyG(NumGLP,kk + 1),mmE(kk + 1)
    real(8) EzR(kk + 1),EzL(kk + 1),EzU(kk + 1),EzD(kk + 1),omega1(Nx,Ny,0:Nphi)


    ! The Lh
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

    ! The Limiter
    real(8) M,beta
    real(8) DeltaUR1(NumEq,1),DeltaUL1(NumEq,1),DeltaUU1(NumEq,1),DeltaUD1(NumEq,1),DeltaU1(NumEq,1),DeltaUmod1(NumEq,1)
    real(8) DeltaUR1mod(NumEq,1),DeltaUL1mod(NumEq,1),DeltaUU1mod(NumEq,1),DeltaUD1mod(NumEq,1)
    real(8) R(NumEq,NumEq),L(NumEq,NumEq)
    real(8) DeltaUR(NumEq,1),DeltaUL(NumEq,1),DeltaU(NumEq,1),DeltaUmod(NumEq,1)
    real(8) Is_trouble_cell(Nx,Ny,0:Nphi)
    real(8) change_all(Nx,Ny)

    ! jump filter
    real(8) densityave,momentum1ave,momentum2ave,Enerave
    real(8) phiGR_ver(2,3),phiGL_ver(2,3), phiGU_ver(2,3), phiGD_ver(2,3)
    real(8) phiGR_ver_derx(2,3),phiGL_ver_derx(2,3), phiGU_ver_derx(2,3), phiGD_ver_derx(2,3)
    real(8) phiGR_ver_dery(2,3),phiGL_ver_dery(2,3), phiGU_ver_dery(2,3), phiGD_ver_dery(2,3)
    real(8) jump(Nx,Ny)
    
    integer bcR,bcL,bcU,bcD,direction
    integer myid,myid1,the_id,the_id2
    integer myidx,myidy,the_idx,the_idy
    integer numprocs, namelen, rc,ierr,status(MPI_STATUS_SIZE),myid0
    character * (MPI_MAX_PROCESSOR_NAME) processor_name

    end module com

    !*****************************************************************************************************


    ! Double Mach reflection
    module init4

    use com

    contains

    function rho(x,y,z)

    real(8) x,y,r,r0,z
    real(8) rho

    if (x < 1d0/6d0 + y/3d0**0.5) then
        rho = 8
    else
        rho = 1.4
    end if

    end function rho


    function v1(x,y,z)

    real(8) x,y,z
    real(8) v1
    parameter(pi = 4*atan(1.0d0))

    if (x < 1d0/6d0 + y/3d0**0.5) then
        v1 = 8.25*cos(pi/6d0)
    else
        v1 = 0
    end if

    end function v1


    function v2(x,y,z)

    real(8) x,y,z
    real(8) v2
    parameter(pi = 4*atan(1.0d0))

    if (x < 1d0/6d0 + y/3d0**0.5) then
        v2 = -8.25*sin(pi/6d0)
    else
        v2 = 0
    end if

    end function v2


    function p(x,y,z)

    real(8) x,y,z
    real(8) p

    if (x < 1d0/6d0 + y/3d0**0.5) then
        p = 116.5
    else
        p = 1
    end if

    end function p

    subroutine mesh

    use com

    xa = 0
    xb = 3.2
    ya = 0
    yb = 1

    bcR = 2
    bcL = 2
    bcU = 3
    bcD = 3

    tend = 0.2

    M = 1
    beta = 1

    end subroutine mesh

    end module init4


    !*****************************************************************************************************
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
        call Euler_Forward
    else if (RKorder == 3) then
        call RK3
    else if (RKorder == 4) then
        call RK4
    end if

    call set_bc

    call save_solution
    ! ????????????
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

    ! The initial value:
    ! 4: Double Mach reflection
    
    use init4

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

    ! L2 Pro for Uh
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

    real uhsave(NumEq),jumpa
    integer the_idx1,the_idy1

    !uh(1:Nx,1:Ny,0:Nphi,1,4) = Is_trouble_cell

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

    !do d = 1,dimPk
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
    !end do

     ! test begin
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
    !end do
 !test over
 
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


        end do

    end do


    !basis 1 
    phiGR_ver(1,1) = 1
    phiGL_ver(1,1) = 1
    phiGU_ver(1,1) = 1
    phiGD_ver(1,1) = 1

    phiGR_ver(2,1) = 1
    phiGL_ver(2,1) = 1
    phiGU_ver(2,1) = 1
    phiGD_ver(2,1) = 1

    phiGR_ver_derx(1,1) = 0
    phiGL_ver_derx(1,1) = 0
    phiGU_ver_derx(1,1) = 0
    phiGD_ver_derx(1,1) = 0

    phiGR_ver_derx(2,1) = 0
    phiGL_ver_derx(2,1) = 0
    phiGU_ver_derx(2,1) = 0
    phiGD_ver_derx(2,1) = 0

    phiGR_ver_dery(1,1) = 0
    phiGL_ver_dery(1,1) = 0
    phiGU_ver_dery(1,1) = 0
    phiGD_ver_dery(1,1) = 0

    phiGR_ver_dery(2,1) = 0
    phiGL_ver_dery(2,1) = 0
    phiGU_ver_dery(2,1) = 0
    phiGD_ver_dery(2,1) = 0

    !basis x
    phiGR_ver(1,2) = 1
    phiGL_ver(1,2) = -1
    phiGU_ver(1,2) = -1
    phiGD_ver(1,2) = -1

    phiGR_ver(2,2) = 1
    phiGL_ver(2,2) = -1
    phiGU_ver(2,2) = 1
    phiGD_ver(2,2) = 1

    phiGR_ver_derx(1,2) = 1
    phiGL_ver_derx(1,2) = 1
    phiGU_ver_derx(1,2) = 1
    phiGD_ver_derx(1,2) = 1

    phiGR_ver_derx(2,2) = 1
    phiGL_ver_derx(2,2) = 1
    phiGU_ver_derx(2,2) = 1
    phiGD_ver_derx(2,2) = 1

    phiGR_ver_dery(1,2) = 0
    phiGL_ver_dery(1,2) =  0
    phiGU_ver_dery(1,2) =  0
    phiGD_ver_dery(1,2) =  0

    phiGR_ver_dery(2,2) = 0
    phiGL_ver_dery(2,2) =  0
    phiGU_ver_dery(2,2) = 0
    phiGD_ver_dery(2,2) = 0

    !basis y
    phiGR_ver(1,3) = 1
    phiGL_ver(1,3) = 1
    phiGU_ver(1,3) = 1
    phiGD_ver(1,3) = -1

    phiGR_ver(2,3) = -1
    phiGL_ver(2,3) = -1
    phiGU_ver(2,3) = 1
    phiGD_ver(2,3) = -1

    phiGR_ver_derx(1,3) = 0
    phiGL_ver_derx(1,3) = 0
    phiGU_ver_derx(1,3) = 0
    phiGD_ver_derx(1,3) = 0

    phiGR_ver_derx(2,3) = 0
    phiGL_ver_derx(2,3) = 0
    phiGU_ver_derx(2,3) = 0
    phiGD_ver_derx(2,3) = 0

    phiGR_ver_dery(1,3) = 1
    phiGL_ver_dery(1,3) = 1
    phiGU_ver_dery(1,3) = 1
    phiGD_ver_dery(1,3) = 1

    phiGR_ver_dery(2,3) = 1
    phiGL_ver_dery(2,3) = 1
    phiGU_ver_dery(2,3) = 1
    phiGD_ver_dery(2,3) = 1

    end subroutine get_basis

    !*****************************************************************************************************

    subroutine calculate_L2_Error

    use com

    use init4

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

    if (myid1 == 1) then
        print *,"The L2 Error:"
        print *,"rho    :",L2(1)
        print *,"rho u1 :",L2(2)
        print *,"rho u2 :",L2(3)
        print *,"E      :",L2(4)
    end if

      ! The value of num solution on the GL points
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
            end do ! ????????do k = 0,Nphi
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

    !*****************************************************************************************************
    subroutine apply_jump_filter_limiter

    use com

    if (limitertypeall == 1) then
        call jumpfilter
        call pp_Limiter
    elseif (limitertypeall == 2) then
        call hybridjumpfilter
        call pp_Limiter
    else
        if (myid1 == 1) then
            print *, 'Unsupported jump filter limiter type:', limitertypeall
        end if
        stop
    end if

    end subroutine apply_jump_filter_limiter

    !*****************************************************************************************************

    subroutine Euler_Forward

    use com

    CFL = 0.01
    t = 0

    call apply_jump_filter_limiter
    call calculate_umax

    if (myid1 == 1) then
        print *,t,umax
    end if

    do while (t < tend)

    call calculate_dt

    if (t + dt > tend) then
        dt = tend - t
        t = tend
    else
        t = t + dt
    end if

    tRK = t
    call Lh

    uh = uh + dt*du

    call apply_jump_filter_limiter

    call calculate_umax

    if (myid1 == 1) then
        print *,t,umax
    end if

    end do

    end subroutine Euler_Forward

    !*****************************************************************************************************

    subroutine RK3

    use com

    CFL = 0.1
    t = 0

    call apply_jump_filter_limiter
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

    call apply_jump_filter_limiter

    ! Stage II
    tRK = tRK + dt
    call Lh

    uII = (3d0/4d0)*uh00 + (1d0/4d0)*uh + (1d0/4d0)*dt*du

    uh = uII

    call apply_jump_filter_limiter

    ! Stage III
    tRK = tRK - 0.5*dt
    call Lh

    uh = (1d0/3d0)*uh00 + (2d0/3d0)*uh + (2d0/3d0)*dt*du

    call apply_jump_filter_limiter

    call calculate_umax

    if (myid1 == 1) then
        print *,t,umax
    end if

    end do

    end subroutine RK3

    !*****************************************************************************************************

    subroutine RK4

    use com

    CFL = 0.75
    t = 0

    !%%% ????
    call apply_jump_filter_limiter
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

    ! Stage I
    do i = 1,5

    call Lh

    uI = uh + (dt/6d0)*du

    tRK = tRK + (dt/6d0)

    uh = uI
    !%%% ????
    call apply_jump_filter_limiter
    end do

    uII = 0.04d0*uII + 0.36d0*uI

    uI = 15*uII - 5*uI

    uh = uI

    tRK = tRK - 0.5*dt

    ! Stage II
    do i = 6,9

    call Lh

    uI = uh + (dt/6d0)*du

    tRK = tRK + dt/6d0

    uh = uI
    !%%% ????
    call apply_jump_filter_limiter
    end do

    call Lh

    uh = uII + 0.6d0*uI + (dt/10d0)*du
    !%%% ????
    call apply_jump_filter_limiter
    call calculate_umax

    if (myid1 == 1) then
        !print *,t,umax,sum(Is_trouble_cell)
        !print *,t,umax 
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

                    ! The Uh
                    ! The Right condition
                    if (i == Nx_process) then

                    if (bcR == 1) then ! periodic
                        the_idx = 1
                        the_idy = j
                        the_id2 = the_idx + Nx_process*(the_idy - 1)
                        if (myid1 == the_id2) then
                            call MPI_SEND(uh(1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,ierr)
                        end if
                        if (myid1 == the_id) then
                            call MPI_RECV(uh(Nx1,0:Ny1,k,d,n),Ny + 2,MPI_REAL8,the_id2 - 1,1,MPI_COMM_WORLD,status,ierr)
                        end if
                    else if (bcR == 2) then ! outflow
                        if (myid1 == the_id) then
                            uh(Nx1,0:Ny1,k,d,n) = uh(Nx,0:Ny1,k,d,n)
                        end if
                    else if (bcR == 5) then!Pure Wall boundary condition
                        if (myid1 == the_id) then
                            do ii = 1,Ny
                                call evenex_y(uh(Nx1,ii,k,:,1),uh(Nx,ii,k,:,1))
                                !call evenex_y(uh(ii,0,k,:,2),uh(ii,1,k,:,2))
                                !call oddex_y(uh(ii,0,k,:,3),uh(ii,1,k,:,3))
                                call oddex_y(uh(Nx1,ii,k,:,2),uh(Nx,ii,k,:,2))
                                call evenex_y(uh(Nx1,ii,k,:,3),uh(Nx,ii,k,:,3))
                                call evenex_y(uh(Nx1,ii,k,:,4),uh(Nx,ii,k,:,4))
                            end do
                        end if       
                    else if (bcR == 4) then ! r-periodic
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



                    ! The Left condition
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
                    else if (bcL == 5) then!Pure Wall boundary condition
                        if (myid1 == the_id) then
                            do ii = 1,Ny
                                call evenex_y(uh(0,ii,k,:,1),uh(1,ii,k,:,1))
                                !call evenex_y(uh(ii,0,k,:,2),uh(ii,1,k,:,2))
                                !call oddex_y(uh(ii,0,k,:,3),uh(ii,1,k,:,3))
                                call oddex_y(uh(0,ii,k,:,2),uh(1,ii,k,:,2))
                                call evenex_y(uh(0,ii,k,:,3),uh(1,ii,k,:,3))
                                call evenex_y(uh(0,ii,k,:,4),uh(1,ii,k,:,4))
                            end do
                        end if

                    else if (bcL == 3) then
                        if (myid1 == the_id) then
                            uh(0,:,k,:,:) = 0
                            do jj = 1,Ny
                                if (Yc(jj) < 0.05 .and. Yc(jj) > -0.05) then
                                    uh(0,jj,k,1,1) = 5
                                    uh(0,jj,k,1,2) = 5*800
                                    uh(0,jj,k,1,3) = 0
                                    uh(0,jj,k,1,4) = 0.4127/gamma1 + 0.5*5*800**2
                                else
                                    uh(0,jj,k,1,1) = 0.5
                                    uh(0,jj,k,1,2) = 0
                                    uh(0,jj,k,1,3) = 0
                                    uh(0,jj,k,1,4) = 0.4127/gamma1
                                end if
                            end do
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



                    ! The Up condition
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
                    else if (bcU == 2) then
                        if (myid1 == the_id) then
                            uh(0:Nx1,Ny1,k,d,n) = uh(0:Nx1,Ny,k,d,n)
                        end if
                    else if (bcU == 5) then ! Pure Wall boundary condition
                        if (myid1 == the_id) then
                            do ii = 1,Nx
                                call evenex_y(uh(ii,Ny1,k,:,1),uh(ii,Ny,k,:,1))
                                call evenex_y(uh(ii,Ny1,k,:,2),uh(ii,Ny,k,:,2))
                                call oddex_y(uh(ii,Ny1,k,:,3),uh(ii,Ny,k,:,3))
                                call evenex_y(uh(ii,Ny1,k,:,4),uh(ii,Ny,k,:,4))

                            end do
                        end if        
                    else if (bcU == 3) then
                        if (myid1 == the_id) then
                            uh(:,Ny1,:,:,:) = 0
                            do ii = 1,Nx
                                if (Xc(ii) < 1d0/6d0 + (1 + 20*tRK)/3d0**0.5) then ! post-shock
                                    uh(ii,Ny1,k,1,1) = 8
                                    uh(ii,Ny1,k,1,2) = 8*8.25*cos(pi/6d0)
                                    uh(ii,Ny1,k,1,3) = -8*8.25*sin(pi/6d0)
                                    uh(ii,Ny1,k,1,4) = 116.5d0/gamma1 + 0.5*8*((8.25*cos(pi/6d0))**2 + (8.25*sin(pi/6d0))**2)
                                else ! pre-shock
                                    uh(ii,Ny1,k,1,1) = 1.4
                                    uh(ii,Ny1,k,1,2) = 0
                                    uh(ii,Ny1,k,1,3) = 0
                                    uh(ii,Ny1,k,1,4) = 1d0/gamma1
                                end if
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



                    ! The Down condition
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
                            do ii = 1,Nx
                                if (Xc(ii) < 1d0/6d0) then
                                    uh(ii,0,:,:,:) = uh(ii,1,:,:,:)
                                else
                                    call evenex_y(uh(ii,0,k,:,1),uh(ii,1,k,:,1))
                                    call evenex_y(uh(ii,0,k,:,2),uh(ii,1,k,:,2))
                                    call oddex_y(uh(ii,0,k,:,3),uh(ii,1,k,:,3))
                                    call evenex_y(uh(ii,0,k,:,4),uh(ii,1,k,:,4))
                                end if
                            end do
                        end if
                    else if (bcD == 4) then ! Pure wall boundary
                        if (myid1 == the_id) then
                            do ii = 1,Nx
                                call evenex_y(uh(ii,0,k,:,1),uh(ii,1,k,:,1))
                                call evenex_y(uh(ii,0,k,:,2),uh(ii,1,k,:,2))
                                call oddex_y(uh(ii,0,k,:,3),uh(ii,1,k,:,3))
                                call evenex_y(uh(ii,0,k,:,4),uh(ii,1,k,:,4))
                            end do
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

    ! calculate the Volume integral
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

    ! calculate the Numerical flux

    ! The x-flux
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

    ! The y-Flux
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

    ! calculate Fx hat
    do j1 = 1,NumGLP
        do k = 0,Nphi
            do j = 1,Ny
                do i = 0,Nx
                    call eigenvalueMm(SRmax,SRmin,UR(i,j,k,j1,1),UR(i,j,k,j1,2),UR(i,j,k,j1,3),UR(i,j,k,j1,4),1,0)
                    call eigenvalueMm(SLmax,SLmin,UL(i + 1,j,k,j1,1),UL(i + 1,j,k,j1,2),UL(i + 1,j,k,j1,3),UL(i + 1,j,k,j1,4),1,0)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
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

    ! calculate Fy hat
    do i1 = 1,NumGLP
        do k = 0,Nphi
            do j = 0,Ny
                do i = 1,Nx
                    call eigenvalueMm(SRmax,SRmin,UU(i,j,k,i1,1),UU(i,j,k,i1,2),UU(i,j,k,i1,3),UU(i,j,k,i1,4),0,1)
                    call eigenvalueMm(SLmax,SLmin,UD(i,j + 1,k,i1,1),UD(i,j + 1,k,i1,2),UD(i,j + 1,k,i1,3),UD(i,j + 1,k,i1,4),0,1)
                    !SR = 0.3*min(SRmax,SLmax)
                    !SL = 0.3*max(SRmin,SLmin)
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

    ! calculate the Surface integral
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

    !Rinv = 0
    !Rmat = 0
    !do i = 1,4
    !    Rinv(i,i) = 1
    !    Rmat(i,i) = 1
    !end do

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

    !*****************************************************************************************************

    subroutine pp_Limiter

    use com

    real(8) eta,epsilon,eta1,pbar,pq,sq(NumEq),norm1,norm2,tq
    real(8) uhGLL(NumGLP,NumGLP,NumEq,2)

    epsilon = 1e-13

    ! Limiting the density
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                eta = 1
                uhGLL = 0
                do d = 1,dimPk
                    do n = 1,NumEq
                        uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,k,d,n)*phiGLL(:,:,d,:)
                    end do
                end do

                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        do d = 1,2
                            ! eta1 = (rhobar - epsilon)/(rhobar - rho(xq))
                            eta1 = abs((uh(i,j,k,1,1) - epsilon)/(uh(i,j,k,1,1) - uhGLL(i1,j1,1,d)))
                            if (eta1 < 1) then
                                eta = eta1
                            end if
                        end do
                    end do
                end do

                if (eta < 1) then
                    uh(i,j,k,2:dimPk,1) = 0.9*eta*uh(i,j,k,2:dimPk,1)
                end if
            end do
        end do
    end do

    ! Limiting the pressure
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                eta = 1
                eta1 = 1
                uhGLL = 0
                do d = 1,dimPk
                    do n = 1,NumEq
                        uhGLL(:,:,n,:) = uhGLL(:,:,n,:) + uh(i,j,k,d,n)*phiGLL(:,:,d,:)
                    end do
                end do
                pbar = pressure(uh(i,j,k,1,1),uh(i,j,k,1,2),uh(i,j,k,1,3),uh(i,j,k,1,4),gamma)
                do i1 = 1,NumGLP
                    do j1 = 1,NumGLP
                        do d = 1,2
                            pq = pressure(uhGLL(i1,j1,1,d),uhGLL(i1,j1,2,d),uhGLL(i1,j1,3,d),uhGLL(i1,j1,4,d),gamma)

                            if (pq < 0) then
                                call calculate_tq(uh(i,j,k,1,:),uhGLL(i1,j1,:,d),tq,gamma)
                                sq = tq*uhGLL(i1,j1,:,d) + (1 - tq)*uh(i,j,k,1,:)
                                call norm(sq - uh(i,j,k,1,:),norm1)
                                call norm(uhGLL(i1,j1,:,d) - uh(i,j,k,1,:),norm2)
                                eta1 = norm1/norm2

                                !eta1 = pbar/(pbar - pq)
                            end if

                            if (eta1 < eta) then
                                eta = eta1
                            end if

                        end do
                    end do
                end do

                if (eta < 1) then
                    eta = 0.9*eta
                end if

                uh(i,j,k,2:dimPk,:) = eta*uh(i,j,k,2:dimPk,:)
            end do
        end do
    end do

    end subroutine pp_Limiter


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

   ! a(2) = b(2)
   ! a(3) = -b(3)
      a(2) = 0d0
       a(3) = 0d0

    end subroutine evenex_y

    !*****************************************************************************************************

    subroutine oddex_y(a,b)

    real a(10),b(10)

    a(1) = -b(1)

    !a(2) = -b(2)
    !a(3) = b(3)
  a(2) = 0d0
       a(3) = 0d0

    end subroutine oddex_y
 
    !*****************************************************************************************************
    !*****************************************************************************************************
    !*****************************************************************************************************



    subroutine writetroubledcells
    ! ????????
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%% ????troubled cells?????
    ! ??troubled cells??????????????????
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



   

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subroutine   jumpfilter 

    use com 

    real(kind=8) rhoLinfty,rhouLinfty,rhovLinfty,EnerLinfty,rhoij,entropyij,enthayij,label_density(Nx,Ny),rhouij,label_rhou(Nx,Ny),rhovij,label_rhov(Nx,Ny),Enerij,label_Ener(Nx,Ny)
    real(kind=8) rhoLinfty1,rhouLinfty1,rhovLinfty1,EnerLinfty1
    real(kind=8) UR_ver(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver(Nx1,Ny,0:Nphi,2,NumEq),UU_ver(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver(Nx,Ny1,0:Nphi,2,NumEq) 
    real(kind=8) UR_ver_derx(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derx(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_derx(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_derx(Nx,Ny1,0:Nphi,2,NumEq)
    real(kind=8) UR_ver_dery(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_dery(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_dery(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_dery(Nx,Ny1,0:Nphi,2,NumEq) 
    integer ploydeg ,ss 
    real(kind=8) scal,u1ij,u2ij,Eij,pressureij,cij,betai,betaj 
    real(kind=8) left_edge_jump_density ,right_edge_jump_density,bottom_edge_jump_density,top_edge_jump_density, deltadensity0 
    real(kind=8) left_edge_jump_rhou,right_edge_jump_rhou,bottom_edge_jump_rhou,top_edge_jump_rhou,deltarhou0
    real(kind=8) left_edge_jump_rhov,right_edge_jump_rhov,bottom_edge_jump_rhov,top_edge_jump_rhov,deltarhov0            
    real(kind=8) left_edge_jump_Ener,right_edge_jump_Ener,bottom_edge_jump_Ener,top_edge_jump_Ener,deltaEner0  
    real(kind=8)  delta0max , delta1max , damping
    real(kind=8)  left_edge_jumpderx_density, right_edge_jumpderx_density, bottom_edge_jumpderx_density, top_edge_jumpderx_density 
    real(kind=8)  left_edge_jumpdery_density, right_edge_jumpdery_density, bottom_edge_jumpdery_density, top_edge_jumpdery_density 
    real(kind=8)  deltadensity1,deltarhou1,deltarhov1,deltaEner1
    real(kind=8) left_edge_jumpderx_rhou,right_edge_jumpderx_rhou,bottom_edge_jumpderx_rhou, top_edge_jumpderx_rhou 
    real(kind=8) left_edge_jumpdery_rhou,right_edge_jumpdery_rhou,bottom_edge_jumpdery_rhou,top_edge_jumpdery_rhou 
    real(kind=8) left_edge_jumpderx_rhov,right_edge_jumpderx_rhov,bottom_edge_jumpderx_rhov, top_edge_jumpderx_rhov 
    real(kind=8) left_edge_jumpdery_rhov,right_edge_jumpdery_rhov,bottom_edge_jumpdery_rhov,top_edge_jumpdery_rhov       
    real(kind=8) left_edge_jumpderx_Ener,right_edge_jumpderx_Ener,bottom_edge_jumpderx_Ener, top_edge_jumpderx_Ener 
    real(kind=8) left_edge_jumpdery_Ener,right_edge_jumpdery_Ener,bottom_edge_jumpdery_Ener,top_edge_jumpdery_Ener 
    real(kind=8) rhofor_x,rhoback_x,rhofor_y,rhoback_y,Linftyrhox,Linftyrhoy
    real(kind=8) rhoufor_x,rhouback_x,rhoufor_y,rhouback_y,Linftyrhoux,Linftyrhouy
    real(kind=8) rhovfor_x,rhovback_x,rhovfor_y,rhovback_y,Linftyrhovx,Linftyrhovy
    real(kind=8) Enerfor_x,Enerback_x,Enerfor_y,Enerback_y,LinftyEnerx,LinftyEnery

    real(kind=8) rhomax,rhomin,rhoumax,rhoumin,rhovmax,rhovmin,Enermax,Enermin
    real(kind=8) Machij
  
    real(kind=8) uhmod(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)


    !%%%%% Step1 : betai,betaj,some solution value, jump filter uses local scaling
    ! ?????????
    call set_bc

    uhmod = uh
    ! UL UR UU UB ?????????????????
    ! The x-direction
    UR_ver = 0
    UL_ver = 0
    UR_ver_derx = 0
    UL_ver_derx = 0
    UR_ver_dery = 0
    UL_ver_dery = 0
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR_ver(i,j,k,:,n) = UR_ver(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver(:,d)
                        UL_ver(i + 1,j,k,:,n) = UL_ver(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver(:,d)
                        !phiGR_ver_derx(2,3),phiGL_ver_derx(2,3), phiGU_ver_derx(2,3), phiGD_ver_derx(2,3)
                        UR_ver_derx(i,j,k,:,n) = UR_ver_derx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derx(:,d)
                        UL_ver_derx(i + 1,j,k,:,n) = UL_ver_derx(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derx(:,d)

                        UR_ver_dery(i,j,k,:,n) = UR_ver_dery(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_dery(:,d)
                        UL_ver_dery(i + 1,j,k,:,n) = UL_ver_dery(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_dery(:,d)

                    end do
                end do
            end do
        end do
    end do
    UR_ver_derx =  UR_ver_derx*2/hx
    UL_ver_derx =  UL_ver_derx*2/hx
    UR_ver_dery =  UR_ver_dery*2/hy
    UL_ver_dery =  UL_ver_dery*2/hy
    ! The y-direction
    UU_ver = 0
    UD_ver = 0
    UU_ver_derx = 0
    UD_ver_derx = 0
    UU_ver_dery = 0
    UD_ver_dery = 0
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
                    end do
                end do
            end do
        end do
    end do
    UU_ver_derx =  UU_ver_derx*2/hx
    UD_ver_derx =  UD_ver_derx*2/hx
    UU_ver_dery =  UU_ver_dery*2/hy
    UD_ver_dery =  UD_ver_dery*2/hy

    !%%%%% Step2 : calculate the jump of each componment 
    ploydeg = 1 
    !scal = 1.0d0 ???0.1
    !% moment0 == ??????????
    outerj: do i = 1,Nx
        outeri : do j = 1,Ny
            outerk :  do k = 0,Nphi
            scal = 0.5d0
            rhoij = uh(i,j,k,1,1)
            u1ij = uh(i,j,k,1,2)/rhoij
            u2ij = uh(i,j,k,1,3)/rhoij
            Eij = uh(i,j,k,1,4)
            pressureij = gamma1*(Eij - 0.5d0*rhoij*(u1ij**2+u2ij**2))
            !entropyij =  pressureij/rhoij**gamma
            enthayij  =   ( Eij + pressureij ) /rhoij
            cij = sqrt(abs(gamma*pressureij/rhoij))
            betai = abs(u1ij) + cij 
            betaj = abs(u2ij) + cij

            Machij  = sqrt(u1ij**2 + u2ij**2)/cij
            !if (Machij<1) then
            !    scal = Machij
            !else
            !    scal = 1/Machij
            !end if
            scal = scal/enthayij
 
        
            ! calculate_jump_of zero moment density
            left_edge_jump_density=0
            right_edge_jump_density =0
            bottom_edge_jump_density=0
            top_edge_jump_density=0
            do ss = 1,2
                ! ???density
                left_edge_jump_density = left_edge_jump_density  + abs(UL_ver(i,j,0,ss,1)-UR_ver(i-1,j,0,ss,1)) 
                right_edge_jump_density = right_edge_jump_density  + abs(UL_ver(i+1,j,0,ss,1)-UR_ver(i,j,0,ss,1)) 
                bottom_edge_jump_density = bottom_edge_jump_density + abs(UU_ver(i,j-1,0,ss,1) - UD_ver(i,j,0,ss,1))
                top_edge_jump_density = top_edge_jump_density + abs(UU_ver(i,j,0,ss,1) - UD_ver(i,j+1,0,ss,1))

            end do

            deltadensity0 =  betai*(left_edge_jump_density+right_edge_jump_density) + &
            betaj*(bottom_edge_jump_density+top_edge_jump_density)  


            !  calculate_jump_of zero moment momentumn1
            left_edge_jump_rhou=0
            right_edge_jump_rhou =0
            bottom_edge_jump_rhou=0
            top_edge_jump_rhou=0
            do ss = 1,2
                ! momentumn 1
                left_edge_jump_rhou = left_edge_jump_rhou  + abs(UL_ver(i,j,0,ss,2)-UR_ver(i-1,j,0,ss,2)) 
                right_edge_jump_rhou = right_edge_jump_rhou  + abs(UL_ver(i+1,j,0,ss,2)-UR_ver(i,j,0,ss,2)) 
                bottom_edge_jump_rhou = bottom_edge_jump_rhou + abs(UU_ver(i,j-1,0,ss,2) - UD_ver(i,j,0,ss,2))
                top_edge_jump_rhou = top_edge_jump_rhou+ abs(UU_ver(i,j,0,ss,2) - UD_ver(i,j+1,0,ss,2))

            end do

            deltarhou0 =  betai*(left_edge_jump_rhou+right_edge_jump_rhou)  + &
            betaj*(bottom_edge_jump_rhou+top_edge_jump_rhou)  

            ! calculate_jump_of zero moment momentumn2
            left_edge_jump_rhov=0
            right_edge_jump_rhov =0
            bottom_edge_jump_rhov=0
            top_edge_jump_rhov=0
            do ss = 1,2
                ! momentumn 2
                left_edge_jump_rhov = left_edge_jump_rhov  + abs(UL_ver(i,j,0,ss,3)- UR_ver(i-1,j,0,ss,3)) 
                right_edge_jump_rhov = right_edge_jump_rhov  + abs(UL_ver(i+1,j,0,ss,3)-UR_ver(i,j,0,ss,3)) 
                bottom_edge_jump_rhov = bottom_edge_jump_rhov + abs(UU_ver(i,j-1,0,ss,3) - UD_ver(i,j,0,ss,3))
                top_edge_jump_rhov = top_edge_jump_rhov + abs(UU_ver(i,j,0,ss,3) - UD_ver(i,j+1,0,ss,3))

            end do

            deltarhov0 =  betai*(left_edge_jump_rhov+right_edge_jump_rhov)  + &
            betaj*(bottom_edge_jump_rhov+top_edge_jump_rhov) 

            ! call calculate_jump_of zero moment Ener
            left_edge_jump_Ener=0
            right_edge_jump_Ener =0
            bottom_edge_jump_Ener=0
            top_edge_jump_Ener=0
            do ss = 1,2
                ! Ener
                left_edge_jump_Ener= left_edge_jump_Ener  + abs(UL_ver(i,j,0,ss,4)-UR_ver(i-1,j,0,ss,4)) 
                right_edge_jump_Ener= right_edge_jump_Ener  + abs(UL_ver(i+1,j,0,ss,4)-UR_ver(i,j,0,ss,4)) 
                bottom_edge_jump_Ener= bottom_edge_jump_Ener + abs(UU_ver(i,j-1,0,ss,4) - UD_ver(i,j,0,ss,4))
                top_edge_jump_Ener = top_edge_jump_Ener+ abs(UU_ver(i,j,0,ss,4) - UD_ver(i,j+1,0,ss,4))

            end do

            deltaEner0 =  betai*(left_edge_jump_Ener+right_edge_jump_Ener) + &
            betaj*(bottom_edge_jump_Ener+top_edge_jump_Ener) 

            !! 0??moment????§³ 
            delta0max = max( deltadensity0 ,  deltarhou0 ,  deltarhov0 ,  deltaEner0)

            !!!!!!!!!! ???moment 
            ! calculate_jump_of one moment density
            left_edge_jumpderx_density=0
            right_edge_jumpderx_density =0
            bottom_edge_jumpderx_density=0
            top_edge_jumpderx_density=0

            left_edge_jumpdery_density=0
            right_edge_jumpdery_density =0
            bottom_edge_jumpdery_density=0
            top_edge_jumpdery_density=0
            do ss = 1,2
                ! ???density   UR_ver_derx
                left_edge_jumpderx_density = left_edge_jumpderx_density  + abs(UL_ver_derx(i,j,0,ss,1)-UR_ver_derx(i-1,j,0,ss,1)) 
                right_edge_jumpderx_density = right_edge_jumpderx_density  + abs(UL_ver_derx(i+1,j,0,ss,1)-UR_ver_derx(i,j,0,ss,1)) 
                bottom_edge_jumpderx_density = bottom_edge_jumpderx_density + abs(UU_ver_derx(i,j-1,0,ss,1) - UD_ver_derx(i,j,0,ss,1))
                top_edge_jumpderx_density = top_edge_jumpderx_density + abs(UU_ver_derx(i,j,0,ss,1) - UD_ver_derx(i,j+1,0,ss,1))
                ! ???density   UR_ver_dery
                left_edge_jumpdery_density = left_edge_jumpdery_density  + abs(UL_ver_dery(i,j,0,ss,1)-UR_ver_dery(i-1,j,0,ss,1)) 
                right_edge_jumpdery_density = right_edge_jumpdery_density  + abs(UL_ver_dery(i+1,j,0,ss,1)-UR_ver_dery(i,j,0,ss,1)) 
                bottom_edge_jumpdery_density = bottom_edge_jumpdery_density + abs(UU_ver_dery(i,j-1,0,ss,1) - UD_ver_dery(i,j,0,ss,1))
                top_edge_jumpdery_density = top_edge_jumpdery_density + abs(UU_ver_dery(i,j,0,ss,1) - UD_ver_dery(i,j+1,0,ss,1))

            end do
deltadensity1 =  betai*(left_edge_jumpderx_density+right_edge_jumpderx_density + left_edge_jumpdery_density+right_edge_jumpdery_density)*hx*2 + &
            betaj*(bottom_edge_jumpderx_density + top_edge_jumpdery_density + bottom_edge_jumpderx_density + top_edge_jumpdery_density )*hy*2  

            ! calculate_jump_of one moment rhou
            left_edge_jumpderx_rhou=0
            right_edge_jumpderx_rhou =0
            bottom_edge_jumpderx_rhou=0
            top_edge_jumpderx_rhou=0

            left_edge_jumpdery_rhou=0
            right_edge_jumpdery_rhou =0
            bottom_edge_jumpdery_rhou=0
            top_edge_jumpdery_rhou=0
            do ss = 1,2
                ! ???? rhou   derx
                left_edge_jumpderx_rhou = left_edge_jumpderx_rhou  + abs(UL_ver_derx(i,j,0,ss,2)-UR_ver_derx(i-1,j,0,ss,2)) 
                right_edge_jumpderx_rhou = right_edge_jumpderx_rhou  + abs(UL_ver_derx(i+1,j,0,ss,2)-UR_ver_derx(i,j,0,ss,2)) 
                bottom_edge_jumpderx_rhou = bottom_edge_jumpderx_rhou + abs(UU_ver_derx(i,j-1,0,ss,2) - UD_ver_derx(i,j,0,ss,2))
                top_edge_jumpderx_rhou = top_edge_jumpderx_rhou + abs(UU_ver_derx(i,j,0,ss,2) - UD_ver_derx(i,j+1,0,ss,2))
                ! ???? rhou dery
                left_edge_jumpdery_rhou = left_edge_jumpdery_rhou  + abs(UL_ver_dery(i,j,0,ss,2)-UR_ver_dery(i-1,j,0,ss,2)) 
                right_edge_jumpdery_rhou = right_edge_jumpdery_rhou  + abs(UL_ver_dery(i+1,j,0,ss,2)-UR_ver_dery(i,j,0,ss,2)) 
                bottom_edge_jumpdery_rhou = bottom_edge_jumpdery_rhou + abs(UU_ver_dery(i,j-1,0,ss,2) - UD_ver_dery(i,j,0,ss,2))
                top_edge_jumpdery_rhou = top_edge_jumpdery_rhou + abs(UU_ver_dery(i,j,0,ss,2) - UD_ver_dery(i,j+1,0,ss,2))

            end do
deltarhou1 =  betai*(left_edge_jumpderx_rhou +right_edge_jumpderx_rhou + left_edge_jumpdery_rhou+right_edge_jumpdery_rhou)*hx*2 + &
            betaj*(bottom_edge_jumpderx_rhou + top_edge_jumpderx_rhou + bottom_edge_jumpdery_rhou + top_edge_jumpdery_rhou )*hy*2

            ! calculate_jump_of one moment rhov
            left_edge_jumpderx_rhov=0
            right_edge_jumpderx_rhov =0
            bottom_edge_jumpderx_rhov=0
            top_edge_jumpderx_rhov=0

            left_edge_jumpdery_rhov=0
            right_edge_jumpdery_rhov =0
            bottom_edge_jumpdery_rhov=0
            top_edge_jumpdery_rhov=0
            do ss = 1,2
                ! ???? rhou   
                left_edge_jumpderx_rhov = left_edge_jumpderx_rhov  + abs(UL_ver_derx(i,j,0,ss,3)-UR_ver_derx(i-1,j,0,ss,3)) 
                right_edge_jumpderx_rhov = right_edge_jumpderx_rhov  + abs(UL_ver_derx(i+1,j,0,ss,3)-UR_ver_derx(i,j,0,ss,3)) 
                bottom_edge_jumpderx_rhov = bottom_edge_jumpderx_rhov + abs(UU_ver_derx(i,j-1,0,ss,3) - UD_ver_derx(i,j,0,ss,3))
                top_edge_jumpderx_rhov = top_edge_jumpderx_rhov + abs(UU_ver_derx(i,j,0,ss,3) - UD_ver_derx(i,j+1,0,ss,3))
                ! ???? rhou
                left_edge_jumpdery_rhov = left_edge_jumpdery_rhov+ abs(UL_ver_dery(i,j,0,ss,3)-UR_ver_dery(i-1,j,0,ss,3)) 
                right_edge_jumpdery_rhov = right_edge_jumpdery_rhov + abs(UL_ver_dery(i+1,j,0,ss,3)-UR_ver_dery(i,j,0,ss,3)) 
                bottom_edge_jumpdery_rhov= bottom_edge_jumpdery_rhov + abs(UU_ver_dery(i,j-1,0,ss,3) - UD_ver_dery(i,j,0,ss,3))
                top_edge_jumpdery_rhov = top_edge_jumpdery_rhov + abs(UU_ver_dery(i,j,0,ss,3) - UD_ver_dery(i,j+1,0,ss,3))

            end do

            deltarhov1 =  betai*(left_edge_jumpderx_rhov +right_edge_jumpderx_rhov + left_edge_jumpdery_rhov+right_edge_jumpdery_rhov)*hx*2 + &
            betaj*(bottom_edge_jumpderx_rhov + top_edge_jumpderx_rhov + bottom_edge_jumpdery_rhov + top_edge_jumpdery_rhov )*hy*2 


            ! calculate_jump_of one moment Ener
            left_edge_jumpderx_Ener=0
            right_edge_jumpderx_Ener =0
            bottom_edge_jumpderx_Ener =0
            top_edge_jumpderx_Ener =0

            left_edge_jumpdery_Ener =0
            right_edge_jumpdery_Ener =0
            bottom_edge_jumpdery_Ener=0
            top_edge_jumpdery_Ener=0
            do ss = 1,2
                ! ???? Ener   
                left_edge_jumpderx_Ener = left_edge_jumpderx_Ener  + abs(UL_ver_derx(i,j,0,ss,4)-UR_ver_derx(i-1,j,0,ss,4)) 
                right_edge_jumpderx_Ener = right_edge_jumpderx_Ener  + abs(UL_ver_derx(i+1,j,0,ss,4)-UR_ver_derx(i,j,0,ss,4)) 
                bottom_edge_jumpderx_Ener = bottom_edge_jumpderx_Ener + abs(UU_ver_derx(i,j-1,0,ss,4) - UD_ver_derx(i,j,0,ss,4))
                top_edge_jumpderx_Ener = top_edge_jumpderx_Ener + abs(UU_ver_derx(i,j,0,ss,4) - UD_ver_derx(i,j+1,0,ss,4))
                ! ???? Ener  
                left_edge_jumpdery_Ener = left_edge_jumpdery_Ener + abs(UL_ver_dery(i,j,0,ss,4)-UR_ver_dery(i-1,j,0,ss,4)) 
                right_edge_jumpdery_Ener = right_edge_jumpdery_Ener + abs(UL_ver_dery(i+1,j,0,ss,4)-UR_ver_dery(i,j,0,ss,4)) 
                bottom_edge_jumpdery_Ener = bottom_edge_jumpdery_Ener + abs(UU_ver_dery(i,j-1,0,ss,4) - UD_ver_dery(i,j,0,ss,4))
                top_edge_jumpdery_Ener = top_edge_jumpdery_Ener + abs(UU_ver_dery(i,j,0,ss,4) - UD_ver_dery(i,j+1,0,ss,4))

            end do

            deltaEner1 =  betai*(left_edge_jumpderx_Ener +right_edge_jumpderx_Ener + left_edge_jumpdery_Ener+right_edge_jumpdery_Ener)*hx*2 + &
            betaj*(bottom_edge_jumpderx_Ener + top_edge_jumpderx_Ener + bottom_edge_jumpdery_Ener + top_edge_jumpdery_Ener)*hy*2 


            !! 1??moment????§³ 
            delta1max = max( deltadensity1 ,  deltarhou1 ,  deltarhov1 ,  deltaEner1)

            damping =   delta0max + delta1max 

            !damping = 0.1d0*hx*hy*damping*(8d0)/hx 

            damping =   scal*hx*damping/hx ! ???????1
            !damping =   scal*hx*damping*(8d0)/hx

            uhmod(i,j,k,2:3,1:4) = exp(-dt*damping)*uh(i,j,k,2:3,1:4) 

            end do outerk
        end do outeri 
    end do outerj 


    uh = uhmod


    end  subroutine jumpfilter

    !!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !! ???????jump filter?????????
    subroutine hybridjumpfilter  
    ! ???????????P_1^lim; ???§Ý??
    use com 

    real(kind=8) rhoLinfty,rhouLinfty,rhovLinfty,EnerLinfty,rhoij,entropyij,enthayij,label_density(Nx,Ny),rhouij,label_rhou(Nx,Ny),rhovij,label_rhov(Nx,Ny),Enerij,label_Ener(Nx,Ny)
    real(kind=8) rhoLinfty1,rhouLinfty1,rhovLinfty1,EnerLinfty1
    real(kind=8) UR_ver(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver(Nx1,Ny,0:Nphi,2,NumEq),UU_ver(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver(Nx,Ny1,0:Nphi,2,NumEq) 
    real(kind=8) UR_ver_derx(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_derx(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_derx(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_derx(Nx,Ny1,0:Nphi,2,NumEq)
    real(kind=8) UR_ver_dery(0:Nx,Ny,0:Nphi,2,NumEq),UL_ver_dery(Nx1,Ny,0:Nphi,2,NumEq),UU_ver_dery(Nx,0:Ny,0:Nphi,2,NumEq),UD_ver_dery(Nx,Ny1,0:Nphi,2,NumEq) 
    integer ploydeg ,ss 
    real(kind=8) scal,u1ij,u2ij,Eij,pressureij,cij,betai,betaj 
    real(kind=8) left_edge_jump_density ,right_edge_jump_density,bottom_edge_jump_density,top_edge_jump_density, deltadensity0 
    real(kind=8) left_edge_jump_rhou,right_edge_jump_rhou,bottom_edge_jump_rhou,top_edge_jump_rhou,deltarhou0
    real(kind=8) left_edge_jump_rhov,right_edge_jump_rhov,bottom_edge_jump_rhov,top_edge_jump_rhov,deltarhov0            
    real(kind=8) left_edge_jump_Ener,right_edge_jump_Ener,bottom_edge_jump_Ener,top_edge_jump_Ener,deltaEner0  
    real(kind=8)  delta0max , delta1max , damping
    real(kind=8)  left_edge_jumpderx_density, right_edge_jumpderx_density, bottom_edge_jumpderx_density, top_edge_jumpderx_density 
    real(kind=8)  left_edge_jumpdery_density, right_edge_jumpdery_density, bottom_edge_jumpdery_density, top_edge_jumpdery_density 
    real(kind=8)  deltadensity1,deltarhou1,deltarhov1,deltaEner1
    real(kind=8) left_edge_jumpderx_rhou,right_edge_jumpderx_rhou,bottom_edge_jumpderx_rhou, top_edge_jumpderx_rhou 
    real(kind=8) left_edge_jumpdery_rhou,right_edge_jumpdery_rhou,bottom_edge_jumpdery_rhou,top_edge_jumpdery_rhou 
    real(kind=8) left_edge_jumpderx_rhov,right_edge_jumpderx_rhov,bottom_edge_jumpderx_rhov, top_edge_jumpderx_rhov 
    real(kind=8) left_edge_jumpdery_rhov,right_edge_jumpdery_rhov,bottom_edge_jumpdery_rhov,top_edge_jumpdery_rhov       
    real(kind=8) left_edge_jumpderx_Ener,right_edge_jumpderx_Ener,bottom_edge_jumpderx_Ener, top_edge_jumpderx_Ener 
    real(kind=8) left_edge_jumpdery_Ener,right_edge_jumpdery_Ener,bottom_edge_jumpdery_Ener,top_edge_jumpdery_Ener 
    real(kind=8) rhofor_x,rhoback_x,rhofor_y,rhoback_y,Linftyrhox,Linftyrhoy
    real(kind=8) rhoufor_x,rhouback_x,rhoufor_y,rhouback_y,Linftyrhoux,Linftyrhouy
    real(kind=8) rhovfor_x,rhovback_x,rhovfor_y,rhovback_y,Linftyrhovx,Linftyrhovy
    real(kind=8) Enerfor_x,Enerback_x,Enerfor_y,Enerback_y,LinftyEnerx,LinftyEnery

    real(kind=8) rhomax,rhomin,rhoumax,rhoumin,rhovmax,rhovmin,Enermax,Enermin
    real(kind=8) Machij
    real(kind=8) uhmod(0:Nx1,0:Ny1,0:Nphi,dimPk,NumEq)
    ! ????????????????
    integer,parameter ::  indicator = 1, ishybrid = 1 ,hybridmethod =1
    integer,parameter ::   alpha1 = 1 ,order1 =4 ! hybridmethod =2
    integer,parameter ::   alpha2 = 1 ,order2 =2 ! hybridmethod =3
    real,parameter ::  C_beta = 5 ,alpha_in = 1 
 
    real(kind=8):: omegaimax,omegaimin,tao_m  
    !integer :: detectmore
    integer :: change_jump(Nx,Ny)
    real(kind=8):: left_edge_jump, right_edge_jump , bottom_edge_jump, top_edge_jump,omega_i(Nx,Ny),jumpmin,jumpmax,jumpmax0!jump(Nx,Ny),
    real(kind=8) :: pressureR(0:Nx,Ny,NumGLP),pressureL(Nx1,Ny,NumGLP),pressureU(Nx,0:Ny,NumGLP),pressureD(Nx,Ny1,NumGLP)
    real(kind=8) :: entropyR(0:Nx,Ny,NumGLP),entropyL(Nx1,Ny,NumGLP),entropyU(Nx,0:Ny,NumGLP),entropyD(Nx,Ny1,NumGLP)
    real(kind=8) :: enthayR(0:Nx,Ny,NumGLP),enthayL(Nx1,Ny,NumGLP),enthayU(Nx,0:Ny,NumGLP),enthayD(Nx,Ny1,NumGLP)
    real(kind=8) :: jumpa
    call set_bc

    uhmod = uh


    !trouble_num = 0
    !trouble_numall = 0 
    jumpmin = 10000000d0
    jumpmax = 0d0
    !change_all = 0 
    omegaimax = -10
    omegaimin = 10
    !detectmore = 0 
    change_all = 0 
    
    change_jump = 1 

    !!!!!! ??????????????????????????????????????jump?????????????????§µ????§¹????????(????????????§»?????),?????????,?????Gauss??------??????????????????????????????????????
    pressureR = gamma1*(UR(:,:,0,:,4)-(UR(:,:,0,:,2)**2+UR(:,:,0,:,3)**2)/(2*UR(:,:,0,:,1))) 
    entropyR = pressureR / UR(:,:,0,:,1)**gamma
    enthayR = (UR(:,:,0,:,4) + pressureR)/UR(:,:,0,:,1)

    pressureL = gamma1*(UL(:,:,0,:,4)-(UL(:,:,0,:,2)**2+UL(:,:,0,:,3)**2)/(2*UL(:,:,0,:,1))) 
    entropyL = pressureL / UL(:,:,0,:,1)**gamma
    enthayL = (UL(:,:,0,:,4) + pressureL)/UL(:,:,0,:,1)

    pressureU = gamma1*(UU(:,:,0,:,4)-(UU(:,:,0,:,2)**2+UU(:,:,0,:,3)**2)/(2*UU(:,:,0,:,1))) 
    entropyU = pressureU / UU(:,:,0,:,1)**gamma
    enthayU = (UU(:,:,0,:,4) + pressureU)/UU(:,:,0,:,1)

    pressureD = gamma1*(UD(:,:,0,:,4)-(UD(:,:,0,:,2)**2+UD(:,:,0,:,3)**2)/(2*UD(:,:,0,:,1))) 
    entropyD = pressureD / UD(:,:,0,:,1)**gamma
    enthayD = (UD(:,:,0,:,4) + pressureD)/UD(:,:,0,:,1)
  
    ! UL UR UU UB ?????????????????
    ! The x-direction
    UR_ver = 0
    UL_ver = 0
    UR_ver_derx = 0
    UL_ver_derx = 0
    UR_ver_dery = 0
    UL_ver_dery = 0
    do n = 1,NumEq
        do d = 1,dimPk
            do k = 0,Nphi
                do j = 1,Ny
                    do i = 0,Nx
                        UR_ver(i,j,k,:,n) = UR_ver(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver(:,d)
                        UL_ver(i + 1,j,k,:,n) = UL_ver(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver(:,d)
                        !phiGR_ver_derx(2,3),phiGL_ver_derx(2,3), phiGU_ver_derx(2,3), phiGD_ver_derx(2,3)
                        UR_ver_derx(i,j,k,:,n) = UR_ver_derx(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_derx(:,d)
                        UL_ver_derx(i + 1,j,k,:,n) = UL_ver_derx(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_derx(:,d)

                        UR_ver_dery(i,j,k,:,n) = UR_ver_dery(i,j,k,:,n) + uh(i,j,k,d,n)*phiGR_ver_dery(:,d)
                        UL_ver_dery(i + 1,j,k,:,n) = UL_ver_dery(i + 1,j,k,:,n) + uh(i + 1,j,k,d,n)*phiGL_ver_dery(:,d)

                    end do
                end do
            end do
        end do
    end do
    UR_ver_derx =  UR_ver_derx*2/hx
    UL_ver_derx =  UL_ver_derx*2/hx
    UR_ver_dery =  UR_ver_dery*2/hy
    UL_ver_dery =  UL_ver_dery*2/hy
    ! The y-direction
    UU_ver = 0
    UD_ver = 0
    UU_ver_derx = 0
    UD_ver_derx = 0
    UU_ver_dery = 0
    UD_ver_dery = 0
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
                    end do
                end do
            end do
        end do
    end do
    UU_ver_derx =  UU_ver_derx*2/hx
    UD_ver_derx =  UD_ver_derx*2/hx
    UU_ver_dery =  UU_ver_dery*2/hy
    UD_ver_dery =  UD_ver_dery*2/hy

    
    ! ???????????????????????????
    ploydeg = 1 
    
    ! Step1 calculate jump
    do i = 1,Nx
        do j = 1,Ny
            do k = 0,Nphi
                
            !scal = 1/(2*ploydeg-1)
            scal = 0.1d0
            change = 0
            rhoij = uh(i,j,k,1,1)
            u1ij = uh(i,j,k,1,2)/rhoij
            u2ij = uh(i,j,k,1,3)/rhoij
            Eij = uh(i,j,k,1,4)
            pressureij = gamma1*(Eij - 0.5d0*rhoij*(u1ij**2+u2ij**2))
            entropyij =  pressureij/rhoij**gamma
            enthayij  =   ( Eij + pressureij ) /rhoij
            cij = sqrt(abs(gamma*pressureij/rhoij))
            betai = abs(u1ij) + cij 
            betaj = abs(u2ij) + cij

            !Machij  = sqrt(u1ij**2 + u2ij**2)/cij
            !if (Machij<1) then
            !    scal = Machij
            !else
            !    scal = 1/Machij
            !end if
            scal = scal/enthayij
            if (indicator == 1)  then  
                ! call calculate_jump_indicator
                ! ?????????????????????paper1?????????????????????
                left_edge_jump=0
                right_edge_jump =0
                bottom_edge_jump=0
                top_edge_jump=0

                !  enthy = (Ener + pressure)/rho ??????

                do ss = 1,NumGLP
                    ! ???
                     left_edge_jump = left_edge_jump  + abs(UL(i,j,k,ss,1)-UR(i-1,j,k,ss,1)) 
                     right_edge_jump = right_edge_jump  + abs(UL(i+1,j,k,ss,1)-UR(i,j,k,ss,1)) 
                     bottom_edge_jump = bottom_edge_jump + abs(UU(i,j-1,k,ss,1) - UD(i,j,k,ss,1))
                     top_edge_jump = top_edge_jump + abs(UU(i,j,k,ss,1) - UD(i,j+1,k,ss,1))
                    ! ??
                     !left_edge_jump = left_edge_jump   + abs(entropyL(i,j,ss) -  entropyR(i-1,j,ss)) 
                     !right_edge_jump = right_edge_jump  + abs(entropyL(i+1,j,ss) - entropyR(i,j,ss))
                     !bottom_edge_jump = bottom_edge_jump + abs(entropyU(i,j-1,ss) - entropyD(i,j,ss))
                     !top_edge_jump = top_edge_jump + abs(entropyU(i,j,ss) - entropyD(i,j+1,ss))
                    ! ??
                    ! left_edge_jump = left_edge_jump   + abs(enthayL(i,j,ss) -  enthayR(i-1,j,ss)) 
                    ! right_edge_jump = right_edge_jump  + abs(enthayL(i+1,j,ss) - enthayR(i,j,ss))
                    ! bottom_edge_jump = bottom_edge_jump + abs(enthayU(i,j-1,ss) - enthayD(i,j,ss))
                    ! top_edge_jump = top_edge_jump + abs(enthayU(i,j,ss) - enthayD(i,j+1,ss))
                end do
               
                jump(i,j) = (left_edge_jump +  right_edge_jump + bottom_edge_jump + top_edge_jump)/(4*NumGLP*hx**alpha_in)/C_beta 
            
                jumpmin = min(jump(i,j),jumpmin)

                jumpmax = max(jump(i,j),jumpmax)
! print *, "jumpmax = ",  left_edge_jump +  right_edge_jump + bottom_edge_jump + top_edge_jump


                if (jump(i,j) <= 1) then
                    change_jump(i,j) = 0

                else
                    change_jump(i,j) = 1

                    !trouble_num = trouble_num + 1 
                end if 
            end if ! ???? if (indicator == 1)  then

            ! ???if?§Ø? ???????filter
            if (change_jump(i,j) == 0) then ! ???????¦Ê????????????? 

            else! ????filter????

            ! ?????????
            change_all(i,j) =  2 !jump indicator ??????????
            ! calculate_jump_of zero moment density
            left_edge_jump_density=0
            right_edge_jump_density =0
            bottom_edge_jump_density=0
            top_edge_jump_density=0
            do ss = 1,2
                ! ???density
                left_edge_jump_density = left_edge_jump_density  + abs(UL_ver(i,j,0,ss,1)-UR_ver(i-1,j,0,ss,1)) 
                right_edge_jump_density = right_edge_jump_density  + abs(UL_ver(i+1,j,0,ss,1)-UR_ver(i,j,0,ss,1)) 
                bottom_edge_jump_density = bottom_edge_jump_density + abs(UU_ver(i,j-1,0,ss,1) - UD_ver(i,j,0,ss,1))
                top_edge_jump_density = top_edge_jump_density + abs(UU_ver(i,j,0,ss,1) - UD_ver(i,j+1,0,ss,1))

            end do

            deltadensity0 =  betai*(left_edge_jump_density+right_edge_jump_density) + &
            betaj*(bottom_edge_jump_density+top_edge_jump_density)  


            !  calculate_jump_of zero moment momentumn1
            left_edge_jump_rhou=0
            right_edge_jump_rhou =0
            bottom_edge_jump_rhou=0
            top_edge_jump_rhou=0
            do ss = 1,2
                ! momentumn 1
                left_edge_jump_rhou = left_edge_jump_rhou  + abs(UL_ver(i,j,0,ss,2)-UR_ver(i-1,j,0,ss,2)) 
                right_edge_jump_rhou = right_edge_jump_rhou  + abs(UL_ver(i+1,j,0,ss,2)-UR_ver(i,j,0,ss,2)) 
                bottom_edge_jump_rhou = bottom_edge_jump_rhou + abs(UU_ver(i,j-1,0,ss,2) - UD_ver(i,j,0,ss,2))
                top_edge_jump_rhou = top_edge_jump_rhou+ abs(UU_ver(i,j,0,ss,2) - UD_ver(i,j+1,0,ss,2))

            end do

            deltarhou0 =  betai*(left_edge_jump_rhou+right_edge_jump_rhou)  + &
            betaj*(bottom_edge_jump_rhou+top_edge_jump_rhou)  

            ! calculate_jump_of zero moment momentumn2
            left_edge_jump_rhov=0
            right_edge_jump_rhov =0
            bottom_edge_jump_rhov=0
            top_edge_jump_rhov=0
            do ss = 1,2
                ! momentumn 2
                left_edge_jump_rhov = left_edge_jump_rhov  + abs(UL_ver(i,j,0,ss,3)- UR_ver(i-1,j,0,ss,3)) 
                right_edge_jump_rhov = right_edge_jump_rhov  + abs(UL_ver(i+1,j,0,ss,3)-UR_ver(i,j,0,ss,3)) 
                bottom_edge_jump_rhov = bottom_edge_jump_rhov + abs(UU_ver(i,j-1,0,ss,3) - UD_ver(i,j,0,ss,3))
                top_edge_jump_rhov = top_edge_jump_rhov + abs(UU_ver(i,j,0,ss,3) - UD_ver(i,j+1,0,ss,3))

            end do

            deltarhov0 =  betai*(left_edge_jump_rhov+right_edge_jump_rhov)  + &
            betaj*(bottom_edge_jump_rhov+top_edge_jump_rhov) 

            ! call calculate_jump_of zero moment Ener
            left_edge_jump_Ener=0
            right_edge_jump_Ener =0
            bottom_edge_jump_Ener=0
            top_edge_jump_Ener=0
            do ss = 1,2
                ! Ener
                left_edge_jump_Ener= left_edge_jump_Ener  + abs(UL_ver(i,j,0,ss,4)-UR_ver(i-1,j,0,ss,4)) 
                right_edge_jump_Ener= right_edge_jump_Ener  + abs(UL_ver(i+1,j,0,ss,4)-UR_ver(i,j,0,ss,4)) 
                bottom_edge_jump_Ener= bottom_edge_jump_Ener + abs(UU_ver(i,j-1,0,ss,4) - UD_ver(i,j,0,ss,4))
                top_edge_jump_Ener = top_edge_jump_Ener+ abs(UU_ver(i,j,0,ss,4) - UD_ver(i,j+1,0,ss,4))

            end do

            deltaEner0 =  betai*(left_edge_jump_Ener+right_edge_jump_Ener) + &
            betaj*(bottom_edge_jump_Ener+top_edge_jump_Ener) 

            !! 0??moment????§³ 
            delta0max = max( deltadensity0 ,  deltarhou0 ,  deltarhov0 ,  deltaEner0)

            !!!!!!!!!! ???moment 
            ! calculate_jump_of one moment density
            left_edge_jumpderx_density=0
            right_edge_jumpderx_density =0
            bottom_edge_jumpderx_density=0
            top_edge_jumpderx_density=0

            left_edge_jumpdery_density=0
            right_edge_jumpdery_density =0
            bottom_edge_jumpdery_density=0
            top_edge_jumpdery_density=0
            do ss = 1,2
                ! ???density   UR_ver_derx
                left_edge_jumpderx_density = left_edge_jumpderx_density  + abs(UL_ver_derx(i,j,0,ss,1)-UR_ver_derx(i-1,j,0,ss,1)) 
                right_edge_jumpderx_density = right_edge_jumpderx_density  + abs(UL_ver_derx(i+1,j,0,ss,1)-UR_ver_derx(i,j,0,ss,1)) 
                bottom_edge_jumpderx_density = bottom_edge_jumpderx_density + abs(UU_ver_derx(i,j-1,0,ss,1) - UD_ver_derx(i,j,0,ss,1))
                top_edge_jumpderx_density = top_edge_jumpderx_density + abs(UU_ver_derx(i,j,0,ss,1) - UD_ver_derx(i,j+1,0,ss,1))
                ! ???density   UR_ver_dery
                left_edge_jumpdery_density = left_edge_jumpdery_density  + abs(UL_ver_dery(i,j,0,ss,1)-UR_ver_dery(i-1,j,0,ss,1)) 
                right_edge_jumpdery_density = right_edge_jumpdery_density  + abs(UL_ver_dery(i+1,j,0,ss,1)-UR_ver_dery(i,j,0,ss,1)) 
                bottom_edge_jumpdery_density = bottom_edge_jumpdery_density + abs(UU_ver_dery(i,j-1,0,ss,1) - UD_ver_dery(i,j,0,ss,1))
                top_edge_jumpdery_density = top_edge_jumpdery_density + abs(UU_ver_dery(i,j,0,ss,1) - UD_ver_dery(i,j+1,0,ss,1))

            end do
deltadensity1 =  betai*(left_edge_jumpderx_density+right_edge_jumpderx_density + left_edge_jumpdery_density+right_edge_jumpdery_density)*hx*2 + &
            betaj*(bottom_edge_jumpderx_density + top_edge_jumpdery_density + bottom_edge_jumpderx_density + top_edge_jumpdery_density )*hy*2  

            ! calculate_jump_of one moment rhou
            left_edge_jumpderx_rhou=0
            right_edge_jumpderx_rhou =0
            bottom_edge_jumpderx_rhou=0
            top_edge_jumpderx_rhou=0

            left_edge_jumpdery_rhou=0
            right_edge_jumpdery_rhou =0
            bottom_edge_jumpdery_rhou=0
            top_edge_jumpdery_rhou=0
            do ss = 1,2
                ! ???? rhou   derx
                left_edge_jumpderx_rhou = left_edge_jumpderx_rhou  + abs(UL_ver_derx(i,j,0,ss,2)-UR_ver_derx(i-1,j,0,ss,2)) 
                right_edge_jumpderx_rhou = right_edge_jumpderx_rhou  + abs(UL_ver_derx(i+1,j,0,ss,2)-UR_ver_derx(i,j,0,ss,2)) 
                bottom_edge_jumpderx_rhou = bottom_edge_jumpderx_rhou + abs(UU_ver_derx(i,j-1,0,ss,2) - UD_ver_derx(i,j,0,ss,2))
                top_edge_jumpderx_rhou = top_edge_jumpderx_rhou + abs(UU_ver_derx(i,j,0,ss,2) - UD_ver_derx(i,j+1,0,ss,2))
                ! ???? rhou dery
                left_edge_jumpdery_rhou = left_edge_jumpdery_rhou  + abs(UL_ver_dery(i,j,0,ss,2)-UR_ver_dery(i-1,j,0,ss,2)) 
                right_edge_jumpdery_rhou = right_edge_jumpdery_rhou  + abs(UL_ver_dery(i+1,j,0,ss,2)-UR_ver_dery(i,j,0,ss,2)) 
                bottom_edge_jumpdery_rhou = bottom_edge_jumpdery_rhou + abs(UU_ver_dery(i,j-1,0,ss,2) - UD_ver_dery(i,j,0,ss,2))
                top_edge_jumpdery_rhou = top_edge_jumpdery_rhou + abs(UU_ver_dery(i,j,0,ss,2) - UD_ver_dery(i,j+1,0,ss,2))

            end do
deltarhou1 =  betai*(left_edge_jumpderx_rhou +right_edge_jumpderx_rhou + left_edge_jumpdery_rhou+right_edge_jumpdery_rhou)*hx*2 + &
            betaj*(bottom_edge_jumpderx_rhou + top_edge_jumpderx_rhou + bottom_edge_jumpdery_rhou + top_edge_jumpdery_rhou )*hy*2

            ! calculate_jump_of one moment rhov
            left_edge_jumpderx_rhov=0
            right_edge_jumpderx_rhov =0
            bottom_edge_jumpderx_rhov=0
            top_edge_jumpderx_rhov=0

            left_edge_jumpdery_rhov=0
            right_edge_jumpdery_rhov =0
            bottom_edge_jumpdery_rhov=0
            top_edge_jumpdery_rhov=0
            do ss = 1,2
                ! ???? rhou   
                left_edge_jumpderx_rhov = left_edge_jumpderx_rhov  + abs(UL_ver_derx(i,j,0,ss,3)-UR_ver_derx(i-1,j,0,ss,3)) 
                right_edge_jumpderx_rhov = right_edge_jumpderx_rhov  + abs(UL_ver_derx(i+1,j,0,ss,3)-UR_ver_derx(i,j,0,ss,3)) 
                bottom_edge_jumpderx_rhov = bottom_edge_jumpderx_rhov + abs(UU_ver_derx(i,j-1,0,ss,3) - UD_ver_derx(i,j,0,ss,3))
                top_edge_jumpderx_rhov = top_edge_jumpderx_rhov + abs(UU_ver_derx(i,j,0,ss,3) - UD_ver_derx(i,j+1,0,ss,3))
                ! ???? rhou
                left_edge_jumpdery_rhov = left_edge_jumpdery_rhov+ abs(UL_ver_dery(i,j,0,ss,3)-UR_ver_dery(i-1,j,0,ss,3)) 
                right_edge_jumpdery_rhov = right_edge_jumpdery_rhov + abs(UL_ver_dery(i+1,j,0,ss,3)-UR_ver_dery(i,j,0,ss,3)) 
                bottom_edge_jumpdery_rhov= bottom_edge_jumpdery_rhov + abs(UU_ver_dery(i,j-1,0,ss,3) - UD_ver_dery(i,j,0,ss,3))
                top_edge_jumpdery_rhov = top_edge_jumpdery_rhov + abs(UU_ver_dery(i,j,0,ss,3) - UD_ver_dery(i,j+1,0,ss,3))

            end do

            deltarhov1 =  betai*(left_edge_jumpderx_rhov +right_edge_jumpderx_rhov + left_edge_jumpdery_rhov+right_edge_jumpdery_rhov)*hx*2 + &
            betaj*(bottom_edge_jumpderx_rhov + top_edge_jumpderx_rhov + bottom_edge_jumpdery_rhov + top_edge_jumpdery_rhov )*hy*2 


            ! calculate_jump_of one moment Ener
            left_edge_jumpderx_Ener=0
            right_edge_jumpderx_Ener =0
            bottom_edge_jumpderx_Ener =0
            top_edge_jumpderx_Ener =0

            left_edge_jumpdery_Ener =0
            right_edge_jumpdery_Ener =0
            bottom_edge_jumpdery_Ener=0
            top_edge_jumpdery_Ener=0
            do ss = 1,2
                ! ???? Ener   
                left_edge_jumpderx_Ener = left_edge_jumpderx_Ener  + abs(UL_ver_derx(i,j,0,ss,4)-UR_ver_derx(i-1,j,0,ss,4)) 
                right_edge_jumpderx_Ener = right_edge_jumpderx_Ener  + abs(UL_ver_derx(i+1,j,0,ss,4)-UR_ver_derx(i,j,0,ss,4)) 
                bottom_edge_jumpderx_Ener = bottom_edge_jumpderx_Ener + abs(UU_ver_derx(i,j-1,0,ss,4) - UD_ver_derx(i,j,0,ss,4))
                top_edge_jumpderx_Ener = top_edge_jumpderx_Ener + abs(UU_ver_derx(i,j,0,ss,4) - UD_ver_derx(i,j+1,0,ss,4))
                ! ???? Ener  
                left_edge_jumpdery_Ener = left_edge_jumpdery_Ener + abs(UL_ver_dery(i,j,0,ss,4)-UR_ver_dery(i-1,j,0,ss,4)) 
                right_edge_jumpdery_Ener = right_edge_jumpdery_Ener + abs(UL_ver_dery(i+1,j,0,ss,4)-UR_ver_dery(i,j,0,ss,4)) 
                bottom_edge_jumpdery_Ener = bottom_edge_jumpdery_Ener + abs(UU_ver_dery(i,j-1,0,ss,4) - UD_ver_dery(i,j,0,ss,4))
                top_edge_jumpdery_Ener = top_edge_jumpdery_Ener + abs(UU_ver_dery(i,j,0,ss,4) - UD_ver_dery(i,j+1,0,ss,4))

            end do

            deltaEner1 =  betai*(left_edge_jumpderx_Ener +right_edge_jumpderx_Ener + left_edge_jumpdery_Ener+right_edge_jumpdery_Ener)*hx*2 + &
            betaj*(bottom_edge_jumpderx_Ener + top_edge_jumpderx_Ener + bottom_edge_jumpdery_Ener + top_edge_jumpdery_Ener)*hy*2 


            !! 1??moment????§³ 
            delta1max = max( deltadensity1 ,  deltarhou1 ,  deltarhov1 ,  deltaEner1)

            damping =   delta0max + delta1max 

            !damping = 0.1d0*hx*hy*damping*(8d0)/hx 

            damping =   scal*hx*damping*(8d0)/hx

            uhmod(i,j,k,2:3,1:4) = exp(-dt*damping)*uh(i,j,k,2:3,1:4) 

            end if 

            end do
        end do
    end do

    
   
   
    
 ! ??????
    if  (ishybrid == 1) then
        ! find the biggest jump
        do the_id = 2,N_process

        if (myid1 == the_id) then
            call MPI_SEND(jumpmax,1,MPI_REAL8,0,1,MPI_COMM_WORLD,ierr)
        end if

        if (myid1 == 1) then
            call MPI_RECV(jumpmax0,1,MPI_REAL8,the_id - 1,1,MPI_COMM_WORLD,status,ierr)
            if (jumpmax0 > jumpmax) then
                jumpmax= jumpmax0
            end if
        end if

        end do

        do the_id = 2,N_process

        if (myid1 == 1) then
            call MPI_SEND(jumpmax,1,MPI_REAL8,the_id - 1,2,MPI_COMM_WORLD,ierr)
        end if

        if (myid1 == the_id) then
            call MPI_RECV(jumpmax,1,MPI_REAL8,0,2,MPI_COMM_WORLD,status,ierr)
        end if

        end do


        if (hybridmethod  == 1) then
            do i = 1,Nx
                do j = 1,Ny
                     do k = 0,Nphi
                    if (change_all(i,j) == 2 ) then
                        omega_i(i,j) =  (jumpmax - jump(i,j))/(jumpmax - 1)
                        do d = 1,4
                            uhmod(i , j ,k, : ,d) = omega_i(i,j)*  uh(i , j ,k, : ,d) + (1 - omega_i(i,j))*  uhmod(i , j ,k, : ,d) 
                        end do
                        ! if   (jump(i,j)>=1d0/8d0* jumpmax) then
                        !        omega_i(i,j) = 0
                        !else  
                        !      omega_i(i,j) =  (1d0/8d0*jumpmax - jump(i,j))/(1d0/8d0*jumpmax - 1)
                        ! end if 
                        ! do d = 1,8
                        !     uhmod(i , j ,0, : ,d) = omega_i(i,j)*  uh(i , j ,0, : ,d) + (1 - omega_i(i,j))*  uhmod(i , j ,0, : ,d) 
                        ! end do
                        !omegaimax = max(omega_i(i,j),omegaimax)
                        !omegaimin = min(omega_i(i,j),omegaimin)
                        ! if (omega_i(i,j)>0.9) then
                        !    detectmore = detectmore + 1    
                        !end if 
                    end if 
                    end do
                end do
            end do 

        end if

        if (hybridmethod  == 2) then
            do i = 1,Nx
                do j = 1,Ny
                     do k = 0,Nphi
                    if (change_all(i,j) == 2 ) then
                        omega_i(i,j) = exp(-alpha1*( jump(i,j)-1)**order1/( (jumpmax-1)**order1- ( jump(i,j)-1)**order1 ))       
                        do d = 1,4
                            uhmod(i , j ,k, : ,d) = omega_i(i,j)*  uh(i , j ,k, : ,d) + (1 - omega_i(i,j))*  uhmod(i , j ,k, : ,d) 
                        end do
                    end if 
                    end do 
                end do
            end do 

        end if

        if (hybridmethod  == 3) then
            do i = 1,Nx
                do j = 1,Ny
                     do k = 0,Nphi
                    if (change_all(i,j) == 2 ) then
                        omega_i(i,j) = 1 - exp(-alpha2/( jump(i,j)-1 )**order2  )       
                        do d = 1,4
                            uhmod(i , j ,k, : ,d) = omega_i(i,j)*  uh(i , j ,k, : ,d) + (1 - omega_i(i,j))*  uhmod(i , j ,k, : ,d) 
                        end do
                    end if 
                    end do
                end do
            end do 

        end if

    end if !if  (ishybrid == 1)????

!print *, "jumpmax = ", jumpmax
     uh = uhmod


    end  subroutine hybridjumpfilter 
