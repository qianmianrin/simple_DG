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