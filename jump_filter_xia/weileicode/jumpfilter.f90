
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
    
    
    
