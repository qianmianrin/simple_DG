
     
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