
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