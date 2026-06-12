

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
