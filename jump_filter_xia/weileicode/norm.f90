    subroutine norm(x,d)

    real x(4),d

    d = 0

    do i = 1,4
        d = d + x(i)**2
    end do

    d = d**0.5

    end subroutine norm
