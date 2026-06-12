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
