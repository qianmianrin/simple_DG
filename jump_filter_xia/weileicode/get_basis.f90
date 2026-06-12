
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