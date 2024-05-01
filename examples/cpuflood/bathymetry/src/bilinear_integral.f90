!@des: This function bilinearintegral integrates the bilinear with values z11,z12,z21,z22 over the area defined by rectangle: xim <= x <= xip, yjm <= y <= yjp

function bilinearintegral(xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy, &
                          z11,z12,z21,z22)

    implicit none

    ! declare variables
    double precision, intent(in) :: xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy
    double precision, intent(in) :: z11,z12,z21,z22
    double precision :: bilinearintegral

    ! declare loacl variables
    double precision :: xlow,ylow,xhi,yhi,area,sumxi,sumeta,a,b,c,d

    ! find limits of integral
    xlow = max(xim,x1)
    ylow = max(yjm,y1)
    xhi = min(xip,x2)
    yhi = min(yjp,y2)

    ! find the area of integration
    area = (xhi-xlow)*(yhi-ylow)
    sumxi = (xhi + xlow - 2.d0*x1)/dxx
    sumeta = (yhi + ylow - 2.d0*y1)/dyy

    ! find the coefficients of the bilinear function a*xi + b*eta + c*xi*eta + d
    a = (z21 - z11)
    b = (z12 - z11)
    c = (z22 - z21 - z12 + z11)
    d = z11    

    ! find the integral
    bilinearintegral = (0.5d0*(a*sumxi + b*sumeta) + 0.25d0*c*sumxi*sumeta + d)*area

    return

end function bilinearintegral


!@des: This function bilinearintegral_s integrates the bilinear with values z11,z12,z21,z22 
!      over the area defined by rectangle: xim <= x <= xip, yjm <= y <= yjp. 
!      The integration is actually done on the surface of a sphere

function bilinearintegral_s(xim, xip, yjm, yjp, x1, x2, y1, y2, &
                            dxx, dyy, z11, z12, z21, z22)

    use geoclaw_module, only: rad2deg, deg2rad, earth_radius

    implicit none

    ! declare variables
    double precision, intent(in) ::  xim,xip,yjm,yjp,x1,x2,y1,y2,dxx,dyy
    double precision, intent(in) :: z11,z12,z21,z22
    double precision :: bilinearintegral_s

    ! declare local variables
    double precision :: xlow,ylow,xhi,yhi,delx,dely,a,b,c,d
    double precision :: xdiffhi,xdifflo,ydiffhi,ydifflo,xdiff2
    double precision :: adsinint,cbsinint

    ! find limits of integral
    xlow = max(xim,x1)
    ylow = max(yjm,y1)
    xhi = min(xip,x2)   
    yhi = min(yjp,y2)

    delx = xhi - xlow
    dely = yhi - ylow

    !  find integration terms
    xdiffhi = xhi - x1
    xdifflo = xlow - x1
    ydiffhi = yhi - y1
    ydifflo = ylow - y1
    xdiff2 = 0.5d0*(xdiffhi**2 - xdifflo**2)

    cbsinint = (rad2deg*cos(deg2rad*yhi) + ydiffhi*sin(deg2rad*yhi)) - &
                (rad2deg*cos(deg2rad*ylow) - ydifflo*sin(deg2rad*ylow))

    adsinint = rad2deg*(sin(deg2rad*yhi) - sin(deg2rad*ylow))

    ! find the coefficients of the bilinear function a*xi + b*eta + c*xi*eta + d
    a = (z21 - z11)/dxx
    b = (z12 - z11)/dyy
    c = (z22 - z21 - z12 + z11)/(dxx*dyy)
    d = z11

    ! find the integral
    bilinearintegral_s = ((a*xdiff2 + d*delx)*adsinint + &
                          rad2deg*(c*xdiff2 + b*delx)*cbsinint)*(earth_radius*deg2rad)**2


    return

end function bilinearintegral_s
