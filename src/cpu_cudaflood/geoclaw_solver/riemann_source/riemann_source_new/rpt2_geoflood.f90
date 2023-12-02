! @description: Riemann solver in the transverse direction using an einfeldt Jacobian

subroutine fc2d_geoclaw_rpt2(ixy,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

    use geoclaw_module, only: g => grav, tol => dry_tolerance
    use geoclaw_module, only: coordinate_system,earth_radius,deg2rad

    implicit none

    ! Input variables
    integer, intent(in) :: ixy,imp,maxm,meqn,mwaves,maux,mbc,mx

    ! Input/Output variables
    double precision, intent(inout) :: ql(meqn,1-mbc:maxm+mbc)
    double precision, intent(inout) :: qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(inout) :: aux1(maux,1-mbc:maxm+mbc)
    double precision, intent(inout) :: aux2(maux,1-mbc:maxm+mbc) 
    double precision, intent(inout) :: aux3(maux,1-mbc:maxm+mbc)
    double precision, intent(inout) :: asdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(inout) :: bmasdq(meqn,1-mbc:maxm+mbc)
    double precision, intent(inout) :: bpasdq(meqn,1-mbc:maxm+mbc)

    ! Local variables
    integer :: i,m,mw,mu,mv
    double precision :: s(3)
    double precision :: r(3,3)
    double precision :: beta(3)
    double precision :: abs_tol
    double precision :: hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
    double precision :: uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
    double precision :: delf1,delf2,delf3,dxdcd,dxdcu
    double precision :: dxdcm,dxdcp,topo1,topo3,eta

    abs_tol = tol

    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    end if

    do i = 2-mbc,mx+mbc

        hl = qr(1,i-1)
        hr = ql(1,i)
        hul = qr(mu,i-1)
        hur = ql(mu,i)
        hvl = qr(mv,i-1)
        hvr = ql(mv,i)

        ! === determine velocity from momentum ===
        ! dry inital left state
        if (hl < abs_tol) then
            hl = 0.0d0
            ul = 0.0d0
            vl = 0.0d0
        else
            ul = hul/hl
            vl = hvl/hl
        endif

        ! dry inital right state
        if (hr < abs_tol) then
            hr = 0.0d0
            ur = 0.0d0
            vr = 0.0d0
        else
            ur = hur/hr
            vr = hvr/hr
        endif

        do mw = 1,mwaves
            s(mw) = 0.0d0
            beta(mw) = 0.0d0
            do m = 1,meqn
                r(m,mw) = 0.0d0
            end do
        end do

        dxdcp = 1.0d0
        dxdcm = 1.0d0

        if (hl <= abs_tol .and. hr <= abs_tol) go to 90 ! if entirely dry initially

        ! check and see if call that transverse waves are going in is high and dry
        if (imp == 1) then
            eta = qr(1,i-1) + aux2(1,i-1)
            topo1 = aux1(1,i-1)
            topo3 = aux3(1,i-1)
        else
            eta = ql(1,i) + aux2(1,i)
            topo1 = aux1(1,i)
            topo3 = aux3(1,i)
        endif

        if (eta < max(topo1,topo3)) go to 90

        if (coordinate_system == 2) then
            if (ixy == 2) then
                dxdcp = earth_radius*deg2rad
                dxdcm = dxdcp
            else
                if (ixy == 1) then
                    dxdcp = earth_radius*cos(aux3(3,i-1))*deg2rad
                    dxdcm = earth_radius*cos(aux1(3,i-1))*deg2rad
                else
                    dxdcp = earth_radius*cos(aux3(3,i))*deg2rad
                    dxdcm = earth_radius*cos(aux1(3,i))*deg2rad
                endif
            endif
        endif
        
        ! === Determine some speeds necessary for the Jacobian ===
        vhat = (vr*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + (vl*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
        uhat = (ur*dsqrt(hr))/(dsqrt(hr)+dsqrt(hl)) + (ul*dsqrt(hl))/(dsqrt(hr)+dsqrt(hl))
        hhat = (hr+hl)/2.0d0

        roe1 = vhat - dsqrt(g*hhat)
        roe3 = vhat + dsqrt(g*hhat)

        s1l = vl - dsqrt(g*hl)
        s3r = vr + dsqrt(g*hr)

        s1 = dmin1(s1l,roe1)
        s3 = dmax1(s3r,roe3)

        s2 = 0.5d0*(s1+s3)

        s(1) = s1
        s(2) = s2
        s(3) = s3

        ! === Determine asdq decomposition (beta) ===
        delf1 = asdq(1,i)
        delf2 = asdq(mu,i)
        delf3 = asdq(mv,i)

        beta(1) = (s3*delf1/(s3-s1)) - (delf3/(s3-s1))
        beta(2) = -s2*delf1 + delf2
        beta(3) = (delf3/(s3-s1)) - (s1*delf1/(s3-s1))

        ! === set-up eigenvectors ===
        r(1,1) = 1.0d0
        r(2,1) = s2
        r(3,1) = s1

        r(1,2) = 0.0d0
        r(2,2) = 1.0d0
        r(3,2) = 0.0d0

        r(1,3) = 1.0d0
        r(2,3) = s2
        r(3,3) = s3

90  continue

        ! === compute fluctuations ===
        bmasdq(1,i) = 0.0d0
        bpasdq(1,i) = 0.0d0
        bmasdq(2,i) = 0.0d0
        bpasdq(2,i) = 0.0d0
        bmasdq(3,i) = 0.0d0
        bpasdq(3,i) = 0.0d0
        do mw = 1,3
            if (s(mw) < 0.d0) then
                bmasdq(1,i) = bmasdq(1,i) + dxdcm*s(mw)*beta(mw)*r(1,mw)
                bmasdq(mu,i) = bmasdq(mu,i) + dxdcm*s(mw)*beta(mw)*r(2,mw)
                bmasdq(mv,i) = bmasdq(mv,i) + dxdcm*s(mw)*beta(mw)*r(3,mw)
            elseif (s(mw) > 0.d0) then
                bpasdq(1,i) = bpasdq(1,i) + dxdcp*s(mw)*beta(mw)*r(1,mw)
                bpasdq(mu,i) = bpasdq(mu,i) + dxdcp*s(mw)*beta(mw)*r(2,mw)
                bpasdq(mv,i) = bpasdq(mv,i) + dxdcp*s(mw)*beta(mw)*r(3,mw)
            endif
        end do

    end do

    return

end subroutine fc2d_geoclaw_rpt2