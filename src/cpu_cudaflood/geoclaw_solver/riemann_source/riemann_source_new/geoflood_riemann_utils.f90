! Augumented Riemann solver for the shallow water equations (AUG-JCP)
! @description: This subroutine solves the Riemann problem for the shallow water equations given single left and right states
!               This solver is described in J. Comput. Phys. (6): 3089-3113, March 2008
!               Augmented Riemann Solvers for the Shallow Equations with Steady States and Inundation
! 
!               To use the original solver call with maxiter=1.
!               This solver allows iteration when maxiter > 1. The iteration seems to help with instbilities taht arise (with any solver) 
!               as flow becomes transcritical over variable topo due to loss of hyperbolicity.
subroutine riemann_aug_JCP(maxiter,meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw)

    implicit none

    ! input
    integer, intent(in) :: maxiter,meqn,mwaves

    double precision, intent(inout) :: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,sE1,sE2
    double precision, intent(inout) :: hvL,hvR,vL,vR
    double precision, intent(in) :: drytol,g

    ! output
    double precision, intent(out) :: fw(meqn,mwaves)
    double precision, intent(out) :: sw(mwaves)

    ! local
    integer :: m,mw,k,iter,max_iter
    double precision :: A(3,3),r(3,3)
    double precision :: lambda(3), d(3),beta(3),del(3)

    double precision :: dh,dhu,dphi,db,delnorm
    double precision :: rare1st,rare2st,sdelta,raremin,raremax
    double precision :: criticaltol,convergencetol,raretol
    double precision :: s1s2bar,s1s2tilde,hbar,hLstar,hRstar,hustar
    double precision :: huRstar,huLstar,uRstar,uLstar,hstarHLL,HstarHLL_
    double precision :: ddh,ddphi
    double precision :: s1m,s2m,hm
    double precision :: det1,det2,det3,determinant

    logical :: rare1,rare2,rarecorrector,rarecorrectortest,sonic

    ! determine jump vectors for the left and right states (Augmented system)
    dh   = hR - hL              ! depth jump
    dhu  = huR - huL            ! momentum jump
    dphi = phiR - phiL          ! momentum flux jump -- phi = hu^2 + g h^2/2
    db   = bR - bL              ! bathymetry jump

    delnorm = dh**2 + dphi**2    

    ! call the Riemann type solver to determine the Riemann structure (wave-type in each family)
    max_iter = 1 ! use the original solver
    call riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,max_iter,drytol,g)

    ! For the solver to handle depth negativity, depth dh is included in the decompostion which 
    ! gives as acess to using the depth positive semidefinite solver (HLLE). This makes the system 
    ! to have 3 waves instead of 2. where the 1st and 3rd are the eigenpairs are related to the flux
    ! Jacobian matrix of the original SWE (since s1<s2<s3, and have been modified by Einfeldt to handle
    ! depth non-negativity) and the 2nd is refered to as the the entropy corrector wave since its introduced 
    ! to correct entropy violating solutions with only 2 waves.

    ! the 1st and 3rd speeds are the eigenvalues of the Jacobian matrix of the original SWE modified by Einfeldt's 
    ! for use with the HLLE solver. 
    lambda(1) = min(sE1,s2m) ! modified by Einfeldt speed; sE1 - flux Jacobian eigen value s2m - Roe speed
    lambda(3) = max(sE2,s1m) ! modified by Einfeldt speed; sE2 - flux Jacobian eigen value s2m - Roe speed

    ! Einfeldt speeds
    sE1 = lambda(1)
    sE2 = lambda(3)

    ! the 2nd speed is the entropy corrector wave speed
    lambda(2) = 0.d0  ! no strong or significant rarefaction wave

    ! determine the middle state in the HLLE solver
    HstarHLL_ = (huL - huR + sE2*hR - sE1*hL)/(sE2 - sE1)
    hstarHLL = max(HstarHLL_,0.d0) ! middle state between the two discontinuities (positive semidefinite depth)

    ! === determine the middle entropy corrector wave ===
    ! rarecorrectortest = .true. provides a more accurate Riemann solution but is more expensive. This is because
    ! a nonlinear Riemann solution with  2 nonlinear waves as a linear Riemann solution 3 (or 2 jump discontionuities 
    ! to approximate 1 smooth nonlinear rarefaction if it's large). When rarecorrectortest = .false. the approximate 
    ! solution has only 2 jump discontinuities instead of 3, so its less accurate but faster.
    rarecorrectortest = .false. 
    rarecorrector     = .false.
    if (rarecorrectortest) then
        sdelta = lambda(3) - lambda(1)
        raremin = 0.5d0 ! indicate a large rarefaction wave but not large 
        raremax = 0.9d0 ! indicate a very large rarefaction wave
        ! i.e (the total speed difference between the fastest and slowest wave in the Riemann solution = 0.5)
        
        if (rare1 .and. sE1*s1m < 0.d0) raremin = 0.2d0
        if (rare2 .and. sE2*s2m < 0.d0) raremax = 0.2d0

        if (rare1 .or. rare2) then
            ! check which rarefaction is the strongest
            rare1st = 3.d0*(sqrt(g*hL) - sqrt(g*hm))
            rare2st = 3.d0*(sqrt(g*hR) - sqrt(g*hm))
            if (max(rare1st,rare2st) < raremin*sdelta .and. min(rare1st,rare2st) < raremax*sdelta) then
                rarecorrector = .true.
                if (rare1st > rare2st) then
                    lambda(2) = s1m
                elseif (rare2st > rare1st) then
                    lambda(2) = s2m
                else
                    lambda(2) = 0.5d0*(s1m + s2m)
                end if
            end if
        end if
        if (hstarHLL < min(hL,hR)/5.d0) rarecorrector = .false.
    end if

    ! determining modified eigen vectors for the 
    do mw=1,mwaves
        r(1,mw) = 1.d0
        r(2,mw) = lambda(mw)
        r(3,mw) = lambda(mw)**2
    end do
    
    ! no strong rarefaction wave
    if (.not.rarecorrector) then
        lambda(2) = 0.5d0*(lambda(1) + lambda(3))
        r(1,2) = 0.d0
        r(2,2) = 0.d0
        r(3,2) = 1.d0
    end if

    ! === Determine the steady state wave ===
    criticaltol = 1.d-6 
    ddh = -dh 
    ddphi = -g*0.5d0*(hR+hL)*db ! some approximation of the source term \int_{x_{l}}^{x_{r}} -g h b_x dx

    ! determine a few quantities needed for steady state wave if iterated
    hLstar = hL
    hRstar = hR
    uLstar = uL
    uRstar = uR
    huLstar = uLstar*hLstar
    huRstar = uRstar*hRstar

    ! iterate to find the steady state wave
    convergencetol = 1.d-6
    do iter = 1,maxiter
        ! determine the steady state wave (this will be subtracted from the delta vectors)
        if (min(hLstar,hRstar) < drytol .and. rarecorrector) then
            rarecorrector = .false.
            hLstar = hL
            hRstar = hR
            uLstar = uL
            uRstar = uR
            huLstar = uLstar*hLstar
            huRstar = uRstar*hRstar
            lambda(2) = 0.5d0*(lambda(1) + lambda(3))
            r(1,2) = 0.d0
            r(2,2) = 0.d0
            r(3,2) = 1.d0
        end if

        ! For any two states; Q_i and Q_i-1, eigen values of SWE must satify: lambda(q_i)*lambda(q_i-1) = u^2 -gh, 
        ! writing this conditon as a function of Q_i and Q_i-1, u and h become averages in lambda(q_i)*lambda(q_i-1) = u^2 -gh
        ! and these averages are denoted by bar and tilde.
        hbar = max(0.5d0*(hLstar + hRstar),0.d0)
        s1s2bar = 0.25d0*(uLstar + uRstar)**2 -g*hbar
        s1s2tilde = max(uLstar*uRstar,0.d0) - g*hbar

        ! Based on the above conditon, smooth staedy state over slopping bathymetry cannot have a sonic point. 
        ! Therefore, for regions with monotonically varying bathymetry, steady-state flow is either entirely 
        ! subsonic (-u^2 +gh > 0) or entirely supersonic. 
        sonic = .false.
        if (abs(s1s2bar) <= criticaltol) sonic = .true.
        if (s1s2bar*s1s2tilde <= criticaltol) sonic = .true.
        if (s1s2bar*sE1*sE2 <= criticaltol) sonic = .true.
        if (min(abs(sE1),abs(sE2)) <= criticaltol) sonic = .true.
        if (sE1 < 0.d0 .and. s1m > 0.d0) sonic = .true.
        if (sE2 > 0.d0 .and. s2m < 0.d0) sonic = .true.
        if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)) < 0.d0) sonic = .true.
        if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)) < 0.d0) sonic = .true.

        ! find jump in h, ddh
        if (sonic) then
            ddh = -dh
        else
            ddh =db*g*hbar/s1s2bar
        end if

        ! find bounds in case of critical state resonance, or negative states
        if (sE1 < -criticaltol .and. sE2 > criticaltol) then
            ddh = min(ddh,hstarHLL*(sE2-sE1/sE2))
            ddh = max(ddh,hstarHLL*(sE1-sE2/sE1))
        else if (sE1 >= criticaltol) then
            ddh = min(ddh,hstarHLL*(sE2-sE1/sE1))
            ddh = max(ddh,-hL)
        else if (sE2 <= -criticaltol) then
            ddh = min(ddh,hR)
            ddh = max(ddh,hstarHLL*(sE1-sE2/sE2))
        endif

        ! find jump in phi, ddphi
        if (sonic) then
            ddphi = -g*hbar*db
        else
            ddphi = -db*g*hbar*s1s2tilde/s1s2bar
        end if

        ! find bounds in case of critical state resonance, or negative states
        ddphi = min(ddphi,g*max(-hLstar*db,-hRstar*db))
        ddphi = max(ddphi,g*min(-hLstar*db,-hRstar*db))

        ! Determine coefficients beta(k) using crammer's rule
        ! first determine the determinant of the eigenvector matrix
        det1=r(1,1)*(r(2,2)*r(3,3)-r(2,3)*r(3,2))
        det2=r(1,2)*(r(2,1)*r(3,3)-r(2,3)*r(3,1))
        det3=r(1,3)*(r(2,1)*r(3,2)-r(2,2)*r(3,1))
        determinant=det1-det2+det3

        ! determine the delta vectors
        del(1) = dh - ddh
        del(2) = dhu
        del(3) = dphi - ddphi

        ! solve for beta(k)
        do k=1,3
            do mw = 1,3
                A(1,mw) = r(1,mw)
                A(2,mw) = r(2,mw)
                A(3,mw) = r(3,mw)
            end do
            A(1,k) = del(1)
            A(2,k) = del(2)
            A(3,k) = del(3)
            det1=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
            det2=A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
            det3=A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
            beta(k)=(det1-det2+det3)/determinant
        end do

        !  exist if things aren't changing
        if (abs(del(1)**2 + del(3)**2 - delnorm) < convergencetol) exit
        delnorm = del(1)**2 + del(3)**2

        ! find new states qLstar and qRstar on either side of interface
        hLstar = hL
        hRstar = hR
        uLstar = uL
        uRstar = uR
        huLstar = uLstar*hLstar
        huRstar = uRstar*hRstar

        ! left state depth and momentum updates
        do mw = 1,mwaves
            if (lambda(mw) < 0.d0) then
                hLstar = hLstar + beta(mw)*r(1,mw)
                huLstar = huLstar + beta(mw)*r(2,mw)
            endif
        end do

        ! right state depth and momentum updates
        do mw = mwaves,1,-1
            if (lambda(mw) > 0.d0) then
                hRstar = hRstar + beta(mw)*r(1,mw)
                huRstar = huRstar + beta(mw)*r(2,mw)
            endif
        end do

        ! left state velocity update
        if (hLstar > drytol) then
            uLstar = huLstar/hLstar
        else ! dry bed
            hLstar = max(hLstar,0.d0)
            uLstar = 0.d0
        end if

        ! right state velocity update
        if (hRstar > drytol) then
            uRstar = huRstar/hRstar
        else ! dry bed
            hRstar = max(hRstar,0.d0)
            uRstar = 0.d0
        end if

    end do ! end of iteration on Riemann problem

    ! === determine the fwaves and speeds===
    do mw = 1,mwaves
        sw(mw) = lambda(mw)
        fw(1,mw) = beta(mw)*r(2,mw)
        fw(2,mw) = beta(mw)*r(3,mw)
        fw(3,mw) = beta(mw)*r(2,mw)
    end do

    ! find transverse components (ie huv jumps)
    fw(3,1) = fw(3,1)*vL
    fw(3,3) = fw(3,3)*vR
    fw(3,2) = hR*uR*vR - hL*uL*vL - fw(3,1) - fw(3,3)

    return
end subroutine riemann_aug_JCP


! fwave-approach
! @description: solve shallow water equations given single left and right states  solutin has two waves
!               flux - source is decomposed
subroutine riemann_fwave(meqn,mwaves,hL,hR,huL,huR,hvL,hvR,bL,bR,uL,uR,vL,vR,phiL,phiR,s1,s2,drytol,g,sw,fw)

    implicit none

    ! input
    integer, intent(in) :: meqn,mwaves

    double precision, intent(in) :: hL,hR,huL,huR,bL,bR,uL,uR,phiL,phiR,s1,s2
    double precision, intent(in) :: hvL,hvR,vL,vR
    double precision, intent(in) :: drytol,g

    ! output
    double precision, intent(out) :: fw(meqn,mwaves)
    double precision, intent(out) :: sw(mwaves)

    ! local
    double precision :: dh,dhu,dphi,db,dh_decomp,dphi_decomp
    double precision :: h_ave,ddphi
    double precision :: beta1,beta2

    ! determine jump vectors for the left and right states (Augmented system)
    dh = hR - hL              ! depth jump
    dhu = huR - huL           ! momentum jump
    dphi = phiR - phiL        ! momentum flux jump -- phi = hu^2 + g h^2/2
    db = bR - bL              ! bathymetry jump

    ! At steady state, u+ = u- = 0,  eta = h + b = 0 => h = -b
    h_ave = 0.5d0*(hL+hR)      ! arithmetic average of the depth
    ddphi = -g*h_ave*db        ! some approximation of the source term \int_{x_{l}}^{x_{r}} -g h b_x dx
    dphi_decomp = dphi - ddphi ! decompose the momentum flux jump in presence of source term

    ! flux decomposition
    beta1  = (s2*dhu - dphi_decomp)/(s2-s1)
    beta2  = (dphi_decomp - s1*dhu)/(s2-s1)

    ! speeds
    sw(1) = s1
    sw(2) = 0.5d0*(s1+s2)
    sw(3) = s2

    !  1st nonlinear wave
    fw(1,1) = beta1
    fw(2,1) = beta1*s1
    fw(3,1) = beta1*vL

    !  2nd nonlinear wave
    fw(1,3) = beta2
    fw(2,3) = beta2*s2
    fw(3,2) = beta2*vR

    ! advection of the transverse wave
    fw(1,2) = 0.d0
    fw(2,2) = 0.d0
    fw(3,3) = hR*uR*vR - hL*uL*vL - fw(3,1) - fw(3,3)
    
end subroutine riemann_fwave


!  Riemann type
!  @description: Determine the Riemann structure (wave-type in each family)
subroutine riemanntype(hL,hR,uL,uR,hm,s1m,s2m,rare1,rare2,maxiter,drytol,g)

    implicit none

    ! input
    double precision, intent(in) :: hL,hR,uL,uR,drytol,g
    integer, intent(in) :: maxiter

    ! output
    double precision, intent(out) :: s1m,s2m,hm
    logical, intent(out) :: rare1,rare2

    ! local
    double precision :: u1m,u2m,um,du, dhu
    double precision :: h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR
    integer :: iter

    ! Test for the Riemann structure
    h_min = min(hL,hR)
    h_max = max(hL,hR)
    du = uR - uL

    ! dry bed
    if (h_min <= drytol) then
        hm = 0.d0
        um = 0.d0
        s1m = uR + uL - 2.d0*sqrt(g*hR) + 2.d0*sqrt(g*hL)
        s2m = uR + uL - 2.d0*sqrt(g*hR) + 2.d0*sqrt(g*hL)

        if (hL <= 0.d0) then ! dry bed on the left
            rare2 = .true.
            rare1 = .false.
        else                 ! dry bed on the right
            rare1 = .true.
            rare2 = .false.
        end if

    else
        F_min = du + 2.d0*(sqrt(g*h_min) - sqrt(g*h_max))
        F_max = du + (h_max - h_min)*sqrt(0.5d0*g*(1.d0/h_max + 1.d0/h_min))

        if (F_min > 0.d0) then  !  2-rarefactions
            hm = (1.d0/16.d0*g)*max(0.d0, -du+2.d0*(sqrt(g*hL) + sqrt(g*hR)))**2
            um = sign(1.d0,hm)*(uL + 2.d0*sqrt(g*hL) - 2.d0*sqrt(g*hm))

            s1m = uL + 2.d0*sqrt(g*hL) - 3.d0*sqrt(g*hm)
            s2m = uR - 2.d0*sqrt(g*hR) + 3.d0*sqrt(g*hm)

            rare1 = .true.
            rare2 = .true.

        elseif (F_max <= 0.d0) then ! 2-shocks
            ! root finding using Newton's method on sqrt(h) 
            h0 = h_max
            do iter = 1,maxiter
                gL = sqrt(0.5d0*g*(1.d0/h0 + 1.d0/hL))
                gR = sqrt(0.5d0*g*(1.d0/h0 + 1.d0/hR))
                F0 = du + (h0 - hL)*gL - (h0 - hR)*gR
                dfdh = gL - g*(h0-hL)/(4.d0*(h0**2)*gL) + gR - g*(h0-hR)/(4.d0*(h0**2)*gR)
                slope = 2.d0*sqrt(h0)*dfdh
                h0 = (sqrt(h0) - F0/slope)**2
            end do
            hm = h0
            u1m = uL - (hm-hL)*sqrt((0.5d0*g)*(1.d0/hm + 1.d0/hL))
            u2m = uR + (hm-hR)*sqrt((0.5d0*g)*(1.d0/hm + 1.d0/hR))
            um = 0.5d0*(u1m + u2m)

            s1m = u1m - sqrt(g*hm)
            s2m = u2m + sqrt(g*hm)

            rare1 = .false.
            rare2 = .false.
        
        else ! 1-shock and 1-rarefaction
            h0 = h_min

            do iter=1,maxiter
                F0 = du + 2.d0*(sqrt(g*h0) - sqrt(g*h_max)) + (h0-h_min)*sqrt(0.5d0*g*(1.d0/h0 + 1.d0/h_min))
                slope = (F_max - F0)/(h_max - h_min)
                h0 = h0 - F0/slope
            end do

            hm = h0
            if (hL > hR) then
                um = uL + 2.d0*sqrt(g*hL) - 2.d0*sqrt(g*hm)
                s1m = uL + 2.d0*sqrt(g*hL) - 3.d0*sqrt(g*hm)
                s2m = uR - 2.d0*sqrt(g*hR) + sqrt(g*hm)
                rare1 = .true.
                rare2 = .false.
            else
                um = uR - 2.d0*sqrt(g*hR) + 2.d0*sqrt(g*hm)
                s1m = uR - 2.d0*sqrt(g*hR) + sqrt(g*hm)
                s2m = uR - 2.d0*sqrt(g*hR) + 3.d0*sqrt(g*hm)
                rare1 = .false.
                rare2 = .true.
            end if
        end if

        return

    end if

end subroutine riemanntype