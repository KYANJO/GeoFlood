! @description: Solves normal Riemann problems for the 2D SHALLOW WATER equations
!               with topography:
!               #        h_t + (hu)_x + (hv)_y = 0                           #
!               #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!               #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
!
!               On input, ql contains the state vector at the left edge of each cell
!               qr contains the state vector at the right edge of each cell
!           
!               This data is along a slice in the x-direction if ixy=1
!               or the y-direction if ixy=2.
!
!               Note that the i'th Riemann problem has left state qr(i-1,:)
!               and right state ql(i,:)
!               From the basic clawpack routines, this routine is called with
!               ql = qr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
 
subroutine fc2d_geoclaw_rpn2(ixy,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,auxl,auxr,fwave,s,amdq,apdq)

    use geoclaw_module, only: g => grav, drytol => dry_tolerance
    use geoclaw_module, only: earth_radius, deg2rad
    use amr_module, only: mcapa

    implicit none

    ! Input variables
    integer, intent(in) :: ixy, maxm, meqn, mwaves, mbc, mx

    ! Input/Output variables
    double precision, intent(inout) :: fwave(meqn,mwaves,1-mbc:maxm+mbc)
    double precision, intent(inout) :: s(mwaves,1-mbc:maxm+mbc)
    double precision, intent(inout) :: ql(meqn,1-mbc:maxm+mbc), qr(meqn,1-mbc:maxm+mbc)
    double precision, intent(inout) :: auxl(maux,1-mbc:maxm+mbc), auxr(maux,1-mbc:maxm+mbc)
    double precision, intent(inout) :: amdq(meqn,1-mbc:maxm+mbc), apdq(meqn,1-mbc:maxm+mbc)

    ! Local variables
    integer :: m,i,mw,maxiter,mu,nv
    double precision :: wall(3)
    double precision :: fw(3,3)
    double precision :: sw(3)

    double precision :: hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL
    double precision :: bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat
    double precision :: s1m,s2m
    double precision :: hstar,hstartest,hstarHLL,sLtest,sRtest
    double precision :: tw,dxdc

    logical :: rare1,rare2

    !  loop through Riemann problems at each grid cell
    do i=2-mbc,mx+mbc

        !  ==== initializing =====
        ! inform of a bad riemann proble from the start
        if((qr(1,i-1).lt.0.d0).or.(ql(1,i) .lt. 0.d0)) then
            write(*,*) 'Negative input: hl,hr,i=',qr(1,i-1),ql(1,i),i
        endif

        ! Initialize Riemann problem for grid interface
        do mw=1,mwaves
            s(mw,i)=0.d0
            fwave(1,mw,i)=0.d0
            fwave(2,mw,i)=0.d0
            fwave(3,mw,i)=0.d0
        enddo

        ! set normal direction
        if(ixy.eq.1) then
            mu=2
            nv=3
        else
            mu=3
            nv=2
        endif

        ! zero (small) negative values if they exist in the left state
        if(qr(1,i-1).lt.0.d0) then
            qr(1,i-1)=0.d0
            qr(2,i-1)=0.d0
            qr(3,i-1)=0.d0
        endif

        ! zero (small) negative values if they exist in the right state
        if(ql(1,i).lt.0.d0) then
            ql(1,i)=0.d0
            ql(2,i)=0.d0
            ql(3,i)=0.d0
        endif

        ! skip problem if in a completely dry carea
        if((qr(1,i-1).le.drytol).and.(ql(1,i).le.drytol)) then
            go to 30
        endif

        ! Riemann problem variables
        hL = qr(1,i-1)
        hR = ql(1,i)
        huL = qr(mu,i-1)
        huR = ql(mu,i)
        bL = auxr(1,i-1)
        bR = auxl(1,i)

        hvL= qr(nv,i-1)
        hvR= ql(nv,i)

        ! check for wet/dry fonts at the left interface 
        if (hR > drytol) then
            uR = huR/hR
            vR = hvR/hR
            phiR = 0.5d0*g*hR**2 + huR**2/hR
        else
            hR = 0.d0
            huR = 0.d0
            hvR = 0.d0
            uR = 0.d0
            vR = 0.d0
            phiR = 0.d0
        endif

        ! check for wet/dry fonts at the right interface
        if (hL > drytol) then
            uL = huL/hL
            vL = hvL/hL
            phiL = 0.5d0*g*hL**2 + huL**2/hL
        else
            hL = 0.d0
            huL = 0.d0
            hvL = 0.d0
            uL = 0.d0
            vL = 0.d0
            phiL = 0.d0
        endif

        


    end do







end subroutine fc2d_geoclaw_rpn2