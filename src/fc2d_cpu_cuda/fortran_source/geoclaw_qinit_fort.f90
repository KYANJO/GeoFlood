
SUBROUTINE fc2d_geoclaw_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use qinit_module
    use geoclaw_module
    use amr_module, only: mcapa

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    double precision :: xhigher,yhigher,xintlow,xinthi,yintlow,yinthi,dq
    double precision :: x,y,xim,xip,yjm,yjp,xc,yc,xipc,ximc,yipc,yjmc,yjpc
   

    ! Locals
    integer :: i,j,m,mf,istart,iend,jstart,jend

    ! Topography integral function
    double precision :: topointegral_geo

    ! Set flat state based on sea_level
    q = 0.d0
    FORALL(i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall

    xhigher = xlower + (mx-0.5)*dx
    yhigher = ylower + (my-0.5)*dy

    ! ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    endif

    if (.false.) then
        open(23, file='fort.aux',status='unknown',form='formatted')
        print *,'Writing out aux arrays'
        print *,' '
        do j=1,my
            do i=1,mx
                write(23,*) i,j,(q(m,i,j),m=1,meqn)
            enddo
        enddo
        close(23)
    endif

  END SUBROUTINE fc2d_geoclaw_qinit
