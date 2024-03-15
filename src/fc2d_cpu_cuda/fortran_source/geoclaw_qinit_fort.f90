
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
    double precision :: topointegral

    ! Set flat state based on sea_level
    q = 0.d0
    FORALL(i=1-mbc:mx+mbc, j=1-mbc:my+mbc)
        q(1,i,j) = max(0.d0, sea_level - aux(1,i,j))
    end forall

    xhigher = xlower + (mx-0.5)*dx
    yhigher = ylower + (my-0.5)*dy


    ! Add perturbation to initial conditions
    if (qinit_type > 0) then
        do mf = 1,mqinitfiles
            if ((xlower.le.xhiqinit(mf).and.xhigher.ge.xlowqinit(mf)).and. &
                (ylower.le.yhiqinit(mf).and.yhigher.ge.ylowqinit(mf))) then
                
                xintlow = dmax1(xlower,xlowqinit(mf))
                xinthi  = dmin1(xhigher,xhiqinit(mf))
                istart  = max(1-mbc,int(0.5 + (xintlow-xlower)/dx))
                iend    = min(mx+mbc,int(1.0 + (xinthi-xlower)/dx))

                yintlow = dmax1(ylower,ylowqinit(mf))
                yinthi  = dmin1(yhigher,yhiqinit(mf))
                jstart  = max(1-mbc,int(0.5 + (yintlow-ylower)/dy))
                jend    = min(my+mbc,int(1.0 + (yinthi-ylower)/dy))

                do i=istart,iend
                x = xlower + (i-0.5d0)*dx
                xim = x - 0.5d0*dx
                xip = x + 0.5d0*dx
                do j=jstart,jend
                    y = ylower + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy

                    if (xip.gt.xlowqinit(mf).and.xim.lt.xhiqinit(mf) &
                        .and.yjp.gt.ylowqinit(mf) .and.yjm.lt.yhiqinit(mf)) then
                        
                            xipc=min(xip,xhiqinit(mf))
                            ximc=max(xim,xlowqinit(mf))
                            xc=0.5d0*(xipc+ximc)

                            yjpc=min(yjp,yhiqinit(mf))
                            yjmc=max(yjm,ylowqinit(mf))
                            yc=0.5d0*(yjmc+yjpc)

                             dq = topointegral(ximc,xc,xipc,yjmc,yc,yjpc, &
                                               xlowqinit(mf),ylowqinit(mf), &
                                               dxqinit(mf),dyqinit(mf), &
                                               mxqinit(mf),myqinit(mf), &
                                               qinitwork(i0qinit(mf):i0qinit(mf)+mqinit(mf)-1), &
                                               1)
                            if (coordinate_system == 2) then
                                dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(mcapa,i,j))
                            else
                                dq = dq / ((xipc-ximc)*(yjpc-yjmc))
                            endif 

                            if (iqinit(mf).le.meqn) then 
                                q(iqinit(mf),i,j) = q(iqinit(mf),i,j) + dq
                            else
                                q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                            end if
                        end if
                    
                    end do

                end do

            end if

        end do
    endif

    ! ! Add perturbation to initial conditions
    ! if (qinit_type > 0) then
    !     call add_perturbation(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    ! endif

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
