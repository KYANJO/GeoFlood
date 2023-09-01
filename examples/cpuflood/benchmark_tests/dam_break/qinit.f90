! qinit routine for the dam_break problem
subroutine dam_break_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none

    ! declare variables
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    real(kind=8), parameter :: eta_upstream = 8.0d0
    real(kind=8), parameter :: eta_downstream = 0.4d0

    !! local variables
    integer :: i,j
    real(kind=8) :: x,y,eta

    do i=1-mbc,mx+mbc
        x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0)*dy
            
            ! eta = eta_upstream - (eta_upstream-eta_downstream)* &(x-xll)/(xupper-xll)
            if (x .le. 0.0d0) then
                eta = eta_upstream
            else
                eta = eta_downstream
            end if

            q(1,i,j) = max(0.0d0,eta-aux(1,i,j))
            q(2,i,j) = 0.0d0
            q(3,i,j) = 0.0d0
        enddo
    enddo

end subroutine dam_break_qinit