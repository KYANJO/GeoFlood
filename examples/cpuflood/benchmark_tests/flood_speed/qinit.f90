! qinit routine for the dam_break problem
subroutine flood_speed_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none

    ! declare variables
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    !! local variables
    integer :: i,j
    real(kind=8) :: x,y,h

    do i=1-mbc,mx+mbc
        x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0)*dy
            if (abs(y - 1000.0d0)<=10.0d0) then
                h = 0.001d0
            else
                h = 0.0d0
            endif
            q(1,i,j) = 0.0d0 !to avoid the domain starting out dry 
            q(2,i,j) = 0.0d0
            q(3,i,j) = 0.0d0
        enddo
    enddo

end subroutine flood_speed_qinit