! qinit routine for the dam_break problem
subroutine flood_speed_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none

    ! declare variables
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    double  precision :: x0 = 500, y0 = 1000
    !! local variables
    integer :: i,j
    double precision :: x,y,h
! stop
    do i=1-mbc,mx+mbc
    
        x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0)*dy
            !q(1,i,j) = 0.1d0 + exp(-100.d0*(((x-x0)/1000)**2 + ((y-y0)/1000)**2)) 
            q(1,i,j) = 0.001d0
            q(2,i,j) = 0.0d0
            q(3,i,j) = 0.0d0
        enddo
    enddo

end subroutine flood_speed_qinit