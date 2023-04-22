subroutine disconnected_water_body_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none

    ! declare variables
    integer, intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
    real(kind=8), intent(inout) :: q(meqn,mbc+mx,mbc+my)
    real(kind=8), intent(in) :: aux(maux,mbc+mx,mbc+my)
    integer :: i,j

    ! problem parameters
    real(kind=8) :: ax = -10.75d0
    real(kind=8) :: bx = 710.75d0
    real(kind=8) :: ay = -18.25d0
    real(kind=8) :: by = 118.25d0

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            q(1,i,j) = 9.7d0
            q(2,i,j) = 0.d0
            q(3,i,j) = 0.d0
        enddo

    enddo

 
end subroutine disconnected_water_body_qinit