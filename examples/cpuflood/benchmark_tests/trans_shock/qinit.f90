! qinit routine for the trans_shock problem
subroutine trans_shock_qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    implicit none

    ! declare variables
    integer, intent(in) :: meqn,mbc,mx,my,maux
    double precision, intent(in) :: xlower,ylower,dx,dy
    double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    !! computational domain
    real(kind=8), parameter :: xll = 0.0d0
    real(kind=8), parameter :: yll = -10.0d0
    real(kind=8), parameter :: xupper = 25.0d0
    real(kind=8), parameter :: yupper = 10.0d0

    !! local variables
    integer :: i,j
    real(kind=8) :: x,y

    do i=1-mbc,mx+mbc
        ! x = xlower + (i-0.5d0)*dx
        do j=1-mbc,my+mbc
            ! y = ylower + (j-0.5d0)*dy
            
            ! check if we are in the dam break region
            ! if (x .ge. xll .and. x .le. xupper .and. &
            !     y .ge. yll .and. y .le. yupper) then
          
                ! if (x == 0.0d0) then
                    q(2,i,j) = 0.18d0
                ! else if (x == 25.0d0) then
                    q(1,i,j) = 0.33d0
                    q(3,i,j) = 0.0d0
                ! end if
            ! end if

        enddo
    enddo

end subroutine trans_shock_qinit