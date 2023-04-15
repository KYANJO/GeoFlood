program NRM
 implicit none

    ! declare variables
    real(kind =8) :: inflow_discharge, inflow_interp              ! intent(in,out)

    ! call the subroutine
    call inflow_interpolation(inflow_interp,inflow_discharge)

end program NRM



subroutine inflow_interpolation(hu1,inflow_interp)

    ! This subroutine linearly interpolates the inflow boundary condition values from a file which are applied along a 20 m line in the middle of the western side of teh floodplain.

    implicit none

    ! declare variables
    ! integer, intent(in) :: meqn, mbc, mx, my
    ! real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
    real(kind=8), dimension(:), allocatable :: tt, inflow
    ! real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: t,inflow_interp
    real(kind=8) :: hu1, h1 = 9.7 , u1= 0                 ! intent(in,out)

    integer :: i,xindex,yindex,n = 5

    ! Open the file for reading
    open(10,file="bc.txt",status='old',action='read')

    ! read in the boundary condition values
    allocate(tt(n),inflow(n))
    do i=1,n
        read(10,*) tt(i), inflow(i)
    end do
    close(10)

    !  read the value of t 
    write(*,*) 'Enter the value of t: '
    read(*,*) t

    !  find the nearest time values
    do i = 1,n-1
        if (t >= tt(i) .and. t <= tt(i+1)) then
            exit
        end if
    end do

    ! linearly interpolate the inflow boundary condition values 
    inflow_interp = inflow(i) + (inflow(i+1) - inflow(i)) / (tt(i+1) - tt(i)) * (t - tt(i))

    ! constraint the inflow boundary condition values to be positive
    if (inflow_interp < 0.0) then
        inflow_interp = 0.0
    end if

    ! Use Newton-Raphson method to compute the inflow depth
    call Riemann_invariants(inflow_interp,hu1,h1,u1)

    WRITE (*,*) 'Inflow_depth_interp =', inflow_interp, 'inflow_discharge = ', hu1, 'T = ', t
    ! free up memory
    deallocate(tt,inflow) 

! end program
end subroutine inflow_interpolation


! NRM  routine to solve Riemann invariants
subroutine Riemann_invariants(h0,hu0,h1,u1)

implicit none

! declare variables
real(kind=8) :: hu0,h0,h1,u1
real(kind=8) :: g,func,dfunc,tol

integer :: i, max_iter

! initialize variables
g = 9.81 ! gravitational acceleration
tol = 1.0e-6 ! tolerance for convergence
max_iter = 100 ! maximum number of iterations
hu0 = 0.1 ! initial guess for the inflow discharge

! solve Riemann invariants
! if (hu0 == 0.0) then
if (h0 == 0.0) then
    hu0 = 0.0
else
    do i = 1,max_iter
        func = hu0/h0 - 2*sqrt(g*h0) - u1 +2*sqrt(g*h1) ! function to be solved

        ! dfunc = -hu0/(h0**2) - sqrt(g/h0)   ! when hu0 is provided
        dfunc = 1.d0/h0 ! when hu0 is not provided, i.e. h0 is provided

        hu0 = hu0 - func/dfunc ! update the flow depth

        if (abs(func) < tol) exit ! check for convergence
    end do
end if

end subroutine Riemann_invariants
