program bc
    implicit none
    double precision :: flow_depth

    ! call inflow_interpolation(flow_depth)

end program bc

! ==================================================================
subroutine fc2d_geoflood_bc2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,t,dt,mthbc)
!  ==================================================================

! standard boundary condition choices

! At each boundary k = 1 (left), 2 (right), 3 (top) or 4 (bottom):
!  if mthbc(k) = 0, user-supplied BC's (must be inserted!)
!              = 1, zero-order extrapolation
!              = 2, periodic BC's
            ! = 3,  solid walls, assuming this can be implemented by reflecting  the data about the boundary and then negating the 2'nd (for k=1,2) 0r 3'rd (for k=3,4) component of q.
! -------------------------------------------------------------------

! Extend the data from the interior cells (1:mx, 1:my) to the ghost cells outside the region:
! (i,1-jbc)  for jbc = 1,mbc, i = 1-mbc, mx+mbc
! (i,my+jbc) for jbc = 1,mbc, i = 1-mbc, mx+mbc
! (1-ibc,j)  for ibc = 1,mbc, j = 1-mbc, my+mbc
! (mx+ibc,j) for ibc = 1,mbc, j = 1-mbc, my+mbc

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, maux, mthbc(4)
    double precision, intent(in) :: xlower, ylower, dx, dy, t, dt

    double precision, dimension(meqn,1-mbc:mx+mbc,1-mbc:my+mbc), intent(inout) :: q
    double precision, dimension(maux,1-mbc:mx+mbc,1-mbc:my+mbc), intent(inout) :: aux

    integer :: m, i, j, ibc, jbc

    double precision :: flow_depth

    ! -------------------------------------------------------------------
    !  left boundary
    ! -------------------------------------------------------------------
    go to (100,110,120,130), mthbc(1)+1
    ! this is how we skip over this side... if (mthbc(1))+1 is not 1,2,3 or 4, then goto above walls through here ...
    goto 199

    100 continue
    ! user-supplied BC's (must be inserted!)
    !  in this case, we are using the inflow_interpolation subroutine to compute the inflow boundary condition values
    do 105 j = 1-mbc,my+mbc
        do 105 ibc=1,mbc
            aux(1,1-ibc,j) = aux(1,1,j)
            do 105 m=1,meqn
                !  apply only at the middle of the western side of the floodplain
                if (j == 1000+mbc) then
                    call inflow_interpolation(flow_depth,t,dx,dy,xlower,ylower,q,meqn,mbc,mx,my)
                    q(m,1-ibc,j) = flow_depth
                else
                    q(m,1-ibc,j) = q(m,1,j)
                end if
105         continue
    end do
    goto 199

    110 continue
    ! zero-order extrapolation
    do 115 j = 1-mbc,my+mbc
        do 115 ibc=1,mbc
            aux(1,1-ibc,j) = aux(1,1,j)
            do 115 m=1,meqn
                q(m,1-ibc,j) = q(m,1,j)
115         continue
    go to 199

    120 continue
    ! periodic BC's: handled by p4est
    goto 199

    130 continue
    









end subroutine fc2d_geoflood_bc2


subroutine inflow_interpolation(flow_depth,t,dx,dy,xlower,ylower,q,meqn,mbc,mx,my)

    ! This subroutine linearly interpolates the inflow boundary condition values from a file which are applied along a 20 m line in the middle of the western side of teh floodplain.

    implicit none

    ! declare variables
    double precision, dimension(:), allocatable :: tt, inflow
    double precision, dimension(meqn,1-mbc:mx+mbc,1-mbc:my+mbc), intent(inout) :: q
    double precision :: dx,dy,slope,x,y,t,inflow_interp, xlower, ylower
    double precision :: flow_depth                        ! intent(in,out)

    integer :: i,xindex,yindex,n = 5

    ! Open the file for reading
    open(10,file="bc.txt",status='old',action='read')

    ! read in the boundary condition values
    allocate(tt(n),inflow(n))
    do i=1,n
        read(10,*) tt(i), inflow(i)
    end do
    close(10)

    ! middle of the western side of the floodplain
    x = 0.0
    y = 1000.0

    ! compute the slope of the channel
    xindex = int((x - xlower) / dx)   ! index of midpoint along x-axis
    yindex = int((y - ylower) / dy)   ! index of midpoint along y-axis
    do i = 1,meqn
        slope = (q(i,xindex+1,yindex+1) - q(i,xindex-1,yindex-1)) / (2*dx)   ! slope at midpoint
    end do

    ! slope = 0.01

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
    call Newton_Raphson(slope,inflow_interp, flow_depth)

    ! free up memory
    deallocate(tt,inflow) 

! end program
end subroutine inflow_interpolation

!  NRM routine to be used in the main program
subroutine  Newton_Raphson(slope,inflow_interp, flow_depth)
    implicit none

    ! declare variables
    double precision, intent(in) :: slope,inflow_interp
    double precision, intent(out) :: flow_depth
    double precision :: Manning_coefficient,base_width,tol,R,func,dfunc
    integer :: i, max_iter

    ! initialize variables
    Manning_coefficient = 0.05 ! Manning's coefficient
    base_width = 20.0          ! base width of the channel
    tol = 1.0e-6               ! tolerance for convergence
    max_iter = 100             ! maximum number of iterations
    
    ! initial guess for the flow depth
    flow_depth = 0.1

    ! Newton-Raphson method
    if (inflow_interp == 0.0) then
        flow_depth = 0.0
    else
        do i = 1,max_iter
            R = (base_width*flow_depth)/(base_width + (2*flow_depth)) ! hydraulic radius
            func = flow_depth - ((inflow_interp*Manning_coefficient)/((R**(2/3))*base_width*sqrt(slope))) ! function to be solved
            dfunc = 1.0 + (((2*Manning_coefficient*inflow_interp)*(R**(1/3)))/(3*base_width*(flow_depth**2)*sqrt(slope))) ! derivative of the function
            flow_depth = flow_depth - func/dfunc ! update the flow depth
            if (abs(func) < tol) exit ! check for convergence
        end do
    end if
    
    ! write(*,*) "The flow depth is ", flow_depth
  
end subroutine Newton_Raphson



