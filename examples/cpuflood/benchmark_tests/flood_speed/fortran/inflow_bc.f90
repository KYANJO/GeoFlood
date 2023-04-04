program interpolate_boundary_condition_values

    ! This subroutine linearly interpolates the inflow boundary condition values from a file which are applied along a 20 m line in the middle of the western side of teh floodplain.

    implicit none

    ! declare variables
    double precision, dimension(:), allocatable :: tt, inflow
    double precision :: dx,dy,slope,x,y,t,inflow_interp,flow_depth

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

    ! create a 20m line segment along the western side of the floodplain
    dx = 0.05
    dy = 0.05
    slope = 0.01

    write(*,*) "input the t value for interpolation"
    ! t = 0.0
    read(*,*) t

    ! compute the slope
    ! xindex = int((x - x0) / dx)   ! index of midpoint along x-axis
    ! yindex = int((y - y0) / dy)   ! index of midpoint along y-axis
    ! slope = (z(xindex+1,yindex) - z(xindex-1,yindex)) / (2*dx)   ! slope at midpoint z is the elevation data

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
end program interpolate_boundary_condition_values

!  NRM routine to be used in the main program
subroutine  Newton_Raphson(slope,inflow_interp, flow_depth)
    implicit none

    ! declare variables
    double precision :: slope,inflow_interp, flow_depth
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
    
    write(*,*) "The flow depth is ", flow_depth
  
end subroutine Newton_Raphson



