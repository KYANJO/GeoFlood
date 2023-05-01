! ==================================================================
subroutine flood_speed_bc2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc)
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
    real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, dt

    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: m, i, j, ibc, jbc

    real(kind=8) :: flow_depth, inflow_interp, y

    ! -------------------------------------------------------------------
    !  left boundary
    ! -------------------------------------------------------------------
    go to (100,110,120,130), mthbc(1)+1
    ! this is how we skip over this side... if (mthbc(1))+1 is not 1,2,3 or 4, then goto above walls through here ...
    goto 199

    100 continue
    ! user-supplied BC's (must be inserted!)
    !  in this case, we are using the inflow_interpolation subroutine to compute the inflow boundary condition values
    do j = 1-mbc,my+mbc
        
            !  apply only at the middle of the western side of the floodplain
            y  = ylower + (j-0.5d0)*dy
            if (abs(y-1000) < 10) then
                call inflow_interpolation(flow_depth,inflow_interp,t,dx,dy,xlower,ylower,q,meqn,mbc,mx,my)
                ! write (*,*) 'flow_depth = ', flow_depth, ' inflow_interp = ', inflow_interp
                do ibc=1,mbc
                    q(1,1-ibc,j) = flow_depth         ! h
                    q(2,1-ibc,j) = inflow_interp/20.0d0 ! hu = flow_interp/base_width
                    q(3,1-ibc,j) = 0.0d0              ! hv vertical velocity = 0
                end do
            else
                q(m,1-ibc,j) = q(m,1,j)
            end if
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
    ! solid wall (assumes 2'nd component is velocity or momentum in x)
         do 135 j = 1-mbc, my+mbc
         do 135 ibc=1,mbc
            aux(1,1-ibc,j) = aux(1,ibc,j)
            do 135 m=1,meqn
               q(m,1-ibc,j) = q(m,ibc,j)
135       continue
! c     # negate the normal velocity:
        do 136 j = 1-mbc, my+mbc
            do 136 ibc=1,mbc
            q(2,1-ibc,j) = -q(2,ibc,j)
136    continue
        go to 199

199 continue
! c
! c-------------------------------------------------------
! c     # right boundary:
! c-------------------------------------------------------
    go to (200,210,220,230) mthbc(2)+1
    goto 299
! c
200 continue
! c     # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
    stop
    go to 299

210 continue
! c     # zero-order extrapolation:
    do 215 j = 1-mbc, my+mbc
        do 215 ibc=1,mbc
        aux(1,mx+ibc,j) = aux(1,mx,j)
        do 215 m=1,meqn
            q(m,mx+ibc,j) = q(m,mx,j)
215       continue
    go to 299

220 continue
! c     # periodic : Handled elsewhere
    go to 299

230 continue
! c     # solid wall (assumes 2'nd component is velocity or momentum in x):
    do 235 j = 1-mbc, my+mbc
        do 235 ibc=1,mbc
        aux(1,mx+ibc,j) = aux(1,mx+1-ibc,j)
        do 235 m=1,meqn
            q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
235       continue
! c     # negate the normal velocity:
    do 236 j = 1-mbc, my+mbc
        do 236 ibc=1,mbc
        q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
236    continue
    go to 299

299 continue
! c
! c-------------------------------------------------------
! c     # bottom boundary:
! c-------------------------------------------------------
    go to (300,310,320,330) mthbc(3)+1
    goto 399
! c
300 continue
! c     # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
    stop
    go to 399
! c
310 continue
! c     # zero-order extrapolation:
    do 315 jbc=1,mbc
        do 315 i = 1-mbc, mx+mbc
        aux(1,i,1-jbc) = aux(1,i,1)
        do 315 m=1,meqn
            q(m,i,1-jbc) = q(m,i,1)
315       continue
    go to 399

320 continue
! c     # periodic: Handled elsewhere
    go to 399

330 continue
! c     # solid wall (assumes 3'rd component is velocity or momentum in y):
    do 335 jbc=1,mbc
        do 335 i = 1-mbc, mx+mbc
        aux(1,i,1-jbc) = aux(1,i,jbc)
        do 335 m=1,meqn
            q(m,i,1-jbc) = q(m,i,jbc)
335       continue
! c     # negate the normal velocity:
    do 336 jbc=1,mbc
        do 336 i = 1-mbc, mx+mbc
        q(3,i,1-jbc) = -q(3,i,jbc)
336    continue
    go to 399

399 continue
! c
! c-------------------------------------------------------
! c     # top boundary:
! c-------------------------------------------------------
    go to (400,410,420,430) mthbc(4)+1
    goto 499

400 continue
    !  # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
    stop
    go to 499

410 continue
!     # zero-order extrapolation:
    do 415 jbc=1,mbc
        do 415 i = 1-mbc, mx+mbc
        aux(1,i,my+jbc) = aux(1,i,my)
        do 415 m=1,meqn
            q(m,i,my+jbc) = q(m,i,my)
415       continue
    go to 499

420 continue
!     # periodic: Handled elsewhere
    go to 499

430 continue
!  solid wall (assumes 3'rd component is velocity or momentum in y):
    do 435 jbc=1,mbc
        do 435 i = 1-mbc, mx+mbc
        aux(1,i,my+jbc) = aux(1,i,my+1-jbc)
        do 435 m=1,meqn
            q(m,i,my+jbc) = q(m,i,my+1-jbc)
435       continue
!  # negate the normal velocity:
    do 436 jbc=1,mbc
        do 436 i = 1-mbc, mx+mbc
        q(3,i,my+jbc) = -q(3,i,my+1-jbc)
436    continue
    go to 499

499 continue

      return
end subroutine flood_speed_bc2


subroutine inflow_interpolation(flow_depth,inflow_interp,t,dx,dy,xlower,ylower,q,meqn,mbc,mx,my)

    ! This subroutine linearly interpolates the inflow boundary condition values from a file which are applied along a 20 m line in the middle of the western side of teh floodplain.

    implicit none

    ! declare variables
    integer, intent(in) :: meqn, mbc, mx, my
    real(kind=8), intent(in) :: xlower, ylower, dx, dy, t
    real(kind=8), dimension(:), allocatable :: tt, inflow
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8) :: slope,x,y,inflow_interp
    real(kind=8) :: flow_depth , h1 =0, u1 =0                       ! intent(in,out)

    integer :: i,xindex,yindex,n = 5

    ! Open the file for reading
    open(10,file="fortran/bc.txt",status='old',action='read')

    ! read in the boundary condition values
    allocate(tt(n),inflow(n))
    do i=1,n
        read(10,*) tt(i), inflow(i)
    end do
    close(10)

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
    call Riemann_invariants(inflow_interp,flow_depth,h1,u1)
    WRITE (*,*) 'flow_depth = ', flow_depth, 'T = ', t
    ! free up memory
    deallocate(tt,inflow) 

! end program
end subroutine inflow_interpolation
  
! NRM  routine to solve Riemann invariants
subroutine Riemann_invariants(hu0,h0,h1,u1)

implicit none

! declare variables
real(kind=8), intent(in) :: hu0
real(kind=8), intent(out) :: h0
real(kind=8) :: h1,u1,g,func,dfunc,tol

integer :: i, max_iter

! initialize variables
g = 9.81 ! gravitational acceleration
tol = 1.0e-6 ! tolerance for convergence
max_iter = 100 ! maximum number of iterations
h0 = 1

! solve Riemann invariants
if (hu0 == 0.0) then
    h0 = 0.0
else
    do i = 1,max_iter
        func = hu0/h0 - 2*sqrt(g*h0) - u1 +2*sqrt(g*h1) ! function to be solved

        dfunc = -hu0/(h0**2) - sqrt(g/h0)  

        if (dfunc == 0.0) then
            write(*,*) "The derivative is zero"
            exit
        end if

        h0 = h0 - func/dfunc ! update the flow depth

        if (abs(func) < tol) exit ! check for convergence
    end do
end if

end subroutine Riemann_invariants