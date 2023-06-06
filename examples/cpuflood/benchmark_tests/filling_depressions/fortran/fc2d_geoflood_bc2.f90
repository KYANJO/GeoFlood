! ==================================================================
subroutine filling_depressions_bc2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc)
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

    use hydrograph_module, only: inflow_interpolate

    implicit none

    integer, intent(in) :: meqn, mbc, mx, my, maux, mthbc(4)
    real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, dt

    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer :: m, i, j, ibc, jbc

    real(kind=8) ::  y
    real(kind=8), dimension(4) :: q0


    real(kind=8) :: h_, hu_ 

    real(kind=8) :: h1 = 0.01d0, u_1=0.01d0

    ! -------------------------------------------------------------------
    !  left boundary
    ! -------------------------------------------------------------------
    go to (100,110,120,130), mthbc(1)+1
    ! this is how we skip over this side... if (mthbc(1))+1 is not 1,2,3 or 4, then goto above walls through here ...
    goto 199

    100 continue
    ! user-supplied BC's (must be inserted!)
    !  in this case, we are using the inflow_interpolation subroutine to compute the inflow boundary condition values
    ! call inflow_interpolate(t,q0)
    ! write(*,*) 't = ', t, ' q0(1) = ', q0(1), ' q0(2) = ', q0(2),'q1(1) = ', q1(1), ' u1 = ', u1
    call read_file_interpolate('fortran/bc.txt', t, h_, hu_, h1, u_1,dx,dy,dt)
    ! write(*,*) 't = ', t, ' h_ = ', h_, ' hu_ = ', hu_, ' h1 = ', h1, ' u_1 = ', u_1
    ! call inflow_interpolate(t,q0)
    do j = 1-mbc,my+mbc
        y = ylower + (j-0.5d0)*dy
        do ibc=1,mbc
            if (abs(y-1900.0d0) <= 100.0d0) then
            
                aux(1,1-ibc,j) = aux(1,1,j)
                q(1,1-ibc,j) = h_  ! h               
                q(2,1-ibc,j) = hu_             
                q(3,1-ibc,j) = 0.0d0             ! hv vertical velocity = 0
            
            else

                aux(1,1-ibc,j) = aux(1,ibc,j)
                do m=1,meqn
                    q(m,1-ibc,j) = q(m,ibc,j)
                enddo
             
            end if
        enddo
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
end subroutine filling_depressions_bc2


subroutine read_file_interpolate(file_name, t, h0, hu0, h1, u1,dx,dy,dt)

    implicit none

    ! declare variables
    character(len=*), intent(in) :: file_name
    real(kind=8), dimension(:), allocatable :: time,z
    real(kind=8) :: t, zinterp, h0, h1 , u1,hu0, h,dx,dt,dy
    character(len=100) :: line
    real(kind=8) :: slope,b,n,Trap,prec,ptrap
    integer :: i,j,num_rows

    ! ----- read time and z from a file -----------------------
    !  open the file for reading
    open(10,file=file_name,status='old')

    ! count the number of rows in the file
    num_rows = 0
    do 
        read(10,*,iostat=i) line
        if (i /= 0) exit
        num_rows = num_rows + 1
    end do

    ! allocate memory for time and z
    allocate(time(num_rows),z(num_rows))

    ! rewind the file
    rewind(10)

    ! read data
    do i = 1,num_rows
        read(10,*) time(i), z(i)
        ! write(*,*) time(i), z(i)
    end do

    ! close the file
    close(10)

    ! ------ Linear interpolation -----------------------------

    ! initialize zinterp to zero
    zinterp = 0.0d0

    ! check if t is within the time range and set the value of zinterp
    if (t < time(1)) then
        zinterp = z(1)
    else if (t > time(size(time))) then
        zinterp = z(size(z))
    else
        do i = 1,size(time)-1
            if (t >= time(i) .and. t <= time(i+1)) then
                zinterp = z(i) + (((z(i+1) - z(i)) / (time(i+1) - time(i))) * (t - time(i)))
                exit
            end if
        end do
    end if

    ! write(*,*) 'The value of zinterp' , zinterp
    ! ----- end of linear interpolation ------------------------
    !
    ! ----- call the Riemann invariant subroutine --------------
    Trap = 100.0d0 + dx*(2.5d0 + 2.5d0)
    ! hu0 = zinterp/(100.0d0 + 2.5d0*dx)
    b = 100.0d0
    prec = b + 2.0d0*(0.5d0+dx)
    ! ptrap = b + dx*(sqrt(1+2.5d0**2) + sqrt(1+2.5d0**2))
    ptrap = b + dx*(sqrt(1 + (2.5d0**2)) )
    ! hu0 = zinterp/((b + Trap)/2.0d0) !<--- works for the first 2 rows
    hu0 = zinterp/prec
    call newton_raphson(h0,hu0,h1,u1)
    ! write(*,*) 'h_0 = ', h0, 'h_u0 = ', hu0
    ! slope = 1/3000.0d0
    ! slope = 0.001d0
    ! b = 100.0d0
    ! n = 0.03d0
    ! call comput_flow(zinterp,h0,hu0,slope,b,n,h1,u1)
    
    
    !  trying the triton approach of computing flow
    ! h0 = 1e-16
    ! h = (hu0 * dt) / (dx * dx)
    ! h0 = h0 + h
    ! write(*,*) 'h0 = ' , h0, 'h = ', h
    ! stop
    ! free up memory
    deallocate(time,z)
end subroutine read_file_interpolate

! subroutine to back computes flow depth based on Manning's equation
! Q : flow rate (m^3/s)
! h0 : flow depth (m)
! slope : bed slope
! b : channel width (m)
! n : Manning's roughness coefficient
subroutine comput_flow(Q,h0,hu0,slope,b,n, h1,u1)

    implicit none

    ! declare variables
    real(kind=8), intent(in) :: Q,slope,b,n,hu0,h1,u1
    real(kind=8), intent(out) :: h0
    real(kind=8) :: R,A,P,func,dfunc
    real(kind=8) :: tol = 1e-6
    real(kind=8) :: coef,nume1,deno1,nume2

    integer :: i,max_iter = 100

    ! Newton-Raphson method
    ! if (Q == 0.0d0) then
    !     call newton_raphson(h0,hu0,h1,u1)
    ! else
        do i = 1,max_iter
            ! Assuming a rectangular channel
            ! Cross-sectional area
            A = b * h0
            ! Wetted perimeter
            P = b + (2.0d0 * h0)
            ! Hydraulic radius
            R = A / P
            
            ! Manning's equation
            func = Q - (1/n) * A * R**(2.0d0/3.0d0) * (slope**(0.5d0))
            ! write(*,*) 'Q = ', Q, 'func = ', func, 'h0 = ', h0
            ! derivative of Manning's equation with respect to h0
            coef = 2*b*h0*sqrt(slope)
            nume1 = (b/(b+(2*h0))) - ((2*b*h0)/(b+(2*h0))**2)
            deno1 = (3*n)*((b*h0)/(b+(2*h0)))**(1/3.0d0)
            nume2 = (b*sqrt(slope))*(((b*h0)/(b+(2*h0)))**(2.0d0/3.0d0))

            dfunc = - coef*(nume1/deno1) - (nume2/n)

            ! check if dfunc is zero
            if (dfunc == 0.0d0) then
                write(*,*) 'dfunc is zero'
                stop
            end if

            ! update h0
            h0 = h0 - (func/dfunc)
            write(*,*) 'h_0 = ', h0, 'func = ', func, 'dfunc = ', dfunc
            stop
            ! check for convergence
            if (abs(func) < tol) then
                write(*,*) 'h0 = ', h0, 'hu0 = ', hu0
                stop
                return
            end if
        enddo
    ! end if
    
end subroutine comput_flow


subroutine newton_raphson(h0,hu0,h1,u1)

    implicit none

    ! declare variables
    real(kind=8) :: h0,h1,u1,x0,xn,tol,hu0
    real(kind=8) :: func,fxn,dfxn,dfunc_hu0,dfunc_h0

    integer :: i, max_iter

    ! initialize variables
    tol = 1.0e-6 ! tolerance for convergence
    max_iter = 100 ! maximum number of iterations
    x0 = 0.001d0 ! initial guess for the inflow discharge

    ! solve Riemann invariants
    ! if (hu0 == 0.0) then
    !     h0 = 0.01
    ! else
        xn = x0
        do i = 1, max_iter
            fxn = func(hu0,xn,h1,u1)
            if (abs(fxn) < tol) then
                h0 = xn
                return 
            end if
            ! dfxn = dfunc_hu0(xn,h0,h1,u1)
            dfxn = dfunc_h0(hu0,xn,h1,u1)
            xn = xn - fxn/dfxn
        end do
        write(*,*) 'Newton-Raphson did not converge'
        xn = 0.0
    ! endif
end subroutine newton_raphson

real(kind=8) function func(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    func = hu0/h0 - 2*sqrt(g*h0) - u1 + 2*sqrt(g*h1)

end function func

!  given hu0
real(kind=8) function dfunc_h0(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    dfunc_h0 = -hu0/(h0**2) - sqrt(g/h0)

end function dfunc_h0

! given h0
real(kind=8) function dfunc_hu0(hu0,h0,h1,u1)
    implicit none
    real(kind=8) :: hu0,h0,h1,u1
    real(kind=8) :: g
    
    g = 9.81d0 ! gravitational acceleration

    dfunc_hu0 = 1/h0

end function dfunc_hu0


