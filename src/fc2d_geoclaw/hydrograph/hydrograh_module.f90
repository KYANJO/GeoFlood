! ========================================================================
! @author: Brian Kyanjo
! @date: 2023-11-7
! @brief: Module containing hydrograph data
! ========================================================================

module hydrograph_module
    
    use geoclaw_module, only: grav,dry_tolerance

    
    implicit none
    
    logical, private :: module_setup = .false.

    ! ========================================================================
    ! Hydrograph data
    ! ========================================================================
    character(len=20) :: use_hydrograph
    character(len=20) :: hydrograph_type
    character(len=100), dimension(4) :: boundary_location
    real(kind=8), dimension(4) :: q0         ! q0 = [h0,hu0,hv0,b0] ghost data
    real(kind=8), dimension(4) :: q1         ! q1 = [h,hu,hv,b] first interior cell data (initial conditions) just inside the boundary
    real(kind=8) :: u1            ! u1 = initial velocity (first cell data) just inside next to the boundary
    real(kind=8) :: froude         ! froude number
    real(kind=8) :: b,x0,y0        ! channel_width, channel center x and y coordinates
    real(kind=8), allocatable :: time(:), eta(:), hu(:)
    integer, parameter :: GEO_PARM_UNIT = 78
    integer :: num_rows

contains

    ! ========================================================================
    ! Reads in hydrograph data from file
    ! ========================================================================
    subroutine read_hydrograph_data(hydrograph_file)

        implicit none

        ! Arguments
        character(len=*), optional, intent(in) :: hydrograph_file
        character(len=20) :: read_file

        ! Local variables
        integer :: i, unit = 127, iostat
        character(len=256) :: line

        if (.not. module_setup) then
            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*)'--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Hydrograph Control Parameters:'
            write(GEO_PARM_UNIT,*) '------------------------------'

            ! if (present(hydrograph_file)) then
            !     call opendatafile(unit,hydrograph_file)
            ! else
            !     call opendatafile(unit,'hydrograph.data')
            ! end if
            ! open the file for reading
            if (present(hydrograph_file)) then
                open(unit, file=hydrograph_file, status='old', action='read', iostat=iostat)
            else
                open(unit, file='hydrograph.data', status='old', action='read', iostat=iostat)
            end if
            
            if (iostat /= 0) then
                write(GEO_PARM_UNIT,*) ' No hydrograph data provided'
                return
            end if

            !  read use hydrograph data
            read(unit,*) use_hydrograph

            !  read the first line (initial conditions)
            read(unit,*) q1(2), q1(4), u1, q1(1)

            ! read the channel width
            read(unit,*) b

            ! read the channel center coordinates
            read(unit,*) x0, y0

            ! read the boundary location
            read(unit,*) boundary_location(1), boundary_location(2), boundary_location(3), boundary_location(4)

            !  asume momentum in y direction is zero
            q1(3) = 0.0d0

            ! compute the momentum in x direction (just right at the cell boundary)
            q1(2) = q1(1)*u1

            !  read the second line (reading from file or setrun.py)
            read(unit,*) read_file

            !  read the third line (hydrograph type)
            read(unit,*) hydrograph_type

            ! read the fourth line (froude number)
            read(unit,*) froude

            !  read the fifth line (number of rows)
            read(unit,*) num_rows

            if (num_rows == 0) then
                write(GEO_PARM_UNIT,*) ' '
                write(GEO_PARM_UNIT,*) 'No hydrograph data provided'
                write(GEO_PARM_UNIT,*) ' '
                return
            end if

            ! Allocate arrays
            allocate(time(num_rows), eta(num_rows), hu(num_rows))
            if (read_file == 'False') then
                if (hydrograph_type == 'discharge' ) then
                    do i = 1,num_rows
                        read(unit,*) time(i), hu(i)
                    end do
                else
                    do i = 1,num_rows
                        read(unit,*) time(i), eta(i)
                    end do
                end if
            else
                if (hydrograph_type == 'discharge' ) then
                    do i = 1,num_rows
                        read(unit,*) time(i), hu(i), eta(i)
                    end do
                else
                    do i = 1,num_rows
                        read(unit,*) time(i), eta(i)
                    end do
                end if  
                
            end if

            !  write out data to parameter file
            ! write(GEO_PARM_UNIT,*) ' initial_conditons:', q1(2), q1(1), u1, q1(4)
            write(GEO_PARM_UNIT,*) ' channel_width:', b
            write(GEO_PARM_UNIT,*) ' channel_center_x: ', x0
            write(GEO_PARM_UNIT,*) ' channel_center_y: ', y0
            write(GEO_PARM_UNIT,*) ' boundary_location:', boundary_location(:)
            write(GEO_PARM_UNIT,*) ' read_file:', read_file
            write(GEO_PARM_UNIT,*) ' hydrograph_type:', hydrograph_type
            write(GEO_PARM_UNIT,*) ' num_rows:', num_rows
            write(GEO_PARM_UNIT,*) ' time, eta, hu:', time, eta, hu
            write(GEO_PARM_UNIT,*) ' '
            
            module_setup = .true.

        end if         

    end subroutine read_hydrograph_data


    subroutine inflow_interpolate(t,q0,q1)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: t
        real(kind=8), dimension(4) :: q0,q1

        if (hydrograph_type == 'discharge') then
            call interpolation(t,hu,q0,q1)
        else
            call interpolation(t,eta,q0,q1)
        end if

        ! call newton_raphson(q0) ! solve for the h or hu at the boundary

    end subroutine inflow_interpolate

    ! ========================================================================
    ! Interpolates the hydrograph data to the current time step
    ! ========================================================================
    subroutine interpolation(t,inflow,q0,q1)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: t
        real(kind=8), dimension(:), intent(in) :: inflow
        real(kind=8) :: interpolated_value
        real(kind=8), dimension(4) :: q0,q1

        ! Local variables
        integer :: i

        ! intialize variables
        interpolated_value = 0.0d0

        ! check if t is within the time range and set the value of interpolated_value
        if (t < time(1)) then
            interpolated_value = inflow(1)
        else if (t > time(size(time))) then
            interpolated_value = inflow(size(inflow))
        else
            do i = 1,size(time)-1
                if (t >= time(i) .and. t <= time(i+1)) then
                    interpolated_value = inflow(i) + (inflow(i+1) - inflow(i)) / (time(i+1) - time(i)) * (t - time(i))
                    exit
                end if
            end do
        end if

         if (hydrograph_type == 'discharge') then
            q0(2) = interpolated_value
            ! if (q0(2) .ne. 0.0) then
            !     call Riemann_invariants(q0,q1)
            !     if (q0(1) > q1(1)) then
            !         call two_shock(q0,q1)
            !     end if
            ! end if
            
        else
            q0(4) = interpolated_value  ! assuming eta is given
            ! call Riemann_invariants(q0,q1)

        end if

    end subroutine interpolation

    ! ========================================================================
    subroutine Riemann_invariants(q0,q1)

        implicit none

        ! declare variables
        real(kind=8) :: tol,u1
        real(kind=8) :: func,fxn,dfxn,dfunc_hu0,dfunc_h0,Fr
        real(kind=8), dimension(4) :: q0,q1

        integer :: i, max_iter

        ! initialize variables
        tol = 1.0e-6    ! tolerance for convergence
        max_iter = 100  ! maximum number of iterations
        u1 = q1(2)/q1(1) ! velocity of the first interior cells
        Fr = u1/sqrt(grav*q1(1)) ! Froude number

        ! solve Riemann invariants
        if (hydrograph_type == 'discharge') then
            q0(1) = (q0(2)/sqrt(grav)*Fr)**(2.0d0/3.0d0)
            do i = 1, max_iter
                fxn = q0(2)/q0(1) - 2*sqrt(grav*q0(1)) - u1 + 2*sqrt(grav*q1(1))
                if (abs(fxn) < tol) then
                    return 
                end if
                dfxn = -q0(2)/(q0(1)**2) - sqrt(grav/q0(1)) !dfunc_h0()
                q0(1) = q0(1) - fxn/dfxn
            end do
             
        else
            if (q0(1) == 0.0) then ! if h == 0 => hu == 0
                q0(2) = 0.0
            else
                q0(2) = (q0(2)/sqrt(grav)*Fr)**(2.0d0/3.0d0)
                do i = 1, max_iter
                    fxn = q0(2)/q0(1) - 2*sqrt(grav*q0(1)) - u1 + 2*sqrt(grav*q1(1))
                    if (abs(fxn) < tol) then
                        return 
                    end if
                    dfxn =1/q0(1) !dfunc_hu0()
                    q0(2) = q0(2) - fxn/dfxn
                end do
                
            end if
        endif

     end subroutine Riemann_invariants


    ! ======================================================================== 
    subroutine two_shock(q0,q1)

        implicit none

        ! declare variables
        real(kind=8) :: tol
        real(kind=8) :: fxn,dfxn,num,deno,Fr
        real(kind=8), dimension(4) :: q0,q1

        integer :: i, max_iter

        ! initialize variables
        tol = 1.0e-6    ! tolerance for convergence
        max_iter = 100  ! maximum number of iterations
        u1 = q1(2)/q1(1) ! velocity of the first interior cells
        Fr = u1/sqrt(grav*q1(1)) ! Froude number

        ! solve Riemann invariants
        if (hydrograph_type == 'discharge') then
            q0(1) = (q0(2)/sqrt(grav)*Fr)**(2.0d0/3.0d0)
            do i = 1, max_iter

                fxn = q0(2)/q0(1) - u1 - (q0(1) - q1(1))*sqrt((grav/2.0d0)*(1.0d0/q0(1) + 1.0d0/q1(1)))

                if (abs(fxn) < tol) then
                    return 
                end if
                num =  sqrt(grav*(1.0d0/q0(1) + 1.0d0/q1(1)))
                deno = 2*sqrt(2.0d0)*(q0(1)**2)*num
                dfxn = -q0(2)/(q0(1)**2) - num/sqrt(2.0d0) + grav*(q0(1) - q1(1))/deno
                q0(1) = q0(1) - fxn/dfxn
            end do
        
        else
            if (q0(1) == 0.0) then ! if h == 0 => hu == 0
                q0(2) = 0.0
            else
                q0(2) = (q0(2)/sqrt(grav)*Fr)**(2.0d0/3.0d0)
                do i = 1, max_iter
                    fxn = q0(2)/q0(1) - 2*sqrt(grav*q0(1)) - u1 + 2*sqrt(grav*q1(1))
                    if (abs(fxn) < tol) then
                        return 
                    end if
                    dfxn =1/q0(1) !dfunc_hu0()
                    q0(2) = q0(2) - fxn/dfxn
                end do
            
            end if
        endif

     end subroutine two_shock
     
end module hydrograph_module