! Module containing hydrograph data
module hydrograph_module
    
    use geoclaw_module, only: grav
    
    implicit none
    
    logical, private :: module_setup = .false.

    ! ========================================================================
    ! Hydrograph data
    ! ========================================================================
    character(len=20) :: hydrograph_type
    real(kind=8), dimension(4) :: q0         ! q0 = [h0,hu0,hv0,b0] ghost data
    real(kind=8), dimension(4) :: q1         ! q1 = [h,hu,hv,b] cell data (initial conditions) just inside the boundary
    real(kind=8) :: u1            ! u1 = initial velocity (cell data) just inside the boundary
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
        integer :: i, unit = 127
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
                open(unit, file=hydrograph_file, status='old', action='read')
            else
                open(unit, file='hydrograph.data', status='old', action='read')
            end if

            !  read the first line (initial conditions)
            read(unit,*) q1(2), q1(4), u1, q1(1)

            !  asume momentum in y direction is zero
            q1(3) = 0.0d0

            ! compute the momentum in x direction (just right at the cell boundary)
            q1(2) = q1(1)*u1

            !  read the second line (reading from file or setrun.py)
            read(unit,*) read_file

            !  read the third line (hydrograph type)
            read(unit,*) hydrograph_type

            !  read the third line (number of rows)
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
            write(GEO_PARM_UNIT,*) ' initial_conditons:', q1(2), q1(1), u1, q1(4)
            write(GEO_PARM_UNIT,*) ' read_file:', read_file
            write(GEO_PARM_UNIT,*) ' hydrograph_type:', hydrograph_type
            write(GEO_PARM_UNIT,*) ' num_rows:', num_rows
            write(GEO_PARM_UNIT,*) ' time, eta, hu:', time, eta, hu
            write(GEO_PARM_UNIT,*) ' '
            
            module_setup = .true.

        end if         

    end subroutine read_hydrograph_data


    subroutine inflow_interpolate(t,q0)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: t
        real(kind=8), dimension(4), intent(out) :: q0

        if (hydrograph_type == 'discharge') then
            call interpolation(t,hu,q0(1))
        else
            call interpolation(t,eta,q0(4))
        end if


    end subroutine inflow_interpolate

    ! ========================================================================
    ! Interpolates the hydrograph data to the current time step
    ! ========================================================================
    subroutine interpolation(t,inflow,interpolated_value)

        implicit none

        ! Arguments
        real(kind=8), intent(in) :: t
        real(kind=8), dimension(:), intent(in) :: inflow
        real(kind=8), intent(out) :: interpolated_value

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
    end subroutine interpolation


    ! subroutine newton_raphson(h0,xn,h1,u1)
    subroutine newton_raphson(q0)

        implicit none

        ! declare variables
        real(kind=8) :: tol, x0
        real(kind=8) :: func,fxn,dfxn,dfunc_hu0,dfunc_h0
        real(kind=8), dimension(4), intent(out) :: q0

        integer :: i, max_iter

        ! initialize variables
        tol = 1.0e-6    ! tolerance for convergence
        max_iter = 100  ! maximum number of iterations
        x0 = 0.001d0    ! initial guess for the inflow discharge

        ! solve Riemann invariants
        if (hydrograph_type == 'discharge') then
            if (q0(2) == 0.0) then
                q0(1) = 0.0
            else
                q0(1) = x0
                do i = 1, max_iter
                    fxn = q0(2)/q0(1) - 2*sqrt(grav*q0(1)) - u1 + 2*sqrt(grav*q1(1))
                    if (abs(fxn) < tol) then
                        return 
                    end if
                    dfxn = -q0(2)/(q0(1)**2) - sqrt(grav/q0(1)) !dfunc_h0()
                    q0(1) = q0(1) - fxn/dfxn
                end do
                write(*,*) 'Newton-Raphson did not converge'
                q0(1) = 0.0  
            end if

        else

            if (q0(1) == 0.0) then
                q0(2) = 0.0
            else
                q0(2) = x0
                do i = 1, max_iter
                    fxn = q0(2)/q0(1) - 2*sqrt(grav*q0(1)) - u1 + 2*sqrt(grav*q1(1))
                    if (abs(fxn) < tol) then
                        return 
                    end if
                    dfxn =1/q0(1) !dfunc_hu0()
                    q0(2) = q0(2) - fxn/dfxn
                end do
                write(*,*) 'Newton-Raphson did not converge'
                q0(2) = 0.0  
            end if
        endif

     end subroutine newton_raphson  

end module hydrograph_module