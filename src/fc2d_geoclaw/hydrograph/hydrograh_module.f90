! Module containing hydrograph data
module hydrograph_module

    implicit none

    logical, private :: module_setup = .false.

    ! ========================================================================
    ! Hydrograph data
    ! ========================================================================
    character(len=20) :: hydrograph_type
    real(kind=8) :: u0,h0,hu0,hv0,num_rows
    real(kind=8), allocatable :: time(:), eta(:), hu(:)
    integer, parameter :: GEO_PARM_UNIT = 78

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
            read(unit,*) hu0, h0, u0

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
            write(GEO_PARM_UNIT,*) ' initial_conditons:', hu0, h0, u0
            write(GEO_PARM_UNIT,*) ' read_file:', read_file
            write(GEO_PARM_UNIT,*) ' hydrograph_type:', hydrograph_type
            write(GEO_PARM_UNIT,*) ' num_rows:', num_rows
            write(GEO_PARM_UNIT,*) ' time, eta, hu:', time, eta, hu
            write(GEO_PARM_UNIT,*) ' '
            
            module_setup = .true.
        end if         

    end subroutine read_hydrograph_data


end module hydrograph_module