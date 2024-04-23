! extract data from refinement.data for refinement flagging
module refinement_data_module

    use geoclaw_module, only: GEO_PARM_UNIT

    implicit none
    save

    ! public :: speed_tolerance
    logical, private :: module_setup = .false.
    

    ! ========================================================================
    !  Refinement Criteria
    ! ========================================================================
    double precision :: wave_tolerance
    double precision :: deep_depth
    integer :: max_level_deep
    double precision :: max_velocity_depth_product
    double precision, allocatable :: speed_tolerance(:)
    ! double precision :: speed_tolerance
    logical :: varRefTime = .FALSE. ! Choose dt refinement automatically

    ! ========================================================================
    !  Flowgrades for refinement
    ! ========================================================================
    integer :: num_flowgrades
    double precision, allocatable :: flowgradevalue(:)
    integer, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, allocatable :: iflowgrademinlevel(:)


contains

    !  read data from refinement.data
    subroutine read_refinement_data(file_name)

        use utility_module, only: get_value_count

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i, iostat
        character(len=128) :: line
        double precision :: value
    
        ! if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Refinement Control Parameters:'
            write(GEO_PARM_UNIT,*) '------------------------------'


            ! open file
            if (present(file_name)) then
                call opendatafile(iunit, file_name)
            else
                call opendatafile(iunit, 'refinement.data')
            endif

            ! read data
            read(iunit, *) wave_tolerance
            read(iunit, *) deep_depth
            read(iunit, *) max_level_deep
            read(iunit, *) max_velocity_depth_product
            read(iunit,'(a)') line
            allocate(speed_tolerance(get_value_count(line)))
            read(line,*) speed_tolerance
            
            read(iunit,*)
            read(iunit,*) varRefTime
            close(iunit)

            ! Write out data to parameter file
            write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
            write(GEO_PARM_UNIT,*) '   deep_depth:',deep_depth
            write(GEO_PARM_UNIT,*) '   max_level_deep:',max_level_deep
            write(GEO_PARM_UNIT,*) '   max_velocity_depth_product:',max_velocity_depth_product
            write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
            write(GEO_PARM_UNIT,*) '   Variable dt Refinement Ratios:',varRefTime
            write(GEO_PARM_UNIT,*) ''

            ! module_setup = .true.
        ! end if

    end subroutine read_refinement_data

    subroutine set_flow_grades(file_name)

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i, iostat
        character(len=128) :: comment_marker = '=:'

        ! if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'SET FLOW GRADES:'
            write(GEO_PARM_UNIT,*) '------------'

            ! Read user parameters from setflowgrades.data
            if (present(file_name)) then
                call opendatafile(iunit, file_name)
            else
                ! call opendatafile(iunit, 'setflowgrades.data')
                open(iunit, file='setflowgrades.data', status='old', action='read', iostat=iostat)
            endif
            
            if (iostat /= 0) then
                write(GEO_PARM_UNIT,*) '  No setflowgrades.data file found'
                return
            endif

            read(iunit,*) num_flowgrades

            if (num_flowgrades == 0) then
                write(GEO_PARM_UNIT,*) '  No flow grades specified'
                return
            endif

            ! Allocate arrays
            allocate(flowgradevalue(num_flowgrades),iflowgradevariable(num_flowgrades))
            allocate(iflowgradetype(num_flowgrades),iflowgrademinlevel(num_flowgrades))

            do i=1,num_flowgrades
                read(iunit,*) flowgradevalue(i),iflowgradevariable(i), &
                    iflowgradetype(i),iflowgrademinlevel(i)
            enddo

            close(iunit)

            write(GEO_PARM_UNIT,*) '   mflowgrades:',  num_flowgrades

            do i=1,num_flowgrades
                write(GEO_PARM_UNIT,"(d12.3,3i4)") flowgradevalue(i), &
                    iflowgradevariable(i),iflowgradetype(i),iflowgrademinlevel(i)

            enddo
            ! module_setup = .true.
        ! end if

    end subroutine set_flow_grades

end module refinement_data_module