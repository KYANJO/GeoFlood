! Module containing refinement flagging criteria
module refinement_module

    use geoclaw_module, only: GEO_PARM_UNIT

    implicit none
    save

    logical, private :: module_setup = .false.

    ! ========================================================================
    !  Refinement Criteria
    ! ========================================================================
    double precision :: wave_tolerance
    double precision, allocatable :: speed_tolerance(:)
    logical :: varRefTime = .FALSE. ! Choose dt refinement automatically
    
    ! ========================================================================
    !  Flowgrades - Not updated yet, use at your own risk
    ! ========================================================================
    integer :: num_flowgrades
    double precision, allocatable :: flowgradevalue(:)
    integer, allocatable :: iflowgradevariable(:), iflowgradetype(:)
    integer, allocatable :: iflowgrademinlevel(:)

contains
    
    ! =========================================================================
    !  Reads in the refinement control parameters
    ! =========================================================================
    subroutine set_refinement(file_name)
        
        use utility_module, only: get_value_count
        
        implicit none
        
        ! Arguments
        character(len=*), optional, intent(in) :: file_name
        
        ! Locals
        integer, parameter :: unit = 127
        character(len=128) :: line

        if (.not.module_setup) then

            write(GEO_PARM_UNIT,*) ' '
            write(GEO_PARM_UNIT,*) '--------------------------------------------'
            write(GEO_PARM_UNIT,*) 'Refinement Control Parameters:'
            write(GEO_PARM_UNIT,*) '------------------------------'

            if (present(file_name)) then
                call opendatafile(unit, file_name)
            else
                call opendatafile(unit, 'refinement.data')
            endif

            ! Basic criteria
            read(unit,*) wave_tolerance
            read(unit,'(a)') line
            allocate(speed_tolerance(get_value_count(line)))
            read(line,*) speed_tolerance
            read(unit,*)
            read(unit,*) varRefTime
            close(unit)
            
            ! Write out data to parameter file
            write(GEO_PARM_UNIT,*) '   wave_tolerance:',wave_tolerance
            write(GEO_PARM_UNIT,*) '   speed_tolerance:',speed_tolerance
            write(GEO_PARM_UNIT,*) '   Variable dt Refinement Ratios:',varRefTime
            write(GEO_PARM_UNIT,*) ''
            
            module_setup = .true.
        end if

    end subroutine set_refinement
    
    
    ! =========================================================================
    ! TODO: This needs to be updated for the new module
    ! =========================================================================
    subroutine set_flow_grades(file_name)

        implicit none

        ! Input arguments
        character(len=*), optional, intent(in) :: file_name

        ! Locals
        integer, parameter :: iunit = 127
        integer :: i, iostat
        character(len=128) :: comment_marker = '=:'

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

    end subroutine set_flow_grades
    
end module refinement_module
