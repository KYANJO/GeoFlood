module qinit_module

    use amr_module, only: rinfinity

    implicit none
    save

    logical :: module_setup = .false.
    
    ! Type of q initialization
    integer, public :: qinit_type
    
    ! Work array
    double precision, private, allocatable :: qinit(:)

    ! Geometry
    double precision :: x_low_qinit
    double precision :: y_low_qinit
    double precision :: t_low_qinit
    double precision :: x_hi_qinit
    double precision :: y_hi_qinit
    double precision :: t_hi_qinit
    double precision :: dx_qinit
    double precision :: dy_qinit
    
    integer, private :: mx_qinit
    integer, private :: my_qinit

    ! for initializing using force_dry to indicate dry regions below sealevel:

    integer :: mx_fdry, my_fdry
    double precision :: xlow_fdry, ylow_fdry, xhi_fdry, yhi_fdry, dx_fdry, dy_fdry
    integer(kind=1), allocatable :: force_dry(:,:)
    logical :: use_force_dry
    double precision :: tend_force_dry  ! always use mask up to this time

    logical :: variable_eta_init

    ! to initialize using different initial eta values in different regions:
    integer :: etain_mx, etain_my
    double precision :: etain_dx, etain_dy
    double precision, allocatable :: etain_x(:), etain_y(:), etain_eta(:,:)

    ! integer :: mqinitfiles
    integer, allocatable :: minlevel_qinit(:), maxlevel_qinit(:)

    ! Work array
      double precision, allocatable :: qinitwork(:)

      ! Topography file data
      character*150, allocatable :: qinitfname(:)
      integer :: mqinitfiles,mqinitsize
      double precision, allocatable :: xlowqinit(:), ylowqinit(:)
      double precision, allocatable :: xhiqinit(:), yhiqinit(:)
      double precision, allocatable :: dxqinit(:), dyqinit(:)

      integer, allocatable ::  mxqinit(:), myqinit(:)
      integer, allocatable :: i0qinit(:), mqinit(:)
      integer, allocatable :: iqinit(:), qinitftype(:)
      integer, allocatable ::  minlevelqinit(:), maxlevelqinit(:)


contains

    subroutine set_qinit(fname)

      use geoclaw_module

      implicit none

      ! Input arguments
      character*25, intent(in), optional :: fname
      character*25 :: fname_force_dry

      ! Locals
      integer, parameter :: iunit = 7
      integer :: i,j,iqinitfile
      character*25 :: file_name
      logical :: found_file

      integer :: num_force_dry

      ! Open and begin parameter file output
      write(GEO_PARM_UNIT,*) ' '
      write(GEO_PARM_UNIT,*) '--------------------------------------------'
      write(GEO_PARM_UNIT,*) 'SETQINIT:'
      write(GEO_PARM_UNIT,*) '---------'

      if (present(fname)) then
         file_name = fname
      else
         file_name  = 'qinit.data'
      endif
      inquire(file=file_name,exist=found_file)
      if (.not. found_file) then
         print *, 'You must provide a file ', file_name
         stop
      endif

      call opendatafile(iunit, file_name)

      read(iunit,*) mqinitfiles

      if (mqinitfiles==0) then
         write(GEO_PARM_UNIT,*) '   mqinitfiles = 0'
         write(GEO_PARM_UNIT,*) '   no initial perturbation = 0'
         write(GEO_PARM_UNIT,*) '   h will be set max(0-b,0)   '
         return
      endif

      write(GEO_PARM_UNIT,*) '   mqinitfiles = ',mqinitfiles

      ! Read and allocate data parameters for each file
      allocate(mxqinit(mqinitfiles),myqinit(mqinitfiles))
      allocate(xlowqinit(mqinitfiles),ylowqinit(mqinitfiles))
      allocate(xhiqinit(mqinitfiles),yhiqinit(mqinitfiles))
      allocate(dxqinit(mqinitfiles),dyqinit(mqinitfiles))
      allocate(minlevelqinit(mqinitfiles),maxlevelqinit(mqinitfiles))
      allocate(qinitfname(mqinitfiles),qinitftype(mqinitfiles))
      allocate(iqinit(mqinitfiles))
      allocate(i0qinit(mqinitfiles),mqinit(mqinitfiles))

      do i=1,mqinitfiles
         read(iunit,*) qinitfname(i)
         read(iunit,*) qinitftype(i),iqinit(i),minlevelqinit(i), maxlevelqinit(i)

         write(GEO_PARM_UNIT,*) '   '
         write(GEO_PARM_UNIT,*) '   ',qinitfname(i)
         write(GEO_PARM_UNIT,*) '  qinitftype = ', qinitftype(i)
         write(GEO_PARM_UNIT,*) '  iqinit = ', iqinit(i)
         write(GEO_PARM_UNIT,*) '  minlevel, maxlevel = ', &
                                  minlevelqinit(i), maxlevelqinit(i)

         call read_qinit_header(qinitfname(i),qinitftype(i),mxqinit(i), &
                myqinit(i),xlowqinit(i),ylowqinit(i),xhiqinit(i),yhiqinit(i), &
                dxqinit(i),dyqinit(i))
            mqinit(i) = mxqinit(i)*myqinit(i)
      enddo

      ! Indexing into work array
      i0qinit(1)=1
      if (mqinitfiles > 1) then
         do i=2,mqinitfiles
            i0qinit(i)=i0qinit(i-1) + mqinit(i-1)
         enddo
      endif

    !  if variable_eta_init then function set_eta_init is called to set
    !  variable eta when interpolating onto newly refined patches
      read(iunit,*) variable_eta_init

      read(iunit,*) num_force_dry
      use_force_dry = (num_force_dry > 0)

      if (num_force_dry > 1) then
         write(GEO_PARM_UNIT,*) '*** num_force_dry > 1 not yet implemented'
         stop
      endif

      if (use_force_dry) then
          read(iunit,*) fname_force_dry
          read(iunit,*) tend_force_dry
          call read_force_dry(trim(fname_force_dry))
      endif

      ! Read and allocate space in work array for each file
      mqinitsize = sum(mqinit)
      allocate(qinitwork(mqinitsize))

      do i=1,mqinitfiles
            call read_qinit(mxqinit(i),myqinit(i),qinitftype(i),qinitfname(i), &
                qinitwork(i0qinit(i):i0qinit(i)+mqinit(i)-1))
      enddo

   end subroutine set_qinit


    subroutine add_perturbation(meqn,mbc,mx,my,xlow_patch,ylow_patch,dx,dy,q,maux,aux)
    
        use geoclaw_module, only: sea_level, coordinate_system
        use amr_module, only: mcapa
    
        implicit none
    
        ! Subroutine arguments
        integer, intent(in) :: meqn,mbc,mx,my,maux
        double precision, intent(in) :: xlow_patch,ylow_patch,dx,dy
        double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
        double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
        
        ! Local
        integer :: i,j
        double precision :: ximc,xim,x,xip,xipc,yjmc,yjm,y,yjp,yjpc,dq
        
        ! Topography integral function
        double precision :: topointegral
        
        if (qinit_type > 0) then
            do i=1-mbc,mx+mbc
                x = xlow_patch + (i-0.5d0)*dx
                xim = x - 0.5d0*dx
                xip = x + 0.5d0*dx
                do j=1-mbc,my+mbc
                    y = ylow_patch + (j-0.5d0)*dy
                    yjm = y - 0.5d0*dy
                    yjp = y + 0.5d0*dy

                    ! Check to see if we are in the qinit region at this grid point
                    if ((xip > x_low_qinit).and.(xim < x_hi_qinit).and.  &
                        (yjp > y_low_qinit).and.(yjm < y_hi_qinit)) then

                        xipc=min(xip,x_hi_qinit)
                        ximc=max(xim,x_low_qinit)

                        yjpc=min(yjp,y_hi_qinit)
                        yjmc=max(yjm,y_low_qinit)

                        dq = topointegral(ximc,xipc,yjmc,yjpc,x_low_qinit, &
                                          y_low_qinit,dx_qinit,dy_qinit,mx_qinit, &
                                          my_qinit,qinit,1)
                        if (coordinate_system == 2) then
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc)*aux(mcapa,i,j))
                        else
                            dq = dq / ((xipc-ximc)*(yjpc-yjmc))
                        endif 

                        if (qinit_type < 4) then 
                            if (aux(1,i,j) <= sea_level) then
                                q(qinit_type,i,j) = q(qinit_type,i,j) + dq
                            endif
                        else if (qinit_type == 4) then
                            q(1,i,j) = max(dq-aux(1,i,j),0.d0)
                        endif
                    endif
                enddo
            enddo
        endif
        
    end subroutine add_perturbation

        
    ! currently only supports one file type:
    ! x,y,z values, one per line in standard order from NW corner to SE
    ! z is perturbation from standard depth h,hu,hv set in qinit_geo,
    ! if iqinit = 1,2, or 3 respectively.
    ! if iqinit = 4, the z column corresponds to the definition of the 
    ! surface elevation eta. The depth is then set as q(i,j,1)=max(eta-b,0)
    subroutine read_qinit(mx,my,filetype,fname,qinit)

        use geoclaw_module

        implicit none

        ! Arguments
        integer, intent(in) :: mx,my,filetype
        character*150, intent(in) :: fname
        double precision, intent(inout) :: qinit(1:mx*my)

        ! Locals
        integer, parameter :: iunit = 19, miss_unit = 17
        double precision, parameter :: qinit_missing = -150.d0
        logical, parameter :: maketype2 = .false.
        integer :: i,j,num_points,missing,status,qinit_start
        double precision :: no_data_value,x,y,z

        print *, ' '
        print *, 'Reading qinit file  ', fname

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(filetype))
            ! ASCII file with x,y,z values on each line.
            ! (progressing from upper left corner across rows, then down)
            ! Assumes a uniform rectangular grid of data values.
            case(1)
                i = 0
                status = 0
                do while (status == 0)
                    i = i + 1
                    read(iunit,fmt=*,iostat=status) x,y,qinit(i)
                enddo

            ! ================================================================
            ! ASCII file with header followed by z data
            ! (progressing from upper left corner across rows, then down)
            ! one value per line if filetype=2 or
            ! mx values per line if filetype=3
            ! ================================================================
            case(2:3)
                ! Read header
                do i=1,5
                    read(iunit,*)
                enddo
                read(iunit,*) no_data_value

                ! Read in data
                missing = 0
                select case(abs(filetype))
                    case(2)
                        do i=1,mx*my
                            read(iunit,*) qinit(i)
                            if (qinit(i) == no_data_value) then
                                missing = missing + 1
                                qinit(i) = qinit_missing
                            endif
                        enddo
                    case(3)
                        do j=1,my
                            read(iunit,*) (qinit((j-1)*mx + i),i=1,mx)
                            do i=1,mx
                                if (qinit((j-1)*mx + i) == no_data_value) then
                                    missing = missing + 1
                                    qinit((j-1)*mx + i) = qinit_missing
                                endif
                            enddo
                        enddo
                end select

                ! Write a warning if we found and missing values
                if (missing > 0)  then
                    print *, '   WARNING...some missing data values this file'
                    print *, '       ',missing,' missing data values'
                    print *, '              (see fort.missing)'
                    print *, '   These values have arbitrarily been set to ',&
                        qinit_missing
                endif
        end select

        close(unit=iunit)

   end subroutine read_qinit

    subroutine read_force_dry(fname)

        use utility_module, only: parse_values
        character(len=*), intent(in) :: fname
        integer :: iunit,i,j,n
        double precision :: values(10), nodata_value
        character(len=80) :: str

        iunit = 8
    
        open(unit=iunit,file=fname,status='old',form='formatted')
        !read(iunit,*) tend_force_dry
        !write(6,*) 'tend_force_dry = ',tend_force_dry
        read(iunit,*) mx_fdry
        read(iunit,*) my_fdry
        read(iunit,*) xlow_fdry
        read(iunit,*) ylow_fdry

        read(iunit,'(a)') str
        call parse_values(str, n, values)
        dx_fdry = values(1)
        if (n == 2) then
            dy_fdry = values(2)
          else
            dy_fdry = dx_fdry
          endif

        read(iunit,*) nodata_value
        allocate(force_dry(mx_fdry,my_fdry))

        xhi_fdry = xlow_fdry + mx_fdry*dx_fdry
        yhi_fdry = ylow_fdry + my_fdry*dy_fdry
        write(6,*) '+++ xlow_fdry, xhi_fdry: ',xlow_fdry, xhi_fdry
        write(6,*) '+++ ylow_fdry, yhi_fdry: ',ylow_fdry, yhi_fdry

        do j=1,my_fdry
            read(iunit, *) (force_dry(i,j), i=1,mx_fdry)
            enddo
    
        close(iunit)
        return
    end subroutine read_force_dry

    
    subroutine read_eta_init(file_name)
        ! To read in file specifying different eta value in at different
        ! locations, then used in qinit function.
        ! Uses etain module variables.
        
        implicit none

        ! Input arguments
        character(len=*), intent(in), optional :: file_name
        
        ! local 
        integer, parameter :: iunit = 7
        integer :: i,j
        double precision :: nodata_value, xllower, yllower

        if (present(file_name)) then
            open(unit=iunit, file=file_name, status='unknown',&
                      form='formatted')
        else
            open(unit=iunit, file='eta_init.data', status='unknown',&
                      form='formatted')
        endif
        
        read(iunit,*) etain_mx
        !write(6,*) '+++ etain_mx = ',etain_mx
        read(iunit,*) etain_my
        !write(6,*) '+++ etain_my = ',etain_my
        read(iunit,*) xllower
        read(iunit,*) yllower
        read(iunit,*) etain_dx
        etain_dy = etain_dx
        !read(iunit,*) etain_dy
        read(iunit,*) nodata_value
        
        allocate(etain_x(etain_mx), etain_y(etain_my))
        allocate(etain_eta(etain_mx, etain_my))
        
        do i=1,etain_mx
            etain_x(i) = xllower + etain_dx*(i-1)
            enddo
            
        do j=1,etain_my
            etain_y(j) = yllower + etain_dy*(etain_my-j+1)
            read(iunit,*) (etain_eta(i,j),i=1,etain_mx)
            enddo

        
        close(unit=iunit)
    end subroutine read_eta_init


    ! ========================================================================
    ! subroutine read_qinit_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)
    ! ========================================================================
    !  Read qinit file header to determine space needed in allocatable array
    !
    !  :Input:
    !   - fname - (char) Name of file
    !   - qinit_type - (int) Type of file format (1 < qinit_type < 3)
    !
    !  :Output:
    !   - mx,my - (int) Number of grid points
    !   - xll,yll,xhi,yhi - (float) Lower and upper coordinates for grid
    !   - dx,dy - (float) Spatial resolution of grid
    ! ========================================================================
    subroutine read_qinit_header(fname,qinit_type,mx,my,xll,yll,xhi,yhi,dx,dy)

        use geoclaw_module

        implicit none

        ! Input and Output
        character*150, intent(in) :: fname
        integer, intent(in) :: qinit_type
        integer, intent(out) :: mx,my
        double precision, intent(out) :: xll,yll,xhi,yhi,dx,dy

        ! Local
        integer, parameter :: iunit = 19
        integer :: qinit_size, status
        double precision :: x,y,z,nodata_value
        logical :: found_file

        inquire(file=fname,exist=found_file)
        if (.not. found_file) then
            print *, 'Missing qinitgraphy file:'
            print *, '   ', fname
            stop
        endif

        open(unit=iunit, file=fname, status='unknown',form='formatted')

        select case(abs(qinit_type))
            ! ASCII file with 3 columns
            ! determine data size
            case(1)
                ! Initial size variables
                qinit_size = 0
                mx = 0

                ! Read in first values, determines xlow and yhi
                read(iunit,*) xll,yhi
                qinit_size = qinit_size + 1
                mx = mx + 1

                ! Go through first row figuring out mx, continue to count
                y = yhi
                do while (yhi == y)
                    read(iunit,*) x,y,z
                    qinit_size = qinit_size + 1
                    mx = mx + 1
                enddo
                mx = mx - 1
                ! Continue to count the rest of the lines
                status = 0
                do while (status == 0)
                    read(iunit,fmt=*,iostat=status) x,y,z
                    qinit_size = qinit_size + 1
                enddo
                if (status > 0) then
                    print *, "IO error occured in ",fname,", aborting!"
                    stop
                endif

                ! Calculate remaining values
                my = qinit_size / mx
                xhi = x
                yll = y
                dx = (xhi-xll) / (mx-1)
                dy = (yhi-yll) / (my-1)

            ! ASCII file with header followed by z data
            case(2:3)
                read(iunit,*) mx
                read(iunit,*) my
                read(iunit,*) xll
                read(iunit,*) yll
                read(iunit,*) dx
                read(iunit,*) nodata_value
                dy = dx
                xhi = xll + (mx-1)*dx
                yhi = yll + (my-1)*dy

            case default
                print *, 'ERROR:  Unrecognized qinit_type'
                print *, '    qinit_file_type = ',qinit_type
                print *, '  for qinit file:'
                print *, '   ', fname
                stop
        end select

        close(iunit)

        write(GEO_PARM_UNIT,*) '  mx = ',mx,'  x = (',xll,',',xhi,')'
        write(GEO_PARM_UNIT,*) '  my = ',my,'  y = (',yll,',',yhi,')'
        write(GEO_PARM_UNIT,*) '  dx, dy (meters/degrees) = ', dx,dy

    end subroutine read_qinit_header


end module qinit_module
