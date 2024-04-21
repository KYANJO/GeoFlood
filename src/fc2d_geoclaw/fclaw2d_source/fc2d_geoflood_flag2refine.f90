! @author : Brian KYANJO
! @date : 2023-08-11
! @brief : This is a modified version of the original flag2refine.f90 file from geoclaw
!          The modifications are to allow for refinement based on the flow grade criteria 
!          and specified regions

integer function fc2d_geoclaw_flag2refine(blockno, mx, my, meqn, maux, q, aux, dx, dy, &
                                          xc, yc, x1, y1, x2, y2, t, level, maxlevel, &
                                          init_flag, is_coarsening, i, j, mbc)

    use geoclaw_module, only : dry_tolerance, sea_level
    use refinement_module, only : wave_tolerance, speed_tolerance, num_flowgrades, &
                                  iflowgradevariable, iflowgradetype, iflowgrademinlevel, flowgradevalue

    implicit none

    integer, intent(in) :: blockno, mx, my, meqn, maux, level, maxlevel, init_flag, mbc, i, j
    double precision, intent(in) :: dx, dy, xc, yc, x1, y1, x2, y2, t
    double precision, intent(in) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc), aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    logical, intent(in) :: is_coarsening

    integer :: m, n, th_factor, iflow
    double precision :: momentum, depth, speed, eta, max_num_speeds
    ! double precision :: flowgradenorm, flowgradegrad, flowgrademeasure
    logical :: refinement_needed, allowflag

    include 'regions.i'

    fc2d_geoclaw_flag2refine = 0  ! Default is not to flag for refinement

    ! Early exit if refinement is forced in any specified region
    do m = 1, mregions
        if (level < minlevelregion(m) .and. t >= tlowregion(m) .and. t <= thiregion(m)) then
            if (x2 > xlowregion(m) .and. x1 < xhiregion(m) .and. &
                y2 > ylowregion(m) .and. y1 < yhiregion(m)) then
                fc2d_geoclaw_flag2refine = 1
                return
            endif
        endif
    enddo

    ! Evaluate conditions for refinement based on flow grade criteria
    if (allowflag(xc, yc, t, level)) then

        max_num_speeds = min(size(speed_tolerance),maxlevel)
        momentum = sqrt(q(2,i,j)**2 + q(3,i,j)**2)
        depth = q(1,i,j)
        eta = depth + aux(1,i,j)

        th_factor = 1.0
        if (is_coarsening) th_factor = 0.5  ! Coarsening factor is half the refinement factor

        ! Iterate over flow grades to determine if refinement is necessary
        do m = 1, num_flowgrades
            refinement_needed = .false.
            select case (iflowgradevariable(m))
                case (1)  ! depth
                    refinement_needed = (depth > flowgradevalue(m))
                case (2)  ! momentum
                    refinement_needed = (momentum > flowgradevalue(m))
                case (3)  ! elevation
                    if (depth > dry_tolerance) then
                        refinement_needed = (abs(eta) > flowgradevalue(m))
                    endif
            end select

            ! Refine if criteria met and within allowable level
            if (refinement_needed .and. level < iflowgrademinlevel(m)) then
                fc2d_geoclaw_flag2refine = 1
                return
            endif
        enddo

        ! Perform a separate check for speed criteria if no specific flow grade triggered refinement
        if (depth > dry_tolerance) then
            speed = momentum / depth  
            do n = 1, max_num_speeds
                if (speed > th_factor * speed_tolerance(n) .and. level <= n) then
                    fc2d_geoclaw_flag2refine = 1
                    return
                endif
            enddo
        endif    

        ! Evaluate additional wave criteria if no specific flow grades or speed creteria are defined
        if (num_flowgrades == 0) then
            if (depth > dry_tolerance) then

                ! Check wave criteria
                if (abs(eta - sea_level) > th_factor*wave_tolerance) then
                    if (level < maxlevel) then
                        fc2d_geoclaw_flag2refine = 1
                        return
                    endif
                endif

            endif
        endif
    endif

end function fc2d_geoclaw_flag2refine
