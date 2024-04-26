! @author : Brian KYANJO
! @date : 2023-08-11
! @brief : This is a modified version of the original flag2refine.f90 file from geoclaw
!          The modifications are to allow for refinement based on the flow grade criteria and specified regions

integer function fc2d_geoclaw_flag2refine(blockno, mx1,my1, meqn,maux,qvec, auxvec, dx,dy, & 
    xc,yc,t,level, maxlevel, init_flag, is_coarsening,i,j)

    USE geoclaw_module, ONLY : dry_tolerance, sea_level
    USE refinement_module, only : wave_tolerance, speed_tolerance
    USE refinement_module, only : num_flowgrades,iflowgradevariable,iflowgradetype, &
                                iflowgrademinlevel,flowgradevalue
    USE topo_module
    USE qinit_module
    
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: mx1, my1
    INTEGER :: blockno, init_flag,level, meqn, maux, maxlevel
    DOUBLE PRECISION :: qvec(meqn),auxvec(maux),xc,yc,dx,dy,t
    logical :: is_coarsening

    INTEGER :: flag_patch, m, max_num_speeds,iflow
    DOUBLE PRECISION :: eta, th_factor, speed, deep_depth,momentum,depth
    DOUBLE PRECISION :: flowgradenorm, flowgradegrad, flowgrademeasure
    LOGICAL :: allowflag

    DOUBLE PRECISION :: x,y,x1,x2,y1,y2,xlow,xhi,ylow,yhi
    INTEGER :: i,j,k, mx, my

    include 'regions.i'

    mx = int(mx1)
    my = int(my1)

    deep_depth = 100 ! meters (set to default value sinece its deprecated)

    ! loop over interior points on this grid block
    y = ylower +  (j-0.5d0)*dy
    y1 = ylower + (j-1)*dy
    y2 = ylower + j*dy
    
    x = xlower +  (i-0.5d0)*dx
    x1 = xlower +  (i-1)*dx
    x2 = xlower +  i*dx

    fc2d_geoclaw_flag2refine = 0 ! default is not to flag

    ! check to see IF refinement is forced in any other region
    DO m=1,mregions
        IF (level .lt. minlevelregion(m) .and. &
            t.ge.tlowregion(m) .and. t.le.thiregion(m)) THEN
            xlow = xlowregion(m)
            xhi = xhiregion(m)
            ylow = ylowregion(m)
            yhi = yhiregion(m)
            IF (x2.gt.xlow.and.x1.lt.xhi.and. &
                y2.gt.ylow.and.y1.lt.yhi) THEN
                fc2d_geoclaw_flag2refine = 1
                return !# flagged, so no need to check anything else
            ENDIF
        ENDIF
    ENDDO

    IF (allowflag(x,y,t,level)) THEN
        
        max_num_speeds = min(size(speed_tolerance),maxlevel)
        momentum = sqrt(qvec(2)**2 + qvec(3)**2)
        depth = qvec(1) 
        speed = momentum / depth
        eta = depth + auxvec(1)

        th_factor = 1
        IF (is_coarsening) THEN
            !! Coarsening factor should be 0.5 refinement factor.
            th_factor = 0.5
        ENDIF

        ! Check flow grade criteria are used (best for overland flows)
        DO iflow=1,num_flowgrades 
            if (iflowgradevariable(iflow) == 1) then
                flowgradenorm = depth
                flowgradegrad = depth
            elseif (iflowgradevariable(iflow) == 2) then
                flowgradenorm = momentum
                flowgradegrad = momentum
            elseif (iflowgradevariable(iflow) == 3) then
                if (depth > dry_tolerance) then
                    flowgradenorm = dabs(eta)
                    flowgradegrad = dabs(eta)
                else
                    flowgradenorm = 0.0d0
                    flowgradegrad = 0.0d0
                endif
            endif

            if (iflowgradetype(iflow) == 1) then
                flowgrademeasure = flowgradenorm
            else
                write(*,*) 'only flowgradetype = 1 is supported'
                stop
                flowgrademeasure = flowgradegrad
            endif

            if (flowgrademeasure .gt. flowgradevalue(iflow) .and. &
                level .lt. iflowgrademinlevel(iflow)) then
                fc2d_geoclaw_flag2refine = 1
                return
            else
                IF (depth > dry_tolerance) THEN
                    
                    ! Check speed criteria
                    DO m = 1,max_num_speeds
                        IF (speed > th_factor*speed_tolerance(m) .AND. level <= m) THEN
                            fc2d_geoclaw_flag2refine = 1
                            return
                        ENDIF
                    ENDDO

                ENDIF
            endif
        ENDDO

        IF (num_flowgrades .eq. 0) THEN
            IF (depth > dry_tolerance) THEN
                IF (abs(eta - sea_level) > th_factor*wave_tolerance) THEN
                    IF (level .lt. maxlevel) THEN
                        ! refine to this level in deep water
                        fc2d_geoclaw_flag2refine = 1
                        ! write(*,*) 'eta = ',eta
                        return
                    ENDIF

                    IF (dabs(auxvec(1)).lt. deep_depth ) THEN
                        ! refine to this level in shallow water (shoreregion or river banks or flood edges) 
                        fc2d_geoclaw_flag2refine = 1
                        ! write(*,*) 'mx,my = ',mx1,my1
                        return
                    ENDIF

                    fc2d_geoclaw_flag2refine = 0
                    return
                ENDIF

                ! don't refine in deep water if already at maxlevel
                if (abs(auxvec(1))> deep_depth .and. speed < speed_tolerance(1)) then
                    fc2d_geoclaw_flag2refine = 0
                    return
                endif

                ! refine at maximum velocity-depth product is 0.5 m/s
                if (speed/qvec(1) > 0.5) then
                    fc2d_geoclaw_flag2refine = 1
                    return
                endif

                ! Check speed criteria
                DO m = 1,max_num_speeds
                    IF (speed > th_factor*speed_tolerance(m) .AND. level <= m) THEN
                        fc2d_geoclaw_flag2refine = 1
                        return
                    ENDIF
                ENDDO
            ENDIF
        ENDIF
    ELSE
        !  It isn't clear what this means;  do we not refine the entire patch
        ! IF a single cell cannot be refined? 
        fc2d_geoclaw_flag2refine = 0
        return
    ENDIF

end function fc2d_geoclaw_flag2refine

