SUBROUTINE fc2d_geoclaw_flag2refine(blockno, meqn,maux,qvec, auxvec, dx,dy, & 
    xc,yc,t,level, maxlevel, init_flag, is_coarsening)

    USE geoclaw_module, ONLY : dry_tolerance, sea_level
    USE refinement_module, only : wave_tolerance, speed_tolerance
    USE topo_module
    USE qinit_module
    
    IMPLICIT NONE

    INTEGER :: blockno, init_flag,level, meqn, maux, maxlevel,mx,my
    DOUBLE PRECISION :: qvec(meqn),auxvec(maux), xc,yc,dx,dy,t
    logical :: is_coarsening

    INTEGER :: flag_patch, m, max_num_speeds
    DOUBLE PRECISION :: eta, th_factor, speed
    LOGICAL :: allowflag

    DOUBLE PRECISION :: x,y,x1,x2,y1,y2,xlow,xhi,ylow,yhi
    INTEGER :: i,j

    include 'regions.i'

    mx = 16
    my = 16
    do 200 j = 1,my
        y = ylower +  (j-0.5d0)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy
        do 100 i = 1,mx
          x = xlower +  (i-0.5d0)*dx
          x1 = xlower +  (i-1)*dx
          x2 = xlower +  i*dx
        
          flag_patch = 0

           do 30 m=1,mtopofiles
            if (level .lt. 1 .and. &
                t.ge.tlowtopo(m) .and. t.le.thitopo(m)) then
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 if (x2.gt.xlow.and.x1.lt.xhi.and. &
                    y2.gt.ylow.and.y1.lt.yhi) then
                    flag_patch = 1
                    go to 100 !# flagged, so no need to check anything else
                    endif
              endif
   30       continue
            do 40 m=1,mregions
                if (level .lt. minlevelregion(m) .and. &
                    t.ge.tlowregion(m) .and. t.le.thiregion(m)) then
                xlow = xlowregion(m)
                xhi = xhiregion(m)
                ylow = ylowregion(m)
                yhi = yhiregion(m)
                    if (x2.gt.xlow.and.x1.lt.xhi.and. &
                        y2.gt.ylow.and.y1.lt.yhi) then
                        flag_patch = 1
                        go to 100 !# flagged, so no need to check anything else
                        endif
                endif
   40       continue

        !  do m=1,qinit_type
        !     if (abs(t).lt.1.d0) then
        !        if (level.lt.1.d0 .and. &
        !             x2.gt.x_low_qinit(m).and.x1.lt.x_hi_qinit(m).and. &
        !             y2.gt.y_low_qinit(m).and.y1.lt.y_hi_qinit(m)) then
        !             flag_patch = 1
        !         go to 100 !# flagged, so no need to check anything else
        !         endif
        !      endif
        !  enddo


    if (allowflag(x,y,t,level)) then
       
        max_num_speeds = min(size(speed_tolerance),maxlevel)

        th_factor = 1
        if (is_coarsening) then
            !! Coarsening factor should be 0.5 refinement factor.
            th_factor = 0.5
        endif

        ! flag_patch = 0
        if (qvec(1) > dry_tolerance) then

            !! Check wave height criteria
            eta = qvec(1) + auxvec(1)
            if (abs(eta - sea_level) > th_factor*wave_tolerance) then
                if (level .lt. maxlevel) then
                    flag_patch = 1
                    ! write(*,*) 'flagging patch: ', eta
                    return
                endif

                if (dabs(auxvec(1)).lt. 100 ) then
                    flag_patch = 1
                    return
                endif
            endif

            ! Check speed criteria
            speed = sqrt(qvec(2)**2 + qvec(3)**2) / qvec(1)
            DO m = 1,max_num_speeds
                IF (speed > th_factor*speed_tolerance(m) .AND. level <= m) THEN
                ! IF (speed > th_factor*speed_tolerance(m)) THEN
                    flag_patch = 1
                    RETURN
                ENDIF
            ENDDO
        ENDIF
    else
        !  It isn't clear what this means;  do we not refine the entire patch
        ! if a single cell cannot be refined? 
        flag_patch = -1
        return
    endif
100    continue  !# end loop on i
200    continue  !# end loop on j
    return
end subroutine fc2d_geoclaw_flag2refine