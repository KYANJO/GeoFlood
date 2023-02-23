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
    DOUBLE PRECISION :: eta, th_factor, speed, deep_depth
    LOGICAL :: allowflag

    DOUBLE PRECISION :: x,y,x1,x2,y1,y2,xlow,xhi,ylow,yhi
    INTEGER :: i,j

    include 'regions.i'

    mx = 16
    my = 16
    deep_depth = 100

    ! loop over interior points on this grid block
    DO 200 j = 1,my
        y = ylower +  (j-0.5d0)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy
        DO 100 i = 1,mx
          x = xlower +  (i-0.5d0)*dx
          x1 = xlower +  (i-1)*dx
          x2 = xlower +  i*dx
        
          flag_patch = 0 ! default is not to flag

           ! check to seeIF refinement is forced in any topo domain
           DO 30 m=1,mtopofiles
            IF (level .lt. 1 .and. &
                t.ge.tlowtopo(m) .and. t.le.thitopo(m)) THEN
              xlow = xlowtopo(m)
              xhi = xhitopo(m)
              ylow = ylowtopo(m)
              yhi = yhitopo(m)
                 IF (x2.gt.xlow.and.x1.lt.xhi.and. &
                    y2.gt.ylow.and.y1.lt.yhi) THEN
                    flag_patch = 1
                    GO TO 100 !# flagged, so no need to check anything else
                    ENDIF
              ENDIF
   30       CONTINUE

              ! check to see IF refinement is forced in any other region
            DO 40 m=1,mregions
                IF (level .lt. minlevelregion(m) .and. &
                    t.ge.tlowregion(m) .and. t.le.thiregion(m)) THEN
                xlow = xlowregion(m)
                xhi = xhiregion(m)
                ylow = ylowregion(m)
                yhi = yhiregion(m)
                    IF (x2.gt.xlow.and.x1.lt.xhi.and. &
                        y2.gt.ylow.and.y1.lt.yhi) THEN
                        flag_patch = 1
                        GO TO 100 !# flagged, so no need to check anything else
                        ENDIF
                ENDIF
   40       CONTINUE

        IF (allowflag(x,y,t,level)) THEN
        
            max_num_speeds = min(size(speed_tolerance),maxlevel)

            th_factor = 1
            IF (is_coarsening) THEN
                !! Coarsening factor should be 0.5 refinement factor.
                th_factor = 0.5
            ENDIF

            ! flag_patch = 0
            IF (qvec(1) > dry_tolerance) THEN

                !! Check wave height criteria
                eta = qvec(1) + auxvec(1)
                IF (abs(eta - sea_level) > th_factor*wave_tolerance) THEN
                    IF (level .lt. maxlevel) THEN
                        ! refine to this level in deep water
                        flag_patch = 1
                        RETURN
                    ENDIF

                    IF (dabs(auxvec(1)).lt. deep_depth ) THEN
                        ! refine to this level in shallow water (shoeregion or river banks or flood edges) 
                        flag_patch = 1
                        RETURN
                    ENDIF
                ENDIF

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
    ELSE
            !  It isn't clear what this means;  do we not refine the entire patch
            ! IF a single cell cannot be refined? 
            flag_patch = -1
            RETURN
        ENDIF
100    CONTINUE  !# end loop on i
200    CONTINUE  !# end loop on j
    RETURN
end subroutine fc2d_geoclaw_flag2refine