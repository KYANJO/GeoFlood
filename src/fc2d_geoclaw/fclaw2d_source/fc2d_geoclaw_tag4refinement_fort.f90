SUBROUTINE fc2d_geoclaw_fort_tag4refinement(mx,my,mbc,meqn,maux,xlower,ylower, &
    dx,dy,t,blockno,q,aux,level,maxlevel,init_flag,tag_patch)

    IMPLICIT NONE

    INTEGER :: mx,my, mbc, meqn, maux, tag_patch, init_flag
    INTEGER :: blockno, level, maxlevel
    DOUBLE PRECISION :: xlower, ylower, dx, dy, t, x1,x2,y1,y2
    DOUBLE PRECISION, intent(in) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION, INTENT(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    !! Local variables
    INTEGER :: i,j, mq,m
    DOUBLE PRECISION :: xc,yc,xupper, yupper, qvec(meqn), auxvec(maux)
    logical :: is_coarsening

!!    INTEGER :: tag_patch_regions, fc2d_geoclaw_refine_using_regions
    INTEGER :: fc2d_geoclaw_flag2refine, flag_patch

!!    xupper = xlower + mx*dx
!!    yupper = ylower + my*dy
!!
!!    tag_patch_regions = fc2d_geoclaw_refine_using_regions(level,xlower,ylower,xupper,yupper,t)
!!
!!    if (tag_patch_regions .ge. 0) then
!!        !! # Tagging based on regions is conclusive
!!        tag_patch = tag_patch_regions
!!        return
!!    endif

    !! We are not refining based on regions and so will use other criteria
    tag_patch = 0 

    is_coarsening = .false.   !! Don't loop over ghost cells.
    DO j = 1,my
        yc = ylower + (j-0.5)*dy
        y1 = ylower + (j-1)*dy
        y2 = ylower + j*dy

        DO i = 1,mx
            xc = xlower + (i-0.5)*dx
            x1 = xlower +  (i-1)*dx
            x2 = xlower +  i*dx

            flag_patch = fc2d_geoclaw_flag2refine( & 
                    blockno,mx,my, meqn,maux, q, aux, dx,dy,xc,yc,x1,y1,x2,y2,t,level, & 
                    maxlevel, init_flag, is_coarsening,i,j,mbc)

!!          # -1 : Not conclusive (possibly ghost cell); don't tag for refinement
!!          # 0  : Does not pass threshold (don't tag for refinement)      
!!          # 1  : Passes threshold (tag for refinement)
            if (flag_patch .gt. 0) then
                tag_patch = 1
                return
            endif
        end do 
    enddo 

END SUBROUTINE fc2d_geoclaw_fort_tag4refinement
