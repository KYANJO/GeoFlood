SUBROUTINE fc2d_geoclaw_step2(maxm,meqn,maux,mbc,mx,my, &
     qold,aux,dx,dy,dt,cflgrid, &
     fm,fp,gm,gp,rpn2,rpt2,block_corner_count)
!     ==========================================================
!
!     # clawpack routine ...  modified for AMRCLAW
!
!     # Take one time step, updating q.
!     # On entry, qold gives
!     #    initial data for this step
!     #    and is unchanged in this version.
!
!     # fm, fp are fluxes to left and right of single cell edge
!     # See the flux2 documentation for more information.
!
!     # modified again for GeoClaw
!------------------step2_geo.f-----------------------
!     This version of step2 relimits the fluxes in order to
!     maintain positivity.
!     to do so set relimit=.true.
!     The only modification is the 101 loop
!------------last modified 12/30/04--------------------------
!

    use geoclaw_module, only: dry_tolerance
    use amr_module, only: mwaves, mcapa

    implicit none

    external rpn2, rpt2

    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    double precision, intent(in) :: dx,dy,dt
    double precision, intent(inout) :: cflgrid
    double precision, intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    double precision, intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Local storage for flux accumulation
    double precision :: faddm(meqn,1-mbc:maxm+mbc)
    double precision :: faddp(meqn,1-mbc:maxm+mbc)
    double precision :: gaddm(meqn,1-mbc:maxm+mbc,2)
    double precision :: gaddp(meqn,1-mbc:maxm+mbc,2)

    ! Scratch storage for Sweeps and Riemann problems
    double precision ::  q1d(meqn,1-mbc:maxm+mbc)
    double precision :: aux1(maux,1-mbc:maxm+mbc)
    double precision :: aux2(maux,1-mbc:maxm+mbc)
    double precision :: aux3(maux,1-mbc:maxm+mbc)
    double precision :: dtdx1d(1-mbc:maxm+mbc)
    double precision :: dtdy1d(1-mbc:maxm+mbc)

 !   double precision ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
 !   double precision ::     s(mwaves, 1-mbc:maxm + mbc)
 !   double precision ::  amdq(meqn,1-mbc:maxm + mbc)
 !   double precision ::  apdq(meqn,1-mbc:maxm + mbc)
!     double precision ::  cqxx(meqn,1-mbc:maxm + mbc)
!     double precision :: bmadq(meqn,1-mbc:maxm + mbc)
!     double precision :: bpadq(meqn,1-mbc:maxm + mbc)

    INTEGER block_corner_count(0:3), ibc,jbc

    ! Looping scalar storage
    integer :: i,j,m,thread_num
    double precision :: dtdx,dtdy,cfl1d,p,phi,cm,dtdxij,dtdyij

    ! Common block storage
    double precision dtcom, dxcom,dycom,tcom
    integer :: icom,jcom
    common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

    ! Parameters
    ! Relimit fluxes to maintain positivity
    logical, parameter :: relimit = .false.

    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

    fm = 0.d0
    fp = 0.d0
    gm = 0.d0
    gp = 0.d0

    ! ==========================================================================
    ! Perform X-Sweeps
    ! do j = 0,my+1
    do  j = 2-mbc,my+mbc-1

        ! Copy old q into 1d slice
       q1d(:,1-mbc:mx+mbc) = qold(:,1-mbc:mx+mbc,j)

       !! FORESTCLAW change to handle cubed sphere
      !  IF (j .EQ. 0) THEN
      !     DO m = 1,meqn
      !        IF (block_corner_count(0) .EQ. 3) THEN
      !           DO ibc = 1,mbc
      !              q1d(m,1-ibc) = qold(m,ibc,0)
      !           ENDDO
      !        ENDIF
      !        IF (block_corner_count(1) .EQ. 3) THEN
      !           DO ibc = 1,mbc
      !              q1d(m,mx+ibc) = qold(m,mx-ibc+1,0)
      !           ENDDO
      !        ENDIF
      !     ENDDO
      !  ELSE IF (j .EQ. my+1) THEN
      !     DO m = 1,meqn
      !        IF (block_corner_count(2) .EQ. 3) THEN
      !           DO ibc = 1,mbc
      !              q1d(m,1-ibc) = qold(m,ibc,my+1)
      !           ENDDO
      !        ENDIF
      !        IF (block_corner_count(3) .EQ. 3) THEN
      !           DO ibc = 1,mbc
      !              q1d(m,mx+ibc) = qold(m,mx-ibc+1,my+1)
      !           ENDDO
      !        ENDIF
      !     ENDDO
      !  ENDIF


        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif

        ! Copy aux array into slices
        if (maux > 0) then
            aux1(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j-1)
            aux2(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j  )
            aux3(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j+1)
        endif

        ! Store value of j along the slice into common block
        ! *** WARNING *** This may not working with threading
        jcom = j

        ! Compute modifications fadd and gadd to fluxes along this slice:
        call fc2d_geoclaw_flux2(1,maxm,meqn,maux,mbc,mx, &
                                q1d,dtdx1d,aux1,aux2,aux3, &
                                faddm,faddp,gaddm,gaddp,cfl1d,rpn2,rpt2)

        cflgrid = dmax1(cflgrid,cfl1d)
        ! write(53,*) 'x-sweep: ',cfl1d,cflgrid

        ! Update fluxes
        ! fm(:,1:mx+1,j) = fm(:,1:mx+1,j) + faddm(:,1:mx+1)
        ! fp(:,1:mx+1,j) = fp(:,1:mx+1,j) + faddp(:,1:mx+1)
        ! gm(:,1:mx+1,j) = gm(:,1:mx+1,j) + gaddm(:,1:mx+1,1)
        ! gp(:,1:mx+1,j) = gp(:,1:mx+1,j) + gaddp(:,1:mx+1,1)
        ! gm(:,1:mx+1,j+1) = gm(:,1:mx+1,j+1) + gaddm(:,1:mx+1,2)
        ! gp(:,1:mx+1,j+1) = gp(:,1:mx+1,j+1) + gaddp(:,1:mx+1,2)
        
        fm(:,2-mbc:mx+mbc,j) = fm(:,2-mbc:mx+mbc,j) + faddm(:,2-mbc:mx+mbc)
        fp(:,2-mbc:mx+mbc,j) = fp(:,2-mbc:mx+mbc,j) + faddp(:,2-mbc:mx+mbc)
        gm(:,2-mbc:mx+mbc,j) = gm(:,2-mbc:mx+mbc,j) + gaddm(:,2-mbc:mx+mbc,1)
        gp(:,2-mbc:mx+mbc,j) = gp(:,2-mbc:mx+mbc,j) + gaddp(:,2-mbc:mx+mbc,1)
        gm(:,2-mbc:mx+mbc,j+1) = gm(:,2-mbc:mx+mbc,j+1) + gaddm(:,2-mbc:mx+mbc,2)
        gp(:,2-mbc:mx+mbc,j+1) = gp(:,2-mbc:mx+mbc,j+1) + gaddp(:,2-mbc:mx+mbc,2)


    enddo

    ! ============================================================================
    !  y-sweeps
    !
    ! do i = 0, mx+1
    do i = 2-mbc,mx+mbc-1

       ! Copy data along a slice into 1d arrays:
       q1d(:,1-mbc:my+mbc) = qold(:,i,1-mbc:my+mbc)

       !! FORESTCLAW change to handle cubed sphere
      !  IF (i .EQ. 0) THEN
      !     DO m = 1,meqn
      !        IF (block_corner_count(0) .EQ. 3) THEN
      !           DO jbc = 1,mbc
      !              q1d(m,1-jbc) = qold(m,0,jbc)
      !           ENDDO
      !        ENDIF
      !        IF (block_corner_count(2) .EQ. 3) THEN
      !           DO jbc = 1,mbc
      !              q1d(m,my+jbc) = qold(m,0,my-jbc+1)
      !           ENDDO
      !        ENDIF
      !     ENDDO
      !  ELSE IF (i .EQ. mx+1) THEN
      !     DO m = 1,meqn
      !        IF (block_corner_count(1) .EQ. 3) THEN
      !           DO jbc = 1,mbc
      !              q1d(m,1-jbc) = qold(m,mx+1,jbc)
      !           ENDDO
      !        ENDIF
      !        IF (block_corner_count(3) .EQ. 3) THEN
      !           DO jbc= 1,mbc
      !              q1d(m,my+jbc) = qold(m,mx+1,my-jbc+1)
      !           ENDDO
      !        ENDIF
      !     ENDDO
      !  ENDIF


        ! Set dt/dy ratio in slice
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! Copy aux slices
        if (maux .gt. 0)  then
            aux1(:,1-mbc:my+mbc) = aux(:,i-1,1-mbc:my+mbc)
            aux2(:,1-mbc:my+mbc) = aux(:,i,1-mbc:my+mbc)
            aux3(:,1-mbc:my+mbc) = aux(:,i+1,1-mbc:my+mbc)
        endif

        ! Store the value of i along this slice in the common block
        ! *** WARNING *** This may not working with threading
        icom = i

        ! Compute modifications fadd and gadd to fluxes along this slice
        call fc2d_geoclaw_flux2(2,maxm,meqn,maux,mbc,my,q1d, & 
                                dtdy1d,aux1,aux2,aux3, &
                                faddm,faddp,gaddm,gaddp,cfl1d,rpn2,rpt2)

        cflgrid = dmax1(cflgrid,cfl1d)
        ! write(53,*) 'y-sweep: ',cfl1d,cflgrid

        ! Update fluxes
        ! gm(:,i,1:my+1) = gm(:,i,1:my+1) + faddm(:,1:my+1)
        ! gp(:,i,1:my+1) = gp(:,i,1:my+1) + faddp(:,1:my+1)
        ! fm(:,i,1:my+1) = fm(:,i,1:my+1) + gaddm(:,1:my+1,1)
        ! fp(:,i,1:my+1) = fp(:,i,1:my+1) + gaddp(:,1:my+1,1)
        ! fm(:,i+1,1:my+1) = fm(:,i+1,1:my+1) + gaddm(:,1:my+1,2)
        ! fp(:,i+1,1:my+1) = fp(:,i+1,1:my+1) + gaddp(:,1:my+1,2)
        gm(:,i,2-mbc:my+mbc) = gm (:,i,2-mbc:my+mbc) + faddm(:,2-mbc:my+mbc)
        gp(:,i,2-mbc:my+mbc) = gp(:,i,2-mbc:my+mbc) + faddp(:,2-mbc:my+mbc)
        fm(:,i,2-mbc:my+mbc) = fm(:,i,2-mbc:my+mbc) + gaddm(:,2-mbc:my+mbc,1)
        fp(:,i,2-mbc:my+mbc) = fp(:,i,2-mbc:my+mbc) + gaddp(:,2-mbc:my+mbc,1)
        fm(:,i+1,2-mbc:my+mbc) = fm(:,i+1,2-mbc:my+mbc) + gaddm(:,2-mbc:my+mbc,2)
        fp(:,i+1,2-mbc:my+mbc) = fp(:,i+1,2-mbc:my+mbc) + gaddp(:,2-mbc:my+mbc,2)

    end do

    ! Relimit correction fluxes if they drive a cell negative
    if (relimit) then
        dtdxij = dtdx
        dtdyij = dtdy
        do i=1,mx
            do j=1,my
                if (mcapa > 0) then
                    dtdxij = dtdx / aux(mcapa,i,j)
                    dtdyij = dtdy / aux(mcapa,i,j)
                endif
                p = max(0.d0,dtdxij*fm(1,i+1,j)) + max(0.d0,dtdyij*gm(1,i,j+1)) &
                  - min(0.d0,dtdxij*fp(1,i,j)) - min(0.d0,dtdyij*gp(1,i,j))
                phi = min(1.d0,abs(qold(1,i,j) / (p+dry_tolerance)))

                if (phi < 1.d0) then
                    do m=1,meqn
                        if (fp(1,i,j) < 0.d0) then
                            cm = fp(m,i,j) - fm(m,i,j)
                            fm(m,i,j) = phi * fm(m,i,j)
                            fp(m,i,j) = fm(m,i,j) + cm
                        endif
                        if (gp(1,i,j) < 0.d0) then
                            cm = gp(m,i,j) - gm(m,i,j)
                            gm(m,i,j) = phi * gm(m,i,j)
                            gp(m,i,j) = gm(m,i,j) + cm
                        endif
                        if (fm(1,i+1,j) > 0.d0) then
                            cm = fp(m,i+1,j) - fm(m,i+1,j)
                            fp(m,i+1,j) = phi * fp(m,i+1,j)
                            fm(m,i+1,j) = fp(m,i+1,j) - cm
                        endif
                        if (gm(1,i,j+1) > 0.d0) then
                            cm = gp(m,i,j+1) - gm(m,i,j+1)
                            gp(m,i,j+1) = phi * gp(m,i,j+1)
                            gm(m,i,j+1) = gp(m,i,j+1) - cm
                        endif
                    end do
                endif
            enddo
        enddo
    endif

  END SUBROUTINE fc2d_geoclaw_step2
