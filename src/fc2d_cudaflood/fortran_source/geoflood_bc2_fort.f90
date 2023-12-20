! author: Brian Kyanjo
! date: 2023-11-11
! purpose: This file contains modified boundary conditions for the 2D shallow water equations with hydrograph support capabilities
!          The boundary conditions are based on the clawpack 4.0 code

! ==================================================================
subroutine fc2d_geoclaw_bc2(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc)
!  ==================================================================
      ! standard boundary condition choices

      ! At each boundary k = 1 (left), 2 (right), 3 (top) or 4 (bottom):
      !  if mthbc(k) = 0, user-supplied BC's (must be inserted!)
      !              = 1, zero-order extrapolation
      !              = 2, periodic BC's
                  ! = 3,  solid walls, assuming this can be implemented by reflecting  the data about the boundary and then negating the 2'nd (for k=1,2) 0r 3'rd (for k=3,4) component of q.
      ! -------------------------------------------------------------------

      ! Extend the data from the interior cells (1:mx, 1:my) to the ghost cells outside the region:
      ! (i,1-jbc)  for jbc = 1,mbc, i = 1-mbc, mx+mbc
      ! (i,my+jbc) for jbc = 1,mbc, i = 1-mbc, mx+mbc
      ! (1-ibc,j)  for ibc = 1,mbc, j = 1-mbc, my+mbc
      ! (mx+ibc,j) for ibc = 1,mbc, j = 1-mbc, my+mbc

      use hydrograph_module, only: inflow_interpolate, Riemann_invariants, two_shock ! contains the interpolation, Riemann invariant, and two-shock ssunroutines
      use hydrograph_module, only: b, x0, y0, boundary_location,use_hydrograph
      use geoclaw_module, only: dry_tolerance, grav

      implicit none

      integer, intent(in) :: meqn, mbc, mx, my, maux, mthbc(4)
      double precision, intent(in) :: xlower, ylower, dx, dy, t, dt

      double precision, intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision, intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer :: m, i, j, ibc, jbc

      double precision ::  x,y
      double precision, dimension(4) :: q0,q1

      ! -------------------------------------------------------------------
      !  left boundary
      ! -------------------------------------------------------------------
      go to (100,110,120,130), mthbc(1)+1
      ! this is how we skip over this side... if (mthbc(1))+1 is not 1,2,3 or 4, then goto above walls through here ...
      goto 199

      100 continue
      ! user-supplied BC's (must be inserted!)
      !  in this case, we are using the inflow_interpolation subroutine to compute the inflow boundary condition values

      ! c     # inflow boundary condition:
      if ((use_hydrograph == 'True') .and. (boundary_location(1)) == 'left') then
            ! apply the hydrograph boundary condition
            do j = 1-mbc,my+mbc
                  y = ylower + (j-0.5d0)*dy
                  q1 = q(:,1,j) ! this is the vector q at the boundary
                  call inflow_interpolate(t,q0,q1)

                  if (abs(y-y0) <= b/2.0d0) then
                        do ibc=1,mbc
                              if (q1(1) < dry_tolerance) then
                                    q(1,1-ibc,j) = max((q0(2)/sqrt(grav))**(2.0d0/3.0d0), 0.001d0) 
                                    q(2,1-ibc,j) = q0(2)
                                    q(3,1-ibc,j) = 0.0d0
                              else 
                                    if (q0(2) .ne. 0.0d0) then
                                          call Riemann_invariants(q0,q1)
                                          if (q0(1) > q1(1)) then
                                                call two_shock(q0,q1)
                                          endif
                                          q(1,1-ibc,j) = q0(1)
                                          q(2,1-ibc,j) = q0(2)
                                          q(3,1-ibc,j) = 0.0d0
                                    else
                                          aux(1,1-ibc,j) = aux(1,ibc,j)
                                          do m=1,meqn
                                                q(m,1-ibc,j) = q(m,ibc,j)
                                          enddo

                                          ! c     # negate the normal velocity:   
                                          q(2,1-ibc,j) = -q(2,ibc,j)
                                    end if
                              endif
                        enddo
                  else
                        do ibc=1,mbc                            
                              aux(1,1-ibc,j) = aux(1,ibc,j)
                              do m=1,meqn
                                    q(m,1-ibc,j) = q(m,ibc,j)
                              enddo
                              ! c     # negate the normal velocity:   
                              q(2,1-ibc,j) = -q(2,ibc,j)
                        enddo
                  endif
            end do
      else
            ! # user-specified boundary conditions go here in place of error output
            write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
            stop
      endif
      goto 199

      110 continue
      ! zero-order extrapolation
      do 115 j = 1-mbc,my+mbc
            do 115 ibc=1,mbc
                  aux(1,1-ibc,j) = aux(1,1,j)
                  do 115 m=1,meqn
                  q(m,1-ibc,j) = q(m,1,j)
      115         continue
      go to 199

      120 continue
      ! periodic BC's: handled by p4est
      goto 199

      130 continue
      ! solid wall (assumes 2'nd component is velocity or momentum in x)
            do 135 j = 1-mbc, my+mbc
            do 135 ibc=1,mbc
                  aux(1,1-ibc,j) = aux(1,ibc,j)
                  do 135 m=1,meqn
                  q(m,1-ibc,j) = q(m,ibc,j)
      135       continue
      ! c     # negate the normal velocity:
            do 136 j = 1-mbc, my+mbc
                  do 136 ibc=1,mbc
                  q(2,1-ibc,j) = -q(2,ibc,j)
      136    continue
            go to 199

      199 continue
      ! c
      ! c-------------------------------------------------------
      ! c     # right boundary:
      ! c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
      goto 299
      ! c
      200 continue
      !  inflow boundary condition:
      if ((use_hydrograph == 'True') .and. (boundary_location(2) == 'right')) then
            ! apply the hydrograph boundary condition
            do j = 1-mbc,my+mbc
                  y = ylower + (j-0.5d0)*dy
                  q1 = q(:,mx,j) ! this is the vector q at the boundary
                  call inflow_interpolate(t,q0,q1)

                  if (abs(y-y0) <= b/2.0d0) then
                        do ibc=1,mbc
                              if (q1(1) < dry_tolerance) then
                                    q(1,mx+ibc,j) = max((q0(2)/sqrt(grav))**(2.0d0/3.0d0), 0.001d0) 
                                    q(2,mx+ibc,j) = q0(2)
                                    q(3,mx+ibc,j) = 0.0d0
                              else 
                                    if (q0(2) .ne. 0.0d0) then
                                          call Riemann_invariants(q0,q1)
                                          if (q0(1) > q1(1)) then
                                                call two_shock(q0,q1)
                                          endif
                                          q(1,mx+ibc,j) = q0(1)
                                          q(2,mx+ibc,j) = q0(2)
                                          q(3,mx+ibc,j) = 0.0d0
                                    else
                                          aux(1,mx+ibc,j) = aux(1,mx+1-ibc,j)
                                          do m=1,meqn
                                                q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
                                          enddo

                                          ! c     # negate the normal velocity:   
                                          q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
                                    end if
                              endif
                        enddo
                  else
                        do ibc=1,mbc                            
                              aux(1,mx+ibc,j) = aux(1,mx+1-ibc,j)
                              do m=1,meqn
                                    q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
                              enddo
                              ! c     # negate the normal velocity:   
                              q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
                        enddo
                  endif
            end do
      else
            ! c     # user-specified boundary conditions go here in place of error output
            write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
            stop
      endif
      go to 299

      210 continue
      ! c     # zero-order extrapolation:
      do 215 j = 1-mbc, my+mbc
            do 215 ibc=1,mbc
            aux(1,mx+ibc,j) = aux(1,mx,j)
            do 215 m=1,meqn
                  q(m,mx+ibc,j) = q(m,mx,j)
      215       continue
      go to 299

      220 continue
      ! c     # periodic : Handled elsewhere
      go to 299

      230 continue
      ! c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 j = 1-mbc, my+mbc
            do 235 ibc=1,mbc
            aux(1,mx+ibc,j) = aux(1,mx+1-ibc,j)
            do 235 m=1,meqn
                  q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
      235       continue
      ! c     # negate the normal velocity:
      do 236 j = 1-mbc, my+mbc
            do 236 ibc=1,mbc
            q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
      236    continue
      go to 299

      299 continue
      ! c
      ! c-------------------------------------------------------
      ! c     # bottom boundary:
      ! c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
      goto 399
      ! c
      300 continue
      ! inflow boundary condition:
      if ((use_hydrograph == 'True') .and. (boundary_location(4) == 'bottom')) then
            ! apply the hydrograph boundary condition
            do i = 1-mbc,mx+mbc
                  x = xlower + (i-0.5d0)*dx
                  q1 = q(:,i,1) ! this is the vector q at the boundary
                  call inflow_interpolate(t,q0,q1)

                  if (abs(x-x0) <= b/2.0d0) then
                        do jbc=1,mbc
                              if (q1(1) < dry_tolerance) then
                                    q(1,i,1-jbc) = max((q0(2)/sqrt(grav))**(2.0d0/3.0d0), 0.001d0) 
                                    q(2,i,1-jbc) = q0(2)
                                    q(3,i,1-jbc) = 0.0d0
                              else 
                                    if (q0(2) .ne. 0.0d0) then
                                          call Riemann_invariants(q0,q1)
                                          if (q0(1) > q1(1)) then
                                                call two_shock(q0,q1)
                                          endif
                                          q(1,i,1-jbc) = q0(1)
                                          q(2,i,1-jbc) = q0(2)
                                          q(3,i,1-jbc) = 0.0d0
                                    else
                                          aux(1,i,1-jbc) = aux(1,i,jbc)
                                          do m=1,meqn
                                                q(m,i,1-jbc) = q(m,i,jbc)
                                          enddo

                                          ! c     # negate the normal velocity:   
                                          q(3,i,1-jbc) = -q(3,i,jbc)
                                    end if
                              endif
                        enddo
                  else
                        do jbc=1,mbc   
                              aux(1,i,1-jbc) = aux(1,i,jbc)
                              do m=1,meqn
                                    q(m,i,1-jbc) = q(m,i,jbc)
                              enddo
                              ! c     # negate the normal velocity:   
                              q(3,i,1-jbc) = -q(3,i,jbc)
                        enddo
                  endif
            end do                         
                             
      else
            ! c     # user-specified boundary conditions go here in place of error output
            write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
            stop
      endif
      go to 399
      ! c
      310 continue
      ! c     # zero-order extrapolation:
      do 315 jbc=1,mbc
            do 315 i = 1-mbc, mx+mbc
            aux(1,i,1-jbc) = aux(1,i,1)
            do 315 m=1,meqn
                  q(m,i,1-jbc) = q(m,i,1)
      315       continue
      go to 399

      320 continue
      ! c     # periodic: Handled elsewhere
      go to 399

      330 continue
      ! c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 jbc=1,mbc
            do 335 i = 1-mbc, mx+mbc
            aux(1,i,1-jbc) = aux(1,i,jbc)
            do 335 m=1,meqn
                  q(m,i,1-jbc) = q(m,i,jbc)
      335       continue
      ! c     # negate the normal velocity:
      do 336 jbc=1,mbc
            do 336 i = 1-mbc, mx+mbc
            q(3,i,1-jbc) = -q(3,i,jbc)
      336    continue
      go to 399

      399 continue
      ! c
      ! c-------------------------------------------------------
      ! c     # top boundary:
      ! c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
      goto 499

      400 continue
      !  inflow boundary condition:
      if ((use_hydrograph == 'True') .and. (boundary_location(3) == 'top')) then
            ! apply the hydrograph boundary condition
            do i = 1-mbc,mx+mbc
                  x = xlower + (i-0.5d0)*dx
                  q1 = q(:,i,my) ! this is the vector q at the boundary
                  call inflow_interpolate(t,q0,q1)

                  if (abs(x-x0) <= b/2.0d0) then
                        do jbc=1,mbc
                              if (q1(1) < dry_tolerance) then
                                    q(1,i,my+jbc) = max((q0(2)/sqrt(grav))**(2.0d0/3.0d0), 0.001d0) 
                                    q(2,i,my+jbc) = q0(2)
                                    q(3,i,my+jbc) = 0.0d0
                              else 
                                    if (q0(2) .ne. 0.0d0) then
                                          call Riemann_invariants(q0,q1)
                                          if (q0(1) > q1(1)) then
                                                call two_shock(q0,q1)
                                          endif
                                          q(1,i,my+jbc) = q0(1)
                                          q(2,i,my+jbc) = q0(2)
                                          q(3,i,my+jbc) = 0.0d0
                                    else
                                          aux(1,i,my+jbc) = aux(1,i,my+1-jbc)
                                          do m=1,meqn
                                                q(m,i,my+jbc) = q(m,i,my+1-jbc)
                                          enddo

                                          ! c     # negate the normal velocity:   
                                          q(3,i,my+jbc) = -q(3,i,my+1-jbc)
                                    end if
                              endif
                        enddo
                  else
                        do jbc=1,mbc   
                              aux(1,i,my+jbc) = aux(1,i,my+1-jbc)
                              do m=1,meqn
                                    q(m,i,my+jbc) = q(m,i,my+1-jbc)
                              enddo
                              ! c     # negate the normal velocity:   
                              q(3,i,my+jbc) = -q(3,i,my+1-jbc)
                        enddo
                  endif
            end do
      else
            !  # user-specified boundary conditions go here in place of error output
            write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
            stop
      endif
      go to 499

      410 continue
      !     # zero-order extrapolation:
      do 415 jbc=1,mbc
            do 415 i = 1-mbc, mx+mbc
            aux(1,i,my+jbc) = aux(1,i,my)
            do 415 m=1,meqn
                  q(m,i,my+jbc) = q(m,i,my)
      415       continue
      go to 499

      420 continue
      !     # periodic: Handled elsewhere
      go to 499

      430 continue
      !  solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
            aux(1,i,my+jbc) = aux(1,i,my+1-jbc)
            do 435 m=1,meqn
                  q(m,i,my+jbc) = q(m,i,my+1-jbc)
      435       continue
      !  # negate the normal velocity:
      do 436 jbc=1,mbc
            do 436 i = 1-mbc, mx+mbc
            q(3,i,my+jbc) = -q(3,i,my+1-jbc)
      436    continue
      go to 499

      499 continue

            return


end subroutine fc2d_geoclaw_bc2
