! @des: topointegral integrates a surface over a rectangular region that is the 
!       intersection with a cartessian grid (of the topography data). 
!       The surface integrated ys defined by a piecewise bilinear through the 
!       cartessian grid.
!       The rectangular intersection has cordinates:
!                                         xim <= x <= xip and yjm <= y <= yjp.
!       The Cartesian grid has coordinates: 
!                                         xxlow <= x <= xxhi and yylow <= y <= yyhi.
!       with grid cell size dxx and dyy and mxx by myy cells.


function topointegral(xim,xip,yjm,yjp,xxlow,yylow,dxx,dyy,mxx &
                     ,myy,zz,intmethod)

    use geoclaw_module

    implicit none

    !  declare input variables
    double precision, intent(in) :: xxlow, yylow, dxx, dyy
    double precision, intent(in) :: zz(1:mxx,1:myy)
    integer, intent(in) :: mxx, myy, intmethod

    !  declare output variables
    double precision :: topointegral, xim, xip, yjm, yjp


    !  declare local variables
    double precision :: theintegral
    double precision :: xxhi, yyhi, dx, dy, djjstart, diistart, diiend
    double precision :: djjend, y1, y2, x1, x2, z11, z12, z21, z22
    integer :: iistart, jjstart, iiend, jjend, jj, ii, jjz1, jjz2

    ! initialize
    theintegral = 0.d0

    xxhi = xxlow + (mxx-1)*dxx
    yyhi = yylow + (myy-1)*dyy

    !  test for small rounding errors
    if ((xim - xxlow) < 0.d0 .or. (xip - xxhi) > 0.d0) then
        xim = max(xim,xxlow)
        xip = min(xip,xxhi)
    end if

    if ((yjm - yylow) < 0.d0 .or. (yjp - yyhi) > 0.d0) then
        yjm = max(yjm,yylow)
        yjp = min(yjp,yyhi)
    end if

    dx = xip - xim
    dy = yjp - yjm

    !  integrate piecewise bilinear over rectangular region
    if (intmethod == 1) then ! use bilinear method
        djjstart=(yjm-yylow)/dyy
        jjstart=idint(djjstart)+1

        diistart=(xim-xxlow)/dxx
        iistart=idint(diistart)+1

        diiend=(xip-xxlow)/dxx
        iiend=ceiling(diiend) + 1

        djjend=(yjp-yylow)/dyy
        jjend=ceiling(djjend)+1

        iistart=max(iistart,1)
        jjstart=max(jjstart,1)
        iiend=min(mxx,iiend)
        jjend=min(myy,jjend)

         do jj=jjstart,jjend-1
            y1=yylow + (jj-1.d0)*dyy
            y2=yylow + (jj)*dyy
            ! the array zz is indexed from north to south: jjz is the actual index
            ! of interest in the array zz
            jjz1= myy-jj+1
            jjz2= jjz1-1

            do ii=iistart,iiend-1
               x1=xxlow + (ii-1.d0)*dxx
               x2=xxlow + (ii)*dxx

               z11 = zz(ii,jjz1)
               z12 = zz(ii,jjz2)
               z21 = zz(ii+1,jjz1)
               z22 = zz(ii+1,jjz2)

               if (coordinate_system.eq.1) then !cartesian rectangle
                    theintegral = theintegral + bilinearintegral(
                                xim, xip, yjm, yjp, &
                                x1, x2, y1, y2, &
                                dxx, dyy, &
                                z11, z12, z21, z22)
               elseif (coordinate_system.eq.2) then !integrate on surface of sphere
                    theintegral = theintegral + bilinearintegral_s(
                                xim, xip, yjm, yjp, &
                                x1, x2, y1, y2, &
                                dxx, dyy, &
                                z11, z12, z21, z22)

                else
                  write(*,*)  'TOPOINTEGRAL: coordinate_system error'
                endif
            enddo
        enddo
    else
        write(*,*) 'TOPOINTEGRAL: only intmethod = 1,2 is supported'
    endif

    topointegral= theintegral
    return

end function topointegral