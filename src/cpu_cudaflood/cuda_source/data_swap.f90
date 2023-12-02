! @author: Brian Kyanjo
! @date: 2023/13/08
! @purpose: This program reads in a (m,i,j) or (i,j,m) array and swaps the indices
!          to the other form. This is useful for exchanging data between geoclaw (m,i,j) 
!          routines and cudaclaw (i,j,m)routines

subroutine cudaclaw_swap_data(mx,my,mbc,meqn,maux,qold_geoclaw,qold_cudaclaw,aux_geoclaw,aux_cudaclaw,geoclaw2cudaclaw)

  implicit none

  ! declare variables
  integer, intent(in) :: mx,my,mbc,meqn,maux,geoclaw2cudaclaw

   ! declare variables
  real(kind=8) :: qold_geoclaw(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
  real(kind=8) :: qold_cudaclaw(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
  real(kind=8) :: aux_geoclaw(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
  real(kind=8) :: aux_cudaclaw(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
  
  integer :: i,j,m

  ! geoclaw2cudaclaw = 1 if qold is (m,i,j) and 0 if qold is (i,j,m)
  if (geoclaw2cudaclaw == 1) then
      ! swap q(m,i,j) to q(i,j,m)
      do m=1,meqn
        do i=1-mbc,mx+mbc
          do j=1-mbc,my+mbc
            qold_cudaclaw(i,j,m) = qold_geoclaw(m,i,j)
          end do
        end do
      end do

      ! swap aux_geoclaw(m,i,j) to aux_geoclaw(i,j,m)
      do m=1,maux
        do i=1-mbc,mx+mbc
          do j=1-mbc,my+mbc
            aux_cudaclaw(i,j,m) = aux_geoclaw(m,i,j)
          end do
        end do
      end do

  else 
      ! swap q(i,j,m) to q(m,i,j)
      do m=1,meqn
        do i=1-mbc,mx+mbc
          do j=1-mbc,my+mbc
            qold_geoclaw(m,i,j) = qold_cudaclaw(i,j,m)
          end do
        end do
      end do
  end if


end subroutine cudaclaw_swap_data

