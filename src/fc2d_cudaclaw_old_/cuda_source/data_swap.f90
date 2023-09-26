! @author: Brian Kyanjo
! @date: 2023/13/08
! @purpose: This program reads in a (m,i,j) or (i,j,m) array and swaps the indices
!          to the other form. This is useful for exchanging data between geoclaw (m,i,j) 
!          routines and cudaclaw (i,j,m)routines

subroutine data_swap(mx,my,mbc,meqn,maux,qold,qold_transpose,aux,aux_transpose,flag2transpose)

  implicit none

  ! declare variables
  integer, intent(in) :: mx,my,mbc,meqn,maux,flag2transpose

  integer :: i,j,m

  ! flag2transpose = 1 if qold is (m,i,j) and 0 if qold is (i,j,m)
  if (flag2transpose == 1) then
    
    ! declare variables
    real(kind=8), intent(in) :: qold(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(out) :: qold_transpose(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    real(kind=8), intent(in) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(out) :: aux_transpose(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
  
    ! swap q(m,i,j) to q(i,j,m)
    do m=1,meqn
      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          qold_transpose(i,j,m) = qold(m,i,j)
        end do
      end do
    end do

    ! swap aux(m,i,j) to aux(i,j,m)
    do m=1,maux
      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          aux_transpose(i,j,m) = aux(m,i,j)
        end do
      end do
    end do

  else 
    ! declare variables
    real(kind=8), intent(in) :: qold(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
    real(kind=8), intent(out) :: qold_transpose(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in) :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)
    real(kind=8), intent(out) :: aux_transpose(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! swap q(i,j,m) to q(m,i,j)
    do m=1,meqn
      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          qold_transpose(m,i,j) = qold(i,j,m)
        end do
      end do
    end do

    ! swap aux(i,j,m) to aux(m,i,j)
    do m=1,maux
      do i=1-mbc,mx+mbc
        do j=1-mbc,my+mbc
          aux_transpose(m,i,j) = aux(i,j,m)
        end do
      end do
    end do
  end if
end subroutine data_swap

