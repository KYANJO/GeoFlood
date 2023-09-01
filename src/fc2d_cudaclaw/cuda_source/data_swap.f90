program test
  implicit none
  real, allocatable :: q(:,:,:)
  integer :: i,j,m

  interface
    subroutine swap_mij_ijm(q)
      implicit none
      real, intent(inout) :: q(:,:,:)
    end subroutine swap_mij_ijm
  end interface

  allocate(q(4,4,4))
  
  ! read qin from file
  open(unit=10,file='qin.txt',status='old',action='read')
  do m=1,4
    do i=1,4
      do j=1,4
        read(10,*) q(i,j,m)
      end do
    end do
  end do
  close(10)
  write (*,*) q(1,1,1)



end program test

! swap q(m,i,j) to q(i,j,m)
subroutine swap_mij_ijm(q)
  implicit none
  real, intent(inout) :: q(:,:,:)
  real, allocatable :: q_tmp(:,:,:)
  integer :: i,j,m

  allocate(q_tmp(size(q,1),size(q,2),size(q,3)))

  do m=1,size(q,3)
    do i=1,size(q,1)
      do j=1,size(q,2)
        q_tmp(i,j,m) = q(m,i,j)
      end do
    end do
  end do

  q = q_tmp

end subroutine swap_mij_ijm

! swap q(i,j,m) to q(m,i,j)
subroutine swap_ijm_mij(q)
  implicit none
  real, intent(inout) :: q(:,:,:)
  real, allocatable :: q_tmp(:,:,:)
  integer :: i,j,m

  allocate(q_tmp(size(q,3),size(q,1),size(q,2)))

  do m=1,size(q,3)
    do i=1,size(q,1)
      do j=1,size(q,2)
        q_tmp(m,i,j) = q(i,j,m)
      end do
    end do
  end do

  q = q_tmp

end subroutine swap_ijm_mij