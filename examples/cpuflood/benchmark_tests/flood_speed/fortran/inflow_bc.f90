PROGRAM linear_interpolation

! This program interploates linearly the inflow values applied as boundary 1 ! ! conditions along 20m line in the middle of the western side of the floodplain.

IMPLICIT NONE

! Declare variables
integer :: n, i
real :: inflow_
real, dimension(:), allocatable :: time, inflow


! Read the inflow values from the file
open(10,file = "../scratch/BC.csv")
read(10,*) n ! number of inflow values
allocate(time(n),inflow(n))
do i = 1, n
    read(10,*) time(i), inflow(i)
end do
close(10)

! linerar interpolation of inflow values
do i = 1, n-1
    inflow_ = inflow(i) + (inflow(i+1)-inflow(i)) * (time(i+2)-time(i)) / (time(i+1)-time(i))
end do
