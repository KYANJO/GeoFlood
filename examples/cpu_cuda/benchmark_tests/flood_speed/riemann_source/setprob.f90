! @desc: Writes geflood parameters to setprob.data file to be read and
!        accessible by the Riemann solvers in the setprob_cuda module
!        in the same order.

subroutine flood_speed_setprob(grav, drytol, earth_radius, coord_system, mcapa)

    implicit none

    double precision :: grav
    double precision :: drytol
    double precision :: earth_radius
    integer :: coord_system
    integer :: mcapa

    ! open setprob.data file
    open(10, file='setprob.data', status='unknown')
    read(10,*) grav
    read(10,*) drytol
    read(10,*) earth_radius
    read(10,*) coord_system
    read(10,*) mcapa
    close(10)

    return

end subroutine flood_speed_setprob