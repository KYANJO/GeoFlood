! Extracts geoclaw parameters to be used in CUDA Riemann solvers
subroutine get_geoclaw_parameters(mcapa_, coord_sys, grav_, dry_tol_, earth_rad_, deg2rad_, &
                                  theta_0_, omega_, coriolis_forcing_, friction_forcing_, &
                                  friction_depth_, variable_friction_, num_manning_, &
                                  friction_index_, manning_coefficient_, manning_break_)
    ! =================== module imports ===================
    USE amr_module, only: mcapa
    USE geoclaw_module, only: coordinate_system, grav, dry_tolerance, earth_radius, deg2rad
    USE geoclaw_module, only: theta_0, omega, coriolis_forcing, friction_forcing, friction_depth
    USE geoclaw_module, only: num_manning, manning_break, manning_coefficient
    USE friction_module, only: variable_friction, friction_index

    implicit none

    ! =================== Declare output scalar variables ===================
    integer, intent(out) :: mcapa_, coord_sys, num_manning_, friction_index_
    double precision, intent(out) :: grav_, dry_tol_, earth_rad_, deg2rad_
    double precision, intent(out) :: theta_0_, omega_, friction_depth_
    logical, intent(out) :: coriolis_forcing_, friction_forcing_, variable_friction_

    ! =================== Declare output array variables ===================
    ! double precision, dimension(:), allocatable :: manning_coefficient_, manning_break_
    double precision, intent(out) :: manning_coefficient_, manning_break_
    
    ! ==================== Assign scalar output variables ====================
    mcapa_ = mcapa
    coord_sys = coordinate_system
    grav_ = grav
    dry_tol_ = dry_tolerance
    earth_rad_ = earth_radius
    deg2rad_ = deg2rad
    theta_0_ = theta_0
    omega_ = omega
    coriolis_forcing_ = coriolis_forcing
    friction_forcing_ = friction_forcing
    friction_depth_ = friction_depth
    variable_friction_ = variable_friction
    num_manning_ = num_manning
    friction_index_ = friction_index

    ! ==================== Assign array output variables ====================
    ! allocate(manning_coefficient_(num_manning))
    ! allocate(manning_break_(num_manning))
    manning_coefficient_ = manning_coefficient(num_manning)
    ! manning_break_ = manning_break
    manning_break_ = manning_break(num_manning)

    ! The deallocation statements have been removed to allow the use of
    ! manning_coefficient_ and manning_break_ outside of this subroutine.
end subroutine get_geoclaw_parameters

! subroutine cleanup_geoclaw_parameters(manning_coefficient, manning_break)
!     implicit none
!     double precision, dimension(:), allocatable :: manning_coefficient, manning_break

!     if (allocated(manning_coefficient)) then
!         deallocate(manning_coefficient)
!     endif

!     if (allocated(manning_break)) then
!         deallocate(manning_break)
!     endif
! end subroutine cleanup_geoclaw_parameters
