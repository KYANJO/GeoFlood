! Extracts geoclaw parameters to be used in CUDA Riemann solvers
subroutine get_geoclaw_parameters(mcapa_,coord_sys,grav_,dry_tol,earth_rad,deg2rad_)

    ! =================== module imports ===================
    USE amr_module, only: mcapa
    USE geoclaw_module, only:  coordinate_system, grav, dry_tolerance, earth_radius, deg2rad

    implicit none

    ! =================== Declare output variables ===================
    integer, intent(out) :: mcapa_, coord_sys
    double precision, intent(out) :: grav_, dry_tol, earth_rad, deg2rad_
    
    ! ==================== Assign scalar output variables ====================
    mcapa_ = mcapa
    coord_sys = coordinate_system
    grav_ = grav
    dry_tol = dry_tolerance
    earth_rad = earth_radius
    deg2rad_ = deg2rad
    
end subroutine get_geoclaw_parameters
