! Extracts b4step2 parameters to be used in C++ code
subroutine get_b4step2_parameters(num_dtopo_val, aux_finalized_val, t0dtopo_val, tfdtopo_val, dt_max_dtopo_val, NEEDS_TO_BE_SET_val)
    USE topo_module, ONLY: num_dtopo, aux_finalized, t0dtopo, tfdtopo, dt_max_dtopo
    USE amr_module, ONLY: NEEDS_TO_BE_SET

    implicit none

    ! Declare output variables
    integer, intent(out) :: num_dtopo_val
    integer, intent(out) :: aux_finalized_val
    real(kind=8), intent(out) :: dt_max_dtopo_val
    real(kind=8), intent(out), allocatable :: t0dtopo_val(:), tfdtopo_val(:)
    real(kind=8), intent(out) :: NEEDS_TO_BE_SET_val

    ! Assign scalar output variables
    num_dtopo_val = num_dtopo
    aux_finalized_val = aux_finalized
    dt_max_dtopo_val = dt_max_dtopo
    NEEDS_TO_BE_SET_val = NEEDS_TO_BE_SET

    ! Allocate and assign array output variables
    if (num_dtopo > 0) then
        allocate(t0dtopo_val(num_dtopo), tfdtopo_val(num_dtopo))
        t0dtopo_val = t0dtopo(1:num_dtopo)
        tfdtopo_val = tfdtopo(1:num_dtopo)
    endif

end subroutine get_b4step2_parameters
