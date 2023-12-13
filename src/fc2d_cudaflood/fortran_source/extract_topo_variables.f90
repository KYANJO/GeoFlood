! Extracts b4step2 parameters to be used in C++ code
subroutine get_b4step2_parameters(num_dtopo_val, aux_finalized_val, t0dtopo_val, tfdtopo_val, &
                                  dt_max_dtopo_val, NEEDS_TO_BE_SET_val, variable_friction_val, &
                                  friction_index_val, xupper_val, yupper_val, xlower_val, ylower_val, &
                                  test_topography_val, mtopofiles_val, topowork_val, &
                                  xlowtopo_val, ylowtopo_val, xhitopo_val, yhitopo_val, dxtopo_val, &
                                    dytopo_val, mxtopo_val, mytopo_val, mtopoorder_val, i0topo_val, &
                                    mtopo_val,mtoposize_val)


    ! =================== module imports ===================
    USE topo_module
    USE amr_module, ONLY: xupper, yupper, xlower, ylower, NEEDS_TO_BE_SET
    USE friction_module, ONLY: variable_friction, friction_index
    USE friction_module, ONLY: set_friction_field

    implicit none

    ! =================== Declare output variables ===================
    integer, intent(out) :: num_dtopo_val, test_topography_val
    integer, intent(out) :: aux_finalized_val, mtopofiles_val, mtoposize_val
    real(kind=8), intent(out) :: dt_max_dtopo_val
    real(kind=8), intent(out), allocatable :: t0dtopo_val(:), tfdtopo_val(:)
    real(kind=8), intent(out) :: NEEDS_TO_BE_SET_val
    logical, intent(out) :: variable_friction_val
    integer, intent(out) :: friction_index_val
    real(kind=8), intent(out) :: xupper_val, yupper_val, xlower_val, ylower_val
    real(kind=8), intent(out), allocatable :: topowork_val(:), xlowtopo_val(:), ylowtopo_val(:) 
     real(kind=8), intent(out), allocatable :: xhitopo_val(:), yhitopo_val(:), dxtopo_val(:), dytopo_val(:)
    integer, intent(out), allocatable :: mxtopo_val(:), mytopo_val(:), mtopoorder_val(:), i0topo_val(:), mtopo_val(:)

    ! ==================== Assign scalar output variables ====================
    ! ====== topo_variables
    num_dtopo_val = num_dtopo
    aux_finalized_val = aux_finalized
    dt_max_dtopo_val = dt_max_dtopo
    test_topography_val = test_topography
    mtopofiles_val = mtopofiles
    mtoposize_val = mtoposize
    
    ! ======= friction_variables
    variable_friction_val = variable_friction
    friction_index_val = friction_index

    ! ======= amr_variables (domain description variables) 
    xupper_val = xupper
    yupper_val = yupper
    xlower_val = xlower
    ylower_val = ylower
    NEEDS_TO_BE_SET_val = NEEDS_TO_BE_SET


    ! ==================== Allocate and assign array output variables ========
    if (num_dtopo > 0) then
        allocate(t0dtopo_val(num_dtopo), tfdtopo_val(num_dtopo))
        t0dtopo_val = t0dtopo(1:num_dtopo)
        tfdtopo_val = tfdtopo(1:num_dtopo)
    endif

    if (mtopofiles > 0) then
        allocate(topowork_val(mtopofiles), xlowtopo_val(mtopofiles), ylowtopo_val(mtopofiles), xhitopo_val(mtopofiles), &
        yhitopo_val(mtopofiles), dxtopo_val(mtopofiles), dytopo_val(mtopofiles))
        allocate(mxtopo_val(mtopofiles), & 
        mytopo_val(mtopofiles), mtopoorder_val(mtopofiles), i0topo_val(mtopofiles), mtopo_val(mtopofiles))
        topowork_val = topowork(1:mtopofiles)
        xlowtopo_val = xlowtopo(1:mtopofiles)
        ylowtopo_val = ylowtopo(1:mtopofiles)
        xhitopo_val = xhitopo(1:mtopofiles)
        yhitopo_val = yhitopo(1:mtopofiles)
        dxtopo_val = dxtopo(1:mtopofiles)
        dytopo_val = dytopo(1:mtopofiles)
        mxtopo_val = mxtopo(1:mtopofiles)
        mytopo_val = mytopo(1:mtopofiles)
        mtopoorder_val = mtopoorder(1:mtopofiles)
        i0topo_val = i0topo(1:mtopofiles)
        mtopo_val = mtopo(1:mtopofiles)
    endif

end subroutine get_b4step2_parameters
