
! Extracts dt from b4step2 to be used in bc2, since there is a bug not yet fixed in the bc2 code for dt
module extract_dt
    implicit none
    real(kind=8) :: dt_extract
end module extract_dt