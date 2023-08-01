subroutine dintob(bmatx ,dbmat ,dmatx ,nevab ,nstre)
    !
    ! CALCULATE D INTO B
    !
    implicit none

    integer :: nevab, nstre, istre, ievab, jstre
    real :: bmatx(4, 18), dbmat(4, 18), dmatx(4,4)

    do istre=1, nstre
        DO ievab=1,nevab
            dbmat(istre, ievab) =0.0
            DO jstre=1,nstre
                dbmat(istre, ievab) = dbmat( istre, ievab)+dmatx(istre, jstre)*bmatx(jstre, ievab)
            end do
        end do
    end do

    return
end subroutine

