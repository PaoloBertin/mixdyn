subroutine address(maxai, mhigh, neqns, nwktl, mkoun)
! **********************************************************************
!
!     EVALUATES ADRESSES OF DIAGONAL ELEMENTS
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: neqns, nwktl, mkoun, neqnn, ieqnn, ieqns
    integer :: maxai(mpoin*mdime), mhigh(mpoin*mdime)

    neqnn = neqns+1
    do ieqnn =1,neqnn
        maxai(ieqnn)=1
    end do

    maxai(2) = 2
    mkoun = 0

    if(neqns .eq. 1) then
        mkoun = mkoun + 1
        nwktl = maxai(neqns + 1) - maxai(1)
        return
    endif

    do ieqns = 2,neqns
        if(mhigh(ieqns) .gt. mkoun) mkoun= mhigh(ieqns)
        maxai(ieqns +1) = maxai(ieqns) + mhigh(ieqns) +1
    end do

    return
end subroutine

