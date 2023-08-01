subroutine colmht(mhigh, nevab, leqns)
! **********************************************************************
!
! EVALUATES THE COLUMN HEIGHT OF STIFFNESS MATRIX
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: maxam, ievab, ieqns, jhigh, nevab
    integer :: leqns(1), mhigh(1)

    maxam=100000

    do ievab=1, nevab
!C        IF(LEQNS(IEVAB)) 110,100,110
!C110     IF(LEQNS(IEVAB)-MAXAM) 120,100,100
!C120     MAXAM=LEQNS(IEVAB)
        if(leqns(ievab) .lt. 0 .or. leqns(ievab) .gt. 0) then
            if(leqns(ievab) - maxam .lt. 0) then
                maxam = leqns(ievab)
            endif
        endif
    end do

    do ievab=1, nevab
        ieqns = leqns(ievab)
        if(ieqns .eq. 0) cycle
        jhigh= ieqns -maxam
        if(jhigh .gt. mhigh(ieqns)) mhigh(ieqns) = jhigh
    end do

    return
end subroutine colmht

