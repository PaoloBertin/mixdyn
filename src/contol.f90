subroutine contol(ndofn, nelem, nmats, npoin)
! **********************************************************************
!
!     READ CONTROL DATA AND CHECK FOR DIMENSIONS
!
! **********************************************************************
    implicit none

    include  'param.inc'

    integer :: ndofn, nelem, nmats, npoin

    read(5, 110) npoin, nelem, ndofn, nmats

    if(nelem .gt. melem .or. npoin .gt. mpoin .or. nmats .gt. mmats) then
        write(66, 120)
        stop
    endif

    return

110 format(4I5)
120 format(/'SET DIMENSION EXCEEDED - CONTOL CHECK '/)

end subroutine
