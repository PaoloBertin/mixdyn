subroutine contol()
! **********************************************************************
!
!     READ CONTROL DATA AND CHECK FOR DIMENSIONS
!
! **********************************************************************
    use model

    implicit none

    !include  'param.inc'


    read(5, 110) npoin, nelem, ndofn, nmats

    if(nelem .gt. melem .or. npoin .gt. mpoin .or. nmats .gt. mmats) then
        write(66, 120)
        stop
    endif

    return

110 format(4I5)
120 format(/'SET DIMENSION EXCEEDED - CONTOL CHECK '/)

end subroutine
