subroutine prevos(force, ndofn, nelem, ngaus, npoin, nprev, strin)
! **********************************************************************
!
! GRAVITY LOADS AND STRESSES
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ndofn, nelem, ngaus, npoin, nprev, nstr1, ngau2, ngash, nposn, ielem, igaus, kgaus, istri

    real :: force(1), strin(4, 1), xgash, ygash

    if(nprev .eq. 0) return

    nstr1 = 4
    ngau2 = ngaus*ngaus

! READ GRAVITY LOADS
    write(6, 920)
    do while(ngash .ne. npoin)
        read(5, 930) ngash, xgash, ygash
        nposn = (ngash-1)*ndofn +1
        force(nposn) = xgash
        nposn = nposn + 1
        force(nposn) = ygash
        write(6, 940) ngash, xgash, ygash
    end do

! READ GRAVITY STRESS
    write(6, 950)
    do ielem =1, NELEM
        do igaus =1, ngau2
            read(5,  930) kgaus, (strin(istri ,kgaus) , istri=1, nstr1)
            write(6, 930) kgaus, (STRIN(ISTRI ,KGAUS) , istri=1, nstr1)
        end do
    end do

    return

! 900   format()
920 format(' NODE', ' GRAVITY X-LOAD', ' GRAVITY Y-LOAD')
930 format( /I5, 4F10.3)
940 format(I5, F10.3, 10X, F10.3)
950 format('GAUSS PT.' , 'GRAVITY X-STRESS', 'GRAVITY Y-STRESS', 'GRAVITY.XY-STRESS', 'GRAVITY Z-STRESS')

end subroutine
