subroutine addban(stiff, maxai, estif, leqns, nevab)
! **********************************************************************
!
! *** aasembly of total stiffness vector
!
!  DA TESTARE !!!!!
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: nevab, ievab, kount, ieqns, kevab, jevab, imaxa, jeqns, ijeqn, isize, jsize
    integer :: maxai(1), leqns(1)

    !real ::
    real :: stiff(1), estif(1)

    kount = 0

    do ievab = 1,nevab
        ieqns = leqns(ievab)
        if(ieqns .lt. 0) then
            kount = kount + nevab - ievab
        else if(ieqns .eq. 0 ) then
            kevab=kevab +nevab - jevab
        else
            imaxa=maxai(ieqns)
        endif

        kevab = ievab
        do jevab=1,nevab
            jeqns = leqns (jevab)
!         if(jeqns) 220,220,110
!110      ijeqn=ieqns-jeqns
!         if(ijeqn) 220,210,210
!210      isize = imaxa+ijeqn
            jsize = kevab
            if(jevab .ge. ievab) jsize = jevab+KOUNT
            stiff(isize)=stiff(isize) + estif(jsize)
!220      kevab = kevab + nevab - jevab
        end do
    end do

    return
end subroutine


