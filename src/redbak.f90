subroutine redbak(stiff, force, maxai, neqns )
    !
    ! To reduce and back-substitute iteration vectors
    !
    implicit none

    integer :: neqns, ieqns, lower, kuper, jeqns, icolm, kmaxa, keqns
    integer :: maxai(1)

    real :: sumcc
    real :: stiff(1), force(1)

    do 400 ieqns = 1, neqns
        lower = maxai(ieqns) + 1
        kuper=maxai( ieqns+1) - 1
        if (kuper-lower) 400,410,410
410     jeqns = ieqns
        sumcc = 0.0
        do 420 icolm=lower, kuper
            jeqns=jeqns-1
420     sumcc = sumcc + stiff(icolm)*force(jeqns)
        force(ieqns) =force( IEQNS) -sumcc
400 continue

    do 480 ieqns=1,neqns
        kmaxa = maxai(ieqns)
480 force( ieqns) =force( ieqns)/stiff(kmaxa)

    if(neqns .eq. 1) return

    jeqns=neqns

    do 500 ieqns=2,neqns
        lower= maxai(jeqns) +1
        kuper = maxai(jeqns+1) -1
        if(kuper-lower) 500,510,510
510     keqns=jeqns
        do 520 icolm=lower, kuper
            keqns=keqns - 1
520     force(keqns) = force(keqns) - stiff( icolm)*force( jeqns)
500 jeqns=jeqns-1

    return

end subroutine
