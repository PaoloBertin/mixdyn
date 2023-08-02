subroutine multpy(final ,amatx ,start ,maxai ,neqns ,nwmtl)
    !
    ! TO EVALUATE PRODUCT OF B TIMES RR AND STORE RESULT IN TT
    !
    implicit none

    include 'param.inc'

    integer :: neqns, nwmtl, ieqns, lower, kuper, jeqns, icolm
    integer :: maxai(1)

    real :: termi, sumaa
    real :: final(melem*mdofn), amatx(melem*mdofn), start(melem*mdofn)

    if(nwmtl .gt. neqns) goto 20
    do ieqns=1,neqns
        final(ieqns) = amatx(ieqns)*start( ieqns)
    end do
    return

20  do ieqns=1,neqns
        final(ieqns) = 0.0
    end do

    do ieqns=1, neqns
        lower = maxai(ieqns)
        kuper = maxai(ieqns+1) -1
        jeqns = ieqns + 1
        termi = start(ieqns)
        do icolm=lower, kuper
            jeqns = jeqns - 1
            final(jeqns)=final(jeqns) + amatx(icolm) * termi
        end do
    end do

    if(neqns .eq. 1) return
    do ieqns=2,neqns
        lower=maxai(ieqns) +1
        kuper=maxai(ieqns+1) -1
        if((kuper-lower) .lt. 0) cycle
        jeqns=ieqns
        sumaa=0.0
        do icolm=lower , kuper
            jeqns=jeqns-1
            sumaa=sumaa+amatx(icolm)*start(jeqns)
        end do
        final (ieqns)=final(ieqns)+sumaa
    end do

    return

end subroutine
