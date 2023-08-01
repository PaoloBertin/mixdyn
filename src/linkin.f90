subroutine linkin(force, ifpre, intgr, leqns, lnods, maxai, maxaj, mhigh, ndofn, nelem, neqns, nnode, npoin, nwktl, nwmtl, xmass,  &
    ymass)
! **********************************************************************
!
! LINKS WITH PROFILE SOLVER
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ndofn, nelem, neqns, nnode, npoin, nwktl, nwmtl, imass, nevab, ipoin, idofn, ielem, meqns, ievab, inode, ident,     &
        ieqns, mkoun, nposm, nposn

    integer :: lnods(nelem,1), ifpre(ndofn, 1), maxai(1), maxaj(1), intgr(1)

    !   real ::
    real :: force(1), emass(171), ymass(1), xmass(1)
    real :: mhigh(1), leqns(18, 1)

    imass=1

    rewind 3

    nevab=nnode*ndofn

    ! number of unknowns
    neqns=0
    do ipoin=1,npoin
        do idofn=1,ndofn
            ! IF(ifpre(idofn,ipoin)) 110,120,110
            if(ifpre(idofn, ipoin) .lt. 0) then
                ifpre(idofn, ipoin) = 0
            else if(ifpre(idofn,ipoin) .GT. 0) then
                ifpre(idofn, ipoin) = 0
            else
                neqns=neqns+1
                ifpre(idofn, ipoin) = neqns
            endif
            ! 120       neqns=neqns+1
            !          ifpre(idofn, ipoin) = neqns
            ! GO TO 150
            ! 110       ifpre(idofn, ipoin)=0
        end do

        write(6,7) ipoin, (ifpre( idofn, ipoin), idofn=1,ndofn)
    end do

    meqns = 1 + neqns

    ! connectivity array leqns
    do ielem=1,nelem
        do ievab=1,nevab
            leqns( ievab, ielem) = 0
        end do
    end do

    do ielem=1,nelem
        ievab=1
        do inode=1,nnode
            ident=lnods(ielem, inode)
            do idofn=1,ndofn
                leqns(ievab, ielem) = ifpre(idofn, ident)
                ievab=ievab+1
            end do
        end do
        ! write(6,6) ielem, (leqns( ievab, ielem) , ievab=1,nevab)
    end do

    ! LOOP OVER ALL ELEMENTS
250 do ielem=1,nelem
        if(intgr(ielem) .ne. imass) cycle
        call colmht(mhigh, nevab, leqns(1, ielem))
    end do

    ! addresses of diagonal elements — maxa array
    call address(maxaj, mhigh, neqns, nwktl, mkoun)
    if(imass .eq. 2) goto 205
    do ieqns=1, meqns
        maxai(ieqns)=maxaj(ieqns)
    end do

    imass=2
    nwmtl=nwktl
    goto 250

205 continue
    write(6,920) neqns, nwmtl, nwktl
    write(6,930) (maxai(ieqns), ieqns=1,meqns)
    write(6,930) (maxaj(ieqns), ieqns=1,meqns)
    if(nwktl .GT. 6000) goto 210
    goto 220
210 write(6,910)
    stop
220 continue

    ! global mass matrix
    DO ielem=1,nelem
        imass=intgr( ielem)
        if(imass .eq. 2) cycle
        read(3) emass
        call addban(xmass,maxai,emass, leqns(1, ielem) , nevab)
    end do

    ! global mass vector
    nposm=0
    do ipoin =1,npoin
        do idofn =1,ndofn
            nposm=nposm+1
            nposn=ifpre(idofn, ipoin)
            IF(nposn.EQ.0) cycle
            ymass(nposn)=ymass(nposm)
            force(nposn)=force(nposm)
        end do
    end do
    return

7   format(4I10)
930 format(5X, 20I5)
920 format(/5X, 'neqns=',I5,5X, 'nwmtl=',I5,5X, 'nwktl=', 15/)
910 format (/'SET DIMENSION EXCEEDED - CHECK LINKIN '/)

END
