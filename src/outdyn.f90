subroutine outdyn(iiter, istep, tdisp, vivel)
!
! Output routine
!
    use model

    implicit none

    integer :: iiter, ngaut
    integer :: istep, nstr1, kount, koutd, nposn, nposm, ipoin, ireqd, ikoun, igaus, istr1, ireqs, nodet, idofn, ngasi, ngasj,     &
        ngask, mgasj, mgask, jpoin, kpoin, igasi, mgasi, igasj, igask, kgaus, ielem, kelgs, jgaus, istre

    real :: ttime, xtime, xgash, xgish, xgesh, xgosh
    real :: strsp(3), vivel(5,1), tdisp(1)

    nstr1 = 4
    kstep = istep
    ngaut = nelem*ngaus*ngaus
    if(istep .eq. 1) write( 10,925)

    ttime = ttime + dtime

    ! Writes displacment history at requested nodal points on tape 10
    ! and stress history at requested gauss points at every noutd steps:
    kount = 0
    koutd = (istep/noutd)*noutd
    if(koutd .ne.istep) goto 510

    do ipoin = 1, npoin
        do ireqd =1, nreqd
            if(ipoin .ne. nprqd(ireqd)) cycle
            nposn=(ipoin - 1) *ndofn+1
            nposm=nposn+1
            kount=kount + 1
            displ(kount) =tdisp( nposn)
            kount=kount+1
            displ(kount) =tdisp( nposm)
        end do
    end do
    write(10, 960) (displ(ikoun) , ikoun=1,kount) , ttime

    do igaus = 1, ngaut
        do ireqs = 1, nreqs
            if(igaus .ne. ngrqs(ireqs)) cycle
            write(11,950) (strsg(istr1,igaus), istr1=1,nstr1)
        end do
    end do

510 koutd=(kstep/noutp ) *noutp
    if(koutd .ne. kstep) return
    xtime = float(kstep) * dtime
    write(6,608) kstep, xtime

    ! Rearrange displacement vector
    nodet = 0
    do ipoin = 1, npoin
        do idofn =1, ndofn
            nodet = nodet +1
            displ(nodet) = tdisp(nodet)
        end do
    end do

    ! Qutput displacements
    write(6, 990)
    do ipoin = 1, npoin, 3
        ngasi = ndofn * ipoin-1
        ngasj = ngasi + ndofn
        ngask = ngasj + ndofn
        ngasi = ngasi + 1
        mgasj = ngasj + 1
        mgask = ngask + 1
        jpoin = ipoin + 1
        kpoin = jpoin + 1

        ! Write displacements on the tape 13 for deformation plot
        write(13,910) ipoin ,(displ(igasi) , igasi=ngasi, mgasi)
        if(jpoin .gt. npoin) cycle
        write(13,910) jpoin ,(displ(igasj), igasj=ngasj, mgasj)
        if(kpoin .gt. npoin) cycle
        write(13,910) kpoin, (displ(igask), igask=ngask, mgask)

        ! WriteS displAcements on output file
        write(6,920) ipoin,displ(ngasi) ,displ(ngasi), jpoin,displ(ngasj) ,displ(mgasj), kpoin, displ( ngask) , displ( mgask)
    end do

    ! Write stresses on output file
    write(6, 900)
    if(ntype .ne. 3) write(6,970)
    if(ntype .eq. 3) write(6,975)

    kgaus=-0
    do ielem = 1, nelem
        kelgs=0
        write(6,930) ielem

        do igaus =1, ngaus
            do jgaus = 1, ngaus
                kgaus = kgaus + 1
                kelgs = kelgs + 1
                xgash = (strsg(1, kgaus) + strsg(2, kgaus)) * 0.5
                xgish = (strsg(1, kgaus) - strsg(2, kgaus)) * 0.5
                xgesh = strsg(3, kgaus)
                xgosh=sqrt(xgish * xgish+xgesh*xgesh )
                strsp(1) = xgash + xgosh
                strsp(2) = xgash - xgosh
                if(xgish .eq. 0.0) xgish = 0.0
                strsp( 3) = atan(xgesh/xgish)*28.647889757

                ! Writes complete stress state on tape 4
                write(4,950) (strsg(istr1, kgaus), istr1 = 1, nstr1), (strsp(istre), istre=1, 3)
                write(6,940) kelgs, (strsg(istr1,kgaus) , istr1=1,nstr1), (strsp(istre) , istre=1, 3) , vivel(5, kgaus)
            end do
        end do
    end do

    return

!980 format(1X,60I2)
960 format(1X, 10E11.4)
950 format (7E10.4)
940 format(I5,2X,6E14.6,F8.3,E14.6)
900 format(/, 10X, 8HSTRESSES ,/)
920 format(3(1X,I5,2E12.5))
910 format(I5,2E15.6)
608 format(//5X,28H Displacements at TIME STEP ,I10,5X,SHTIME ,E20.11)
925 format(5X,' Displacements')
990 format( /3(1X, 'NNODE' , 3X, 'X-DISP',6X, 'Y-DISP', 3X) / )
970 format(1X,'G.P.',6X, 'XX-STRESS' , 5X, 'YY-STRESS' , 5X, 'XY-STRESS', 5X, 'ZZ-STRESS', 6X, 'MAX P.S.', 6X, 'MIN P.S.',3X,      &
        'ANGLE', 3X, 'P.S.')
975 format('G.P.', 6X, 'RR-STRESS' , 5X, 'ZZ-STRESS' ,5X, 'RZ-STRESS', 5X, 'SHTT-STRESS', 6X, 'MAX P.S.',6X,'MIN P.S.',3X,'ANGLE', &
        3X,'P.S.')
930 format(5X, 'ELEMENT NO. =',I5)

end subroutine

