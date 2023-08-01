subroutine gstiff(coord, epstn, intgr, istep, kstep, leqns, lnods, matno, maxai, maxaj, ncrit, ndime, ndofn, nelem, ngaus, nlaps,  &
    nmats, nnode, npoin ,nstre ,ntype ,nwmtl ,nwktl ,posgp, props ,stiff ,stifi ,strsg ,tdisp, weigp)
! **********************************************************************
!
! EVALUATES GEOMETRIcallY NONLINEAR stiffNESS MATRIX
! FOR 2-D PLANE stresS/STRAIN 2-D ELEMENT
!
! **********************************************************************
    implicit none

    INCLUDE 'param.inc'

    integer :: istep, kstep, ncrit, ndime, ndofn, nelem, ngaus, nlaps, nmats, nnode,                                               &
        npoin, nstre, ntype, nwmtl, nwktl, kount, lstep, kgaus, nevab, iwktl, ielem, iposn, inode, idime, lnode, nposn,            &
        ievab, jevab, kgasp, nstr1, lprop, isize, igaus, jgaus, istr1, istre, jstre

    integer :: lnods(nelem, 1), maxai(1), intgr(1), matno(1), leqns(18, 1), maxaj(1)

    real :: djacm, deriv, shape, twopi, dispt, young, thick, hards, frict,                                                         &
        exisp, etasp, djacb, telem, dvolu, sint3, steff, theta, varj2, yield, abeta, poiss

    real :: dlcod(2,9), elcod(2, 9), estif(171), stres(4), dmatx(4,4), gpcod(2,9), bmatx(4, 18), stiff(1), posgp(1),               &
        props(mmats, mprop), stifi(1), weigp(1), avect(4), cartd(2,9), dvect(4), dbmat(4,18), devia(4),                            &
        coord(npoin,1), epstn(1), strsg(4, 1), tdisp(1)


    if(istep .eq. 1) goto 200
    kount=(lstep/kstep)*kstep
    if(kount .ne. istep) return
200 continue
    twopi=6.283185307179586
    kgaus=0

    ! Loop over each element
    nstr1=0
    nevab=ndofn*nnode
    do iwktl=1, nwktl
        stiff(iwktl) = 0.0
        stifi(iwktl) = 0.0
    end do

    do 70 ielem=1,nelem
        lprop=matno(ielem)

        ! EVALUATE THE coordINATES OF THE ELEMENT NODAL POINTS
        iposn=0
        do inode=1,nnode
            lnode=lnods( ielem, inode)
            do idime=1,ndime
                iposn=iposn+1
                nposn=leqns( iposn, ielem)
                if(nposn .eq. 0) dispt=0.
                if(nposn .ne. 0) dispt=tdisp(nposn)
                dlcod(idime, inode) = coord(lnode, idime)+dispt
                elcod(idime, inode) = coord(lnode, idime)
            end do
        end do
        young=props(lprop, 1)
        poiss=props(lprop, 2)
        thick=props(lprop, 3)
        hards=props(lprop, 7)
        frict=props(lprop, 8)

        ! Initialize the element stiffness matrix 171=nevab*(nevab+1)/2
        do isize=1,171
            estif(isize)=0.0
        end do
        kgasp=0

        ! Enter loops for area numerical integration
        do 50 igaus=1,ngaus
            exisp=posgp( igaus)
            do 55 jgaus=1,ngaus
                etasp=posgp( jgaus)
                kgasp=kgasp+1
                kgaus=kgaus+1
                call modps(dmatx, lprop, nmats, nstre, ntype, props)
                call sfr2(deriv,nnode, shape, exisp, etasp)
                call jacob2(cartd, deriv, djacb, elcod, gpcod, telem, kgasp, nnode, shape)
                call jacobd(cartd, dlcod, djacm, ndime, nlaps, nnode)
                dvolu=djacb*weigp(igaus)*weigp(jgaus)
                if(ntype .eq. 3) dvolu=dvolu*twopi*gpcod(1,kgasp)
                if(ntype .eq. 1) dvolu=dvolu*thick

                ! Evaluate the B and DB matrices
                call blarge(bmatx, cartd, djacm, dlcod, gpcod, kgasp, nlaps, nnode, ntype, shape)
                if(nlaps .eq. 2 .OR. nlaps .eq. 0) GO TO 80
                if(istep .eq. 1) GO TO 80
                if(epstn(kgaus) .eq.0.0) GO TO 80
                do istr1=1,nstr1
                    stres(istr1)=strsg(istr1,kgaus)
                end do
                call invar(devia, lprop, ncrit, nmats, props, sint3, steff, stres, theta, varj2, yield)
                call yieldf(avect, devia, frict, ncrit, sint3, steff, theta, varj2)
                call flowpl(avect, abeta, dvect, hards, ntype, poiss, young)
                do istre=1, nstre
                    do jstre=1,nstre
                        dmatx(istre, jstre)=dmatx(istre,jstre)- abeta*dvect(istre)*dvect(jstre)
                    end do
                end do
80              continue
                call dintob(bmatx, dbmat, dmatx, nevab, nstre)

                ! Evaluate geometric stiffness terms
                if(nlaps.LT.2) goto 85
                call geomst(cartd, dvolu, estif, kgaus, ndofn, nnode, strsg, shape, ntype, gpcod, kgasp)

                ! Calculate the element stiffnesses
85              kount=0
                do ievab=1,nevab
                    do jevab=ievab, nevab
                        kount=kount+1
                        do istre=1,nstre
                            estif(kount)=estif(kount)+bmatx(istre,ievab)*dbmat(istre,jevab)*dvolu
                        end do
                    end do
                end do
55          continue
50      continue

        ! GENERATES GLOBAL stiffNSS MATRIX IN COMPACTED COLUMN FORM
        if(intgr(ielem) .eq. 2) GO TO 210
        call addban(stifi, maxai, estif, leqns(1, ielem), nevab)
210     call addban(stiff, maxaj, estif, leqns(1, ielem), nevab)
70  continue
    write(6,900) (stifi(iwktl) ,iwktl=1,nwmtl)

    return

900 format(10E12.4)

end subroutine

