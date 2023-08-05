subroutine resepl (eload, iiter, istep)
    !
    ! Evaluates residual forces
    !
    use model

    implicit none

    integer :: iiter, istep, idofn, nevab, ntotv, nstr1, ielem, ievab, kgaus, lprop, itotv, iposn, inode, lnode, idime, nposn,     &
        kgasp, igaus, jgaus, istr1, mstep, jstep, istre, mgash, lmveb

    real :: frict, twopi, young, poiss, thick, uniax, hards, exisp, etasp, djacb, dvolu, preys, sint3, steff, theta, varj2, yield, &
        espre, escur, rfact, astep, reduc, abeta, agash, dlamd, bgash, curys, bring, dispv

    real :: deriv(2, 9), dmatx(4, 4), avect(4), dlcod(2, 9), eload(nelem, 1), cartd(2, 9), shape(9), desig(4), dvect(4), stran(4), &
        devia(4), bmatx(4, 18), gpcod(2, 9), djacm(2, 2), stres(4), elcod(2,9), sigma(4), sgtot(4), eldis(2, 9)

    twopi = 6.283185307179586
    nevab = nnode * ndofn
    ntotv = npoin * ndofn
    nstr1 = 4
    do ielem = 1, nelem
        if(intgr(ielem) .eq. 2 .and. iiter .gt. 1 .and. ipred .eq. 1) cycle
        do ievab = 1, nevab
            eload(ielem, ievab) = 0.0
        end do
    end do

    do itotv =1, ntotv
        resid(itotv)=-0.0
    end do

    kgaus = 0
    do ielem =1, nelem
        if(intgr(ielem) .eq. 2 .and. iiter .gt. 1 .AND. ipred .eq.1) cycle
        lprop = matno(ielem)
        young=props(lprop, 1)
        poiss=props(lprop, 2)
        thick=props(lprop, 3)
        uniax=props(lprop, 6)
        hards=props(lprop, 7)
        frict=props(lprop, 8)
        frict=frict*0.017453292
        if(ncrit .eq. 3) uniax = uniax*cos(frict)
        if(ncrit .eq. 4) uniax = 6.0 * uniax * cos(frict)/(1.73205080757*(3.0-sin(frict) ))

        ! Compute coordinate and incremental displacements of the element nodal points
        iposn = 0
        do inode = 1, nnode
            lnode = lnods(ielem, inode)
            do idime = 1, ndime
                iposn = iposn + 1
                nposn = leqns(iposn, ielem)
                if(nposn .eq. 0) dispv = 0.0
                if(nposn .ne. 0) dispv = displ(nposn)
                dlcod(idime, inode) = coord(lnode, idime) + dispv
                elcod(idime, inode) = coord(lnode, idime)
                eldis(idime, inode) = dispv
            end do
        end do

        ! matrice legame costitutivo
        call modps(dmatx,lprop, nmats ,nstre, ntype, props)
        kgasp=0
        do igaus=1,ngaus
            do jgaus=1,ngaus
                exisp=posgp(igaus)
                etasp=posgp(jgaus)
                kgaus=kgaus+1
                kgasp=kgasp+1
                call sfr2(deriv, nnode, shape, exisp ,etasp)
                call jacob2(cartd, deriv, djacb, elcod, gpcod, ielem, kgasp, nnode, shape)
                call jacobd(cartd,dlcod , DJACM, ndime , nlaps, nnode)
                dvolu=djacb*weigp(igaus)*weigp(jgaus)
                if(ntype .eq. 3) dvolu=dvolu*twopi*gpcod(1,kgasp)
                if(ntype .eq. 1) dvolu=dvolu*thick
                call blarge (bmatx, cartd, djacb, dlcod ,GPCOD, kgasp , nlaps, nnode, ntype, shape)
                call lingnl(cartd, djacm, dmatx, eldis, gpcod, kgasp, kgaus , ndofn, nlaps, nnode, nstre, ntype, poiss , shape,    &
                    stran, stres, strag)
                do istr1 =1, nstr1
                    strag(istr1, kgaus) = strag(istr1 ,kgaus) + stran(istr1)
                end do

                if(istep .gt. 1 .and. iiter .gt.1) goto 160
                do istr1 = 1,nstr1
                    stres(istr1) = stres(istr1) + strin(istr1, kgaus)
                end do
160             continue
                preys = uniax + epstn(kgaus) * hards
                do istr1=1, nstr1
                    desig(istr1) = stres(istr1)
                    sigma(istr1) = strsg(istr1,kgaus) + stres(istr1)
                end do

                if(nlaps .eq. 2 .or. nlaps .eq. 0) goto 60
                call invar(devia, lprop, ncrit, nmats, props, sint3, steff, sigma, theta, varj2, yield)
                espre = effst(kgaus) - preys
                if(espre .ge. 0.0) goto 50
                escur = yield-preys
                IF(escur .le. 0.0) goto 60
                rfact=escur/(yield-effst(kgaus))
                goto 70
50              escur=yield-effst(kgaus)
                if(escur .le. 0.0) goto 60
                rfact=1.0
70              mstep = escur*8.0/uniax+1.0
                if(mstep .gt. 10) mstep=10
                astep=mstep
                reduc=1.0-rfact
                do istr1=1,nstr1
                    sgtot(istr1)=strsg(istr1 , kgaus) + reduc*stres(istr1)
                    stres(istr1)=rfact*stres(istr1)/astep
                end do
                DO jstep=1,mstep
                    call invar(devia, lprop, ncrit, nmats, props, sint3, steff, sgtot , theta, varj2, yield)
                    call yieldf(avect, devia,frict,ncrit,sint3,steff, theta, varj2)
                    call flowpl(avect, abeta, dvect, hards, ntype, poiss, young)
                    agash=0.0
                    DO istr1=1,nstr1
                        agash=agash+avect(istr1)*stres(istr1)
                    end do
                    dlamd = agash*abeta
                    if(dlamd .lt. 0.0) dlamd=0.0
                    bgash=0.0
                    do istr1=1,nstr1
                        bgash=bgash+avect(istr1) *sgtot( istr1)
                        sgtot(istr1)=sgtot ( istr1)+stres(istr1) -dlamd*dvect(istr1)
                    end do
                    epstn(kgaus) = epstn(kgaus) + dlamd*bgash/ yield
                end do

                call invar(devia,lprop,ncrit,nmats, props, sint3,steff, sgtot , theta, varj2, yield)

                curys = uniax + epstn(kgaus)*hards
                bring=1.0
                if(yield .gt. curys) bring = curys/yield
                DO istr1=1,nstr1
                    strsg(istr1,kgaus)=bring*sgtot (istr1)
                end do
                effst(kgaus) =bring*yield

                ! Alternative location of stress reduction loop termination card
                ! 90 CONTINUE
                goto 190
60              DO istr1 = 1, nstr1
                    strsg(istr1,kgaus)=strsg(istr1,kgaus)+DESIG(istr1)
                end do
                effst(kgaus)= yield

                ! Calculate the equivalent nodal forces and associate with the element nodes
190             mgash = 0
                do inode=1,nnode
                    do idofn=1,ndofn
                        mgash=mgash+1
                        do istre =1, nstre
                            eload(ielem, mgash) = eload(ielem, mgash) + bmatx(ISTRE, mgash)*strsg(istre , kgaus)*dvolu
                        end do
                    end do
                end do
            end do
        end do
    end do

    do ielem =1, nelem
        do ievab =1, nevab
            lmveb = leqns(ievab, ielem)
            if(lmveb .ne. 0) then
                resid(lmveb) = resid(lmveb) + eload(ielem, ievab)
            endif
        end do
    end do

    return

end subroutine


