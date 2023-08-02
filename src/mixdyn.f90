program mixdyn
! **********************************************************************
!
! TIME INTEGRATION IMPLICIT-EXPLICITY ALGORITHM
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: ndofn, nelem, nmats, npoin, nconm, ncrit, ndime, ngaum, ngaus, nlaps, nnode, nprev, ntype, nstre, ifixd, ifunc,     &
        kstep, miter, noutd, noutp, nreqd, nreqs, nstep, ipred, neqns, nwktl, nwmtl, istep, iiter, nchek

    integer :: ifpre(mdofn*mpoin), lnods(melem, mnode), matno(melem), maxai(mpoin*mdime), maxaj(mpoin*mdime), mhigh(mpoin*mdime),  &
        intgr(melem), ngrqs(mpoin), nprqd(mpoin), leqns(mevab, mpoin), niter(2000)

    real :: aalfa, afact, azero, beeta, bzero, delta, dtime, dtend, gaama, omega, toler, consd, consf

    real :: coord(mpoin, mdime), props(mmats, mprop), dispi(mpoin*mdime), veloi(mpoin*mdime), accei(mpoin*mdime),                  &
        displ(mpoin*mdime), ymass(mpoin*mdime), accek(mpoin*mdime), accej(mpoin*mdime), resid(mpoin*mdime), dispt(mpoin*mdime),    &
        dispq(mpoin*mdime), tempe(mpoin), force(mpoin*mdime), acceh(mpoin*mdime), velot(mpoin*mdime), velol(mpoin*mdime),          &
        accel(mpoin*mdime), accev(mpoin*mdime), stiff(maxma), stifs(maxma), stifi(maxma), xmass(maxma), dampi(maxma),              &
        dampg(maxma), rload(melem, mnode*mdime), strin(4, 450), strag(4, 450), strsg(4, 450), epstn(450), effst(450), posgp(4),    &
        weigp(4)

    common stiff, xmass, dampg, stifi, stifs, dampi

!  Open the file I/O
    open(5, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/input.dat',  status='OLD')
    open(6, file='/home/paolo-bertin/Documenti/FortranProjets/mixdyn/data/result.dat', status='OLD')

    call contol(ndofn, nelem, nmats, npoin)

    call inputd(coord, ifpre, lnods, matno, nconm, ncrit, ndime, ndofn, nelem, ngaum, ngaus, nlaps, nmats, nnode, npoin, nprev,    &
        nstre, ntype, posgp, props, weigp)

    call intime(aalfa, acceh, accev, afact, azero, beeta, bzero, delta, dtime, dtend, gaama, ifixd, ifunc, intgr, kstep, miter,    &
        ndofn, nelem, ngrqs, noutd, noutp, npoin, nprqd, nreqd, nreqs, nstep, omega, dispi, toler, veloi, ipred)

    call prevos(force, ndofn, nelem, ngaus, npoin, nprev, strin)

    call loadpl(coord, force, lnods, matno, ndime, ndofn, nelem, ngaus, nmats, nnode, npoin, nstre, ntype, posgp, props, rload,    &
        strin, tempe, weigp)

    call lumass(coord, intgr, lnods, matno, NCONM, ndime, ndofn, nelem, NGAUM, nmats, nnode, npoin, ntype, props, ymass)

    call linkin(force, ifpre, intgr, leqns, lnods, maxai, maxaj, MHIGH, ndofn, nelem, neqns, nnode, npoin, nwktl, nwmtl, xmass,    &
        ymass)

    do istep=1, nstep
        do iiter=1, miter
            call gstiff(coord, epstn, intgr, istep, kstep, leqns, lnods, matno, maxai, maxaj, ncrit, ndime, ndofn, nelem, ngaus,   &
                nlaps, nmats, nnode, npoin, nstre, ntype, nwmtl, nwktl, posgp, props, stiff, stifi, strsg, dispt, weigp)

            call impexp(aalfa, acceh, accei, accej, accek, accel, accev, afact, azero, beeta, bzero, consd, consf, dampi, dampg,   &
                dispt, dtend, dtime, gaama, ifixd, ifpre, delta, dispi, displ, ifunc, iiter, istep, kstep, maxai, maxaj, ndofn,    &
                neqns, npoin, nwktl, nwmtl, omega, force, stiff, stifi, stifs, veloi, velol, velot, xmass, ymass, ipred)

            call resepl(coord, dispt, effst, rload, epstn, iiter, intgr, leqns, lnods, matno ,ncrit, ndime, ndofn, nelem, ngaus,   &
                nlaps, nmats, nnode, npoin, nstre, ntype, posgp, props, resid, STRAG, strin, strsg, weigp, ipred, istep)

            call itrate(accei, accel, consd, consf, xmass, dispi, displ, dispt, maxai, nchek, neqns, nwmtl, resid, stifs, toler,   &
                veloi, velol, velot, iiter ,miter)

            if(nchek.EQ.1) goto 520
        end do

520     call outdyn(dispq, dtime, epstn, ifpre, iiter, istep, ndofn, nelem, ngaus, ngrqs, niter, noutd, noutp, npoin, nprqd, nreqd,&
            nreqs, ntype, strsg, dispi)
    end do

!   close file I/O
    close(5)
    close(6)

    stop
end
