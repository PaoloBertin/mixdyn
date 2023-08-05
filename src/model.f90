module model
    implicit none

    integer, parameter :: mpoin=200, mdime=2, melem=50, mnode=9, mdofn=2, mmats=10, mstre=4, mprop=13, mgaus=3, mgau2=4, mgau3=9,  &
        mevab=mnode*mdofn, maxma=6000

    integer :: numel, npoin, ndime, ndofn, nelem, nnode, ngaus, ngaum, nmats, nprop, ntype, nvfix, nstre, ncrit, nprev, nconm,     &
        nlaps, nrads, neqns

    integer :: nstep, noutd, noutp, nreqd, nreqs, nacce, ifunc, ifixd, miter, kstep, ipred, nwktl, nwmtl

    integer :: lnods(melem, mnode), matno(melem), ifpre(mdofn, mpoin)
    integer :: intgr(melem), nprqd(mnode), ngrqs(melem*mgaus), leqns(mnode*mdofn, melem)
    integer :: maxai(mpoin*mdime), maxaj(mpoin*mdime), mhigh(mpoin*mdime)

    real :: coord(mpoin, mdime), props(mmats, mprop)
    real :: posgp(4), weigp(4)
    real :: dtime, dtend, dtrec, aalfa, beeta, delta, gaama, azero, bzero, omega, toler, afact
    real :: force(mpoin*mdofn), dispt(mpoin*mdime), dispi(mpoin*mdime), veloi(mpoin*mdime), accei(mpoin*mdime), accej(mpoin*mdime)
    real :: tempe(mpoin)
    real :: accek(mpoin*mdime), acceh(600), accev(600), strin(mgau2, melem),rload(melem, mnode*mdime), ymass(mpoin*mdime)
    real :: accel(mpoin*mdime), velol(mpoin*mdime), displ(mpoin*mdime), stifs(maxma), velot(mpoin*mdime)
    real :: stiff(maxma), stifi(maxma), xmass(maxma), dampi(maxma), dampg(maxma), strsg(mstre, melem*mgaus)
    real :: epstn(melem*mgaus), effst(melem*mgaus), resid(mpoin*mevab), strag(mstre, melem*mgau3)
end module model
