module model
    implicit none
    include 'param.inc'

    integer :: npoin, ndime, ndofn, nelem, nnode, ngaus, ngaum, nprop, ntype, nvfix, nstre, ncrit, nprev, nconm, nlaps, nrads
    integer :: lnods(melem, mnode), matno(melem), ifpre(mdofn, mpoin)
    real :: coord(mpoin, mdime), props(mmats, mprop)

end module model