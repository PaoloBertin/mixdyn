subroutine blarge(bmatx, cartd, djacm, dlcod, gpcod, kgasp, nlaps, nnode, ntype, shape)
    !
    ! Large displacement B matrix
    !

    implicit none

    include 'param.inc'

    integer :: ngash, inode, jnode, nnode, ntype, nlaps, kgasp

    real :: fmult
    real :: cartd(2, 9), bmatx(4, 18), djacm(2, 2), dlcod(2, 9), gpcod(2, 9), shape(9)

    ngash=0
    do inode=1,nnode
        nnode=ngash+1
        ngash=nnode+1
        bmatx(1, nnode) = cartd(1, inode) * djacm(1, 1)
        bmatx(1, ngash) = cartd(1, inode) * djacm(2, 1)
        bmatx(2, nnode) = cartd(2, inode) * djacm(1, 2)
        bmatx(2, ngash) = cartd(2, inode) * djacm(2, 2)
        bmatx(3, nnode) = cartd(2, inode) * djacm(1, 1) + cartd(1, inode) * djacm(1, 2)
        bmatx(3, ngash) = cartd(1, inode) * djacm(2, 2) + cartd(2, inode) * djacm(2, 1)
    end do

    if(ntype .ne. 3) return
    fmult=1.
    if(nlaps .lt. 2) goto 40
    fmult=0.0

    do jnode=1,nnode
        fmult=fmult+DLCOD(1, jnode)*shape(jnode)
    end do 
    fmult=fmult/GPCOD(1, KGASP)

40  ngash=0
    do inode=1,nnode
        nnode=ngash+1
        ngash=nnode+1
        bmatx (4, nnode) =shape( inode) *fmult/GPCOD( 1, KGASP)
        bmatx(4,ngash)=0.0
    end do

    return
end subroutine



