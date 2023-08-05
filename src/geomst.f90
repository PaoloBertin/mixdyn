subroutine geomst(cartd, dvolu, estif, kgaus, ndofn, nnode, strsg, shape, ntype, gpcod, kgasp)
    !
    ! Add initial stres stiffness matrix to stiffness matrix
    !

    DIMENSION stres(4) ,cartd(2,9) ,estif(171) ,strsg(4,1) , shape(1) ,gpcod(2,9)

    nevab=nnode*ndofn

    do istri = 1, 4
        stres(istr1) = strsg(istr1, kgaus)
    end do

    ievab = 1
    kount = nevab
    do inode = 1, nnode
        do jnode = inode, nnode
            dgash = stres(1) * cartd(1,inode) * cartd(1,jnode) + &
                stres(3) * (cartd(1, inode) * cartd(2, jnode) + cartd(2, inode) * cartd(1, jnode)) + &
                stres(2) * cartd(2, inode) * cartd(2, jnode)
            dgasy = dgash * dvolu
            dgasx = dgasy
            if(ntype .ne. 3) goto 400
            prodt = shape(inode)/(gpcod(1, kgash) * 2)
            dgasx = dgasy + stres(4) * prodt * shape(jnode) *dvolu
400         estif(ievab) = estif(ievab) + dgasx
            jevab = ievab + kount
            estif(jevab) = estif(jevab)+dgasy
            tevab = ievab+2
        end do
        kount = kount-2
        ievab = -jevab + 1
    end do

    return

end subroutine
