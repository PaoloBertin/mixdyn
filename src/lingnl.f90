subroutine lingnl(cartd, djacm, dmatx, eldis, gpcod, kgasp, kgaus, ndofn, nlaps, nnode, nstre, ntype, poiss, shape, stran, stres,  &
    strin)
    !
    ! Elastic strain and stresses
    !
    real :: cartd(2,9) ,stran(4), dmatx(4,4) ,strin(4,1), eldis(2,9), stres(4), djacm(2,2), agash(2,2), gpcod(2,9) ,shape(9)

    ! Calculate strains from deformation jacobian

    if(nlaps .lt. 2) goto 15
    stran(1)=0.5 * (djacm(1, 1) * djacm(1, 1) + djacm(2, 1) * djacm(2, 1) - 1.0)
    stran(2)=0.5 * (djacm(1,2) * djacm(1,2) + djacm(2,2) * djacm(2,2) -1.0)
    stran(3) = djacm(1, 1) * djacm( 1,2) + djacm(2,1) * djacm(2, 2)

    ! For small displacements
    goto 25
15  continue
    do idofn=1, ndofn
        do jdofn=1, ndofn
            bgash=0.0
            do  inode=1, nnode
                bgash = bgash + cartd(jdofn, inode)*eldis(idofn,inode)
            end do
        end do
        agash(idofn, jdofn) =bgash
    end do

    stran(1)=agash(1, 1)
    stran(2)=agash(2,2)
    stran(3)=agash(1,2)+agash(2, 1)
25  continue
    if(ntype.LT.3) GO TO 90
    stran(4)=0.0
    do inode=1,nnode
        stran(4)=stran(4)+eldis(1,inode)*shape(inode)/gpcod(1,kgasp)
    end do
    extra=0.0
    do inode=1,nnode
        extra= extra+eldis(1,inode)*shape(inode)/gpcod(1,kgasp)
    end do
    stran(4)=stran(4)+0.5*extra*extra
90  do istre=1,4
        stran(istre)=stran( istre)-strin(istre,kgaus)
    end do

    ! and the corresponding stresses
    do  istre=1,nstre
        stres(istre)=0.0
        do  jstre=1,nstre
            stres(istre) = stres(istre)+dmatx(istre, jstre) *stran(jstre)
        end do
    end do

    if(ntype .eq. 1) stres(4) = 0.0
    if(ntype .eq. 2) stres(4) = poiss * (stres(1)+stres(2))

    return

end subroutine lingnl

