subroutine lumass()
! **********************************************************************
!
! CALCULATES LUMPED MASS FOR 4, 8 AND 9 NODED ELEMENT
!
! **********************************************************************
    use model

    implicit none

    integer :: nevab, ntotv, itotv, ielem, ievab, imass, kgasp, lprop, inode, lnode, idime, igaus, jgaus, kount, jnode,jevab, &
        iposn, iconm, ipoin, nposn, idofn

    real :: twopi, tarea, thick, rhoel, exisp, etasp, djacb, dvolu, shapi, shapt, dmass, sumas, xcmas, ycmas

    real :: elcod(2, 9), diagm(9), cartd(2, 9), shape(9), gpcod(2,9), emass(171), deriv(2,9)

    rewind 3

    twopi=6.283185307179586
    nevab=nnode*ndofn
    ntotv=npoin*ndofn
    do itotv =1, ntotv
        ymass(itotv)=0.0
    end do

    call gaussq(ngaum, posgp, weigp)
    do  ielem=1, nelem
        do ievab=1,171
            emass(ievab)=0.0
        end do
        imass=intgr(ielem)
        kgasp=0
        tarea=0.0
        lprop=matno(ielem)
        thick=props(lprop, 3)
        rhoel=props(lprop, 4)
        do inode =1,nnode
            diagm(inode)=0.0
            lnode=lnods(ielem, inode)
            do idime=1, ndime
                elcod(idime, inode)=coord(lnode, idime)
            end do
        end do
        do igaus=1,ngaum
            exisp=POSGP(IGAUS)
            do jgaus=1,ngaum
                kgasp=kgasp+1
                etasp=posgp(jgaus)
                call sfr2(deriv, nnode, shape, exisp, etasp)
                call jacob2(cartd, deriv, djacb, elcod, gpcod, ielem, kgasp, nnode, shape)
                dvolu=djacb*weigp(igaus)*weigp(jgaus)
                if(ntype .eq. 1) dvolu=dvolu*thick
                IF(ntype .eq. 3) dvolu=dvolu*twopi*gpcod(1, kgasp)
                IF(imass .eq. 1) goto 210
                DO inode =1, nnode
                    shapt = shape(inode)
                    diagm(inode)=diagm(inode) + shapt*shapi*dvolu
                end do
                tarea=tarea+dvolu
210             if(imass .eq. 2) cycle
                dvolu=dvolu*rhoel
                ievab=1
                kount=nevab
                do inode=1, nnode
                    shapi=shape(inode)
                    do jnode=inode, nnode
                        dmass=dvolu*shapi*shape(jnode)
                        emass(ievab)=emass(ievab) + dmass
                        jevab = ievab + kount
                        emass(jevab)=emass(jevab)+dmass
                        ievab=ievab+2
                    end do
                    kount=kount-2
                    ievab=jevab+1
                end do
            end do
        end do

        ! Writes consistent mass matrix on tape 3
        if(imass .eq. 2) goto 200
        write(3) emass
        ! write(6, 90) (emass(i), i = 1, 171)
200     if(imass .eq. 1) GO TO 100
        ! Generates lumped mass matrix proportional to diagonal
        sumas=0.0
        do inode=1, nnode
            sumas=sumas + diagm(inode)
        end do
        tarea=tarea*rhoel
        sumas=tarea/sumas
        do inode=1, nnode
            lnode = lnods(ielem, inode)
            iposn = (lnode-1)*ndofn
            do idofn =1, ndofn
                IPOSN=IPOSN+1
                YMASS(IPOSN)=YMASS(IPOSN)+DIAGM(INODE)*SUMAS
            end do
        end do
!90      format(ex, 9e12.3)
    end do
100 continue

    ! Cncentrated masses
    if(nconm .eq. 0) return
    write(6, 900)
    do iconm=1, nconm
        read(5,910) ipoin, xcmas, ycmas
        write(6,910) ipoin, xcmas, ycmas
        nposn=(ipoin-1)*ndofn+1
        ymass(nposn)=ymass(nposn)+xcmas
        nposn = nposn+1
        ymass(nposn)=ymass(nposn)+ycmas
    end do

!    write(6, 90) (ymass(i), i= 1, ntotv)

    return

900 format(5X, 'CONCENTRATED MASSES')
910 format(I5,2F10.3)

end

