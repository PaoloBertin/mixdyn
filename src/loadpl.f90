subroutine loadpl()
! **********************************************************************
!
! TIME INTEGRATION IMPLICIT-EXPLICITY ALGORITHM
!
! **********************************************************************
    use model

    implicit none

    integer :: nevab, ielem, ievab, iplod, igrav, iedge, itemp, lodpt, idofn, inode, nloca, ngash, lprop, lnode, idime, kgasp,     &
        igaus, jgaus, mgash, nedge, ncode, neass, knode, iodeg, jnode, kount, istre, kgast, jstre, ipoin, nodpt, kevab, nposn, nodeg

    integer :: noprs(4)

    real :: twopi, theta, gravy, thick, dense, gxcom, gycom, exisp, etasp, dvolu, djacb, radus, pxcom, pycom, eigen, young, poiss, &
        therm, extra, alpha

    real :: gpcod(2, 9), shape(9), deriv(2, 9), point(2), elcod(2,9), cartd(2,9), press(3,2), pgash(2), &
        dgash(2), stran(4), stres(4), dmatx(4,4)

    character*80 title

    twopi = 6.283185307179586
    nevab = nnode*ndofn
    do ielem=1, nelem
        do ievab=1, nevab
            rload(ielem, ievab) = 0.0
        end do
    end do

    read(5, 900)
    read(5, 901) title
    write(6, 903) title

    ! Read data controlling loading types to be inputed
    read(5, 919) 
    write(6, 990)
    write(6, 991) iplod, igrav, iedge, itemp

    ! Read nodal point loads
    if(iplod .ne. 0) then
        write(6, 998)
        do while(lodpt .lt. npoin)
            read (5, 931) lodpt, (point(idofn), idofn = 1, ndofn)
            write(6, 933) lodpt, (point(idofn), idofn = 1, ndofn)

            ! Associate the nodal point loads with an element
            do ielem = 1, nelem
                do inode=1,nnode
                    nloca= iabs(lnods(ielem, inode) )
                    if(lodpt .eq. nloca) goto 40
                end do
            end do

40          do idofn=1, ndofn
                ngash = (inode-1) * ndofn + idofn
                rload(ielem, ngash) = point(idofn)
            end do
        end do
    endif

    if(igrav .ne. 0) then

        ! Read gravity angle and gravitazional constant
        read(5,906) theta, gravy
        write(6, 911) theta, gravy

        theta = theta/57.295779514
        do ielem=1, nelem
            ! set up preliminary constants
            lprop = matno(ielem)
            thick = props(lprop, 3)
            dense = props(lprop, 4)
            if(dense .eq. 0.0) cycle
            gxcom = dense*gravy*sin(theta)
            gycom=-dense*gravy*cos(theta)

            ! Compute coordinates of the element nodal point
            do inode = 1, nnode
                lnode = iabs(lnods(ielem, inode) )
                do idime=1,ndime
                    elcod(idime, inode)=coord(lnode, idime)
                end do
            end do

            ! Enter loops for area numerical integration
            kgasp = 0
            do igaus=1,ngaus
                do jgaus=1,ngaus
                    kgasp = kgasp +1
                    exisp=posgp(igaus)
                    etasp=posgp(jgaus)

                    ! Compute the shape functions at  the sampling points and elemental volume
                    call sfr2(deriv, nnode, shape, exisp, etasp)
                    call jacob2(cartd, deriv, djacb, elcod, gpcod, ielem, kgasp, nnode, shape)
                    dvolu = djacb*weigp(igaus)*weigp(jgaus)
                    if(ntype .eq. 1) dvolu = dvolu*thick
                    if(ntype .eq. 3) dvolu = dvolu*twopi*gpcod(1,kgasp)

                    ! Calculate loads and associate with element nodal point
                    do inode=1, nnode
                        ngash = (inode-1)*ndofn + 1
                        mgash = (inode-1)*ndofn + 2
                        rload(ielem, ngash) = rload(ielem, ngash) + gxcom*shape(inode)*dvolu
                        rload(ielem, mgash) = rload(ielem, mgash) + gycom*shape(inode)*dvolu
                    end do
                end do
            end do
        end do
    endif

    if(iedge .ne. 0) then
        ! Distribuited edge loads section
        read(5,932) nedge
        write(6,912) nedge
        write(6, 915)
        nodeg = 3
        ncode = nnode
        if(nnode .eq. 4) nodeg = 2
        if(nnode .eq. 9) ncode = 8

        ! Loop over each loaded edge
        do iedge = 1, nedge
            ! Read data locating the loaded edge and applied load
            read(5, 902) neass, (noprs(iodeg), iodeg = 1, nodeg)
            write(6,913) neass, (noprs(iodeg), iodeg = 1, nodeg)
            read(5, 914) ((press(iodeg, idofn), iodeg = 1, nodeg), idofn = 1, ndofn)
            write(6,914) ((press(iodeg, idofn), iodeg = 1, nodeg), idofn = 1, ndofn)

            !cycle
            etasp = -1.0

            ! Calculate the coordinates of the nodes  of the the element edge
            do  iodeg = 1,nodeg
                lnode = noprs(iodeg)
                do idime = 1, ndime
                    elcod(idime, iodeg) = coord(lnode, idime)
                end do
            end do

            ! Enter loop for linear numerical integration
            do igaus = 1,ngaus
                exisp = posgp(igaus)

                ! Evaluate the shape functions at the sampling points
                call sfr2(deriv, nnode, shape, exisp,etasp)

                ! Calculate components of the equivalent nodal loads
                do idofn = 1,ndofn
                    pgash(idofn) = 0.0
                    dgash(idofn) = 0.0
                    do iodeg = 1,nodeg
                        pgash(idofn)=pgash(idofn) + press(iodeg, idofn)*shape(iodeg)
                        pgash(idofn)=dgash(IDOFN) + elcod(idofn, iodeg)*deriv(1, iodeg)
                    end do
                end do
            end do
            dvolu = weigp(igaus)
            pxcom = dgash(1)*pgash(2)-dgash(2)*pgash(1)
            pycom = dgash(1)*pgash(1)+dgash(2)*pgash(2)
            if(ntype .ne. 3) goto 115
            radus = 0.0
            do iodeg = 1, nodeg
                radus = radus+shape(iodeg)*elcod(1, iodeg)
                dvolu = dvolu*twopi*radus
            end do
115         continue

            ! Associate the equivalent nodal edge loads with an element
            do inode = 1,nnode
                nloca = iabs(lnods(neass, inode) )
                if(nloca .eq. noprs(1)) goto 130
            end do
130         jnode = inode + nodeg - 1
            kount = 0
            do knode = inode, jnode
                kount = kount + 1
                ngash = (knode - 1) * ndofn + 1
                mgash = (knode - 1) * ndofn + 2
                if(knode .gt. ncode) ngash = 1
                if(knode .gt. ncode) mgash = 2
                rload(neass, ngash) = rload(neass, ngash) + shape(kount)*pxcom*dvolu
                rload(neass, mgash) = rload(neass, mgash) + shape(kount)*pycom*dvolu
            end do
        end do
    endif

    if(itemp .ne. 0) then
        ! Initialize and input the nodal temperatures temperatures
        do ipoin=1, npoin
            tempe(ipoin)=0.0
        end do

        write(6,917)
        do while(nodpt .lt. npoin)
            read (5,916) nodpt, tempe(nodpt)
            write(6,916) nodpt, tempe(nodpt)
        end do

        kgast=0

        ! Loop over each element
        do ielem=1, nelem
            lprop = matno(ielem)
            do inode =1,nnode
                lnode = iabs(lnods(ielem, inode))

                ! Identity the coordinates and temperature of each element node point
                do idime=1,ndime
                    elcod(idime, inode) = coord(lnode, idime)
                end do
                elcod(2, inode) = tempe(lnode)
            end do

            ! Set up materials properties
            call modps(dmatx, lprop, nmats, nstre, ntype, props)
            young = props(lprop, 1)
            poiss = props(lprop, 2)
            thick = props(lprop, 3)
            alpha = props(lprop, 5)

            ! Enter loops for area numeriacal integration
            kgasp = 0
            do igaus=1, ngaus
                do jgaus=1, ngaus
                    kgast = kgast + 1
                    kgasp = kgasp + 1
                    exisp = posgp(igaus)
                    etasp = posgp(jgaus)

                    ! Evaluate the shape functions and temperature at the sampling point elemental volume and cartesian derivates
                    call sfr2(deriv, nnode, shape, exisp, etasp)
                    call jacob2(cartd, deriv, djacb, elcod, gpcod, ielem, kgasp, nnode, shape)
                    therm = 0.0
                    do inode = 1, nnode
                        therm = therm + elcod(2,inode) * shape(inode)
                    end do
                    dvolu = djacb * weigp(igaus) * weigp(jgaus)
                    if(ntype .eq. 1) dvolu = dvolu*thick
                    if(ntype .eq. 3) dvolu = dvolu*twopi*gpcod(1, kgasp)

                    ! Evaluate the initial thermal strain
                    eigen = therm*alpha
                    if(ntype .eq. 2) goto 220
                    stran(1) = eigen
                    stran(2) = eigen
                    stran(3) = 0.0
                    goto 230
220                 stran(1) = -(1.0 + poiss) * eigen
                    stran(2) = -(1.0 + poiss) * eigen
                    stran(3) = 0.0

                    ! and the corresponding initial stresses
230                 do istre = 1, nstre
                        stres(istre) = 0.0
                        do jstre=1, nstre
                            stres(istre) = stres(istre) + dmatx(istre, jstre) * stran(jstre)
                            strin(istre, kgast) = stres(istre)
                        end do
                    end do
                    if(ntype .eq. 2) strin(4, kgast) =-young*eigen
                    if(ntype .eq. 1) strin(4, kgast) = 0.0

                    ! calculate the equivalent nodal forces and associate with the element nodes
                    extra = 0.0
                    do inode = 1, nnode
                        if(ntype .eq. 3) then
                            extra = dvolu * shape(inode) * stres(4)/gpcod(1, kgasp)
                        endif
                        ngash = (inode - 1) * ndofn + 1
                        mgash = (inode - 1) * ndofn + 2
                        rload(ielem, ngash) = rload(ielem, ngash) + extra - (cartd(1, inode) * stres(1) + &
                            cartd(2, inode)*stres(3))*dvolu
                        rload(ielem, mgash) = rload(ielem, mgash) - (cartd(1, inode)*stres(3) + &
                            cartd(2, inode)*stres(2))*dvolu
                        rload(ielem, ngash) = rload(ielem, ngash) + extra - (cartd(1, inode)* stres(1) + &
                            cartd(2, inode)*stres(3))*dvolu
                        rload(ielem, mgash) = rload(ielem, mgash) + (cartd(1, inode)*stres(3) + &
                            cartd(2, inode) * stres(2)) * dvolu
                    end do
                end do
            end do
        end do
    endif

    write(6, 907)
    do ielem=1,nelem
        write(6,905) ielem, (rload(ielem, ievab) , ievab=1,nevab)
    end do

    do ielem=1,nelem
        kevab=0
        do inode = 1,nnode
            lnode = lnods(ielem , inode)
            nposn = (lnode-1) * ndofn
            do idofn = 1, ndofn
                kevab=kevab + 1
                nposn=nposn + 1
                force(nposn) = force(nposn) + rload(ielem, kevab)
            end do
        end do
    end do

    return

906 format (2F10.3)
900 format()
901 format(A80)
903 format(//'LOAD CASE TITLE: ', A80)
911 format('GRAVITY ANGLE =', F10.3, 'GRAVITY CONSTANT',F10.3)
913 format(i10, 5X, 4i5)
919 format(4I5)
990 format(//'LOAD INPUT PARAMETERS')
991 format('POINT LOADS = ' ,I5/, 'GRAVITY     = ', I5/ , 'HEDGE LOADS = ', I5/, 'TEMPERATURE = ', I5)
998 format(/5X ,' NODE' , '       PX' , '      PY'/)
931 format(I5, 2F10.3)
933 format(I5, 2F10.3)
932 format(I5)
912 format(5X, 'NO. OF LOADED EDGES = ', I5)
915 format('LIST OF LOADED EDGES AND APPLIED LOADS')
902 format(4I5)
914 format(6F10.3)
907 format('TOTAL NODAL FORCES FOR EACH ELEMENT')
917 format(5X, 'PRESCRIBED NODAL TEMPERATURES')
916 format(I5, F10.3)
905 format(1X,I4,5X,8E12.4/(10X, 8E12.4) )

end subroutine

