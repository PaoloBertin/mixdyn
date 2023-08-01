subroutine inputd(ifpre, lnods, matno, nconm, ncrit, ndime, ndofn, nelem, ngaum, ngaus, nlaps, nmats, nnode, npoin, nprev, &
    nstre, ntype, posgp, props, weigp)
! **********************************************************************
!
!  MIXDYN INPUT ROUTINE
!
! **********************************************************************
    implicit none

    include 'param.inc'

    integer :: nconm, ncrit, ndime, ndofn, nelem, ngaum, ngaus, nlaps, nmats, nnode, npoin, nprev, nstre, ntype, nvfix, nprop,     &
        nrads, ielem, numel, inode, ipoin, idime, idofn, ivfix, imats, imate, iprop
    integer :: lnods(melem,mnode), ifpre(mdofn,mpoin), matno(melem)

    real :: coord(mpoin,mdime), props(mmats, mprop), posgp(4), weigp(4)

    character*80 title

    read(5,  913) title
    write(6, 914) title

! READ THE FIRST DATA CARD, AND ECHO 1T IMMEDIATELY
    read(5, 900) nvfix, ntype, nnode, nprop, ngaus, ndime, nstre, ncrit, nprev, nconm, nlaps, ngaum, nrads
    write(6, 901) npoin, nelem, nvfix, ntype, nnode, ndofn, nmats, nprop, ngaus, ndime, nstre, ncrit, nprev, nconm, nlaps, ngaum,  &
        nrads

! READ THE ELEMENT NODAL CONNECTIONS, AND THE PROPERTY NUMBERS.
    read(5, 911)
    write(6, 902)
    do ielem = 1,nelem
        read(5, 900) numel, matno(numel), (lnods(numel, inode), inode=1, nnode)
        write(6, 903) numel, matno(numel), (lnods(numel, inode), inode=1, nnode)
    end do

! ZERO ALL THE NODAL COORD1NATES, PRIOR TO READING SOME OF THEM
    do IPOIN=1, MPOIN
        do IDIME=1, MDIME
            coord(ipoin, idime) = 0.0
        end do
    end do

! READ SOME COORDINATES, FINISHING WITH THE LAST NODE OF ALL
    read(5, 911)
    do while(ipoin .ne. npoin)
        read(5, 905) ipoin, (coord(ipoin, idime), idime=1, ndime)
        write(6, 905)ipoin, (coord(ipoin, idime), idime=1, ndime)
    end do

! INTERPOLATE COORDINATES OF MID-SIDE NODES
    call nodxyr(coord, lnods, nelem, nnode, npoin, nrads, ntype)

! PRINT NODE COORDINATES
    write(6, 904)
    do ipoin =1, npoin
        write(6, 905) ipoin, (coord(ipoin, idime), idime=1, ndime)
    end do

! READ THE FIXED VALUES
    write(6, 907)
    do IPOIN = 1, npoin
        do IDOFN=1, ndofn
            ifpre(idofn, ipoin)=0
        end do
    end do

    read(5, 911)
    do ivfix=1, nvfix
        read(5,  908) ipoin, (ifpre(idofn, ipoin), idofn=1, ndofn)
        write(6, 909) ipoin, (ifpre(idofn, ipoin), idofn=1, ndofn)
    end do

! READ THE AVAILABLE SELECTION OF ELEMENT PROPERTIES.
    write(6, 910)
    read (5, 911)
    do imats=1, nmats
        read(5, 917)  imate, (props(imate, iprop), iprop=1, nprop)
        write(6, 912) imate, (props(imate, iprop), iprop=1, npoin)
    end do

! SET UP GAUSSIAN 1NTEGRATION CONSTANTS
    call gaussq(ngaus, posgp, weigp)

    return

913 format(A80)
914 format(A80/ )
900 format(16I5)
901 format('CONTROL PARAMETERS '/, &
        'NPOIN =', I5, 5X, 'NELEM =', I5, 5X, 'NVF1X =', I5/ &
        'NTYPE =', I5, 5X, 'NNODE =', I5, 5X, 'NDOFN =', I5/ &
        'NMATS =', I5, 5X, 'NPROP =', I5, 5X, 'NGAUS =', I5/ &
        'NDIME =', I5, 5X, 'NSTRE =', I5, 5X, 'NCR1T =', I5/ &
        'NPREV =', I5, 5X, 'NCONM =', I5, 5X, 'NLAPS =', I5/ &
        'NGAUM =', I5, 5X, 'NRADS =', I5)
911 format(I5)
902 format(//3X, 'ELEMENT', 2X, 'PROPERTY', 15X, 'NODE NUMBERS')
903 format(2I10, 5X, 9I5)
904 format(//' NODE', 9X, 'X', 9X, 'Y')
905 format(I5, 6F10.5)
907 format(//' NODE', 1X, 'CODE')
908 format(I5, 3X, 2I1)
909 format(I5, 3X ,2I1)
910 format(//'MATERIAL PROPERTIES')
917 format(I5, 11E10.3)
912 format('MATERIAL NO ',I5/, &
        'YOUNG MODULUS = ', E10.3,/ &
        'POISSON RATIO = ', E10.3,/ &
        'THICKNESS     = ', E10.3,/ &
        'MASS DENSITY  = ', E10.3,/ &
        'ALPHA TEMPR   = ', E10.3,/ &
        'REFERENCE FO  = ', E10.3,/ &
        'HARNING PAR   = ', E10.3,/ &
        'FRICT ANGLE   = ', E10.3,/ &
        'FLUIDITY PAR  = ', E10.3,/ &
        'EXP delta     = ', E10.3,/ &
        'FLOW CODE     = ', E10.3)

end subroutine
