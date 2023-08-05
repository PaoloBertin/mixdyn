subroutine inputd()
! **********************************************************************
!
!  MIXDYN INPUT ROUTINE
!
! **********************************************************************
    use model

    implicit none

    ! include 'param.inc'

    integer :: ielem, inode, ipoin, idime, idofn, ivfix, imats, imate, iprop

    character*80 title

    read(5,  913) title
    write(6, 914) title

    ! Read the first data card, and echo it immediately
    read(5, 900) nvfix, ntype, nnode, nprop, ngaus, ndime, nstre, ncrit, nprev, nconm, nlaps, ngaum, nrads
    write(6, 901) npoin, nelem, nvfix, ntype, nnode, ndofn, nmats, nprop, ngaus, ndime, nstre, ncrit, nprev, nconm, nlaps, ngaum,  &
        nrads

    ! Read the element nodal connections, and the property numbers.
    read(5, 911)
    write(6, 902)
    do ielem = 1,nelem
        read(5, 900) numel, matno(numel), (lnods(numel, inode), inode=1, nnode)
        write(6, 903) numel, matno(numel), (lnods(numel, inode), inode=1, nnode)
    end do

    ! Zero all the nodal coordinates, prior to reading some of them
    do ipoin=1, mpoin
        do idime=1, mdime
            coord(ipoin, idime) = 0.0
        end do
    end do

    ! Read some coordinates, finishing with the last node of all
    read(5, 911)
    do while(ipoin .ne. npoin)
        read(5, 905) ipoin, (coord(ipoin, idime), idime=1, ndime)
        write(6, 905)ipoin, (coord(ipoin, idime), idime=1, ndime)
    end do

    ! Interpolate coordinates of mid-side nodes
    call nodxyr(coord, lnods, nelem, nnode, npoin, nrads, ntype)

    ! Print node coordinates
    write(6, 904)
    do ipoin =1, npoin
        write(6, 905) ipoin, (coord(ipoin, idime), idime=1, ndime)
    end do

    ! Read the fixed values
    write(6, 907)
    do ipoin = 1, npoin
        do IDOFN=1, ndofn
            ifpre(idofn, ipoin)=0
        end do
    end do

    read(5, 911)
    do ivfix=1, nvfix
        read(5,  908) ipoin, (ifpre(idofn, ipoin), idofn=1, ndofn)
        write(6, 909) ipoin, (ifpre(idofn, ipoin), idofn=1, ndofn)
    end do

    ! Read the  available selection of the element properties.
    write(6, 910)
    read (5, 911)
    do imats=1, nmats
        read(5, 917)  imate, (props(imate, iprop), iprop=1, nprop)
        write(6, 912) imate, (props(imate, iprop), iprop=1, nprop)
    end do

    ! Set up gaussian integration constants
    call gaussq()

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
